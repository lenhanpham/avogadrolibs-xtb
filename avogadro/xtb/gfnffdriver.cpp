/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_eg.f90 (gfnff_eg, 0d non-periodic path),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#include "gfnffdriver.h"

#include "constants.h"
#include "environment.h"
#include "gffgraph.h"
#include "gfnffdata.h"
#include "gfnffdisp.h"
#include "gfnffegbond.h"
#include "gfnffegnonbond.h"
#include "gfnffhb.h"
#include "gfnffparams.h"
#include "gfnffsetup.h"
#include "gfnfftopo.h"

#include <cmath>
#include <utility>
#include <vector>

namespace Avogadro {
namespace Xtb {

bool gfnffSinglePoint(const DriverInput& input, const GffData& param,
                      const HbParams& hbpar, DriverResult& result,
                      Environment& env)
{
  int n = static_cast<int>(input.numbers.size());
  if (n <= 0 || static_cast<int>(input.xyz.size()) != 3 * n ||
      static_cast<int>(input.qa.size()) != n ||
      static_cast<int>(input.chieeq.size()) != n) {
    env.error("empty driver setup", "gfnffSinglePoint");
    return false;
  }
  std::vector<double> g(3 * n, 0.0);
  double sigma[3][3] = { { 0 } };
  result.ees = result.edisp = result.erep = result.ebond = 0.0;
  result.eangl = result.etors = result.ehb = result.exb = result.ebatm = 0.0;

  double dispThr, cnThr, repThr, hbThr1, hbThr2;
  gffThresholds(input.accuracy, dispThr, cnThr, repThr, hbThr1, hbThr2);

  // Packed distances + D3 list (verbatim driver loop, incl. diagonal init).
  int npair = n * (n + 1) / 2;
  std::vector<double> sqrab(npair, 0.0), srab(npair, 0.0);
  std::vector<std::pair<int, int>> d3pairs;
  for (int i = 0; i < n; ++i) {
    int ijbase = i * (i + 1) / 2;
    for (int j = 0; j < i; ++j) {
      int k = ijbase + j;
      double dx = input.xyz[3 * i] - input.xyz[3 * j];
      double dy = input.xyz[3 * i + 1] - input.xyz[3 * j + 1];
      double dz = input.xyz[3 * i + 2] - input.xyz[3 * j + 2];
      double r2 = dx * dx + dy * dy + dz * dz;
      sqrab[k] = r2;
      if (r2 < dispThr)
        d3pairs.emplace_back(i, j);
      srab[k] = std::sqrt(r2);
    }
    sqrab[ijbase + i] = 0.0;
    srab[ijbase + i] = 0.0;
  }

  // erf CN (0d dlogcoord).
  std::vector<double> rcov(103);
  for (int z = 0; z < 103; ++z)
    rcov[z] = covalentRadiusD3(z + 1);
  std::vector<double> cn, dlogCn;
  gffCoordinationNumber(n, input.numbers, input.xyz, srab, rcov, 4.4, cnThr,
                        cn, dlogCn);

  // EEQ charges over the setup fragments (single fragment by default).
  std::vector<int> frags(n, 1);
  std::vector<double> qfrag(1, 0.0);
  int nfrag = 1;
  if (static_cast<int>(input.fragOfAtom.size()) == n) {
    frags = input.fragOfAtom;
    nfrag = 0;
    for (int f : frags) {
      if (f > nfrag)
        nfrag = f;
    }
    if (nfrag < 1)
      nfrag = 1;
    qfrag.assign(nfrag, 0.0);
    for (size_t f = 0;
         f < input.fragCharges.size() && f < qfrag.size(); ++f)
      qfrag[f] = input.fragCharges[f];
  }
  std::vector<double> charges, gamma, erf;
  double ees = 0.0;
  if (!eeqCharges(n, input.numbers, srab, input.chieeq, input.gameeq,
                  input.alpeeq, cn, param.cnf, nfrag, frags, qfrag,
                  input.totalCharge, charges, gamma, erf, ees, env))
    return false;
  result.ees = ees;

  // D3(BJ) 0d.
  DispModel model;
  if (!buildDispersionModel(input.numbers, model, env))
    return false;
  std::vector<double> zetac6;
  if (!buildZetaC6(input.numbers, input.qa, zetac6, env))
    return false;
  {
    std::vector<double> gd(3 * n, 0.0);
    double edisp = 0.0;
    if (!d3Gradient(n, input.numbers, input.xyz, d3pairs, model, zetac6,
                    param.d3r0, param.dispScale, cn, dlogCn, edisp, gd, env))
      return false;
    result.edisp = edisp;
    for (int k = 0; k < 3 * n; ++k)
      g[k] += gd[k];
  }

  // ES gradient.
  {
    std::vector<double> ge(3 * n, 0.0);
    if (!esGradient(n, input.xyz, srab, sqrab, gamma, erf, charges, cn,
                    param.cnf, input.numbers, dlogCn, ge, env))
      return false;
    for (int k = 0; k < 3 * n; ++k)
      g[k] += ge[k];
  }

  // Non-bonded repulsion (0d, single cell).
  {
    std::vector<double> gr(3 * n, 0.0);
    double erep = 0.0;
    double sr[3][3] = { { 0 } };
    if (!repulsionEnergyGradient(
          n, input.numbers, input.xyz, input.nonbonded.alphanb,
          input.nonbonded.bpair, param.repz, param.repScaleN, repThr,
          input.mcfNrep, erep, gr, sr, env))
      return false;
    result.erep = erep;
    for (int k = 0; k < 3 * n; ++k)
      g[k] += gr[k];
    for (int r = 0; r < 3; ++r) {
      for (int c = 0; c < 3; ++c)
        sigma[r][c] += sr[r][c];
    }
  }

  // Bonded part: rab estimates, bond loop, bonded repulsion.
  int nbond = static_cast<int>(input.bonds.size());
  if (nbond > 0) {
    std::vector<std::pair<int, int>> blist;
    std::vector<double> rabIn, vb2, vb3;
    for (const auto& b : input.bonds) {
      blist.emplace_back(b.first, b.second);
      rabIn.push_back(b.shift);
      vb2.push_back(b.steepness);
      vb3.push_back(b.prefactor);
    }
    std::vector<double> rabOut, grabOut, rabdcnOut;
    if (!bondLengthEstimates(n, input.numbers, cn, dlogCn, blist, rabIn,
                             rabOut, grabOut, rabdcnOut, env))
      return false;
    std::vector<double> dEdcn(n, 0.0);
    double ebond = 0.0;
    // HB erf CN, gated on the per-bond HB counts like dncoord_erf.
    std::vector<double> hbCn, hbDcn;
    int nrHbSum = 0;
    for (int v : input.nrHb)
      nrHbSum += v;
    if (nrHbSum > 0) {
      if (!hbErfCoordination(n, input.numbers, input.xyz, rcov, input.bondHb,
                             hbCn, hbDcn, env))
        return false;
    }
    std::vector<char> consideredAbh(
      static_cast<size_t>(input.bondHb.mapNab) * input.bondHb.mapNab *
        input.bondHb.mapNh,
      0);
    for (int k = 0; k < nbond; ++k) {
      int jat = blist[k].first, iat = blist[k].second;
      double d[3] = { input.xyz[3 * iat] - input.xyz[3 * jat],
                      input.xyz[3 * iat + 1] - input.xyz[3 * jat + 1],
                      input.xyz[3 * iat + 2] - input.xyz[3 * jat + 2] };
      double rab = norm2(d);
      double rij = rabOut[k];
      std::vector<double> drij(3 * n);
      for (int m = 0; m < n; ++m) {
        for (int c = 0; c < 3; ++c)
          drij[3 * m + c] = grabOut[(k * n + m) * 3 + c];
      }
      double drijdcn[2] = { rabdcnOut[2 * k], rabdcnOut[2 * k + 1] };
      // Bonds in A-H...B groups take the CN-weakened egbond_hb path.
      if (k < static_cast<int>(input.nrHb.size()) && input.nrHb[k] >= 1) {
        if (!egbondHb(k, iat, jat, rab, rij, drij, drijdcn, input.numbers,
                      hbCn, hbDcn, input.xyz, param.vbondScale, vb2[k], vb3[k],
                      input.bondHb, consideredAbh, dEdcn, ebond, g, env))
          return false;
        continue;
      }
      bondEnergyGradient(k, iat, jat, rab, rij, drij, drijdcn, n, input.xyz,
                         rabIn, vb2, vb3, ebond, g, dEdcn, sigma, env);
    }
    result.ebond = ebond;
    // NOTE: the dEdcn->sigma contraction reads uninitialized dcndL on the
    // 0d reference path; it is skipped (gradient/energy-neutral).
    {
      std::vector<double> grb(3 * n, 0.0);
      double erepb = 0.0;
      double srb[3][3] = { { 0 } };
      if (!bondedRepulsion(n, input.numbers, input.xyz, blist, param.repa,
                           param.repz, param.repScaleB, erepb, grb, srb, env))
        return false;
      result.erep += erepb;
      for (int k = 0; k < 3 * n; ++k)
        g[k] += grb[k];
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c)
          sigma[r][c] += srb[r][c];
      }
    }
  }

  // Bends.
  {
    double eangl = 0.0;
    for (const auto& bd : input.bends) {
      double ag[3][3] = { { 0 } };
      double ds[3][3] = { { 0 } };
      double e = 0.0;
      bendEnergyGradient(bd.center, bd.first, bd.second, bd.equilibrium,
                         bd.forceConstant, input.numbers, input.xyz, rcov,
                         param.angleCutA, e, ag, ds);
      for (int c = 0; c < 3; ++c) {
        g[3 * bd.center + c] += ag[0][c];
        g[3 * bd.first + c] += ag[1][c];
        g[3 * bd.second + c] += ag[2][c];
      }
      eangl += e;
    }
    result.eangl = eangl;
  }

  // Torsions.
  {
    double etors = 0.0;
    for (const auto& t : input.torsions) {
      double ag[4][3] = { { 0 } };
      double ds[3][3] = { { 0 } };
      double e = 0.0;
      torsionEnergyGradient(t.i, t.j, t.k, t.l, t.multiplicity, t.phase,
                            t.forceConstant, input.numbers, input.xyz, rcov,
                            param.angleCutT, false, e, ag, ds);
      const int ids[4] = { t.i, t.j, t.k, t.l };
      for (int row = 0; row < 4; ++row) {
        for (int c = 0; c < 3; ++c)
          g[3 * ids[row] + c] += ag[row][c];
      }
      etors += e;
    }
    result.etors = etors;
  }

  // Bonded ATM triples.
  {
    double ebatm = 0.0;
    for (const auto& t : input.triples) {
      double ag[3][3] = { { 0 } };
      double e = 0.0;
      batmEnergyGradient(t.iat, t.jat, t.kat, input.numbers, input.xyz,
                         charges, param.zb3atm, e, ag, env);
      const int ids[3] = { t.iat, t.jat, t.kat };
      for (int row = 0; row < 3; ++row) {
        for (int c = 0; c < 3; ++c)
          g[3 * ids[row] + c] += ag[row][c];
      }
      ebatm += e;
    }
    result.ebatm = ebatm;
  }

  // HB (nhb1): A-H...B without neighbour orientation.
  {
    double ehb = 0.0;
    for (const auto& e1 : input.hb1) {
      double g3[3][3] = { { 0 } };
      double etmp = 0.0;
      double s1[3][3] = { { 0 } };
      if (!abhEg1(e1.a, e1.b, e1.h, input.numbers, input.xyz, input.qa,
                  input.hbBas, input.hbAci, param.rad, hbpar, input.mcfEhb,
                  etmp, g3, s1, env))
        return false;
      for (int c = 0; c < 3; ++c) {
        g[3 * e1.a + c] += g3[c][0] * input.mcfEhb;
        g[3 * e1.b + c] += g3[c][1] * input.mcfEhb;
        g[3 * e1.h + c] += g3[c][2] * input.mcfEhb;
      }
      ehb += etmp * input.mcfEhb;
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c)
          sigma[r][c] += s1[r][c];
      }
    }
    result.ehb = ehb;
  }

  // HB (nhb2) with carbonyl/nitro/rnr dispatch.
  {
    double ehb = result.ehb;
    for (const auto& e2 : input.hb2) {
      double etmp = 0.0;
      std::vector<double> g5(3 * n, 0.0);
      double s2[3][3] = { { 0 } };
      int b = e2.b;
      int nbb = static_cast<int>(input.neighbours[b].size());
      bool done = false;
      if (input.numbers[b] == 8 && nbb == 1) {
        int nbk = input.neighbours[b][0];
        int atnb = input.numbers[nbk];
        int nbnbk = static_cast<int>(input.neighbours[nbk].size());
        if (atnb == 6 && nbnbk > 1) {
          std::vector<int> cnbrs;
          for (int nb : input.neighbours[nbk]) {
            if (nb != b)
              cnbrs.push_back(nb);
          }
          // C's neighbour list incl. B like the reference tlist setup.
          std::vector<int> full = input.neighbours[nbk];
          int npair = n * (n + 1) / 2;
          std::vector<double> sq(npair, 0.0), sr(npair, 0.0);
          if (!abhEg3(e2.a, b, e2.h, nbk, full, input.numbers, input.xyz,
                       input.qa, input.hbBas, input.hbAci, sq, sr, param.rad,
                       hbpar, input.mcfEhb, etmp, g5, s2, env))
            return false;
          done = true;
        } else if (atnb == 7 && nbnbk > 1) {
          std::vector<int> full = input.neighbours[nbk];
          int npair = n * (n + 1) / 2;
          std::vector<double> sq(npair, 0.0), sr(npair, 0.0);
          if (!abhEg3(e2.a, b, e2.h, nbk, full, input.numbers, input.xyz,
                       input.qa, input.hbBas, input.hbAci, sq, sr, param.rad,
                       hbpar, input.mcfEhb, etmp, g5, s2, env))
            return false;
          done = true;
        }
      } else if (input.numbers[b] == 7 && nbb == 2) {
        if (!abhEg2rnr(e2.a, b, e2.h, input.neighbours[b], input.numbers,
                       input.xyz, input.qa, input.hbBas, input.hbAci,
                       param.rad, param.repz, hbpar, input.mcfEhb, etmp, g5,
                       s2, env))
          return false;
        done = true;
      }
      if (!done) {
        if (!abhEg2new(e2.a, b, e2.h, input.neighbours[b], input.numbers,
                       input.xyz, input.qa, input.hbBas, input.hbAci,
                       param.rad, hbpar, input.mcfEhb, etmp, g5, s2, env))
          return false;
      }
      for (int k = 0; k < 3 * n; ++k)
        g[k] += g5[k] * input.mcfEhb;
      ehb += etmp * input.mcfEhb;
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c)
          sigma[r][c] += s2[r][c];
      }
    }
    result.ehb = ehb;
  }

  // XB.
  {
    double exb = 0.0;
    for (const auto& xb : input.xb) {
      double g3[3][3] = { { 0 } };
      double etmp = 0.0;
      double s3[3][3] = { { 0 } };
      double cx = param.xbAci[input.numbers[xb.x] - 1];
      if (!rbxEg(xb.a, xb.b, xb.x, 1.0, cx, input.numbers, input.xyz,
                 input.qa, param.rad, hbpar, etmp, g3, s3, env))
        return false;
      for (int c = 0; c < 3; ++c) {
        g[3 * xb.a + c] += g3[c][0];
        g[3 * xb.b + c] += g3[c][1];
        g[3 * xb.x + c] += g3[c][2];
      }
      exb += etmp;
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c)
          sigma[r][c] += s3[r][c];
      }
    }
    result.exb = exb;
  }

  result.etot = result.ees + result.edisp + result.erep + result.ebond +
                result.eangl + result.etors + result.ehb + result.exb +
                result.ebatm;
  result.gradient = g;
  result.charges = charges;
  result.cn = cn;
  result.dlogCn = dlogCn;
  result.dipole[0] = result.dipole[1] = result.dipole[2] = 0.0;
  for (int i = 0; i < n; ++i) {
    result.dipole[0] += input.xyz[3 * i] * charges[i];
    result.dipole[1] += input.xyz[3 * i + 1] * charges[i];
    result.dipole[2] += input.xyz[3 * i + 2] * charges[i];
  }
  double g2 = 0.0;
  for (int k = 0; k < 3 * n; ++k)
    g2 += g[k] * g[k];
  result.gnorm = std::sqrt(g2);
  (void)sigma;
  (void)repThr;
  (void)hbThr1;
  (void)hbThr2;
  return true;
}

} // namespace Xtb
} // namespace Avogadro
