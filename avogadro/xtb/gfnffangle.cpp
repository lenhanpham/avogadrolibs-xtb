/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_ini.f90 (alist/vangl setup), src/constr.f90 (banglPBC,
  mode 1 only) and src/gfnff/gfnff_ini2.F90 (ringsbend),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#include "gfnffangle.h"

#include "constants.h"
#include "environment.h"
#include "gffgraph.h"
#include "gfnffbonds.h"
#include "gfnffdata.h"
#include "gfnfftopo.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace Avogadro {
namespace Xtb {

double bondAnglePbc(const std::vector<double>& xyz, int end1, int central,
                    int end2)
{
  double dx1 = xyz[3 * end1] - xyz[3 * central];
  double dy1 = xyz[3 * end1 + 1] - xyz[3 * central + 1];
  double dz1 = xyz[3 * end1 + 2] - xyz[3 * central + 2];
  double dx2 = xyz[3 * central] - xyz[3 * end2];
  double dy2 = xyz[3 * central + 1] - xyz[3 * end2 + 1];
  double dz2 = xyz[3 * central + 2] - xyz[3 * end2 + 2];
  double dx3 = xyz[3 * end1] - xyz[3 * end2];
  double dy3 = xyz[3 * end1 + 1] - xyz[3 * end2 + 1];
  double dz3 = xyz[3 * end1 + 2] - xyz[3 * end2 + 2];
  // Left-associated sums like sum() over the component squares.
  double d2ij = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
  double d2jk = dx2 * dx2 + dy2 * dy2 + dz2 * dz2;
  double d2ik = dx3 * dx3 + dy3 * dy3 + dz3 * dz3;
  double xy = std::sqrt(d2ij * d2jk);
  double temp = 0.5 * (d2ij + d2jk - d2ik) / xy;
  if (temp > 1.0)
    temp = 1.0;
  if (temp < -1.0)
    temp = -1.0;
  return std::acos(temp);
}

int smallestRingAngle(const std::vector<Ring>& ringsI,
                      const std::vector<Ring>& ringsJ,
                      const std::vector<Ring>& ringsK, int i, int j, int k)
{
  // Zero-padded tables like c(10,20,n)/s(20,n): s[19] is the ring count.
  int si[20] = { 0 }, sj[20] = { 0 }, sk[20] = { 0 };
  int ci[10][20] = { { 0 } }, cj[10][20] = { { 0 } },
      ck[10][20] = { { 0 } };
  const std::vector<Ring>* tabs[3] = { &ringsI, &ringsJ, &ringsK };
  int* ss[3] = { si, sj, sk };
  int(*cs[3])[20] = { ci, cj, ck };
  for (int t = 0; t < 3; ++t) {
    int nr = std::min<int>(static_cast<int>(tabs[t]->size()), 19);
    ss[t][19] = nr;
    for (int m = 0; m < nr; ++m) {
      const std::vector<int>& mem = (*tabs[t])[m].members;
      ss[t][m] = static_cast<int>(mem.size());
      for (size_t l = 0; l < mem.size() && l < 10; ++l)
        cs[t][l][m] = mem[l];
    }
  }
  if (si[19] == 0 || sj[19] == 0 || sk[19] == 0)
    return 0;
  int rings1 = 99, rings2 = 99, rings3 = 99;
  for (int m = 0; m < si[19]; ++m) {
    int itest = 0;
    for (int l = 0; l < si[m]; ++l) {
      if (ci[l][m] == j || ci[l][m] == k)
        ++itest;
    }
    if (itest == 2 && si[m] < rings1)
      rings1 = si[m];
  }
  for (int m = 0; m < sj[19]; ++m) {
    int itest = 0;
    for (int l = 0; l < sj[m]; ++l) {
      if (cj[l][m] == i || cj[l][m] == k)
        ++itest;
    }
    if (itest == 2 && sj[m] < rings2)
      rings2 = sj[m];
  }
  for (int m = 0; m < sk[19]; ++m) {
    int itest = 0;
    for (int l = 0; l < sk[m]; ++l) {
      if (ck[l][m] == i || ck[l][m] == j)
        ++itest;
    }
    // NOTE: s(m,j) instead of s(m,k), kept verbatim from the reference.
    if (itest == 2 && sj[m] < rings3)
      rings3 = sk[m];
  }
  int rings = std::min(rings1, std::min(rings2, rings3));
  if (rings == 99)
    rings = 0;
  return rings;
}

bool buildAngleList(int n, const std::vector<std::vector<int>>& neighbours,
                    const std::vector<int>& numbers,
                    const std::vector<double>& xyz,
                    const std::vector<double>& angl,
                    const std::vector<double>& angl2, double fcThr,
                    const std::vector<int>& metalTable,
                    std::vector<Angle>& angles, Environment& env)
{
  if (n <= 0) {
    env.error("empty angle setup", "buildAngleList");
    return false;
  }
  angles.clear();
  for (int i = 0; i < n; ++i) {
    int nn = static_cast<int>(neighbours[i].size());
    if (nn <= 1)
      continue;
    if (nn > 6)
      continue;
    int ati = numbers[i];
    for (int j = 0; j < nn; ++j) {
      for (int k = 0; k < j; ++k) {
        int jj = neighbours[i][j];
        int kk = neighbours[i][k];
        int atj = numbers[jj], atk = numbers[kk];
        double fijk =
          angl[ati - 1] * angl2[atj - 1] * angl2[atk - 1];
        if (fijk < fcThr)
          continue;
        double phi = bondAnglePbc(xyz, jj, i, kk);
        // Single-precision 180./60. like the reference.
        if (metalTable[ati - 1] > 0 && (phi * 180.0f) / pi < 60.0f)
          continue; // eta cases
        Angle a;
        a.center = i;
        a.first = jj;
        a.second = kk;
        a.phi = phi;
        angles.push_back(a);
      }
    }
  }
  return true;
}

bool buildAngleTerms(const std::vector<Angle>& angles,
                     const std::vector<AngleAtom>& atoms,
                     const std::vector<std::vector<int>>& neighbours,
                     const std::vector<std::vector<Ring>>& ringsPerAtom,
                     const std::vector<double>& pboPacked,
                     const std::vector<int>& group,
                     const std::vector<int>& metalTable,
                     const std::vector<double>& angl,
                     const std::vector<double>& angl2, const GffData& param,
                     const GffGenerator& gen, std::vector<AngleTerm>& terms,
                     Environment& env)
{
  (void)param;
  (void)angl;
  (void)angl2;
  int n = static_cast<int>(atoms.size());
  if (n <= 0) {
    env.error("empty angle setup", "buildAngleTerms");
    return false;
  }
  // Layout adapter for isAmideNitrogen (numctr = 1): its
  // [(cell * n + atom) * maxNeighbors + slot] layout.
  std::vector<int> nbLists(n * maxNeighbors, 0), nbCounts(n, 0),
    numbers(n), hyb(n), piAtoms(n);
  for (int i = 0; i < n; ++i) {
    numbers[i] = atoms[i].element;
    hyb[i] = atoms[i].hyb;
    piAtoms[i] = atoms[i].pi;
    nbCounts[i] = static_cast<int>(neighbours[i].size());
    for (size_t s = 0; s < neighbours[i].size(); ++s)
      nbLists[i * maxNeighbors + s] = neighbours[i][s];
  }
  terms.assign(angles.size(), AngleTerm());
  for (size_t a = 0; a < angles.size(); ++a) {
    int ii = angles[a].center, jj = angles[a].first, kk = angles[a].second;
    double phi = angles[a].phi;
    int ati = atoms[ii].element;
    int atj = atoms[jj].element, atk = atoms[kk].element;
    int nn = static_cast<int>(neighbours[ii].size());
    double feta = 1.0;
    if (atoms[ii].imetal == 2 && atoms[jj].itag == -1 && atoms[jj].pi > 0)
      feta = 0.3;
    if (atoms[ii].imetal == 2 && atoms[kk].itag == -1 && atoms[kk].pi > 0)
      feta = feta * 0.3;
    int nh = 0, nnn = 0, no = 0, nheav = 0, nsi = 0, nc = 0, nmet = 0,
        npi = 0;
    if (atj == 1)
      ++nh;
    if (atk == 1)
      ++nh;
    if (atj == 7)
      ++nnn;
    if (atk == 7)
      ++nnn;
    if (atj == 8)
      ++no;
    if (atk == 8)
      ++no;
    if (atj > 14)
      ++nheav;
    if (atk > 14)
      ++nheav;
    if (atj == 14)
      ++nsi;
    if (atk == 14)
      ++nsi;
    if (atj == 6)
      ++nc;
    if (atk == 6)
      ++nc;
    if (metalTable[atj - 1] != 0)
      ++nmet;
    if (metalTable[atk - 1] != 0)
      ++nmet;
    if (atoms[jj].pi != 0)
      ++npi;
    if (atoms[kk].pi != 0)
      ++npi;
    int rings = smallestRingAngle(ringsPerAtom[ii], ringsPerAtom[jj],
                                  ringsPerAtom[kk], ii, jj, kk);
    // Kept verbatim (redundant second hyb(ii) test like the reference).
    bool triple = (atoms[ii].hyb == 1 || atoms[jj].hyb == 1) ||
                  (atoms[ii].hyb == 1 || atoms[kk].hyb == 1);
    double fqq = 0.0;
    if (atoms[ii].imetal == 0 && atoms[jj].imetal == 0 &&
        atoms[kk].imetal == 0)
      fqq = 1.0 - (atoms[ii].charge * atoms[jj].charge +
                    atoms[ii].charge * atoms[kk].charge) *
                     gen.qFacBen;
    else
      fqq = 1.0 - (atoms[ii].charge * atoms[jj].charge +
                    atoms[ii].charge * atoms[kk].charge) *
                     gen.qFacBen * 2.5;
    double f2 = 1.0, fn = 1.0;
    // Unsuffixed literals below are single precision, widened on use.
    double r0 = 100.0f;
    if (atoms[ii].hyb == 1)
      r0 = 180.0f;
    if (atoms[ii].hyb == 2)
      r0 = 120.0f;
    if (atoms[ii].hyb == 3)
      r0 = 109.5f;
    if (atoms[ii].hyb == 3 && atoms[ii].element > 10) {
      if (nn <= 3)
        r0 = gen.aheavy3;
      if (nn >= 4)
        r0 = gen.aheavy4;
      if (nn == 4 && group[ati - 1] == 5)
        r0 = 109.5f;
      if (nn == 4 && group[ati - 1] == 4 && ati > 49)
        r0 = 109.5f;
      if (group[ati - 1] == 4)
        r0 = r0 - static_cast<double>(nh * 5.0f);
      if (group[ati - 1] == 5)
        r0 = r0 - static_cast<double>(nh * 5.0f);
      if (group[ati - 1] == 6)
        r0 = r0 - static_cast<double>(nh * 5.0f);
    }
    if (atoms[ii].hyb == 5) {
      r0 = 90.0f;
      f2 = 0.11f;
      if ((phi * 180.0f) / pi > gen.linThr)
        r0 = 180.0f;
    }
    if (ati == 5) {
      if (atoms[ii].hyb == 3)
        r0 = 115.0f;
      if (atoms[ii].hyb == 2)
        r0 = 115.0f;
    }
    if (ati == 6) {
      if (atoms[ii].hyb == 3 && nh == 2)
        r0 = 108.6f;
      if (atoms[ii].hyb == 3 && no == 1)
        r0 = 108.5f;
      if (atoms[ii].hyb == 2 && no == 2)
        r0 = 122.0f;
      if (atoms[ii].hyb == 2 && no == 1)
        f2 = 0.7f;
      if (atoms[ii].hyb == 1 && no == 2) {
        triple = false;
        f2 = 2.0f;
      }
      if (atoms[ii].hyb == 3 && nn > 4) {
        if ((phi * 180.0f) / pi > gen.linThr)
          r0 = 180.0f;
      }
    }
    if (ati == 8 && nn == 2) {
      r0 = 104.5f;
      if (nh == 2) {
        r0 = 100.0f;
        f2 = 1.20f;
      }
      r0 = r0 + static_cast<double>(7.0f * nsi);
      r0 = r0 + static_cast<double>(14.0f * nmet);
      if (npi == 2)
        r0 = 109.0f;
      if (nmet > 0 && (phi * 180.0f) / pi > gen.linThr) {
        r0 = 180.0f;
        f2 = 0.3f;
      }
    }
    if (ati == 7 && nn == 2) {
      f2 = 1.4f;
      r0 = 115.0f;
      if (rings != 0)
        r0 = 105.0f;
      if (atoms[kk].element == 8 || atoms[jj].element == 8)
        r0 = 103.0f;
      if (atoms[kk].element == 9 || atoms[jj].element == 9)
        r0 = 102.0f;
      if (atoms[ii].hyb == 1)
        r0 = 180.0f;
      if (atoms[jj].imetal == 2 && atoms[ii].hyb == 1 && atk == 7)
        r0 = 135.0f;
      if (atoms[kk].imetal == 2 && atoms[ii].hyb == 1 && atj == 7)
        r0 = 135.0f;
    }
    if (ati == 7 && atoms[ii].hyb == 3) {
      if (npi > 0) {
        if (isAmideNitrogen(n, numbers, hyb, nbLists, nbCounts, 1, piAtoms,
                            ii)) {
          r0 = 115.0f;
          f2 = 1.2;
        } else {
          double sumppi = pboPacked[packedIndex(ii, jj)] +
                          pboPacked[packedIndex(ii, kk)];
          r0 = 113.0f;
          f2 = 1.0 - sumppi * 0.7;
        }
      } else {
        r0 = 104.0f;
        f2 = 0.40f;
        f2 = f2 + static_cast<double>(nh * 0.19f);
        f2 = f2 + static_cast<double>(no * 0.25f);
        f2 = f2 + static_cast<double>(nc * 0.01f);
      }
    }
    if (rings == 3)
      r0 = 82.0f;
    if (rings == 4)
      r0 = 96.0f;
    if (rings == 5 && ati == 6)
      r0 = 109.0f;
    if (rings == 0) {
      if (atoms[ii].ring == 3) {
        int ringsj = atoms[jj].ring, ringsk = atoms[kk].ring;
        if (ringsj + ringsk == 102)
          r0 = r0 + 4.0;
      }
    }
    if (triple) {
      f2 = 0.60;
      if (atj == 7 || atk == 7)
        f2 = 1.00;
      if ((atoms[jj].imetal == 2 || atoms[kk].imetal == 2) &&
          (phi * 180.0f) / pi > gen.linThr) {
        if (ati == 6 && atj == 6)
          f2 = 3.0f;
        if (ati == 6 && atk == 6)
          f2 = 3.0f;
        if (ati == 6 && atj == 7)
          f2 = 3.0f;
        if (ati == 6 && atk == 7)
          f2 = 3.0f;
        if (ati == 6 && group[atj - 1] == 6)
          f2 = 14.0f;
        if (ati == 6 && group[atk - 1] == 6)
          f2 = 14.0f;
        if (ati == 7 && atj == 7)
          f2 = 10.0f;
        if (ati == 7 && atj == 6)
          f2 = 10.0f;
        if (ati == 7 && atk == 6)
          f2 = 10.0f;
        if (ati == 7 && atj == 8) {
          r0 = 180.0f;
          f2 = 12.0f;
        }
        if (ati == 7 && atk == 8) {
          r0 = 180.0f;
          f2 = 12.0f;
        }
      }
    }
    if (group[ati - 1] == 4 && nn == 2 && atoms[ii].itag == 1) {
      if (ati == 6)
        r0 = 145.0f;
      if (ati > 6)
        r0 = 90.0f;
    }
    if (group[ati - 1] == 6 && nn == 4 && no >= 1)
      r0 = 115.0f;
    if (group[ati - 1] == 7 && atoms[ii].hyb == 1) {
      if (ati == 9)
        r0 = 90.0f;
      if (ati == 17)
        r0 = 90.0f;
      if (ati == 35)
        r0 = 90.0f;
      if (ati == 53)
        r0 = 90.0f;
      if (ati > 9 && (phi * 180.0f) / pi > gen.linThr)
        r0 = 180.0f;
      f2 = 0.6f / std::pow(static_cast<double>(ati),
                           static_cast<double>(0.15f));
    }
    if (atoms[ii].hyb == 3 && group[ati - 1] == 4 && ati > 32 &&
        atoms[ii].charge > 0.4) {
      if ((phi * 180.0f) / pi > 140.0f)
        r0 = 180.0f;
      if ((phi * 180.0f) / pi < 100.0f)
        r0 = 90.0f;
      f2 = 1.0;
    }
    if (atoms[ii].imetal > 0) {
      if (atoms[ii].hyb == 0) {
        r0 = 90.0f;
        f2 = 1.35f;
      }
      if (atoms[ii].hyb == 1)
        r0 = 180.0f;
      if (atoms[ii].hyb == 2)
        r0 = 120.0f;
      if (atoms[ii].hyb == 3)
        r0 = 109.5f;
      if ((phi * 180.0f) / pi > gen.linThr)
        r0 = 180.0f;
    }
    double dnn = static_cast<double>(nn);
    fn = 1.0 - 2.36 / (dnn * dnn);
    AngleTerm term;
    term.equilibrium = r0 * pi / 180.0f;
    double diff = term.equilibrium - pi;
    double fbsmall =
      1.0 - gen.fbs1 * std::exp(static_cast<double>(-0.64f) * (diff * diff));
    term.forceConstant = fijkOf(angles[a], atoms, angl, angl2) * fqq * f2 *
                         fn * fbsmall * feta;
    terms[a] = term;
  }
  return true;
}

// angl-table product for one angle (fijk in the reference loop).
double fijkOf(const Angle& angle, const std::vector<AngleAtom>& atoms,
              const std::vector<double>& angl,
              const std::vector<double>& angl2)
{
  return angl[atoms[angle.center].element - 1] *
         angl2[atoms[angle.first].element - 1] *
         angl2[atoms[angle.second].element - 1];
}

} // namespace Xtb
} // namespace Avogadro
