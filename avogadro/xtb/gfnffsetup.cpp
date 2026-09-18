/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_param.f90 (newGFNFFGenerator, gfnff_set_param,
  gfnff_thresholds, gfnff_load_param, loadGFNFFAngewChem2020),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#include "gfnffsetup.h"

#include "constants.h"
#include "coordination.h"
#include "eeq.h"
#include "environment.h"
#include "gffgraph.h"
#include "gfnffangle.h"
#include "gfnffbonds.h"
#include "gfnffdisp.h"
#include "gfnffegbond.h"
#include "gfnffegnonbond.h"
#include "gfnffhb.h"
#include "gfnffhuckel.h"
#include "gfnffhyb.h"
#include "gfnffparams.h"
#include "gfnffrab.h"
#include "gfnfftopo.h"
#include "gfnfftorsion.h"
#include "gfnffvbond.h"

#include <cmath>

namespace Avogadro {
namespace Xtb {

void GffData::init(int ndim)
{
  en.assign(ndim, 0.0);
  rad.assign(ndim, 0.0);
  rcov.assign(ndim, 0.0);
  metal.assign(ndim, 0);
  group.assign(ndim, 0);
  normCn.assign(ndim, 0);
  repa.assign(ndim, 0.0);
  repan.assign(ndim, 0.0);
  repz.assign(ndim, 0.0);
  zb3atm.assign(ndim, 0.0);
  xhAci.assign(ndim, 0.0);
  xhBas.assign(ndim, 0.0);
  xbAci.assign(ndim, 0.0);
  chi.assign(ndim, 0.0);
  gam.assign(ndim, 0.0);
  cnf.assign(ndim, 0.0);
  alp.assign(ndim, 0.0);
  bond.assign(ndim, 0.0);
  angl.assign(ndim, 0.0);
  angl2.assign(ndim, 0.0);
  tors.assign(ndim, 0.0);
  tors2.assign(ndim, 0.0);
  d3r0.assign(ndim * (1 + ndim) / 2, 0.0);
}

void makeDefaultGenerator(GffGenerator& gen)
{
  gen = GffGenerator();
  gen.cnMax = 4.4f;
  gen.linThr = 160.0f;
  gen.fcThr = 1.0e-3;
  gen.tdistThr = 12.0f;
  gen.rThr = 1.25f;
  gen.rThr2 = 1.00f;
  gen.rShrink = 0.23f;
  gen.hQaThr = 0.01f;
  gen.qaBThr = 0.10f;
  gen.srb1 = 0.3731f;
  gen.srb2 = 0.3171f;
  gen.srb3 = 0.2538f;
  gen.qRepScal = 0.3480f;
  gen.nRepScal = -0.1270f;
  gen.hhFac = 0.6290f;
  gen.hh13Rep = 1.4580f;
  gen.hh14Rep = 0.7080f;
  gen.bstren = { 1.00, 1.24, 1.98, 1.22, 1.00, 0.78, 3.40, 3.40, 0.0 };
  gen.qFacBen = -0.54f;
  gen.qFacTor = 12.0;
  gen.fr3 = 0.3f;
  gen.fr4 = 1.0f;
  gen.fr5 = 1.5f;
  gen.fr6 = 5.7f;
  // torsf(4) is never assigned by the reference; it stays 0.0 here while
  // Fortran leaves it indeterminate.
  gen.torsf = { 1.00f, 1.18f, 1.05f, 0.0f, 0.50f, -0.90f, 0.70f, -2.00f };
  gen.fbs1 = 0.50f;
  gen.batmScal = 0.30;
  gen.mchiShift = -0.09;
  gen.rabShift = -0.110f;
  gen.rabShiftH = -0.050f;
  gen.hyperShift = 0.03f;
  gen.hShift3 = -0.11f;
  gen.hShift4 = -0.11f;
  gen.hShift5 = -0.06f;
  gen.metal1Shift = 0.2f;
  gen.metal2Shift = 0.15f;
  gen.metal3Shift = 0.05f;
  gen.etaShift = 0.040f;
  gen.qfacbm = { 1.0, -0.2, -0.2, 0.70, 0.50 };
  gen.qfacbm0 = 0.047f;
  gen.rfgoed1 = 1.175f;
  gen.hTriple = 1.45;
  gen.hueckelP2 = 1.00;
  gen.hueckelP3 = -0.24;
  // hdiag/hoffdiag use 1-based indices in xtb; unset entries stay 0.0.
  gen.hdiag[4] = -0.5;
  gen.hdiag[5] = 0.00;
  gen.hdiag[6] = 0.14;
  gen.hdiag[7] = -0.38;
  gen.hdiag[8] = -0.29;
  gen.hdiag[15] = -0.30;
  gen.hdiag[16] = -0.30;
  gen.hoffdiag[4] = 0.5;
  gen.hoffdiag[5] = 1.00;
  gen.hoffdiag[6] = 0.66;
  gen.hoffdiag[7] = 1.10;
  gen.hoffdiag[8] = 0.23;
  gen.hoffdiag[15] = 0.60;
  gen.hoffdiag[16] = 1.00;
  gen.hIter = 0.700;
  gen.hueckelP = 0.340;
  gen.bzRef = 0.370;
  gen.bzRef2 = 0.315;
  gen.pilpf = 0.530;
  gen.maxHIter = 5.0;
  gen.d3a1 = 0.58;
  gen.d3a2 = 4.80;
  gen.split0 = 0.670;
  gen.fringbo = 0.020;
  gen.aheavy3 = 89.0f;
  gen.aheavy4 = 100.0f;
  gen.split1 = 1.0 - gen.split0;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j)
      gen.bsmat[i][j] = -999.0f;
  }
  gen.bsmat[0][0] = gen.bstren[0];
  gen.bsmat[3][0] = gen.bstren[0];
  gen.bsmat[3][3] = gen.bstren[0];
  gen.bsmat[2][2] = gen.bstren[1];
  gen.bsmat[1][1] = gen.bstren[2];
  gen.bsmat[1][0] = gen.split0 * gen.bstren[0] + gen.split1 * gen.bstren[2];
  gen.bsmat[3][1] = gen.split0 * gen.bstren[0] + gen.split1 * gen.bstren[2];
  gen.bsmat[2][1] = gen.split0 * gen.bstren[1] + gen.split1 * gen.bstren[2];
  gen.bsmat[2][0] = gen.split0 * gen.bstren[0] + gen.split1 * gen.bstren[1];
  gen.bsmat[3][2] = gen.split0 * gen.bstren[0] + gen.split1 * gen.bstren[1];
  gen.bstren[8] = 0.5 * (gen.bstren[6] + gen.bstren[7]);
}

void gffThresholds(double accuracy, double& dispThr, double& cnThr,
                   double& repThr, double& hbThr1, double& hbThr2)
{
  dispThr = 1500.0 - std::log10(accuracy) * 1000.0;
  cnThr = 100.0 - std::log10(accuracy) * 50.0;
  repThr = 400.0 - std::log10(accuracy) * 100.0;
  hbThr1 = 200.0 - std::log10(accuracy) * 50.0;
  hbThr2 = 400.0 - std::log10(accuracy) * 50.0;
}

static void setGffParams(GffData& param, GffGenerator& gen)
{
  param.cnMax = 4.4f;
  param.angleCutA = 0.595;
  param.angleCutT = 0.505;
  param.angleCutANci = 0.395;
  param.angleCutTNci = 0.305;
  param.repScaleB = 1.7583f;
  param.repScaleN = 0.4270;
  param.dispScale = 1.0; // set%dispscale default (gfnff_setup copies it)
  param.hbAngleCut = 49.0f;
  param.hbShortCut = 22.0f;
  param.xbAngleCut = 70.0f;
  param.xbShortCut = 5.0f;
  param.hbSf = 1.0f;
  param.hbSt = 15.0f;
  param.xbSf = 0.03f;
  param.xbSt = 15.0f;
  param.hbAlp = 6.0f;
  param.hbLongCut = 85.0f;
  param.hbLongCutXb = 70.0f;
  param.hbAbMix = 0.80f;
  param.hbNbCut = 11.20f;
  param.torsHb = 0.94f;
  param.bendHb = 0.20f;
  param.vbondScale = 0.9f;
  param.xhAciGlobAbH = 0.268f;
  param.xhAciCoh = 0.350f;
  param.xhAciGlob = 1.50f;
  // Element indices below are 0-based (xtb uses 1-based element numbers).
  param.xhBas[5] = 0.80;
  param.xhBas[6] = 1.68;
  param.xhBas[7] = 0.67;
  param.xhBas[8] = 0.52;
  param.xhBas[13] = 4.0;
  param.xhBas[14] = 3.5;
  param.xhBas[15] = 2.0;
  param.xhBas[16] = 1.5;
  param.xhBas[34] = 1.5;
  param.xhBas[52] = 1.9;
  param.xhBas[32] = param.xhBas[14];
  param.xhBas[33] = param.xhBas[15];
  param.xhBas[50] = param.xhBas[14];
  param.xhBas[51] = param.xhBas[15];
  param.xhAci[5] = 0.75f;
  param.xhAci[6] = param.xhAciGlob + 0.1f;
  param.xhAci[7] = param.xhAciGlob;
  param.xhAci[8] = param.xhAciGlob;
  param.xhAci[14] = param.xhAciGlob;
  param.xhAci[15] = param.xhAciGlob;
  param.xhAci[16] = param.xhAciGlob + 1.0f;
  param.xhAci[34] = param.xhAciGlob + 1.0f;
  param.xhAci[52] = param.xhAciGlob + 1.0f;
  param.xbAci[14] = 1.0;
  param.xbAci[15] = 1.0;
  param.xbAci[16] = 0.5;
  param.xbAci[32] = 1.2;
  param.xbAci[33] = 1.2;
  param.xbAci[34] = 0.9;
  param.xbAci[50] = 1.2;
  param.xbAci[51] = 1.2;
  param.xbAci[52] = 1.2;

  // 3-atom bond prefactors and packed D3 R0^2 table.
  double cubeRoot = std::cbrt(gen.batmScal);
  int k = 0;
  for (int i = 0; i < gffElements; ++i) {
    param.zb3atm[i] = -(i + 1) * cubeRoot;
    for (int j = 0; j <= i; ++j) {
      double dum = sqrtZr4r2Value(i + 1) * sqrtZr4r2Value(j + 1) * 3.0;
      param.d3r0[k++] =
        (gen.d3a1 * std::sqrt(dum) + gen.d3a2) * (gen.d3a1 * std::sqrt(dum) + gen.d3a2);
    }
  }
  param.zb3atm[0] = -0.25 * cubeRoot;
}

bool loadGffParams(int version, GffData& param, GffGenerator& gen)
{
  if (version != GffAngewChem2020 && version != GffAngewChem2020_1 &&
      version != GffAngewChem2020_2 && version != GffHarmonic2020 &&
      version != GffMcGfnFF2023) {
    return false;
  }
  param.init(gffElements);
  for (int i = 0; i < gffElements; ++i) {
    param.en[i] = gffEn[i];
    param.rad[i] = gffRad[i];
    // Mirrors covalentRadD3: ([...] * aatoau) * 4/3 with xtb's own
    // Angstrom definition (0.52917726), not the CODATA one.
    param.rcov[i] =
      ((covalentRadD3Angstrom[i] * xtbAngstromToBohr) * 4.0) / 3.0;
    param.metal[i] = gffMetal[i];
    param.group[i] = gffGroup[i];
    param.normCn[i] = gffNormCn[i];
    param.repz[i] = gffRepz[i];
    param.chi[i] = gffChi[i];
    param.gam[i] = gffGam[i];
    param.cnf[i] = gffCnf[i];
    param.alp[i] = gffAlp[i];
    param.bond[i] = gffBond[i];
    param.repa[i] = gffRepa[i];
    param.repan[i] = gffRepan[i];
    param.angl[i] = gffAngl[i];
    param.angl2[i] = gffAngl2[i];
    param.tors[i] = gffTors[i];
    param.tors2[i] = gffTors2[i];
  }
  makeDefaultGenerator(gen);
  setGffParams(param, gen);
  return true;
}

HbParams makeHbParams(const GffData& param)
{
  HbParams hp;
  hp.hbacut = param.hbAngleCut;
  hp.hblongcut = param.hbLongCut;
  hp.hbscut = param.hbShortCut;
  hp.hbalp = param.hbAlp;
  hp.hbst = param.hbSt;
  hp.hbsf = param.hbSf;
  hp.hbabmix = param.hbAbMix;
  hp.hbnbcut = param.hbNbCut;
  hp.xbacut = param.xbAngleCut;
  hp.xbscut = param.xbShortCut;
  hp.xbst = param.xbSt;
  hp.xbsf = param.xbSf;
  hp.hblongcutXb = param.hbLongCutXb;
  hp.xhaciGlobAbh = param.xhAciGlobAbH;
  hp.xhaciCoh = param.xhAciCoh;
  hp.torsHb = param.torsHb;
  hp.bendHb = param.bendHb;
  return hp;
}

bool gfnffSetup0d(int version, const std::vector<int>& numbers,
                  const std::vector<double>& xyz, double totalCharge,
                  double accuracy, DriverInput& in, GffData& param,
                  GffGenerator& gen, Environment& env)
{
  int n = static_cast<int>(numbers.size());
  if (n <= 0 || static_cast<int>(xyz.size()) != 3 * n) {
    env.error("empty setup geometry", "gfnffSetup0d");
    return false;
  }
  if (!loadGffParams(version, param, gen))
    return false;
  std::vector<int> metal(gffElements), group(gffElements);
  for (int i = 0; i < gffElements; ++i) {
    metal[i] = gffMetal[i];
    group[i] = gffGroup[i];
  }
  double dispThr, cnThr, repThr, hbThr1, hbThr2;
  gffThresholds(accuracy, dispThr, cnThr, repThr, hbThr1, hbThr2);

  // Neighbour pipeline (hyb-test pattern, qa = 0 first pass).
  std::vector<double> cn0(n);
  for (int i = 0; i < n; ++i)
    cn0[i] = param.normCn[numbers[i] - 1];
  std::vector<double> rtmp;
  gfnffBondGuesses(n, numbers, cn0, rtmp);
  std::vector<double> qa0(n, 0.0);
  scaleBondGuesses(n, numbers, qa0, metal, gen.rShrink, rtmp);
  std::vector<double> dist(n * n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double dx = xyz[3 * j] - xyz[3 * i];
      double dy = xyz[3 * j + 1] - xyz[3 * i + 1];
      double dz = xyz[3 * j + 2] - xyz[3 * i + 2];
      dist[i * n + j] = std::sqrt(dx * dx + dy * dy + dz * dz);
    }
  }
  std::vector<double> mch0(n, 0.0);
  std::vector<int> nbf, nbfc, nb, nbc, nbm, nbmc, full(n, 0);
  fillNeighborList(n, numbers, rtmp, dist, mch0.data(), 1, 1.25, 1.0,
                   gffMetal.data(), gffGroup.data(), gffNormCn.data(),
                   full.data(), 1, nbf, nbfc);
  for (int i = 0; i < n; ++i)
    full[i] = nbfc[i];
  fillNeighborList(n, numbers, rtmp, dist, mch0.data(), 2, 1.25, 1.0,
                   gffMetal.data(), gffGroup.data(), gffNormCn.data(),
                   full.data(), 1, nb, nbc);
  fillNeighborList(n, numbers, rtmp, dist, mch0.data(), 3, 1.25, 1.0,
                   gffMetal.data(), gffGroup.data(), gffNormCn.data(),
                   full.data(), 1, nbm, nbmc);
  std::vector<int> hyb, itag;
  if (!assignHybridization(n, numbers, xyz, nbf, nbfc, nb, nbc, nbm, nbmc,
                           qa0, metal.data(), group.data(), 160.0, 1, true,
                           hyb, itag, env))
    return false;
  std::vector<std::vector<int>> neighbours(n);
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < nbfc[i]; ++k)
      neighbours[i].push_back(nbf[i * maxNeighbors + k]);
  }
  std::vector<int> imetal =
    effectiveMetals(n, numbers, nbc, 1, metal.data(), group.data());
  PiSystem pi = buildPiSystem(n, numbers, hyb, nb, nbc, 1);
  std::vector<int> piFlags(n, 0);
  for (int p : pi.piIndex) {
    if (p >= 0 && p < n)
      piFlags[p] = 1;
  }
  std::vector<int> counts(n);
  for (int i = 0; i < n; ++i)
    counts[i] = static_cast<int>(neighbours[i].size());
  std::vector<int> firstNb(n, -1);
  for (int i = 0; i < n; ++i) {
    if (!neighbours[i].empty())
      firstNb[i] = neighbours[i][0];
  }

  // EEQ xi + initial parameters.
  std::vector<double> dxi, chi0, gam0, alp0;
  if (!eeqXiCorrections(n, numbers, itag, imetal, piFlags, neighbours,
                        counts, group.data(), dxi, env))
    return false;
  if (!eeqInitialParams(n, numbers, param.chi, param.gam, param.alp,
                        param.cnf, imetal, counts, gen.cnMax, gen.mchiShift,
                        dxi, chi0, gam0, alp0, env))
    return false;

  // Coordination numbers + topology (goedeckera) charges for the setup.
  int npair = n * (n + 1) / 2;
  std::vector<double> srab(npair, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j)
      srab[packedIndex(i, j)] = dist[i * n + j];
  }
  std::vector<double> rcov(gffElements);
  for (int z = 0; z < gffElements; ++z)
    rcov[z] = covalentRadiusD3(z + 1);
  std::vector<double> cn, dlogCn;
  gffCoordinationNumber(n, numbers, xyz, srab, rcov, gen.cnMax, cnThr, cn,
                        dlogCn);
  std::vector<float> rabdF;
  std::vector<double> rtmpP;
  estimateBondLengths(n, numbers, nbc, nb, param.rad, gen.rfgoed1,
                      gen.tdistThr, rabdF, rtmpP);
  std::vector<int> fragOfAtom;
  int nfrag = findFragments(neighbours, fragOfAtom);
  if (nfrag < 1)
    return false;
  std::vector<double> fragCharges(nfrag, 0.0);
  fragCharges[0] = totalCharge;
  std::vector<double> qaEst;
  double esEst = 0.0;
  if (!topologyCharges(n, rtmpP, chi0, gam0, alp0, nfrag, fragOfAtom,
                       fragCharges, env, qaEst, esEst))
    return false;

  // Gamma + final EEQ parameters from the setup charges.
  std::vector<double> dgam, chieeq, gameeq, alpeeq;
  if (!eeqGammaCorrections(n, numbers, hyb, imetal, piFlags, neighbours,
                           group.data(), qaEst, dgam, env))
    return false;
  if (!eeqFinalParams(n, numbers, hyb, imetal, piFlags, neighbours,
                      param.chi, param.gam, param.alp, dxi, dgam, qaEst,
                      group.data(), chieeq, gameeq, alpeeq, env))
    return false;

  // HB strengths, bpair flags, bonds + types.
  std::vector<double> hbBas, hbAci;
  if (!hbBasicity(n, numbers, itag, neighbours, counts, param.xhBas, hbBas,
                  env))
    return false;
  if (!hbAcidity(n, numbers, neighbours, hyb, piFlags, param.xhAci, hbAci,
                 env))
    return false;
  std::vector<int> bpair;
  if (!bondPairFlags(n, neighbours, bpair, env))
    return false;
  std::vector<Bond> bonds = buildBondList(n, nbf, nbfc, 1);
  std::vector<std::pair<int, int>> abonds;
  for (const auto& b : bonds)
    abonds.emplace_back(b.first, b.second);

  // Ring lists over the icase-3 lists + per-bond ring fields, mirroring
  // the nbrngs/getring36 setup and the ringsbond/ringsatom queries.
  std::vector<std::vector<int>> nbmNeighbours(n);
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < nbmc[i]; ++k)
      nbmNeighbours[i].push_back(nbm[i * maxNeighbors + k]);
  }
  std::vector<std::vector<Ring>> ringsPerAtom(n);
  for (int i = 0; i < n; ++i) {
    if (!findRingsThrough(n, numbers, nbmNeighbours, i, ringsPerAtom[i],
                          env))
      return false;
  }

  // Iterative Hueckel pi bond orders (post-Hueckel piadr below).
  int npairPi = n * (n + 1) / 2;
  std::vector<double> pibo(abonds.size(), 0.0), pbo(npairPi, 0.0);
  std::vector<int> piadrOut(n, 0);
  {
    std::vector<double> hdiag(gen.hdiag.begin(), gen.hdiag.end());
    std::vector<double> hoffdiag(gen.hoffdiag.begin(), gen.hoffdiag.end());
    if (!huckelPiBondOrders(n, numbers, hyb, itag, qaEst, xyz, abonds,
                            pi.piIndex, pi.fragOfPi, pi.fragments, hdiag,
                            gen.hueckelP3, gen.pilpf, hoffdiag, gen.hIter,
                            gen.hTriple, gen.maxHIter, pibo, pbo, piadrOut,
                            env))
      return false;
  }
  std::vector<int> piPost;
  for (int i = 0; i < n; ++i) {
    if (piadrOut[i] != 0)
      piPost.push_back(i);
  }
  std::vector<int> btypes =
    assignBondTypes(bonds, numbers, hyb, itag, piPost, imetal,
                    group.data());

  // HB perception + triplets + maps + geometry lists.
  HbDonorLists donors;
  if (!hbDonorLists(n, numbers, hyb, piFlags, qaEst, neighbours, firstNb,
                    group.data(), param.xhBas, gen.hQaThr, gen.qaBThr,
                    bpair, hbBas, hbAci, donors, env))
    return false;
  std::vector<HbBondTriplet> triplets;
  if (!hbBondTriplets(n, donors.hatAB, donors.hatH, xyz, hbThr1, bpair,
                      triplets, env))
    return false;
  HbBondMaps bondHb;
  if (!hbAhbMaps(n, numbers, abonds, triplets, bondHb, env))
    return false;
  std::vector<HbTriple> hb1, hb2;
  std::vector<XbTriple> xb;
  if (!hbTripletLists(n, xyz, donors.hatAB, donors.hatH, donors, bpair,
                      hbThr1, hbThr2, hb1, hb2, xb, env))
    return false;

  // Non-bonded tables from the setup charges.
  std::vector<double> alphanb, zetac6dummy;
  if (!buildNonbondedTables(n, numbers, qaEst, counts, param.repan, metal,
                            gen.nRepScal, gen.qRepScal, gen.hhFac,
                            gen.hh13Rep, gen.hh14Rep, bpair, alphanb,
                            zetac6dummy, env))
    return false;

  // vbond terms (pibo + post-Hueckel pi + ring fields from the
  // ringsbond/ringsatom queries).
  std::vector<double> mchar(n);
  metallicCharacter(n, numbers, param.en, cn, dlogCn, mchar);
  std::vector<VbondBond> vb;
  for (size_t b = 0; b < bonds.size(); ++b) {
    VbondBond bb;
    bb.first = bonds[b].first;
    bb.second = bonds[b].second;
    bb.guess = rtmpP[packedIndex(bb.first, bb.second)];
    bb.pibo = pibo[b];
    bb.ring = smallestRingBond(ringsPerAtom[bb.first],
                               ringsPerAtom[bb.second], bb.first, bb.second);
    bb.ringFirst = smallestRingThrough(ringsPerAtom[bb.first]);
    bb.ringSecond = smallestRingThrough(ringsPerAtom[bb.second]);
    bb.coordFirst = static_cast<int>(neighbours[bb.first].size());
    bb.coordSecond = static_cast<int>(neighbours[bb.second].size());
    vb.push_back(bb);
  }
  std::vector<VbondAtom> va(n);
  for (int i = 0; i < n; ++i) {
    va[i].element = numbers[i];
    va[i].hyb = hyb[i];
    va[i].charge = qaEst[i];
    va[i].itag = itag[i];
    va[i].imetal = imetal[i];
    va[i].pi = piadrOut[i];
    va[i].mchar = mchar[i];
  }
  std::vector<int> row6(n);
  for (int i = 0; i < n; ++i)
    row6[i] = elementRow6(numbers[i]);
  std::vector<VbondTerm> vterms;
  std::vector<int> outBtypes;
  if (!buildVbondTerms(vb, va, neighbours, numbers, group, metal, param.en,
                       param.bond, row6, param, gen, vterms, outBtypes, env))
    return false;

  // Bends (rings from perception; pbo from Hueckel).
  std::vector<Angle> angles;
  if (!buildAngleList(n, neighbours, numbers, xyz, param.angl, param.angl2,
                      gen.fcThr, param.metal, angles, env))
    return false;
  std::vector<AngleAtom> aatoms(n);
  for (int i = 0; i < n; ++i) {
    aatoms[i].element = numbers[i];
    aatoms[i].hyb = hyb[i];
    aatoms[i].itag = itag[i];
    aatoms[i].imetal = imetal[i];
  }
  std::vector<AngleTerm> aterms;
  if (!buildAngleTerms(angles, aatoms, neighbours, ringsPerAtom, pbo, group,
                       metal, param.angl, param.angl2, param, gen, aterms,
                       env))
    return false;

  // Torsions (Hueckel pibo, post-Hueckel pi, setup charges in fqq).
  std::vector<TorsionBond> tb;
  for (size_t b = 0; b < bonds.size(); ++b) {
    TorsionBond t;
    t.first = bonds[b].first;
    t.second = bonds[b].second;
    t.btype = btypes[b];
    t.pibo = pibo[b];
    tb.push_back(t);
  }
  std::vector<TorsionAtom> tatoms(n);
  for (int i = 0; i < n; ++i) {
    tatoms[i].element = numbers[i];
    tatoms[i].hyb = hyb[i];
    tatoms[i].charge = qaEst[i];
    tatoms[i].imetal = imetal[i];
    tatoms[i].pi = piadrOut[i];
  }
  std::vector<Torsion> torsions;
  std::vector<TorsionTerm> tterms;
  if (!buildTorsions(tb, tatoms, neighbours, xyz, ringsPerAtom, group,
                     metal, param.tors, param.tors2, param, gen, torsions,
                     tterms, env))
    return false;

  // Bonded-ATM triples over bpair == 3 (ini b3list rule, 0d).
  struct Triple
  {
    int i, j, k;
  };
  std::vector<Triple> b3;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      if (bpair[j * n + i] != 3)
        continue;
      for (int k : neighbours[j]) {
        if (i == k || j == k)
          continue;
        b3.push_back({ i, j, k });
      }
      for (int k : neighbours[i]) {
        if (i == k || j == k)
          continue;
        b3.push_back({ i, j, k });
      }
    }
  }

  // Assemble the driver input.
  in = DriverInput();
  in.numbers = numbers;
  in.xyz = xyz;
  in.accuracy = accuracy;
  in.totalCharge = totalCharge;
  for (size_t b = 0; b < bonds.size(); ++b) {
    DriverBond db;
    db.first = bonds[b].first;
    db.second = bonds[b].second;
    db.shift = vterms[b].shift;
    db.steepness = vterms[b].steepness;
    db.prefactor = vterms[b].prefactor;
    in.bonds.push_back(db);
  }
  for (size_t a = 0; a < angles.size(); ++a) {
    DriverBend db;
    db.center = angles[a].center;
    db.first = angles[a].first;
    db.second = angles[a].second;
    db.equilibrium = aterms[a].equilibrium;
    db.forceConstant = aterms[a].forceConstant;
    in.bends.push_back(db);
  }
  for (size_t t = 0; t < torsions.size(); ++t) {
    DriverTorsion dt;
    dt.i = torsions[t].outer1;
    dt.j = torsions[t].center1;
    dt.k = torsions[t].center2;
    dt.l = torsions[t].outer2;
    dt.multiplicity = torsions[t].multiplicity;
    dt.phase = tterms[t].phase;
    dt.forceConstant = tterms[t].forceConstant;
    in.torsions.push_back(dt);
  }
  for (const auto& t : b3)
    in.triples.push_back({ t.i, t.j, t.k });
  for (const auto& t : hb1)
    in.hb1.push_back({ t.a, t.b, t.h });
  for (const auto& t : hb2)
    in.hb2.push_back({ t.a, t.b, t.h });
  for (const auto& t : xb)
    in.xb.push_back({ t.a, t.b, t.x });
  in.neighbours = neighbours;
  in.chieeq = chieeq;
  in.gameeq = gameeq;
  in.alpeeq = alpeeq;
  in.qa = qaEst;
  in.hbBas = hbBas;
  in.hbAci = hbAci;
  in.nrHb = bondHb.nrHb;
  in.bondHb = bondHb;
  in.fragOfAtom = fragOfAtom;
  in.fragCharges = fragCharges;
  in.nonbonded.alphanb = alphanb;
  in.nonbonded.bpair = bpair;
  return true;
}

} // namespace Xtb
} // namespace Avogadro
