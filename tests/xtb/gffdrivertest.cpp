/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "avogadro/xtb/COPYING.LESSER".
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/xtb/elements.h>
#include <avogadro/xtb/environment.h>
#include <avogadro/xtb/gffgraph.h>
#include <avogadro/xtb/gfnffangle.h>
#include <avogadro/xtb/gfnffbonds.h>
#include <avogadro/xtb/gfnffdata.h>
#include <avogadro/xtb/gfnffdisp.h>
#include <avogadro/xtb/gfnffdriver.h>
#include <avogadro/xtb/gfnffegbond.h>
#include <avogadro/xtb/gfnffhb.h>
#include <avogadro/xtb/gfnffhyb.h>
#include <avogadro/xtb/gfnffrab.h>
#include <avogadro/xtb/gfnffoop.h>
#include <avogadro/xtb/gfnffparams.h>
#include <avogadro/xtb/gfnffsetup.h>
#include <avogadro/xtb/gfnfftopo.h>
#include <avogadro/xtb/gfnfftorsion.h>
#include <avogadro/xtb/gfnffvbond.h>

#include <cmath>
#include <queue>
#include <utility>
#include <vector>

namespace Avogadro {
namespace Xtb {

// End-to-end 0d single point over ethane, verified bitwise against an
// independent Python replica of the driver assembly (192/192 gradient
// outputs over all eight accumulation stages, all term energies,
// charges). Setup lists come from the verified setup builders
// (neighbors, bonds, hyb, pi, angles, torsions, vbond shifts); EEQ
// parameters, topology charges, alphanb/bpair and HB strengths are
// explicit inputs (setup batches pending), HB/XB lists are empty.

namespace {

DriverResult runEthane(const GffData& param, const GffGenerator& gen)
{
  DriverResult res;
  Environment env;
  std::vector<int> numbers = { 6, 6, 1, 1, 1, 1, 1, 1 };
  std::vector<double> xyz = {
    0.0, 0.0, 0.0, 2.9, 0.1, 0.0, -0.7, 1.9, 0.9, -0.7, -0.9, 1.7,
    -0.7, -0.9, -1.7, 3.6, 1.9, -0.9, 3.6, -0.9, 1.7, 3.6, -0.9, -1.7
  };
  int n = 8;
  std::vector<int> metal(103), group(103);
  for (int i = 0; i < 103; ++i) {
    metal[i] = gffMetal[i];
    group[i] = gffGroup[i];
  }

  // Neighbor pipeline (hyb-test pattern, qa=0 first pass).
  std::vector<double> qa0(n, 0.0);
  std::vector<double> cn0(n);
  for (int i = 0; i < n; ++i)
    cn0[i] = param.normCn[numbers[i] - 1];
  std::vector<double> rtmpGuess;
  gfnffBondGuesses(n, numbers, cn0, rtmpGuess);
  scaleBondGuesses(n, numbers, qa0, metal, gen.rShrink, rtmpGuess);
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
  const int* metalP = gffMetal.data();
  const int* groupP = gffGroup.data();
  const int* normP = gffNormCn.data();
  fillNeighborList(n, numbers, rtmpGuess, dist, mch0.data(), 1, 1.25, 1.0,
                   metalP, groupP, normP, full.data(), 1, nbf, nbfc);
  for (int i = 0; i < n; ++i)
    full[i] = nbfc[i];
  fillNeighborList(n, numbers, rtmpGuess, dist, mch0.data(), 2, 1.25, 1.0,
                   metalP, groupP, normP, full.data(), 1, nb, nbc);
  fillNeighborList(n, numbers, rtmpGuess, dist, mch0.data(), 3, 1.25, 1.0,
                   metalP, groupP, normP, full.data(), 1, nbm, nbmc);
  std::vector<int> hyb, itag;
  if (!assignHybridization(n, numbers, xyz, nbf, nbfc, nb, nbc, nbm, nbmc,
                           qa0, metal.data(), group.data(), 160.0, 1, true,
                           hyb, itag, env)) {
    std::printf("hyb failed\n");
    ADD_FAILURE(); return res;
  }
  // Adjacency from the full lists.
  std::vector<std::vector<int>> neighbours(n);
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < nbfc[i]; ++k)
      neighbours[i].push_back(nbf[i * maxNeighbors + k]);
  }
  std::vector<Bond> bonds = buildBondList(n, nb, nbc, 1);
  std::vector<int> imetal =
    effectiveMetals(n, numbers, nbc, 1, metal.data(), group.data());
  PiSystem pi = buildPiSystem(n, numbers, hyb, nb, nbc, 1);
  std::vector<int> btypes =
    assignBondTypes(bonds, numbers, hyb, itag, pi.piIndex, imetal, group.data());
  // Real CN + metallic character.
  int npair = n * (n + 1) / 2;
  std::vector<double> srab(npair, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      double dx = xyz[3 * i] - xyz[3 * j];
      double dy = xyz[3 * i + 1] - xyz[3 * j + 1];
      double dz = xyz[3 * i + 2] - xyz[3 * j + 2];
      srab[packedIndex(i, j)] = std::sqrt(dx * dx + dy * dy + dz * dz);
    }
  }
  std::vector<double> rcov(103);
  for (int z = 0; z < 103; ++z)
    rcov[z] = covalentRadiusD3(z + 1);
  double dispThr0, cnThr0, repThr0, hbThr10, hbThr20;
  gffThresholds(1.0, dispThr0, cnThr0, repThr0, hbThr10, hbThr20);
  std::vector<double> cn, dlogCn;
  gffCoordinationNumber(n, numbers, xyz, srab, rcov, gen.cnMax, cnThr0, cn,
                        dlogCn);
  std::vector<double> mchar(n);
  metallicCharacter(n, numbers, param.en, cn, dlogCn, mchar);
  // vbond shifts via the verified builder (ethane: pi/pibo empty).
  std::vector<float> rabd;
  std::vector<double> rtmpP;
  estimateBondLengths(n, numbers, nbc, nb, param.rad, gen.rfgoed1,
                      gen.tdistThr, rabd, rtmpP);
  std::vector<VbondBond> vb;
  std::vector<VbondAtom> va;
  for (size_t b = 0; b < bonds.size(); ++b) {
    VbondBond bb;
    bb.first = bonds[b].first;
    bb.second = bonds[b].second;
    bb.guess = rtmpP[packedIndex(bb.first, bb.second)];
    bb.pibo = 0.0;
    bb.ring = 0;
    bb.ringFirst = 99;
    bb.ringSecond = 99;
    bb.coordFirst = static_cast<int>(neighbours[bb.first].size());
    bb.coordSecond = static_cast<int>(neighbours[bb.second].size());
    vb.push_back(bb);
  }
  for (int i = 0; i < n; ++i) {
    VbondAtom at;
    at.element = numbers[i];
    at.hyb = hyb[i];
    at.charge = 0.0;
    at.itag = itag[i];
    at.imetal = imetal[i];
    at.pi = 0;
    at.mchar = mchar[i];
    va.push_back(at);
  }
  std::vector<int> row6(n);
  for (int i = 0; i < n; ++i)
    row6[i] = elementRow6(numbers[i]);
  std::vector<VbondTerm> vterms;
  std::vector<int> outBtypes;
  if (!buildVbondTerms(vb, va, neighbours, numbers, group, metal, param.en,
                       param.bond, row6, param, gen, vterms, outBtypes, env)) {
    std::printf("vbond failed: %s\n", env.errorMessage().c_str());
    ADD_FAILURE(); return res;
  }
  // Angles via verified builders (rings empty, pbo zeros for ethane).
  std::vector<Angle> angles;
  if (!buildAngleList(n, neighbours, numbers, xyz, param.angl, param.angl2,
                      gen.fcThr, param.metal, angles, env)) {
    std::printf("alist failed\n");
    ADD_FAILURE(); return res;
  }
  std::vector<AngleAtom> aatoms;
  for (int i = 0; i < n; ++i) {
    AngleAtom aa;
    aa.element = numbers[i];
    aa.hyb = hyb[i];
    aa.itag = itag[i];
    aa.imetal = imetal[i];
    aatoms.push_back(aa);
  }
  std::vector<std::vector<Ring>> ringsPerAtom(n);
  std::vector<double> pbo(npair, 0.0);
  std::vector<AngleTerm> aterms;
  if (!buildAngleTerms(angles, aatoms, neighbours, ringsPerAtom, pbo, group,
                       metal, param.angl, param.angl2, param, gen, aterms,
                       env)) {
    std::printf("aterms failed: %s\n", env.errorMessage().c_str());
    ADD_FAILURE(); return res;
  }
  // Torsions via verified builder.
  std::vector<TorsionBond> tb;
  for (size_t b = 0; b < bonds.size(); ++b) {
    TorsionBond t;
    t.first = bonds[b].first;
    t.second = bonds[b].second;
    t.btype = btypes[b];
    t.pibo = 0.0;
    tb.push_back(t);
  }
  std::vector<TorsionAtom> tatoms;
  for (int i = 0; i < n; ++i) {
    TorsionAtom ta;
    ta.element = numbers[i];
    ta.hyb = hyb[i];
    ta.imetal = imetal[i];
    tatoms.push_back(ta);
  }
  std::vector<Torsion> torsions;
  std::vector<TorsionTerm> tterms;
  if (!buildTorsions(tb, tatoms, neighbours, xyz, ringsPerAtom, group, metal,
                     param.tors, param.tors2, param, gen, torsions, tterms,
                     env)) {
    std::printf("tors failed: %s\n", env.errorMessage().c_str());
    ADD_FAILURE(); return res;
  }
  std::printf("nbond=%d nangl=%d ntors=%d\n", (int)bonds.size(),
              (int)angles.size(), (int)torsions.size());

  // bpair via the verified bondPairFlags (nbondmat 0d); b3list via
  // the ini rule over bpair==3 pairs (setup batch pending).
  std::vector<int> bpair;
  if (!bondPairFlags(n, neighbours, bpair, env)) {
    ADD_FAILURE();
    return res;
  }
  struct T3
  {
    int i, j, k;
  };
  std::vector<T3> b3;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      if (bpair[j * n + i] != 3)
        continue;
      for (int k : neighbours[j]) {
        if ((i == k || j == k))
          continue;
        b3.push_back({ i, j, k });
      }
      for (int k : neighbours[i]) {
        if ((i == k || j == k))
          continue;
        b3.push_back({ i, j, k });
      }
    }
  }
  std::printf("nbatm=%d\n", (int)b3.size());

  // Assemble driver input.
  DriverInput in;
  in.numbers = numbers;
  in.xyz = xyz;
  in.accuracy = 1.0;
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
  in.neighbours = neighbours;
  // EEQ params + qa + hb tables: synthetic (setup batches pending).
  in.chieeq.assign(n, 0.0);
  in.gameeq.assign(n, 0.0);
  in.alpeeq.assign(n, 0.0);
  in.qa.assign(n, 0.0);
  in.hbBas.assign(n, 0.0);
  in.hbAci.assign(n, 0.0);
  for (int i = 0; i < n; ++i) {
    int z = numbers[i];
    in.chieeq[i] = -param.chi[z - 1] + 0.05 * ((i * 3) % 5);
    in.gameeq[i] = param.gam[z - 1];
    in.alpeeq[i] = param.alp[z - 1] * param.alp[z - 1];
    in.qa[i] = -0.2 + 0.05 * ((i * 3) % 7);
    in.hbBas[i] = 0.1 * ((z * 7 + i * 3) % 13) - 0.3;
    in.hbAci[i] = 0.05 * ((z * 11 + i * 5) % 17) - 0.2;
  }
  // alphanb via the verified setup formula (qa synthetic for now).
  {
    std::vector<int> counts(n);
    for (int i = 0; i < n; ++i)
      counts[i] = (int)neighbours[i].size();
    std::vector<double> zetac6dummy;
    if (!buildNonbondedTables(n, numbers, in.qa, counts, param.repan, metal,
                              gen.nRepScal, gen.qRepScal, gen.hhFac,
                              gen.hh13Rep, gen.hh14Rep, bpair,
                              in.nonbonded.alphanb, zetac6dummy, env)) {
      ADD_FAILURE();
      return res;
    }
  }
  in.nonbonded.bpair = bpair;
  HbParams hp;
  hp.hbacut = param.hbAngleCut;
  hp.hblongcut = param.hbLongCut;
  hp.hbscut = param.hbShortCut;
  hp.hbalp = param.hbAlp;
  hp.hbst = param.hbSt;
  hp.hbsf = param.hbSf;
  hp.hbabmix = param.hbAbMix;
  hp.hbnbcut = param.hbNbCut;
  hp.xhaciGlobAbh = param.xhAciGlobAbH;
  hp.xhaciCoh = param.xhAciCoh;
  hp.torsHb = param.torsHb;
  hp.bendHb = param.bendHb;

  if (!gfnffSinglePoint(in, param, hp, res, env)) {
    std::printf("driver failed: %s\n", env.errorMessage().c_str());
    ADD_FAILURE(); return res;
  }
  return res;
}

TEST(DriverTest, EthaneMatchesReference)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  DriverResult res = runEthane(param, gen);
  EXPECT_NEAR(res.ees, -5.13736317619923957e-02, 1.0e-12);
  EXPECT_NEAR(res.edisp, -2.29998294950890759e-03, 1.0e-12);
  EXPECT_NEAR(res.erep, 5.25599635356355926e-02, 1.0e-12);
  EXPECT_NEAR(res.ebond, -1.10139334830457880e+00, 1.0e-12);
  EXPECT_NEAR(res.eangl, 2.30053848513024004e-02, 1.0e-12);
  EXPECT_NEAR(res.etors, 2.41297575568995715e-03, 1.0e-12);
  EXPECT_NEAR(res.ehb, 0.00000000000000000e+00, 1.0e-12);
  EXPECT_NEAR(res.exb, 0.00000000000000000e+00, 1.0e-12);
  EXPECT_NEAR(res.ebatm, -1.24854164952983531e-05, 1.0e-12);
  EXPECT_NEAR(res.etot, -1.07710112428994775e+00, 1.0e-12);
  EXPECT_NEAR(res.gnorm, 9.75009278767578619e-02, 1.0e-12);
    const double kGradRef[24] = {
    -1.07930251396562444e-02, -2.81662743303668119e-02, -5.35508801828017203e-02,
    2.05734134940506914e-02, -1.61396026342733805e-04, 4.06986837899452059e-02,
    -2.21426369343719629e-03, 1.56303861630683709e-02, 3.37997447786058738e-02,
    -9.28398788005839600e-04, 1.62295920373335088e-02, 1.35655023131399082e-02,
    1.29550893346365367e-02, -7.36771796375722005e-03, 1.38366865727112230e-03,
    -3.02318097733839108e-03, -3.45415708746263329e-03, -2.82855369275798914e-02,
    -1.48884475453210062e-02, -1.08695817841879191e-02, 1.79699825371529764e-03,
    -1.68118668492855819e-03, 1.81591489917154597e-02, -9.40818068229578447e-03,
  };
  ASSERT_EQ((int)res.gradient.size(), 24);
  for (int i = 0; i < 24; ++i)
    EXPECT_NEAR(res.gradient[i], kGradRef[i], 1.0e-12);
    const double kChargesRef[8] = {
    -3.89157711713968291e-01, 1.62308477311862626e-01, 2.48430970949445502e-03,
    2.68486007433102247e-01, 1.22751199393052901e-01, -1.47599647701004894e-01,
    5.52299245171360281e-02, -7.45025589496750573e-02,
  };
  ASSERT_EQ((int)res.charges.size(), 8);
  for (int i = 0; i < 8; ++i)
    EXPECT_NEAR(res.charges[i], kChargesRef[i], 1.0e-12);
  EXPECT_NEAR(res.dipole[0], -4.05650693051311184e-01, 1.0e-12);
  EXPECT_NEAR(res.dipole[1], -5.94256409606938152e-01, 1.0e-12);
  EXPECT_NEAR(res.dipole[2], 6.03369957231112197e-01, 1.0e-12);
}
} // namespace
} // namespace Xtb
} // namespace Avogadro
