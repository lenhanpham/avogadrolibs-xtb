/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "avogadro/xtb/COPYING.LESSER".
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/xtb/environment.h>
#include <avogadro/xtb/gfnffdata.h>
#include <avogadro/xtb/gfnffhb.h>
#include <avogadro/xtb/gfnffparams.h>
#include <avogadro/xtb/gfnffsetup.h>

#include <utility>
#include <vector>

namespace Avogadro {
namespace Xtb {

// Reference values verified bitwise: perception, triplet and map builders,
// HB erf coordination, geometry lists and the bonded HB correction against
// verbatim copies of the ini HAB/HV blocks, the real ini2 bond_hbset(0),
// bond_hb_AHB_set(0/1), hbonds, xatom, amideH and nbondmat_pbc routines and
// the real eg egbond_hb/dncoord_erf, compiled with gfortran
// -ffp-contract=off (64 sections, 310 integer + 1664 float values over a
// water dimer, a formamide dimer, a bifurcated water trimer and a
// CH3Cl/H2O halogen-bonded pair, 0 mismatches; dumped bpair also matches
// real nbondmat_pbc on all inputs).

namespace {

struct HbFixture
{
  std::vector<int> numbers;
  std::vector<double> xyz;
  std::vector<double> qa;
  std::vector<int> hyb, itag, imetal, pi, counts, firstNb;
  std::vector<std::vector<int>> neighbours;
  std::vector<int> bpair; // pair(a,b) at b*n+a
  std::vector<std::pair<int, int>> blist; // (jj,ii) reference order
};

HbFixture makeWaterDimer()
{
  HbFixture f;
  f.numbers = { 8, 1, 1, 8, 1, 1 };
  f.xyz = { 0.0, 0.0, 0.0, 1.75, 0.0, 0.0, -0.5, 1.68, 0.0,
            5.0, 0.2, 0.0, 6.2, 1.55, 0.0, 5.8, -1.1, 0.0 };
  f.qa = { -2.02557504182036718e-01, 1.23846691413292243e-01,
           9.84319142931487645e-02, -2.44111700848554630e-01,
           9.39810468292040374e-02, 1.30409552494946318e-01 };
  f.hyb = { 3, 0, 0, 3, 0, 0 };
  f.itag = { 0, 0, 0, 0, 0, 0 };
  f.imetal = { 0, 0, 0, 0, 0, 0 };
  f.pi = { 0, 0, 0, 0, 0, 0 };
  f.neighbours = { { 1, 2 }, { 0 }, { 0 }, { 4, 5 }, { 3 }, { 3 } };
  f.counts = { 2, 1, 1, 2, 1, 1 };
  f.firstNb = { 1, 0, 0, 4, 3, 3 };
  f.bpair = { 0, 1, 1, 5, 5, 5, 1, 0, 2, 5, 5, 5, 1, 2, 0, 5, 5, 5,
              5, 5, 5, 0, 1, 1, 5, 5, 5, 1, 0, 2, 5, 5, 5, 1, 2, 0 };
  f.blist = { { 1, 0 }, { 2, 0 }, { 4, 3 }, { 5, 3 } };
  return f;
}

void groupTable(std::vector<int>& group)
{
  group.assign(103, 0);
  for (int i = 0; i < 103; ++i)
    group[i] = gffGroup[i];
}

TEST(HbSetTest, WaterDimerPerception)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  HbFixture f = makeWaterDimer();
  std::vector<int> group;
  groupTable(group);
  Environment env;
  std::vector<double> hbBas, hbAci;
  ASSERT_TRUE(hbBasicity(6, f.numbers, f.itag, f.neighbours, f.counts,
                         param.xhBas, hbBas, env));
  ASSERT_TRUE(hbAcidity(6, f.numbers, f.neighbours, f.hyb, f.pi, param.xhAci,
                        hbAci, env));
  EXPECT_NEAR(hbBas[0], 6.70000000000000040e-01, 0.0);
  EXPECT_NEAR(hbBas[1], 0.0, 0.0);
  EXPECT_NEAR(hbAci[0], 1.50000000000000000e+00, 0.0);
  HbDonorLists lists;
  ASSERT_TRUE(hbDonorLists(6, f.numbers, f.hyb, f.pi, f.qa, f.neighbours,
                           f.firstNb, group.data(), param.xhBas, gen.hQaThr,
                           gen.qaBThr, f.bpair, hbBas, hbAci, lists, env));
  ASSERT_EQ(lists.hatH.size(), 4u);
  EXPECT_EQ(lists.hatH[0], 1);
  EXPECT_EQ(lists.hatH[1], 2);
  EXPECT_EQ(lists.hatH[2], 4);
  EXPECT_EQ(lists.hatH[3], 5);
  ASSERT_EQ(lists.hatAB.size(), 1u);
  EXPECT_EQ(lists.hatAB[0].first, 3);
  EXPECT_EQ(lists.hatAB[0].second, 0);
  EXPECT_TRUE(lists.xbatA.empty());
}

TEST(HbSetTest, WaterDimerTripletsAndMaps)
{
  HbFixture f = makeWaterDimer();
  Environment env;
  std::vector<HbBondTriplet> triplets;
  ASSERT_TRUE(hbBondTriplets(6, { { 3, 0 } }, { 1, 2, 4, 5 }, f.xyz, 200.0,
                             f.bpair, triplets, env));
  ASSERT_EQ(triplets.size(), 4u);
  EXPECT_EQ(triplets[0].a, 0);
  EXPECT_EQ(triplets[0].b, 3);
  EXPECT_EQ(triplets[0].h, 1);
  EXPECT_EQ(triplets[3].a, 3);
  EXPECT_EQ(triplets[3].b, 0);
  EXPECT_EQ(triplets[3].h, 5);
  HbBondMaps maps;
  ASSERT_TRUE(hbAhbMaps(6, f.numbers, f.blist, triplets, maps, env));
  ASSERT_EQ(maps.ahA.size(), 4u);
  EXPECT_EQ(maps.ahA[0], 0);
  EXPECT_EQ(maps.ahH[0], 1);
  EXPECT_EQ(maps.ahA[3], 3);
  EXPECT_EQ(maps.ahH[3], 5);
  ASSERT_EQ(maps.bAtoms.size(), 4u);
  for (const auto& b : maps.bAtoms)
    ASSERT_EQ(b.size(), 1u);
  EXPECT_EQ(maps.bAtoms[0][0], 3);
  EXPECT_EQ(maps.bAtoms[3][0], 0);
  EXPECT_EQ(maps.bMax, 1);
  EXPECT_EQ(maps.mapNab, 2);
  EXPECT_EQ(maps.mapNh, 4);
  ASSERT_EQ(maps.nrHb.size(), 4u);
  for (int v : maps.nrHb)
    EXPECT_EQ(v, 1);
  for (char v : maps.isAbh)
    EXPECT_EQ((int)v, 1);
  ASSERT_EQ(maps.mapAbh.size(), 6u);
  EXPECT_EQ(maps.mapAbh[0], 0);
  EXPECT_EQ(maps.mapAbh[1], 0);
  EXPECT_EQ(maps.mapAbh[2], 1);
  EXPECT_EQ(maps.mapAbh[3], 1);
  EXPECT_EQ(maps.mapAbh[4], 2);
  EXPECT_EQ(maps.mapAbh[5], 3);
}

TEST(HbSetTest, BifurcatedGrouping)
{
  Environment env;
  // Bifurcated water trimer: donor H(1) shared by two acceptors.
  std::vector<int> numbers = { 8, 1, 1, 8, 1, 1, 8, 1, 1 };
  std::vector<std::pair<int, int>> bonds = { { 1, 0 }, { 2, 0 }, { 4, 3 },
                                             { 5, 3 }, { 5, 4 }, { 7, 5 },
                                             { 7, 6 }, { 8, 6 }, { 8, 7 } };
  std::vector<HbBondTriplet> triplets = { { 0, 3, 1 }, { 0, 3, 2 },
                                          { 3, 0, 5 }, { 0, 6, 1 },
                                          { 0, 6, 2 }, { 6, 0, 7 },
                                          { 3, 6, 5 }, { 6, 3, 7 } };
  HbBondMaps maps;
  ASSERT_TRUE(hbAhbMaps(9, numbers, bonds, triplets, maps, env));
  EXPECT_EQ(maps.bMax, 2);
  ASSERT_EQ(maps.ahA.size(), 4u);
  // Groups follow blist order: (0,1) first with B = {3, 6}.
  EXPECT_EQ(maps.ahA[0], 0);
  EXPECT_EQ(maps.ahH[0], 1);
  ASSERT_EQ(maps.bAtoms[0].size(), 2u);
  EXPECT_EQ(maps.bAtoms[0][0], 3);
  EXPECT_EQ(maps.bAtoms[0][1], 6);
  EXPECT_EQ(maps.ahA[2], 3);
  EXPECT_EQ(maps.ahH[2], 5);
  ASSERT_EQ(maps.bAtoms[2].size(), 2u);
  EXPECT_EQ(maps.bAtoms[2][0], 0);
  EXPECT_EQ(maps.bAtoms[2][1], 6);
  ASSERT_EQ(maps.nrHb.size(), 9u);
  std::vector<int> expect = { 2, 2, 2, 2, 0, 0, 2, 2, 0 };
  EXPECT_EQ(maps.nrHb, expect);
}

TEST(HbSetTest, HalogenTriple)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  std::vector<int> group;
  groupTable(group);
  Environment env;
  // CH3Cl/H2O: Cl(1) on C(0), O(5) acceptor.
  std::vector<int> numbers = { 6, 17, 1, 1, 1, 8, 1, 1 };
  std::vector<std::vector<int>> neighbours = { { 1, 2, 3, 4 }, { 0 }, { 0 },
                                               { 0 },           { 0 },
                                               { 6, 7 },        { 5, 7 },
                                               { 5, 6 } };
  std::vector<int> bpair = { 0, 1, 1, 1, 1, 5, 5, 5, 1, 0, 2, 2, 2, 5, 5, 5,
                             1, 2, 0, 2, 2, 5, 5, 5, 1, 2, 2, 0, 2, 5, 5, 5,
                             1, 2, 2, 2, 0, 5, 5, 5, 5, 5, 5, 5, 5, 0, 1, 1,
                             5, 5, 5, 5, 5, 1, 0, 1, 5, 5, 5, 5, 5, 1, 1, 0 };
  std::vector<int> hyb(8, 0), pi(8, 0), firstNb(8, 0);
  std::vector<double> qa(8, 0.0), hbBas(8, 0.0), hbAci(8, 0.0);
  HbDonorLists lists;
  ASSERT_TRUE(hbDonorLists(8, numbers, hyb, pi, qa, neighbours, firstNb,
                           group.data(), param.xhBas, gen.hQaThr, gen.qaBThr,
                           bpair, hbBas, hbAci, lists, env));
  ASSERT_EQ(lists.xbatA.size(), 1u);
  EXPECT_EQ(lists.xbatA[0], 0);
  EXPECT_EQ(lists.xbatB[0], 5);
  EXPECT_EQ(lists.xbatX[0], 1);
  EXPECT_TRUE(isXAtom(17));
  EXPECT_TRUE(isXAtom(15));
  EXPECT_FALSE(isXAtom(8));
  // Geometry list keeps the A-X...B entry within hbThr2.
  std::vector<double> xyz = { 0.0, 0.0, 0.0, 3.4, 0.0, 0.0, -1.1, 1.1, 1.1,
                              -1.1, -1.1, 1.1, -1.1, 0.0, -1.55, 9.2, 0.6,
                              0.0, 10.4, 1.5, 0.0, 10.4, -0.3, 0.0 };
  std::vector<HbTriple> hb1, hb2;
  std::vector<XbTriple> xb;
  ASSERT_TRUE(hbTripletLists(8, xyz, {}, {}, lists, bpair, 200.0, 400.0, hb1,
                             hb2, xb, env));
  ASSERT_EQ(xb.size(), 1u);
  EXPECT_EQ(xb[0].a, 0);
  EXPECT_EQ(xb[0].b, 5);
  EXPECT_EQ(xb[0].x, 1);
}

TEST(HbSetTest, WaterListsAndErf)
{
  HbFixture f = makeWaterDimer();
  Environment env;
  HbDonorLists donors;
  donors.hatH = { 1, 2, 4, 5 };
  donors.hatAB = { { 3, 0 } };
  std::vector<HbTriple> hb1, hb2;
  std::vector<XbTriple> xb;
  ASSERT_TRUE(hbTripletLists(6, f.xyz, donors.hatAB, donors.hatH, donors,
                             f.bpair, 200.0, 400.0, hb1, hb2, xb, env));
  EXPECT_TRUE(hb1.empty());
  ASSERT_EQ(hb2.size(), 4u);
  EXPECT_EQ(hb2[0].a, 0);
  EXPECT_EQ(hb2[0].b, 3);
  EXPECT_EQ(hb2[0].h, 1);
  EXPECT_EQ(hb2[3].a, 3);
  EXPECT_EQ(hb2[3].b, 0);
  EXPECT_EQ(hb2[3].h, 5);
  EXPECT_TRUE(xb.empty());
  std::vector<double> rcov(103);
  for (int z = 0; z < 103; ++z)
    rcov[z] = covalentRadiusD3(z + 1);
  HbBondMaps maps;
  std::vector<HbBondTriplet> triplets = { { 0, 3, 1 }, { 0, 3, 2 },
                                          { 3, 0, 4 }, { 3, 0, 5 } };
  ASSERT_TRUE(hbAhbMaps(6, f.numbers, f.blist, triplets, maps, env));
  std::vector<double> hbCn, hbDcn;
  ASSERT_TRUE(hbErfCoordination(6, f.numbers, f.xyz, rcov, maps, hbCn, hbDcn,
                               env));
  ASSERT_EQ(hbCn.size(), 6u);
  EXPECT_NEAR(hbCn[0], 0.0, 0.0);
  EXPECT_NEAR(hbCn[1], 1.00000000000000000e+00, 0.0);
  EXPECT_NEAR(hbCn[3], 1.00000000000000000e+00, 0.0);
  ASSERT_EQ(hbDcn.size(), 108u);
  EXPECT_NEAR(hbDcn[0], 5.41626386608350461e-49, 0.0);
  EXPECT_NEAR(hbDcn[1], -1.02722245736066480e-49, 0.0);
  EXPECT_NEAR(hbDcn[2], 0.0, 0.0);
}

TEST(HbSetTest, EgbondHbWaterBond0)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  HbFixture f = makeWaterDimer();
  Environment env;
  std::vector<double> rcov(103);
  for (int z = 0; z < 103; ++z)
    rcov[z] = covalentRadiusD3(z + 1);
  HbBondMaps maps;
  std::vector<HbBondTriplet> triplets = { { 0, 3, 1 }, { 0, 3, 2 },
                                          { 3, 0, 4 }, { 3, 0, 5 } };
  ASSERT_TRUE(hbAhbMaps(6, f.numbers, f.blist, triplets, maps, env));
  std::vector<double> hbCn, hbDcn;
  ASSERT_TRUE(hbErfCoordination(6, f.numbers, f.xyz, rcov, maps, hbCn, hbDcn,
                               env));
  // Bond 0 (H-O donor): dumped rab setup.
  std::vector<double> drij = { 9.58593083116912603e-03,
                               7.17879645201998078e-03,
                               0.00000000000000000e+00,
                               -1.15460595226169523e-02,
                               1.08563308457233006e-05,
                               0.00000000000000000e+00,
                               2.13654653375332671e-03,
                               -7.17879635338054389e-03,
                               0.00000000000000000e+00,
                               -1.76417842305503261e-04,
                               -1.08564294851598191e-05,
                               0.00000000000000000e+00,
                               -3.09155527700356057e-69,
                               -7.72888819250890143e-70,
                               0.00000000000000000e+00,
                               -1.24932345898831080e-53,
                               2.36940656020695968e-54,
                               0.00000000000000000e+00 };
  double dd[2] = { 2.60330908952114304e-01, 1.53426951480975071e-01 };
  std::vector<double> g(18, 0.0), dEdcn(6, 0.0);
  std::vector<char> seen(16, 0);
  double e = 0.0;
  ASSERT_TRUE(egbondHb(0, 0, 1, 1.75000000000000000e+00,
                       1.62580462201213471e+00, drij, dd, f.numbers, hbCn,
                       hbDcn, f.xyz, param.vbondScale,
                       6.49706271153752790e-01, -1.38077009963469860e-01,
                       maps, seen, dEdcn, e, g, env));
  EXPECT_NEAR(e, -1.36837258114995586e-01, 0.0);
  EXPECT_NEAR(g[0], -2.00651605518891545e-02, 0.0);
  EXPECT_NEAR(g[1], -1.42676020911390106e-04, 0.0);
  EXPECT_NEAR(g[2], 0.0, 0.0);
  EXPECT_NEAR(g[3], 2.01041174110283133e-02, 0.0);
  EXPECT_NEAR(dEdcn[0], -5.17398403169401847e-03, 0.0);
}

TEST(HbSetTest, StaleNrHbWithoutGroups)
{
  Environment env;
  // Formamide-dimer topology: no H qualifies, so no triplets or groups,
  // but N/O-H bonds keep stale running counts (verbatim quirk).
  std::vector<int> numbers = { 6, 8, 7, 1, 1, 1, 6, 8, 7, 1, 1, 1 };
  std::vector<std::pair<int, int>> bonds = { { 1, 0 }, { 2, 0 }, { 3, 0 },
                                             { 4, 2 }, { 5, 2 }, { 5, 4 },
                                             { 7, 6 }, { 8, 6 }, { 9, 6 },
                                             { 10, 8 }, { 11, 8 }, { 11, 10 } };
  HbBondMaps maps;
  ASSERT_TRUE(hbAhbMaps(12, numbers, bonds, {}, maps, env));
  EXPECT_TRUE(maps.ahA.empty());
  EXPECT_EQ(maps.bMax, 1);
  ASSERT_EQ(maps.nrHb.size(), 12u);
  std::vector<int> expect = { 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0 };
  EXPECT_EQ(maps.nrHb, expect);
  // Participation flags are still set for every N/O-H bond.
  EXPECT_EQ((int)maps.isAbh[2], 1);
  EXPECT_EQ((int)maps.isAbh[4], 1);
  EXPECT_EQ((int)maps.isAbh[0], 0);
}

} // namespace

} // namespace Xtb
} // namespace Avogadro
