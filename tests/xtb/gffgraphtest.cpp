/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "avogadro/xtb/COPYING.LESSER".
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/xtb/environment.h>
#include <avogadro/xtb/gfnffdata.h>
#include <avogadro/xtb/gfnffparams.h>
#include <avogadro/xtb/gfnffsetup.h>
#include <avogadro/xtb/gffgraph.h>

#include <cmath>
#include <vector>

namespace Avogadro {
namespace Xtb {

// Reference values produced by compiling the unmodified xtb sources
// (mrec.f90, ini2.F90 getring36/chkrng, gfnff_eg.f90 gfnff_dlogcoord) with
// gfortran. Integers match exactly; doubles use 1e-9 tolerances
// (linear-solver headroom where applicable).

TEST(GffGraphTest, fragments)
{
  std::vector<int> frag;
  // Two separated H2 molecules.
  EXPECT_EQ(findFragments({ { 1 }, { 0 }, { 3 }, { 2 } }, frag), 2);
  EXPECT_EQ(frag, std::vector<int>({ 1, 1, 2, 2 }));
  // Water is connected.
  EXPECT_EQ(findFragments({ { 1, 2 }, { 0 }, { 0 } }, frag), 1);
  // Periodic variant agrees on molecular input.
  EXPECT_EQ(findFragmentsPbc(
              { { { 1, 0 }, { 2, 0 } }, { { 0, 0 } }, { { 0, 0 } } },
              1, frag),
            1);
  EXPECT_EQ(findFragmentsPbc(
              { { { 1, 0 } }, { { 0, 0 } }, { { 3, 0 } }, { { 2, 0 } } },
              1, frag),
            2);
  EXPECT_EQ(frag, std::vector<int>({ 1, 1, 2, 2 }));
}

TEST(GffGraphTest, benzeneRings)
{
  // Benzene C ring with H leaves, 0-based.
  std::vector<int> numbers(12, 1);
  for (int i = 0; i < 6; ++i)
    numbers[i] = 6;
  std::vector<std::vector<int>> nb(12);
  for (int i = 0; i < 6; ++i)
    nb[i] = { (i + 1) % 6, (i + 5) % 6, 6 + i };
  for (int i = 6; i < 12; ++i)
    nb[i] = { i - 6 };
  Environment env;
  const std::vector<std::vector<int>> expected = {
    { 1, 2, 3, 4, 5, 0 }, { 2, 3, 4, 5, 0, 1 }, { 3, 4, 5, 0, 1, 2 },
    { 4, 5, 0, 1, 2, 3 }, { 5, 0, 1, 2, 3, 4 }, { 4, 3, 2, 1, 0, 5 }
  };
  for (int a0 = 0; a0 < 6; ++a0) {
    std::vector<Ring> rings;
    ASSERT_TRUE(findRingsThrough(12, numbers, nb, a0, rings, env));
    ASSERT_EQ(rings.size(), 1);
    EXPECT_EQ(rings[0].members, expected[a0]);
    EXPECT_EQ(rings[0].hetero, 0);
    EXPECT_EQ(smallestRingThrough(rings), 6);
  }
  std::vector<Ring> hydrogens;
  ASSERT_TRUE(findRingsThrough(12, numbers, nb, 6, hydrogens, env));
  EXPECT_TRUE(hydrogens.empty());
  EXPECT_EQ(smallestRingThrough(hydrogens), 99);
  // C-C bond of the ring is in a 6-ring; C-H bond in none.
  std::vector<Ring> r0, r6;
  ASSERT_TRUE(findRingsThrough(12, numbers, nb, 0, r0, env));
  ASSERT_TRUE(findRingsThrough(12, numbers, nb, 6, r6, env));
  EXPECT_EQ(smallestRingBond(r0, r6, 0, 6), 0);
  std::vector<Ring> r1;
  ASSERT_TRUE(findRingsThrough(12, numbers, nb, 1, r1, env));
  EXPECT_EQ(smallestRingBond(r0, r1, 0, 1), 6);
}

TEST(GffGraphTest, cyclopropaneDedup)
{
  // Triangle C0-C1-C2 with H leaves, 0-based. The triangle is found in
  // both rotations; dedup must merge them into exactly one ring per
  // carbon (verified against real getring36).
  std::vector<int> numbers = { 6, 1, 1, 6, 1, 1, 6, 1, 1 };
  std::vector<std::vector<int>> nb = { { 1, 2, 3, 6 }, { 0 }, { 0 },
                                       { 0, 4, 5, 6 }, { 3 }, { 3 },
                                       { 0, 3, 7, 8 }, { 6 }, { 6 } };
  Environment env;
  for (int a0 : { 0, 3, 6 }) {
    std::vector<Ring> rings;
    ASSERT_TRUE(findRingsThrough(9, numbers, nb, a0, rings, env));
    ASSERT_EQ(rings.size(), 1u);
    EXPECT_EQ(smallestRingThrough(rings), 3);
  }
  std::vector<Ring> r0, r3, rh;
  ASSERT_TRUE(findRingsThrough(9, numbers, nb, 0, r0, env));
  EXPECT_EQ(r0[0].members, std::vector<int>({ 3, 6, 0 }));
  ASSERT_TRUE(findRingsThrough(9, numbers, nb, 3, r3, env));
  EXPECT_EQ(smallestRingBond(r0, r3, 0, 3), 3);
  ASSERT_TRUE(findRingsThrough(9, numbers, nb, 1, rh, env));
  EXPECT_TRUE(rh.empty());
  EXPECT_EQ(smallestRingBond(r0, rh, 0, 1), 0);
}

TEST(GffGraphTest, waterCoordination)
{
  const int n = 3;
  std::vector<int> numbers = { 8, 1, 1 };
  std::vector<double> xyz = { 0.0, 0.0, 0.0, 1.4305226400, 1.1071846211,
                              0.0, -1.4305226400, 1.1071846211, 0.0 };
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  std::vector<double> rab(n * (n + 1) / 2, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      double dx = xyz[3 * j] - xyz[3 * i];
      double dy = xyz[3 * j + 1] - xyz[3 * i + 1];
      double dz = xyz[3 * j + 2] - xyz[3 * i + 2];
      rab[j + i * (i + 1) / 2] = std::sqrt(dx * dx + dy * dy + dz * dz);
    }
  }
  std::vector<double> rcov(103);
  for (int i = 0; i < 103; ++i)
    rcov[i] = ((covalentRadD3Angstrom[i] * xtbAngstromToBohr) * 4.0) / 3.0;
  std::vector<double> logCn, dlogCn, mchar;
  gffCoordinationNumber(n, numbers, xyz, rab, rcov, param.cnMax, 40.0, logCn,
                        dlogCn);
  ASSERT_EQ(logCn.size(), 3);
  EXPECT_NEAR(logCn[0], 1.916588025003, 1e-9);
  EXPECT_NEAR(logCn[1], 0.974742821016, 1e-9);
  EXPECT_NEAR(logCn[2], 0.974742821016, 1e-9);
  EXPECT_NEAR(dlogCn[(0 * n + 1) * 3 + 1], 0.036501495507, 1e-9);
  std::vector<double> en(103);
  for (int i = 0; i < 103; ++i)
    en[i] = gffEn[i];
  metallicCharacter(n, numbers, en, logCn, dlogCn, mchar);
  ASSERT_EQ(mchar.size(), 3);
  EXPECT_NEAR(mchar[0], 1.637869231316e-44, 1e-54);
  EXPECT_NEAR(mchar[1], 0.003885248757, 1e-9);
  EXPECT_NEAR(mchar[2], 0.003885248757, 1e-9);
}

} // namespace Xtb
} // namespace Avogadro
