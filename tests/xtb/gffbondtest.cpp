/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "avogadro/xtb/COPYING.LESSER".
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/xtb/environment.h>
#include <avogadro/xtb/gfnffbonds.h>
#include <avogadro/xtb/gfnffdata.h>
#include <avogadro/xtb/gfnffhyb.h>
#include <avogadro/xtb/gfnffparams.h>
#include <avogadro/xtb/gfnffrab.h>
#include <avogadro/xtb/gfnffsetup.h>
#include <avogadro/xtb/gfnfftopo.h>

#include <cmath>
#include <vector>

namespace Avogadro {
namespace Xtb {

// Reference values cross-checked against an independent Python
// implementation of the same Fortran sources, with amide() additionally
// verified against the compiled xtb routine. Neighbor lists come from the
// verified fillNeighborList; hybridizations from the verified
// assignHybridization.

struct BondCase
{
  std::vector<int> numbers;
  std::vector<double> xyz;
  std::vector<std::pair<int, int>> bonds; // 0-based, first-seen order
  std::vector<int> piAtoms; // 0-based pi atoms
  std::vector<int> types;
  std::vector<int> amides; // 0-based amide N atoms
};

static BondCase runCase(const std::vector<int>& numbers,
                        const std::vector<double>& xyz)
{
  BondCase out;
  int n = static_cast<int>(numbers.size());
  GffData param;
  GffGenerator gen;
  EXPECT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  std::vector<double> cn(n);
  for (int i = 0; i < n; ++i)
    cn[i] = param.normCn[numbers[i] - 1];
  std::vector<double> rtmp;
  gfnffBondGuesses(n, numbers, cn, rtmp);
  std::vector<double> qa(n, 0.0);
  std::vector<int> metal(103), group(103);
  for (int i = 0; i < 103; ++i) {
    metal[i] = gffMetal[i];
    group[i] = gffGroup[i];
  }
  scaleBondGuesses(n, numbers, qa, metal, gen.rShrink, rtmp);
  std::vector<double> dist(n * n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double dx = xyz[3 * j] - xyz[3 * i];
      double dy = xyz[3 * j + 1] - xyz[3 * i + 1];
      double dz = xyz[3 * j + 2] - xyz[3 * i + 2];
      dist[i * n + j] = std::sqrt(dx * dx + dy * dy + dz * dz);
    }
  }
  std::vector<double> mchar(n, 0.0);
  std::vector<int> nbf, nbfc, nb, nbc, nbm, nbmc, full(n, 0);
  fillNeighborList(n, numbers, rtmp, dist, mchar.data(), 1, 1.25, 1.0,
                   gffMetal.data(), gffGroup.data(), gffNormCn.data(),
                   full.data(), 1, nbf, nbfc);
  for (int i = 0; i < n; ++i)
    full[i] = nbfc[i];
  fillNeighborList(n, numbers, rtmp, dist, mchar.data(), 2, 1.25, 1.0,
                   gffMetal.data(), gffGroup.data(), gffNormCn.data(),
                   full.data(), 1, nb, nbc);
  fillNeighborList(n, numbers, rtmp, dist, mchar.data(), 3, 1.25, 1.0,
                   gffMetal.data(), gffGroup.data(), gffNormCn.data(),
                   full.data(), 1, nbm, nbmc);
  Environment env;
  std::vector<int> hyb, itag;
  ASSERT_TRUE(assignHybridization(n, numbers, xyz, nbf, nbfc, nb, nbc, nbm,
                                  nbmc, qa, metal.data(), group.data(),
                                  160.0, 1, true, hyb, itag, env));
  // Flat nb view for the helpers below.
  std::vector<std::vector<int>> nbv(n);
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < nbc[i]; ++k)
      nbv[i].push_back(nb[i * maxNeighbors + k]);
  }
  std::vector<Bond> bonds = buildBondList(n, nb, nbc, 1);
  for (const Bond& b : bonds)
    out.bonds.emplace_back(b.first, b.second);
  PiSystem pi = buildPiSystem(n, numbers, hyb, nb, nbc, 1);
  for (int p : pi.piIndex)
    out.piAtoms.push_back(p);
  std::vector<int> imetal = effectiveMetals(n, numbers, nbc, 1, metal.data(),
                                            group.data());
  out.types = assignBondTypes(bonds, numbers, hyb, itag, pi.atomPi, imetal,
                              group.data());
  for (int i = 0; i < n; ++i) {
    if (isAmideNitrogen(n, numbers, hyb, nb, nbc, 1, pi.atomPi, i))
      out.amides.push_back(i);
  }
  return out;
}

TEST(GffBondTest, water)
{
  BondCase c = runCase({ 8, 1, 1 },
                       { 0.0, 0.0, 0.0, 1.4305226400, 1.1071846211, 0.0,
                         -1.4305226400, 1.1071846211, 0.0 });
  EXPECT_EQ(c.bonds, std::vector<std::pair<int, int>>({ { 1, 0 }, { 2, 0 } }));
  EXPECT_TRUE(c.piAtoms.empty());
  EXPECT_EQ(c.types, std::vector<int>({ 1, 1 }));
  EXPECT_TRUE(c.amides.empty());
}

TEST(GffBondTest, ethene)
{
  BondCase c = runCase({ 6, 6, 1, 1, 1, 1 },
                       { 0.0, 0.0, 0.0, 2.53, 0.0, 0.0, -1.025, 1.775, 0.0,
                         -1.025, -1.775, 0.0, 3.555, 1.775, 0.0, 3.555,
                         -1.775, 0.0 });
  EXPECT_EQ(c.bonds.size(), 5);
  EXPECT_EQ(c.types, std::vector<int>({ 2, 1, 1, 1, 1 }));
  EXPECT_EQ(c.piAtoms, std::vector<int>({ 0, 1 }));
}

TEST(GffBondTest, benzene)
{
  BondCase c = runCase({ 6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1 },
                       { 2.627, 0.0, 0.0, 1.3135, 2.275, 0.0, -1.3135, 2.275,
                         0.0, -2.627, 0.0, 0.0, -1.3135, -2.275, 0.0, 1.3135,
                         -2.275, 0.0, 4.687, 0.0, 0.0, 2.3435, 4.059, 0.0,
                         -2.3435, 4.059, 0.0, -4.687, 0.0, 0.0, -2.3435,
                         -4.059, 0.0, 2.3435, -4.059, 0.0 });
  EXPECT_EQ(c.bonds.size(), 12);
  EXPECT_EQ(c.types, std::vector<int>({ 2, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1 }));
  EXPECT_EQ(c.piAtoms, std::vector<int>({ 0, 1, 2, 3, 4, 5 }));
}

TEST(GffBondTest, formamide)
{
  BondCase c = runCase({ 6, 8, 7, 1, 1, 1 },
                       { 0.0, 0.0, 0.0, 2.28, 0.0, 0.0, -1.30, 2.25, 0.0,
                         -1.00, -1.80, 0.0, -0.337, 3.918, 0.0, -2.84, 1.09,
                         0.0 });
  EXPECT_EQ(c.bonds.size(), 5);
  EXPECT_EQ(c.types, std::vector<int>({ 2, 2, 1, 1, 1 }));
  EXPECT_EQ(c.amides, std::vector<int>({ 2 }));
}

} // namespace Xtb
} // namespace Avogadro
