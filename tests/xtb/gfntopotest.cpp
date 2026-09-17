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
#include <avogadro/xtb/gfnfftopo.h>

#include <cmath>
#include <vector>

namespace Avogadro {
namespace Xtb {

// Reference values produced by compiling the unmodified xtb sources
// (neighbor.f90 fillnb via TNeigh%get_nb, ini2.F90 goedeckera) with
// gfortran. Neighbor lists match exactly (integers); charges/energies use
// 1e-9 tolerances (linear-solver rounding headroom).

struct Water
{
  std::vector<int> numbers = { 8, 1, 1 };
  std::vector<double> xyz = { 0.0, 0.0, 0.0, 1.4305226400, 1.1071846211,
                              0.0, -1.4305226400, 1.1071846211, 0.0 };
};

static void buildInputs(int n, const std::vector<int>& numbers,
                        const std::vector<double>& xyz,
                        std::vector<double>& dist, std::vector<double>& rad)
{
  dist.assign(n * n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double dx = xyz[3 * j] - xyz[3 * i];
      double dy = xyz[3 * j + 1] - xyz[3 * i + 1];
      double dz = xyz[3 * j + 2] - xyz[3 * i + 2];
      dist[i * n + j] = std::sqrt(dx * dx + dy * dy + dz * dz);
    }
  }
  rad.assign(n * (n + 1) / 2, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j)
      rad[packedIndex(j, i)] =
        covalentRadiusD3(numbers[i]) + covalentRadiusD3(numbers[j]);
  }
}

TEST(GfnTopoTest, packedIndex)
{
  EXPECT_EQ(packedIndex(0, 0), 0);
  EXPECT_EQ(packedIndex(0, 1), 1);
  EXPECT_EQ(packedIndex(1, 0), 1);
  EXPECT_EQ(packedIndex(1, 1), 2);
  EXPECT_EQ(packedIndex(0, 2), 3);
  EXPECT_EQ(packedIndex(2, 2), 5);
  // involution spot check: lin(j,i) == lin(i,j)
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j)
      EXPECT_EQ(packedIndex(i, j), packedIndex(j, i));
  }
}

TEST(GfnTopoTest, waterNeighborLists)
{
  Water mol;
  const int n = 3;
  std::vector<double> dist, rad;
  buildInputs(n, mol.numbers, mol.xyz, dist, rad);
  std::vector<double> mchar(n, 0.0);
  std::vector<int> full(n, 0), list, counts;
  // icase 1: O has 2 neighbours, each H has 1.
  fillNeighborList(n, mol.numbers, rad, dist, mchar.data(), 1, 1.25, 1.0,
                   gffMetal.data(), gffGroup.data(), gffNormCn.data(),
                   full.data(), 1, list, counts);
  ASSERT_EQ(counts, std::vector<int>({ 2, 1, 1 }));
  EXPECT_EQ(list[0 * maxNeighbors + 0], 1);
  EXPECT_EQ(list[0 * maxNeighbors + 1], 2);
  EXPECT_EQ(list[1 * maxNeighbors + 0], 0);
  EXPECT_EQ(list[2 * maxNeighbors + 0], 0);
  for (int i = 0; i < n; ++i)
    full[i] = counts[i];
  // icase 2 and 3 agree here (no metals, nothing hypervalent).
  fillNeighborList(n, mol.numbers, rad, dist, mchar.data(), 2, 1.25, 1.0,
                   gffMetal.data(), gffGroup.data(), gffNormCn.data(),
                   full.data(), 1, list, counts);
  EXPECT_EQ(counts, std::vector<int>({ 2, 1, 1 }));
  fillNeighborList(n, mol.numbers, rad, dist, mchar.data(), 3, 1.25, 1.0,
                   gffMetal.data(), gffGroup.data(), gffNormCn.data(),
                   full.data(), 1, list, counts);
  EXPECT_EQ(counts, std::vector<int>({ 2, 1, 1 }));
}

TEST(GfnTopoTest, sif6Filters)
{
  // SiF6 octahedron: Si-F 3.0 Bohr. The full list bonds Si six times, but
  // the icase-3 list drops everything (Si hypervalent: 6 > normcn 4).
  const int n = 7;
  std::vector<int> numbers = { 14, 9, 9, 9, 9, 9, 9 };
  std::vector<double> xyz = { 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, -3.0, 0.0, 0.0,
                              0.0, 3.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 3.0,
                              0.0, 0.0, -3.0 };
  std::vector<double> dist, rad;
  buildInputs(n, numbers, xyz, dist, rad);
  std::vector<double> mchar(n, 0.0);
  std::vector<int> full(n, 0), list, counts;
  fillNeighborList(n, numbers, rad, dist, mchar.data(), 1, 1.25, 1.0,
                   gffMetal.data(), gffGroup.data(), gffNormCn.data(),
                   full.data(), 1, list, counts);
  EXPECT_EQ(counts[0], 6);
  for (int i = 0; i < n; ++i)
    full[i] = counts[i];
  fillNeighborList(n, numbers, rad, dist, mchar.data(), 3, 1.25, 1.0,
                   gffMetal.data(), gffGroup.data(), gffNormCn.data(),
                   full.data(), 1, list, counts);
  for (int c : counts)
    EXPECT_EQ(c, 0);
}

TEST(GfnTopoTest, naclMetalScaling)
{
  // NaCl at 4.5 Bohr bonds only through the icase-1 metal radius scaling.
  const int n = 2;
  std::vector<int> numbers = { 11, 17 };
  std::vector<double> xyz = { 0.0, 0.0, 0.0, 4.5, 0.0, 0.0 };
  std::vector<double> dist, rad;
  buildInputs(n, numbers, xyz, dist, rad);
  std::vector<double> mchar(n, 0.0);
  std::vector<int> full(n, 0), list, counts;
  fillNeighborList(n, numbers, rad, dist, mchar.data(), 1, 1.25, 1.0,
                   gffMetal.data(), gffGroup.data(), gffNormCn.data(),
                   full.data(), 1, list, counts);
  EXPECT_EQ(counts, std::vector<int>({ 1, 1 }));
  // ... while the icase-3 list excludes the metal entirely.
  for (int i = 0; i < n; ++i)
    full[i] = counts[i];
  fillNeighborList(n, numbers, rad, dist, mchar.data(), 3, 1.25, 1.0,
                   gffMetal.data(), gffGroup.data(), gffNormCn.data(),
                   full.data(), 1, list, counts);
  EXPECT_EQ(counts, std::vector<int>({ 0, 0 }));
}

TEST(GfnTopoTest, waterTopologyCharges)
{
  Water mol;
  const int n = 3;
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  // Plain chi/gam/alp2 inputs with a single fragment, real pair distances.
  std::vector<double> pair(n * (n + 1) / 2, 0.0), chi(n), gam(n), alp(n);
  for (int i = 0; i < n; ++i) {
    chi[i] = -param.chi[mol.numbers[i] - 1];
    gam[i] = param.gam[mol.numbers[i] - 1];
    double a = param.alp[mol.numbers[i] - 1];
    alp[i] = a * a;
    for (int j = 0; j <= i; ++j) {
      double dx = mol.xyz[3 * j] - mol.xyz[3 * i];
      double dy = mol.xyz[3 * j + 1] - mol.xyz[3 * i + 1];
      double dz = mol.xyz[3 * j + 2] - mol.xyz[3 * i + 2];
      pair[packedIndex(j, i)] = std::sqrt(dx * dx + dy * dy + dz * dz);
    }
  }
  Environment env;
  std::vector<int> frag(n, 1);
  std::vector<double> fragCharges = { 0.0 };
  std::vector<double> charges;
  double energy = 0.0;
  ASSERT_TRUE(topologyCharges(n, pair, chi, gam, alp, 1, frag, fragCharges,
                              env, charges, energy))
    << env.errorMessage();
  ASSERT_EQ(charges.size(), 3);
  EXPECT_NEAR(charges[0], -1.171027151812, 1e-9);
  EXPECT_NEAR(charges[1], 0.585513575906, 1e-9);
  EXPECT_NEAR(charges[2], 0.585513575906, 1e-9);
  EXPECT_NEAR(energy, -0.016683893842, 1e-9);
}

TEST(GfnTopoTest, bondLengthEstimator)
{
  // Cross-checked against an independent numpy implementation on five
  // geometries (see port notes); spot values below are the C++ outputs.
  Water mol;
  const int n = 3;
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  std::vector<double> dist, rad;
  buildInputs(n, mol.numbers, mol.xyz, dist, rad);
  std::vector<double> mchar(n, 0.0);
  std::vector<int> full(n, 0), list, counts;
  fillNeighborList(n, mol.numbers, rad, dist, mchar.data(), 2, 1.25, 1.0,
                   gffMetal.data(), gffGroup.data(), gffNormCn.data(),
                   full.data(), 1, list, counts);
  std::vector<float> rabd;
  std::vector<double> rtmp;
  estimateBondLengths(n, mol.numbers, counts, list, param.rad, gen.rfgoed1,
                      12.0, rabd, rtmp);
  ASSERT_EQ(rtmp.size(), 6);
  // O-H estimate: rfgoed1 * (0.64 + 0.32) / 0.52917726 Bohr.
  EXPECT_NEAR(rtmp[packedIndex(0, 1)], 2.1316107367, 1e-9);
  // H-H 1,4 estimate through O.
  EXPECT_NEAR(rtmp[packedIndex(1, 2)], 4.2632214733, 1e-9);
  EXPECT_NEAR(rabd[0 * n + 1], 0.96f, 1e-6);
}

} // namespace Xtb
} // namespace Avogadro
