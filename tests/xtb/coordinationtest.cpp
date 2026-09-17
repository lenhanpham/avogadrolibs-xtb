/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "avogadro/xtb/COPYING.LESSER".
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/xtb/coordination.h>

#include <cmath>
#include <vector>

namespace Avogadro {
namespace Xtb {

// Parity harness: reference values below were produced by compiling the
// unmodified xtb sources (src/disp/ncoord.f90 + mctc accuracy/constants and
// the pbc helpers) with gfortran and calling ncoord_d3/erf/d4/gfn plus the
// dncoord_* variants on the same geometries (coordinates in Bohr).
// The harness program is kept next to this test for re-generation.

struct Water
{
  std::vector<int> numbers = { 8, 1, 1 };
  // O-H 0.9584 Angstrom, H-O-H 104.45 deg, in Bohr.
  std::vector<double> xyz = { 0.0, 0.0, 0.0, 1.4305226400, 1.1071846211, 0.0,
                              -1.4305226400, 1.1071846211, 0.0 };
};

// D3 reference from ncoord_d3 on Water.
TEST(CoordinationTest, d3Water)
{
  Water mol;
  std::vector<double> cn;
  coordinationD3(mol.numbers, mol.xyz, cn);
  ASSERT_EQ(cn.size(), 3);
  EXPECT_NEAR(cn[0], 1.988715307819, 1e-9);
  EXPECT_NEAR(cn[1], 0.995285286192, 1e-9);
  EXPECT_NEAR(cn[2], 0.995285286192, 1e-9);
}

// GFN reference from ncoord_gfn on Water.
TEST(CoordinationTest, gfnWater)
{
  Water mol;
  std::vector<double> cn;
  coordinationGfn(mol.numbers, mol.xyz, cn);
  ASSERT_EQ(cn.size(), 3);
  EXPECT_NEAR(cn[0], 1.924069048663, 1e-9);
  EXPECT_NEAR(cn[1], 0.974540067672, 1e-9);
  EXPECT_NEAR(cn[2], 0.974540067672, 1e-9);
}

// D4 reference from ncoord_d4 on Water.
TEST(CoordinationTest, d4Water)
{
  Water mol;
  std::vector<double> cn;
  coordinationD4(mol.numbers, mol.xyz, cn);
  ASSERT_EQ(cn.size(), 3);
  EXPECT_NEAR(cn[0], 1.611256544646, 1e-9);
  EXPECT_NEAR(cn[1], 0.805628272323, 1e-9);
  EXPECT_NEAR(cn[2], 0.805628272323, 1e-9);
}

// Erf reference from ncoord_erf on Water.
TEST(CoordinationTest, erfWater)
{
  Water mol;
  std::vector<double> cn;
  coordinationErf(mol.numbers, mol.xyz, cn);
  ASSERT_EQ(cn.size(), 3);
  EXPECT_NEAR(cn[0], 1.990429022622, 1e-9);
  EXPECT_NEAR(cn[1], 0.995214511311, 1e-9);
  EXPECT_NEAR(cn[2], 0.995214511311, 1e-9);
}

// Derivative parity vs the dncoord_* reference routines. d(cn[0])/d(xyz[j]).
// Note the upstream sign quirk documented in coordination.h: the erf
// diagonal and d4 off-diagonal entries reproduce xtb exactly but differ in
// sign from finite differences, so they are pinned to reference values here
// instead of checked against finite differences below.
TEST(CoordinationTest, derivativeParityWater)
{
  Water mol;
  std::vector<double> cn, dcn;
  coordinationErfDerivative(mol.numbers, mol.xyz, cn, dcn);
  const double erfg[9] = { 0.0, -0.075427728650, 0.0, -0.048727678953,
                           -0.037713864325, 0.0, 0.048727678953,
                           -0.037713864325, 0.0 };
  for (int k = 0; k < 9; ++k)
    EXPECT_NEAR(dcn[k], erfg[k], 1e-9);
  coordinationD4Derivative(mol.numbers, mol.xyz, cn, dcn);
  const double d4g[9] = { 0.0, 0.061058907428, 0.0, 0.039445160177,
                          0.030529453714, 0.0, -0.039445160177,
                          0.030529453714, 0.0 };
  for (int k = 0; k < 9; ++k)
    EXPECT_NEAR(dcn[k], d4g[k], 1e-9);
}

// Analytic derivatives must match central finite differences of the values.
// (Checks mathematical correctness independently of the Fortran parity.)
TEST(CoordinationTest, derivativesMatchFiniteDifferences)
{
  Water mol;
  const double h = 1e-6;
  const double tol = 1e-6;
  std::vector<double> cn, dcn, cnPlus, cnMinus, xyzPlus, xyzMinus;
  auto check =
    [&](void (*values)(const std::vector<int>&, const std::vector<double>&,
                       std::vector<double>&, double),
        void (*derivs)(const std::vector<int>&, const std::vector<double>&,
                       std::vector<double>&, std::vector<double>&, double)) {
      derivs(mol.numbers, mol.xyz, cn, dcn, coordinationThreshold);
      for (int j = 0; j < 3; ++j) {
        for (int c = 0; c < 3; ++c) {
          xyzPlus = mol.xyz;
          xyzMinus = mol.xyz;
          xyzPlus[3 * j + c] += h;
          xyzMinus[3 * j + c] -= h;
          values(mol.numbers, xyzPlus, cnPlus, coordinationThreshold);
          values(mol.numbers, xyzMinus, cnMinus, coordinationThreshold);
          for (int i = 0; i < 3; ++i) {
            double fd = (cnPlus[i] - cnMinus[i]) / (2.0 * h);
            EXPECT_NEAR(dcn[(i * 3 + j) * 3 + c], fd, tol);
          }
        }
      }
    };
  check(coordinationD3, coordinationD3Derivative);
  check(coordinationGfn, coordinationGfnDerivative);
}

TEST(CoordinationTest, radiiAndElectronegativities)
{
  // rcov(H) = 4/3 * 0.32 / 0.52917726 Bohr, with xtb's single-precision
  // table rounding (see coordination.cpp).
  EXPECT_NEAR(covalentRadius(1), 0.80628305368, 1e-9);
  // EN table entries are single-rounded like rcov (1e-6 tolerance covers
  // the float(2.20) = 2.2000000477 representation xtb itself carries).
  EXPECT_NEAR(paulingElectronegativity(1), 2.20, 1e-6);
  EXPECT_NEAR(paulingElectronegativity(8), 3.44, 1e-6);
  EXPECT_EQ(covalentRadius(0), 0.0);
  EXPECT_EQ(covalentRadius(119), 0.0);
}

TEST(CoordinationTest, logCutoff)
{
  EXPECT_NEAR(logCoordinationCutoff(0.0), 0.0, 1e-12);
  EXPECT_NEAR(dLogCoordinationCutoff(0.0, 4.5), 0.9890130574, 1e-9);
  // Fixed point behaviour: cutoff(x) -> x for small x is smooth.
  std::vector<double> cn = { 1.0 };
  applyLogCoordinationCutoff(cn);
  EXPECT_NEAR(cn[0], logCoordinationCutoff(1.0), 1e-12);
}

} // namespace Xtb
} // namespace Avogadro
