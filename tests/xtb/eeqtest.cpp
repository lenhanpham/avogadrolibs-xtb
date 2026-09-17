/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "avogadro/xtb/COPYING.LESSER".
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/xtb/coordination.h>
#include <avogadro/xtb/eeq.h>
#include <avogadro/xtb/environment.h>
#include <avogadro/xtb/gfn0params.h>

#include <vector>

namespace Avogadro {
namespace Xtb {

// Reference values produced by compiling the unmodified xtb sources
// (src/eeq_model.f90 goedecker_chrgeq with naive BLAS/LAPACK-equivalent
// shims) with gfortran. The C++ port reproduces them to ~1e-16; the 1e-9
// tolerances below leave headroom for production linear solvers
// (LAPACK dsysv / Eigen), which agree up to solver rounding.

struct Water
{
  std::vector<int> numbers = { 8, 1, 1 };
  std::vector<double> xyz = { 0.0, 0.0, 0.0, 1.4305226400, 1.1071846211,
                              0.0, -1.4305226400, 1.1071846211, 0.0 };
};

struct Ammonium
{
  std::vector<int> numbers = { 7, 1, 1, 1, 1 };
  std::vector<double> xyz = { 0.0, 0.0, 0.0, 1.09696552, 1.09696552,
                              1.09696552, 1.09696552, -1.09696552,
                              -1.09696552, -1.09696552, 1.09696552,
                              -1.09696552, -1.09696552, -1.09696552,
                              1.09696552 };
};

static void solve(const std::vector<int>& numbers,
                  const std::vector<double>& xyz, double totalCharge,
                  std::vector<double>& charges, double& energy,
                  std::vector<double>& gradient,
                  std::vector<double>& derivatives)
{
  Environment env;
  ChargeParams params;
  ASSERT_TRUE(buildChargeParams(numbers, params, env));
  std::vector<double> cn, dcn;
  coordinationErf(numbers, xyz, cn);
  coordinationErfDerivative(numbers, xyz, cn, dcn);
  energy = 0.0;
  ASSERT_TRUE(eeqCharges(numbers, xyz, totalCharge, params, cn, &dcn, env,
                         charges, energy, &gradient, &derivatives))
    << env.errorMessage();
}

TEST(EeqTest, waterChargesAndEnergy)
{
  Water mol;
  std::vector<double> charges, gradient, derivatives;
  double energy = 0.0;
  solve(mol.numbers, mol.xyz, 0.0, charges, energy, gradient, derivatives);
  ASSERT_EQ(charges.size(), 3);
  EXPECT_NEAR(charges[0], -0.598247097451, 1e-9);
  EXPECT_NEAR(charges[1], 0.299123548726, 1e-9);
  EXPECT_NEAR(charges[2], 0.299123548726, 1e-9);
  EXPECT_NEAR(charges[0] + charges[1] + charges[2], 0.0, 1e-12);
  EXPECT_NEAR(energy, -0.064561627370, 1e-9);
}

TEST(EeqTest, waterGradient)
{
  Water mol;
  std::vector<double> charges, gradient, derivatives;
  double energy = 0.0;
  solve(mol.numbers, mol.xyz, 0.0, charges, energy, gradient, derivatives);
  ASSERT_EQ(gradient.size(), 9);
  const double ref[9] = { -0.0, -0.048693893400, 0.0, 0.017392751072,
                          0.021921582392, 0.0, -0.017392751073, 0.021921582392,
                          0.0 };
  for (int i = 0; i < 9; ++i)
    EXPECT_NEAR(gradient[i], ref[i], 1e-9);
}

TEST(EeqTest, waterChargeDerivatives)
{
  Water mol;
  std::vector<double> charges, gradient, derivatives;
  double energy = 0.0;
  solve(mol.numbers, mol.xyz, 0.0, charges, energy, gradient, derivatives);
  ASSERT_EQ(derivatives.size(), 27);
  const double ref[27] = {
    -0.0, -0.439974919604, 0.0, 0.135407517445, 0.219987459802, 0.0,
    -0.135407517445, 0.219987459802, 0.0, 0.168425716972, 0.208750388984,
    0.0, -0.150308991475, -0.155531757736, 0.0, -0.018116725498,
    -0.053218631210, 0.0, -0.168425716972, 0.208750388984, 0.0,
    0.018116725498, -0.053218631210, 0.0, 0.150308991475, -0.155531757736,
    0.0
  };
  for (int i = 0; i < 27; ++i)
    EXPECT_NEAR(derivatives[i], ref[i], 1e-9);
}

TEST(EeqTest, ammoniumCation)
{
  Ammonium mol;
  std::vector<double> charges, gradient, derivatives;
  double energy = 0.0;
  solve(mol.numbers, mol.xyz, 1.0, charges, energy, gradient, derivatives);
  ASSERT_EQ(charges.size(), 5);
  EXPECT_NEAR(charges[0], -0.890242825922, 1e-9);
  for (int i = 1; i < 5; ++i)
    EXPECT_NEAR(charges[i], 0.472560706481, 1e-9);
  double sum = 0.0;
  for (double q : charges)
    sum += q;
  EXPECT_NEAR(sum, 1.0, 1e-12);
  EXPECT_NEAR(energy, 1.372981323735, 1e-9);
  // Tetrahedral symmetry: nitrogen gradient vanishes, hydrogen gradients
  // are equal in magnitude along their N-H directions.
  EXPECT_NEAR(gradient[0], 0.0, 1e-9);
  EXPECT_NEAR(gradient[1], 0.0, 1e-9);
  EXPECT_NEAR(gradient[2], 0.0, 1e-9);
  EXPECT_NEAR(gradient[3], 0.012499673415, 1e-9);
  EXPECT_NEAR(gradient[4], 0.012499673415, 1e-9);
  EXPECT_NEAR(gradient[5], 0.012499673415, 1e-9);
}

TEST(EeqTest, ammoniumChargeDerivatives)
{
  Ammonium mol;
  std::vector<double> charges, gradient, derivatives;
  double energy = 0.0;
  solve(mol.numbers, mol.xyz, 1.0, charges, energy, gradient, derivatives);
  ASSERT_EQ(derivatives.size(), 75);
  // Spot checks against the reference (full tensor verified bitwise during
  // porting); d(q0)/d(x0) vanishes by symmetry.
  EXPECT_NEAR(derivatives[0], 0.0, 1e-9);
  EXPECT_NEAR(derivatives[1], 0.0, 1e-9);
  EXPECT_NEAR(derivatives[2], 0.0, 1e-9);
  EXPECT_NEAR(derivatives[(0 * 5 + 1) * 3 + 0], 0.125854140533, 1e-9);
  EXPECT_NEAR(derivatives[(1 * 5 + 0) * 3 + 0], 0.220798061645, 1e-9);
  EXPECT_NEAR(derivatives[(1 * 5 + 1) * 3 + 0], -0.103954542035, 1e-9);
}

TEST(EeqTest, chargeParams)
{
  Environment env;
  ChargeParams params;
  EXPECT_TRUE(buildChargeParams({ 8, 1, 1 }, params, env));
  EXPECT_NEAR(params.en[0], 1.56866440, 1e-12);
  EXPECT_NEAR(params.gam[0], 0.03151644, 1e-12);
  EXPECT_NEAR(params.alpha[1], 0.55159092, 1e-12);
  EXPECT_FALSE(buildChargeParams({ 92 }, params, env));
  EXPECT_TRUE(env.failed());
}

TEST(EeqTest, coulombMatrix)
{
  Water mol;
  Environment env;
  ChargeParams params;
  ASSERT_TRUE(buildChargeParams(mol.numbers, params, env));
  std::vector<double> amat;
  coulombMatrix(mol.numbers, mol.xyz, params, amat);
  ASSERT_EQ(amat.size(), 9);
  // Symmetric with screened off-diagonals below 1/r.
  EXPECT_NEAR(amat[0 * 3 + 1], amat[1 * 3 + 0], 1e-15);
  double r = 1.8089369278; // O-H distance in Bohr
  EXPECT_LT(amat[0 * 3 + 1], 1.0 / r);
  EXPECT_GT(amat[0 * 3 + 1], 0.0);
  EXPECT_GT(amat[0], 0.0);
}

} // namespace Xtb
} // namespace Avogadro
