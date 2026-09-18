/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "avogadro/xtb/COPYING.LESSER".
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/xtb/eeq.h>
#include <avogadro/xtb/environment.h>
#include <avogadro/xtb/gfnffdata.h>
#include <avogadro/xtb/gfnffparams.h>
#include <avogadro/xtb/gfnffsetup.h>

#include <vector>

namespace Avogadro {
namespace Xtb {

// Reference values verified bitwise: the EEQ xi/initial/gamma/final
// kernels against verbatim copies of the gfnff_ini.f90 blocks (dxi loop,
// initial/final parameter blocks, dgam loop) and the ini2 amide/amideH
// predicates, compiled with gfortran -ffp-contract=off (80 sections over
// ethane, water, formamide, nitromethane, BH3, N-methylacetamide,
// perchlorate and CH3Cl/Na+ inputs, 0 mismatches).

namespace {

struct EeqCase
{
  std::vector<int> numbers;
  std::vector<int> hyb, itag, imetal, pi;
  std::vector<std::vector<int>> neighbours;
  std::vector<int> counts;
  std::vector<double> qa;
};

EeqCase makeWater()
{
  EeqCase c;
  c.numbers = { 8, 1, 1 };
  c.hyb = { 3, 0, 0 };
  c.itag = { 0, 0, 0 };
  c.imetal = { 0, 0, 0 };
  c.pi = { 0, 0, 0 };
  c.neighbours = { { 1, 2 }, { 0 }, { 0 } };
  c.counts = { 2, 1, 1 };
  c.qa = { -0.4, 0.2, 0.2 };
  return c;
}

EeqCase makeFormamide()
{
  EeqCase c;
  c.numbers = { 6, 8, 7, 1, 1, 1 };
  c.hyb = { 2, 2, 3, 0, 1, 1 };
  c.itag = { 0, 0, 0, 0, 0, 0 };
  c.imetal = { 0, 0, 0, 0, 0, 0 };
  c.pi = { 1, 1, 1, 0, 0, 0 };
  c.neighbours = { { 1, 2, 3 }, { 0 }, { 0, 4, 5 }, { 0 }, { 2, 5 }, { 2, 4 } };
  c.counts = { 3, 1, 3, 1, 2, 2 };
  c.qa = { 0.1, -0.3, -0.2, 0.1, 0.1, 0.1 };
  return c;
}

EeqCase makeNma()
{
  EeqCase c;
  c.numbers = { 6, 6, 8, 7, 1, 6, 1, 1, 1, 1, 1, 1 };
  c.hyb = { 3, 2, 2, 3, 0, 3, 0, 0, 0, 0, 0, 0 };
  c.itag = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  c.imetal = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  c.pi = { 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0 };
  c.neighbours = { { 1, 6, 7, 8 }, { 0, 2, 3 }, { 1 }, { 1, 4, 5 }, { 3 },
                   { 3, 9, 10, 11 }, { 0 }, { 0 }, { 0 }, { 5 }, { 5 }, { 5 } };
  c.counts = { 4, 3, 1, 3, 1, 4, 1, 1, 1, 1, 1, 1 };
  c.qa = { 0.0, 0.2, -0.3, -0.2, 0.15, 0.0,
           0.05, 0.05, 0.05, 0.05, 0.05, 0.05 };
  return c;
}

void runPipeline(const EeqCase& c, const GffData& param, const GffGenerator& gen,
                 std::vector<double>& dxi, std::vector<double>& chi0,
                 std::vector<double>& gam0, std::vector<double>& alp0,
                 std::vector<double>& dgam, std::vector<double>& chi,
                 std::vector<double>& gam, std::vector<double>& alp)
{
  Environment env;
  int n = static_cast<int>(c.numbers.size());
  std::vector<int> group(103), metal(103);
  for (int i = 0; i < 103; ++i) {
    group[i] = gffGroup[i];
    metal[i] = gffMetal[i];
  }
  ASSERT_TRUE(eeqXiCorrections(n, c.numbers, c.itag, c.imetal, c.pi,
                               c.neighbours, c.counts, group.data(), dxi,
                               env));
  ASSERT_TRUE(eeqInitialParams(n, c.numbers, param.chi, param.gam, param.alp,
                               param.cnf, c.imetal, c.counts, gen.cnMax,
                               gen.mchiShift, dxi, chi0, gam0, alp0, env));
  ASSERT_TRUE(eeqGammaCorrections(n, c.numbers, c.hyb, c.imetal, c.pi,
                                  c.neighbours, group.data(), c.qa, dgam,
                                  env));
  ASSERT_TRUE(eeqFinalParams(n, c.numbers, c.hyb, c.imetal, c.pi, c.neighbours,
                             param.chi, param.gam, param.alp, dxi, dgam, c.qa,
                             group.data(), chi, gam, alp, env));
}

} // namespace

TEST(EeqSetupTest, WaterH2OCorrection)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  EeqCase c = makeWater();
  std::vector<double> dxi, chi0, gam0, alp0, dgam, chi, gam, alp;
  runPipeline(c, param, gen, dxi, chi0, gam0, alp0, dgam, chi, gam, alp);
  // H2O rule (-0.02 float) plus the group-6 oxygen rule (-nh*0.005 float).
  EXPECT_NEAR(dxi[0], -2.99999993294477463e-02, 0.0);
  EXPECT_NEAR(dxi[1], 0.0, 0.0);
  EXPECT_NEAR(dxi[2], 0.0, 0.0);
  EXPECT_NEAR(chi0[0], -1.49867025264935405e+00, 0.0);
  EXPECT_NEAR(chi0[1], -1.21815000000000007e+00, 0.0);
  EXPECT_NEAR(gam0[0], -3.12359999999999999e-02, 0.0);
  EXPECT_NEAR(alp0[0], 8.19653190408999976e-01, 0.0);
  EXPECT_NEAR(dgam[0], 6.00000023841857910e-02, 0.0);
  EXPECT_NEAR(dgam[1], -1.59999996423721320e-02, 0.0);
  EXPECT_NEAR(chi[0], -1.72120099932944770e+00, 0.0);
  EXPECT_NEAR(gam[0], 2.87640023841857911e-02, 0.0);
  EXPECT_NEAR(alp[0], 8.41525517916896781e-01, 0.0);
  EXPECT_NEAR(alp[1], 3.42305734760999958e-01, 0.0);
  EXPECT_FALSE(isAmide(3, c.numbers, c.hyb, c.neighbours, c.pi, 0));
  EXPECT_FALSE(isAmideHydrogen(3, c.numbers, c.hyb, c.neighbours, c.pi, 1));
}

TEST(EeqSetupTest, FormamideAmideBranch)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  EeqCase c = makeFormamide();
  EXPECT_TRUE(isAmide(6, c.numbers, c.hyb, c.neighbours, c.pi, 2));
  EXPECT_FALSE(isAmide(6, c.numbers, c.hyb, c.neighbours, c.pi, 0));
  std::vector<double> dxi, chi0, gam0, alp0, dgam, chi, gam, alp;
  runPipeline(c, param, gen, dxi, chi0, gam0, alp0, dgam, chi, gam, alp);
  // Amide N takes ff=-0.16 (not the pi-N -0.14).
  EXPECT_NEAR(dgam[2], 3.19999992847442641e-02, 0.0);
  EXPECT_NEAR(chi[2], -1.52848500000000009e+00, 0.0);
  EXPECT_NEAR(gam[2], 4.71999928474426539e-03, 0.0);
  EXPECT_NEAR(alp[2], 1.74342446708114496e+00, 0.0);
  EXPECT_NEAR(alp[0], 8.32528505552620590e-01, 0.0);
}

TEST(EeqSetupTest, MethylAcetamideAmideH)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  EeqCase c = makeNma();
  EXPECT_TRUE(isAmide(12, c.numbers, c.hyb, c.neighbours, c.pi, 3));
  EXPECT_TRUE(isAmideHydrogen(12, c.numbers, c.hyb, c.neighbours, c.pi, 4));
  EXPECT_FALSE(isAmideHydrogen(12, c.numbers, c.hyb, c.neighbours, c.pi, 6));
  std::vector<double> dxi, chi0, gam0, alp0, dgam, chi, gam, alp;
  runPipeline(c, param, gen, dxi, chi0, gam0, alp0, dgam, chi, gam, alp);
  // Amide-H chi carries the extra -0.02 reduction.
  EXPECT_NEAR(chi[4], -1.24705399955296525e+00, 0.0);
  EXPECT_NEAR(dgam[3], 3.19999992847442641e-02, 0.0);
  EXPECT_NEAR(dgam[4], -1.19999997317790982e-02, 0.0);
}

TEST(EeqSetupTest, PerchlorateGroup7Branch)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  EeqCase c;
  c.numbers = { 17, 8, 8, 8, 8 };
  c.hyb = { 5, 2, 2, 2, 2 };
  c.itag = { 0, 0, 0, 0, 0 };
  c.imetal = { 0, 0, 0, 0, 0 };
  c.pi = { 1, 1, 1, 1, 1 };
  c.neighbours = { { 1, 2, 3, 4 }, { 0 }, { 0 }, { 0 }, { 0 } };
  c.counts = { 4, 1, 1, 1, 1 };
  c.qa = { 0.6, -0.4, -0.4, -0.4, -0.4 };
  std::vector<double> dxi, chi0, gam0, alp0, dgam, chi, gam, alp;
  runPipeline(c, param, gen, dxi, chi0, gam0, alp0, dgam, chi, gam, alp);
  // Polyvalent Cl without metal neighbours: -nn*0.021 (float).
  EXPECT_NEAR(dxi[0], -8.39999988675117493e-02, 0.0);
  EXPECT_NEAR(dxi[1], 0.0, 0.0);
}

TEST(EeqSetupTest, ChlorideSodiumMetalBranch)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  EeqCase c;
  c.numbers = { 6, 17, 11, 1, 1, 1 };
  c.hyb = { 3, 1, 0, 1, 1, 0 };
  c.itag = { 0, 0, 0, 0, 0, 0 };
  c.imetal = { 0, 0, 1, 0, 0, 0 };
  c.pi = { 0, 1, 0, 0, 0, 0 };
  c.neighbours = { { 1, 3, 4, 5 }, { 0, 2 }, { 1 }, { 0, 4 }, { 0, 3 }, { 0 } };
  c.counts = { 4, 2, 1, 2, 2, 1 };
  c.qa = { 0.0, -0.3, 0.3, 0.0, 0.0, 0.0 };
  std::vector<double> dxi, chi0, gam0, alp0, dgam, chi, gam, alp;
  runPipeline(c, param, gen, dxi, chi0, gam0, alp0, dgam, chi, gam, alp);
  // Cl with a metal neighbour: +nn*0.05 (float).
  EXPECT_NEAR(dxi[1], 1.00000001490116119e-01, 0.0);
}

TEST(EeqSetupTest, NitromethaneNitroOxygen)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  EeqCase c;
  c.numbers = { 6, 7, 8, 8, 1, 1, 1 };
  c.hyb = { 3, 2, 2, 2, 1, 1, 0 };
  c.itag = { 0, 1, 0, 0, 0, 0, 0 };
  c.imetal = { 0, 0, 0, 0, 0, 0, 0 };
  c.pi = { 0, 1, 1, 1, 0, 0, 0 };
  c.neighbours = { { 1, 4, 5, 6 }, { 0, 2, 3 }, { 1 },
                   { 1 }, { 0, 5 }, { 0, 4 }, { 0 } };
  c.counts = { 4, 3, 1, 1, 2, 2, 1 };
  c.qa = { 0.0, 0.3, -0.2, -0.2, 0.05, 0.05, 0.0 };
  std::vector<double> dxi, chi0, gam0, alp0, dgam, chi, gam, alp;
  runPipeline(c, param, gen, dxi, chi0, gam0, alp0, dgam, chi, gam, alp);
  // Nitro oxygens take dxi=+0.05 (single precision).
  EXPECT_NEAR(dxi[2], 5.00000007450580597e-02, 0.0);
  EXPECT_NEAR(dxi[3], 5.00000007450580597e-02, 0.0);
  EXPECT_NEAR(chi[2], -1.64120099925494189e+00, 0.0);
}

TEST(EeqSetupTest, BoraneBoronRule)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  EeqCase c;
  c.numbers = { 5, 1, 1, 1 };
  c.hyb = { 2, 0, 0, 0 };
  c.itag = { 0, 0, 0, 0 };
  c.imetal = { 0, 0, 0, 0 };
  c.pi = { 1, 0, 0, 0 };
  c.neighbours = { { 1, 2, 3 }, { 0 }, { 0 }, { 0 } };
  c.counts = { 3, 1, 1, 1 };
  c.qa = { 0.0, 0.0, 0.0, 0.0 };
  std::vector<double> dxi, chi0, gam0, alp0, dgam, chi, gam, alp;
  runPipeline(c, param, gen, dxi, chi0, gam0, alp0, dgam, chi, gam, alp);
  // Boron: +nh*0.015 (float, nh=3).
  EXPECT_NEAR(dxi[0], 4.49999980628490448e-02, 0.0);
  EXPECT_NEAR(chi[0], -1.14149900193715093e+00, 0.0);
  EXPECT_NEAR(alp[0], 1.36264298632900016e+00, 0.0);
}

} // namespace Xtb
} // namespace Avogadro
