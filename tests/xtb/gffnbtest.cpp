/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "avogadro/xtb/COPYING.LESSER".
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/xtb/elements.h>
#include <avogadro/xtb/environment.h>
#include <avogadro/xtb/gfnffdata.h>
#include <avogadro/xtb/gfnffdisp.h>
#include <avogadro/xtb/gfnffparams.h>
#include <avogadro/xtb/gfnffsetup.h>
#include <avogadro/xtb/gfnfftopo.h>

#include <vector>

namespace Avogadro {
namespace Xtb {

// Reference values verified bitwise against the unmodified nbondmat_pbc
// + countf (src/gfnff/gfnff_ini2.F90) and the ini alphanb loop
// (src/gfnff/gfnff_ini.f90), compiled with gfortran -ffp-contract=off
// (172 bpair/alphanb/zetac6 outputs over ethane and water-sodium inputs;
// the latter covers the unpaired-bond restore and the M...H rule).

namespace {

struct NbFixture
{
  std::vector<int> numbers;
  std::vector<std::vector<int>> neighbours;
  std::vector<double> qa;
};

NbFixture ethaneFixture()
{
  NbFixture f;
  f.numbers = { 6, 6, 1, 1, 1, 1, 1, 1 };
  f.neighbours = { { 1, 2, 3, 4 }, { 0, 5, 6, 7 }, { 0 }, { 0 },
                   { 0 },          { 1 },          { 1 }, { 1 } };
  f.qa = { -0.2, -0.05, 0.1, 0.1, 0.1, 0.05, 0.05, 0.05 };
  return f;
}

NbFixture waterSodiumFixture()
{
  NbFixture f;
  f.numbers = { 8, 1, 1, 11 };
  f.neighbours = { { 1, 2 }, { 0 }, { 0 }, { 1 } };
  f.qa = { -0.4, 0.2, 0.2, 0.3 };
  return f;
}

TEST(NbTest, PairFlagsMatchReference)
{
  Environment env;
  {
    NbFixture f = ethaneFixture();
    int n = static_cast<int>(f.numbers.size());
    std::vector<int> flags;
    ASSERT_TRUE(bondPairFlags(n, f.neighbours, flags, env));
    const int ref[64] = { 0, 1, 1, 1, 1, 2, 2, 2, 1, 0, 2, 2, 2, 1, 1, 1, 1, 2, 0, 2, 2, 3, 3, 3, 1, 2, 2, 0, 2, 3, 3, 3, 1, 2, 2, 2, 0, 3, 3, 3, 2, 1, 3, 3, 3, 0, 2, 2, 2, 1, 3, 3, 3, 2, 0, 2, 2, 1, 3, 3, 3, 2, 2, 0 };
    ASSERT_EQ((int)flags.size(), 64);
    for (int i = 0; i < 64; ++i)
      EXPECT_EQ(flags[i], ref[i]);
  }
  {
    NbFixture f = waterSodiumFixture();
    int n = static_cast<int>(f.numbers.size());
    std::vector<int> flags;
    ASSERT_TRUE(bondPairFlags(n, f.neighbours, flags, env));
    const int ref[16] = { 0, 1, 1, 5, 1, 0, 2, 1, 1, 2, 0, 5, 5, 1, 5, 0 };
    ASSERT_EQ((int)flags.size(), 16);
    for (int i = 0; i < 16; ++i)
      EXPECT_EQ(flags[i], ref[i]);
  }
}

TEST(NbTest, AlphanbMatchesReference)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  std::vector<int> metal(103), group(103);
  for (int i = 0; i < 103; ++i) {
    metal[i] = gffMetal[i];
    group[i] = gffGroup[i];
  }
  Environment env;
  {
    NbFixture f = ethaneFixture();
    int n = static_cast<int>(f.numbers.size());
    std::vector<int> flags;
    ASSERT_TRUE(bondPairFlags(n, f.neighbours, flags, env));
    std::vector<int> counts(n);
    for (int i = 0; i < n; ++i)
      counts[i] = (int)f.neighbours[i].size();
    std::vector<double> alphanb, zetac6;
    ASSERT_TRUE(buildNonbondedTables(
      n, f.numbers, f.qa, counts, param.repan, metal, gen.nRepScal,
      gen.qRepScal, gen.hhFac, gen.hh13Rep, gen.hh14Rep, flags, alphanb,
      zetac6, env));
    const double kAlphaRef[36] = { 5.13788758517586874e-01, 5.28005111102807900e-01, 5.42614825118143518e-01,
    6.64646330859556267e-01, 6.83036858926451496e-01, 6.53077116718755679e-01,
    6.64646330859556267e-01, 6.83036858926451496e-01, 9.52186400052219106e-01,
    6.53077116718755679e-01, 6.64646330859556267e-01, 6.83036858926451496e-01,
    9.52186400052219106e-01, 9.52186400052219106e-01, 6.53077116718755679e-01,
    6.59034679246797039e-01, 6.77269934905353366e-01, 4.58474709243083756e-01,
    4.58474709243083756e-01, 4.58474709243083756e-01, 6.42095727547011608e-01,
    6.59034679246797039e-01, 6.77269934905353366e-01, 4.58474709243083756e-01,
    4.58474709243083756e-01, 4.58474709243083756e-01, 9.36175535247231183e-01,
    6.42095727547011608e-01, 6.59034679246797039e-01, 6.77269934905353366e-01,
    4.58474709243083756e-01, 4.58474709243083756e-01, 4.58474709243083756e-01,
    9.36175535247231183e-01, 9.36175535247231183e-01, 6.42095727547011608e-01, };
    const double kZetaRef[36] = { 1.09053543641771400e+00, 1.05543506155971434e+00, 1.02146444028333061e+00,
    9.15428172955583097e-01, 8.85963864916346600e-01, 7.68437880930821660e-01,
    9.15428172955583097e-01, 8.85963864916346600e-01, 7.68437880930821660e-01,
    7.68437880930821660e-01, 9.15428172955583097e-01, 8.85963864916346600e-01,
    7.68437880930821660e-01, 7.68437880930821660e-01, 7.68437880930821660e-01,
    9.75364265446414747e-01, 9.43970832278642646e-01, 8.18750035685931077e-01,
    8.18750035685931077e-01, 8.18750035685931077e-01, 8.72356292643597153e-01,
    9.75364265446414747e-01, 9.43970832278642646e-01, 8.18750035685931077e-01,
    8.18750035685931077e-01, 8.18750035685931077e-01, 8.72356292643597153e-01,
    8.72356292643597153e-01, 9.75364265446414747e-01, 9.43970832278642646e-01,
    8.18750035685931077e-01, 8.18750035685931077e-01, 8.18750035685931077e-01,
    8.72356292643597153e-01, 8.72356292643597153e-01, 8.72356292643597153e-01, };
    ASSERT_EQ((int)alphanb.size(), 64);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j <= i; ++j)
        EXPECT_NEAR(alphanb[i * n + j], kAlphaRef[i * (i + 1) / 2 + j],
                    1.0e-12);
    }
    ASSERT_EQ((int)zetac6.size(), 36);
    for (int i = 0; i < 36; ++i)
      EXPECT_NEAR(zetac6[i], kZetaRef[i], 1.0e-12);
  }
  {
    NbFixture f = waterSodiumFixture();
    int n = static_cast<int>(f.numbers.size());
    std::vector<int> flags;
    ASSERT_TRUE(bondPairFlags(n, f.neighbours, flags, env));
    std::vector<int> counts(n);
    for (int i = 0; i < n; ++i)
      counts[i] = (int)f.neighbours[i].size();
    std::vector<double> alphanb, zetac6;
    ASSERT_TRUE(buildNonbondedTables(
      n, f.numbers, f.qa, counts, param.repan, metal, gen.nRepScal,
      gen.qRepScal, gen.hhFac, gen.hh13Rep, gen.hh14Rep, flags, alphanb,
      zetac6, env));
    const double kAlphaRef[10] = { 7.10781546110741402e-01, 9.08324001259021019e-01, 6.75039895062243489e-01,
    9.08324001259021019e-01, 9.84208129662194731e-01, 6.75039895062243489e-01,
    4.86913751476177603e-01, 5.08560119962368051e-01, 5.08560119962368051e-01,
    3.33555369682693503e-01, };
    const double kZetaRef[10] = { 1.20022533080658800e+00, 8.56760380303890234e-01, 6.11583783825967564e-01,
    8.56760380303890234e-01, 6.11583783825967564e-01, 6.11583783825967564e-01,
    1.07994789058615948e+00, 7.70902381159695671e-01, 7.70902381159695671e-01,
    9.71723739239626738e-01, };
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j <= i; ++j)
        EXPECT_NEAR(alphanb[i * n + j], kAlphaRef[i * (i + 1) / 2 + j],
                    1.0e-12);
    }
    ASSERT_EQ((int)zetac6.size(), 10);
    for (int i = 0; i < 10; ++i)
      EXPECT_NEAR(zetac6[i], kZetaRef[i], 1.0e-12);
  }
}

} // namespace
} // namespace Xtb
} // namespace Avogadro
