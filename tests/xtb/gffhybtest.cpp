/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "avogadro/xtb/COPYING.LESSER".
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/xtb/environment.h>
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

// Reference values produced by compiling the unmodified xtb sources
// (gfnff_rab.f gfnffrab, ini2.F90 gfnff_neigh) with gfortran. Integers match
// exactly; the two rab spot values are bitwise identical.

// Full qloop-0 pipeline: bond guesses, charge scaling, neighbor lists,
// hybridization. qa is zero (first production pass).
static void runPipeline(int n, const std::vector<int>& numbers,
                        const std::vector<double>& xyz,
                        const std::vector<double>& qa, const GffData& param,
                        const GffGenerator& gen, std::vector<int>& hyb,
                        std::vector<int>& itag)
{
  std::vector<double> cn(n);
  for (int i = 0; i < n; ++i)
    cn[i] = param.normCn[numbers[i] - 1];
  std::vector<double> rtmp;
  gfnffBondGuesses(n, numbers, cn, rtmp);
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
  ASSERT_TRUE(assignHybridization(n, numbers, xyz, nbf, nbfc, nb, nbc, nbm,
                                  nbmc, qa, metal.data(), group.data(), 160.0,
                                  1, true, hyb, itag, env))
    << env.errorMessage();
}

TEST(GffHybTest, water)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  std::vector<int> hyb, itag;
  runPipeline(3, { 8, 1, 1 },
              { 0.0, 0.0, 0.0, 1.4305226400, 1.1071846211, 0.0, -1.4305226400,
                1.1071846211, 0.0 },
              { 0.0, 0.0, 0.0 }, param, gen, hyb, itag);
  EXPECT_EQ(hyb, std::vector<int>({ 3, 0, 0 }));
  EXPECT_EQ(itag, std::vector<int>({ 0, 0, 0 }));
}

TEST(GffHybTest, hydrocarbons)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  std::vector<int> hyb, itag;
  // Ethane: sp3.
  runPipeline(8, { 6, 6, 1, 1, 1, 1, 1, 1 },
              { -1.45, 0.0, 0.0, 1.45, 0.0, 0.0, -2.136, 1.942, 0.0, -2.136,
                -0.971, 1.681, -2.136, -0.971, -1.681, 2.136, 1.942, 0.0,
                2.136, -0.971, 1.681, 2.136, -0.971, -1.681 },
              std::vector<double>(8, 0.0), param, gen, hyb, itag);
  EXPECT_EQ(hyb, std::vector<int>({ 3, 3, 0, 0, 0, 0, 0, 0 }));
  // Ethene: sp2.
  runPipeline(6, { 6, 6, 1, 1, 1, 1 },
              { 0.0, 0.0, 0.0, 2.53, 0.0, 0.0, -1.025, 1.775, 0.0, -1.025,
                -1.775, 0.0, 3.555, 1.775, 0.0, 3.555, -1.775, 0.0 },
              std::vector<double>(6, 0.0), param, gen, hyb, itag);
  EXPECT_EQ(hyb, std::vector<int>({ 2, 2, 0, 0, 0, 0 }));
  // Ethyne: linear sp via the geometry-dependent branch.
  runPipeline(4, { 6, 6, 1, 1 },
              { 0.0, 0.0, 0.0, 2.27, 0.0, 0.0, -1.99, 0.0, 0.0, 4.26, 0.0,
                0.0 },
              std::vector<double>(4, 0.0), param, gen, hyb, itag);
  EXPECT_EQ(hyb, std::vector<int>({ 1, 1, 0, 0 }));
  EXPECT_EQ(itag, std::vector<int>({ 0, 0, 0, 0 }));
  // Benzene: aromatic sp2.
  runPipeline(12, { 6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1 },
              { 2.627, 0.0, 0.0, 1.3135, 2.275, 0.0, -1.3135, 2.275, 0.0,
                -2.627, 0.0, 0.0, -1.3135, -2.275, 0.0, 1.3135, -2.275, 0.0,
                4.687, 0.0, 0.0, 2.3435, 4.059, 0.0, -2.3435, 4.059, 0.0,
                -4.687, 0.0, 0.0, -2.3435, -4.059, 0.0, 2.3435, -4.059, 0.0 },
              std::vector<double>(12, 0.0), param, gen, hyb, itag);
  EXPECT_EQ(hyb, std::vector<int>({ 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0 }));
}

TEST(GffHybTest, heteroDiatomics)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  std::vector<int> hyb, itag;
  runPipeline(3, { 1, 6, 7 }, { 0.0, 0.0, 0.0, 2.01, 0.0, 0.0, 4.189, 0.0,
                               0.0 },
              std::vector<double>(3, 0.0), param, gen, hyb, itag);
  EXPECT_EQ(hyb, std::vector<int>({ 0, 1, 1 }));
  runPipeline(2, { 6, 8 }, { 0.0, 0.0, 0.0, 2.14, 0.0, 0.0 },
              std::vector<double>(2, 0.0), param, gen, hyb, itag);
  EXPECT_EQ(hyb, std::vector<int>({ 1, 1 }));
  runPipeline(2, { 7, 7 }, { 0.0, 0.0, 0.0, 2.07, 0.0, 0.0 },
              std::vector<double>(2, 0.0), param, gen, hyb, itag);
  EXPECT_EQ(hyb, std::vector<int>({ 1, 1 }));
  runPipeline(2, { 8, 8 }, { 0.0, 0.0, 0.0, 2.29, 0.0, 0.0 },
              std::vector<double>(2, 0.0), param, gen, hyb, itag);
  EXPECT_EQ(hyb, std::vector<int>({ 1, 1 }));
}

TEST(GffHybTest, nitrogenCases)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  std::vector<int> hyb, itag;
  // Ammonia: pyramidal sp3.
  runPipeline(4, { 7, 1, 1, 1 },
              { 0.0, 0.0, 0.0, 1.806, 0.0, -0.633, -0.903, 1.564, -0.633,
                -0.903, -1.564, -0.633 },
              std::vector<double>(4, 0.0), param, gen, hyb, itag);
  EXPECT_EQ(hyb, std::vector<int>({ 3, 0, 0, 0 }));
  // NO2: bent sp2 with the Hueckel tag on N.
  runPipeline(3, { 7, 8, 8 },
              { 0.0, 0.0, 0.0, 0.879, 2.071, 0.0, 0.879, -2.071, 0.0 },
              std::vector<double>(3, 0.0), param, gen, hyb, itag);
  EXPECT_EQ(hyb, std::vector<int>({ 2, 2, 2 }));
  // Pyridine: aromatic N stays sp2.
  runPipeline(11, { 7, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1 },
              { 0.0, 2.6, 0.0, -2.25, 1.3, 0.0, -2.25, -1.3, 0.0, 0.0, -2.6,
                0.0, 2.25, -1.3, 0.0, 2.25, 1.3, 0.0, -4.023, 2.323, 0.0,
                -4.023, -2.323, 0.0, 0.0, -4.649, 0.0, 4.023, -2.323, 0.0,
                4.023, 2.323, 0.0 },
              std::vector<double>(11, 0.0), param, gen, hyb, itag);
  EXPECT_EQ(hyb, std::vector<int>({ 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0 }));
  // Formaldehyde: C=O double bonds.
  runPipeline(4, { 6, 8, 1, 1 },
              { 0.0, 0.0, 0.0, 2.28, 0.0, 0.0, -1.025, 1.775, 0.0, -1.025,
                -1.775, 0.0 },
              std::vector<double>(4, 0.0), param, gen, hyb, itag);
  EXPECT_EQ(hyb, std::vector<int>({ 2, 2, 0, 0 }));
}

TEST(GffHybTest, metalCarbonyl)
{
  // Fe(CO)5: TM rules give Fe hyb 0; the icase-2 HC filter isolates Fe
  // (group -8 takes crit 4 < 5 neighbours), so each CO unit sees a
  // singly-bonded C and the terminal O takes the CO branch (hyb 1).
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  std::vector<int> hyb, itag;
  runPipeline(11, { 26, 6, 6, 6, 6, 6, 8, 8, 8, 8, 8 },
              { 0.0, 0.0, 0.0, 0.0, 0.0, 3.4, 0.0, 0.0, -3.4, 3.4, 0.0, 0.0,
                -1.7, 2.945, 0.0, -1.7, -2.945, 0.0, 0.0, 0.0, 5.56, 0.0, 0.0,
                -5.56, 5.56, 0.0, 0.0, -2.78, 4.816, 0.0, -2.78, -4.816,
                0.0 },
              std::vector<double>(11, 0.0), param, gen, hyb, itag);
  EXPECT_EQ(hyb, std::vector<int>({ 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }));
  EXPECT_EQ(itag, std::vector<int>(11, 0));
}

TEST(GffHybTest, chargeOverride)
{
  // Carbon suboxide: central C is linear sp (hyb 1) at zero charge but
  // takes the qa < -0.4 override (hyb 2, itag reset) when charged.
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  std::vector<int> at = { 6, 6, 6, 8, 8 };
  std::vector<double> xyz = { 0.0, 0.0, 0.0, 2.44, 0.0, 0.0, 4.88, 0.0, 0.0,
                              7.08, 0.0, 0.0, -2.19, 0.0, 0.0 };
  std::vector<int> hyb, itag;
  runPipeline(5, at, xyz, std::vector<double>(5, 0.0), param, gen, hyb, itag);
  EXPECT_EQ(hyb, std::vector<int>({ 1, 1, 1, 2, 2 }));
  EXPECT_EQ(itag, std::vector<int>(5, 0));
  runPipeline(5, at, xyz, { 0.0, -0.5, 0.0, 0.0, 0.0 }, param, gen, hyb,
              itag);
  EXPECT_EQ(hyb, std::vector<int>({ 1, 2, 1, 2, 2 }));
  EXPECT_EQ(itag, std::vector<int>(5, 0));
}

TEST(GffHybTest, bondGuesses)
{
  // Direct gfnffrab values (normcn inputs, water).
  std::vector<int> at = { 8, 1, 1 };
  std::vector<double> rab;
  gfnffBondGuesses(3, at, { 4.0, 1.0, 1.0 }, rab);
  ASSERT_EQ(rab.size(), 6);
  EXPECT_NEAR(rab[packedIndex(0, 1)], 2.325926574183, 1e-9);
  EXPECT_NEAR(rab[packedIndex(1, 2)], 1.472800672054, 1e-9);
  EXPECT_EQ(elementRow6(1), 1);
  EXPECT_EQ(elementRow6(10), 2);
  EXPECT_EQ(elementRow6(18), 3);
  EXPECT_EQ(elementRow6(36), 4);
  EXPECT_EQ(elementRow6(54), 5);
  EXPECT_EQ(elementRow6(86), 6);
  EXPECT_EQ(elementRow6(118), 6);
  EXPECT_EQ(elementRow6(0), 0);
}

} // namespace Xtb
} // namespace Avogadro
