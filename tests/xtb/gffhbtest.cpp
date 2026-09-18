/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "avogadro/xtb/COPYING.LESSER".
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/xtb/environment.h>
#include <avogadro/xtb/gfnffdata.h>
#include <avogadro/xtb/gfnffhb.h>
#include <avogadro/xtb/gfnffsetup.h>

#include <vector>

namespace Avogadro {
namespace Xtb {

// Reference values verified bitwise against the unmodified abhgfnff_eg1,
// eg2new, eg2_rnr, eg3 and rbxgfnff_eg in xtb src/gfnff/gfnff_eg.f90
// (plus helpers from constr.f90/basic_geo.f90/math.f), compiled with
// gfortran -ffp-contract=off (119 energy/gradient/virial outputs over
// five 0d fixtures: water dimer, pyridine-water, formaldehyde-water and
// methylchloride-pyridine).

namespace {

HbParams testParams(const GffData& p)
{
  HbParams h;
  h.hbacut = p.hbAngleCut;
  h.hblongcut = p.hbLongCut;
  h.hbscut = p.hbShortCut;
  h.hbalp = p.hbAlp;
  h.hbst = p.hbSt;
  h.hbsf = p.hbSf;
  h.hbabmix = p.hbAbMix;
  h.hbnbcut = p.hbNbCut;
  h.xbacut = p.xbAngleCut;
  h.xbscut = p.xbShortCut;
  h.xbst = p.xbSt;
  h.xbsf = p.xbSf;
  h.hblongcutXb = p.hbLongCutXb;
  h.xhaciGlobAbh = p.xhAciGlobAbH;
  h.xhaciCoh = p.xhAciCoh;
  h.torsHb = p.torsHb;
  h.bendHb = p.bendHb;
  return h;
}

struct HbFixture
{
  std::vector<int> numbers;
  std::vector<double> xyz;
  std::vector<double> qa, hbBas, hbAci;
};

void fillTables(HbFixture& f)
{
  int n = static_cast<int>(f.numbers.size());
  f.qa.assign(n, 0.0);
  f.hbBas.assign(n, 0.0);
  f.hbAci.assign(n, 0.0);
  for (int i = 0; i < n; ++i) {
    int z = f.numbers[i];
    f.qa[i] = -0.3 + 0.07 * ((i * 5) % 11);
    f.hbBas[i] = 0.1 * ((z * 7 + i * 3) % 13) - 0.3;
    f.hbAci[i] = 0.05 * ((z * 11 + i * 5) % 17) - 0.2;
  }
}

HbFixture waterDimer()
{
  HbFixture f;
  f.numbers = { 8, 1, 1, 8, 1, 1 };
  f.xyz = { 0.0, 0.0, 0.0, 0.96, 0.0, 0.0, -0.2485, 0.9269, 0.0,
            2.6, 0.1, 0.0, 3.3, 0.5, 0.0, 2.9, -0.7, 0.0 };
  fillTables(f);
  return f;
}

TEST(HbTest, Eg1MatchesReference)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  HbFixture f = waterDimer();
  Environment env;
  double energy = 0.0;
  double gdr[3][3] = { { 0 } };
  double sigma[3][3] = { { 0 } };
  ASSERT_TRUE(abhEg1(0, 3, 1, f.numbers, f.xyz, f.qa, f.hbBas, f.hbAci,
                     param.rad, testParams(param), 0.727406, energy, gdr,
                     sigma, env));
  EXPECT_NEAR(energy, -6.42518774503854861e-08, 1.0e-12);
  const double kG1Ref[9] = {
    2.16154627851853136e-07, 2.68632088832852068e-08, 0.00000000000000000e+00,
    -2.26099538049840912e-07, 1.93824800181079451e-09, 0.00000000000000000e+00,
    9.94491019798774936e-09, -2.88014568850960013e-08, 0.00000000000000000e+00,
  };
  for (int col = 0; col < 3; ++col) {
    for (int c = 0; c < 3; ++c)
      EXPECT_NEAR(gdr[c][col], kG1Ref[col * 3 + c], 1.0e-12);
  }
  const double kS1Ref[9] = {
    -4.20667389640596379e-07, -1.64466160574682622e-08, 0.00000000000000000e+00,
    -1.64466160574682589e-08, 1.40989322600518292e-10, 0.00000000000000000e+00,
    0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00,
  };
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c)
      EXPECT_NEAR(sigma[r][c], kS1Ref[r * 3 + c], 1.0e-12);
  }
}

TEST(HbTest, Eg2newMatchesReference)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  HbFixture f = waterDimer();
  int n = static_cast<int>(f.numbers.size());
  Environment env;
  double energy = 0.0;
  std::vector<double> gdr(3 * n, 0.0);
  double sigma[3][3] = { { 0 } };
  const std::vector<int> nbrs = { 4, 5 };
  ASSERT_TRUE(abhEg2new(0, 3, 1, nbrs, f.numbers, f.xyz, f.qa, f.hbBas,
                        f.hbAci, param.rad, testParams(param), 0.727406,
                        energy, gdr, sigma, env));
  EXPECT_NEAR(energy, -1.02567019914852543e-07, 1.0e-12);
  const double kG2Ref[18] = {
    4.79499296553245829e-07, 4.51916038295307198e-08, 0.00000000000000000e+00,
    -2.09073295526299839e-07, -5.91358202666922308e-08, 0.00000000000000000e+00,
    0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00,
    -2.48944756439302996e-07, 3.77405614000344796e-09, 0.00000000000000000e+00,
    -7.17038452995970919e-09, -2.49422309568689691e-09, 0.00000000000000000e+00,
    -1.43108600576831405e-08, 1.26643833928449604e-08, 0.00000000000000000e+00,
  };
  for (int i = 0; i < 3 * n; ++i)
    EXPECT_NEAR(gdr[i], kG2Ref[i], 1.0e-12);
  const double kS2Ref[9] = {
    -6.64216599801112220e-07, -1.34294174851653427e-08, 0.00000000000000000e+00,
    -1.34294174851653559e-08, -7.08113324089212369e-09, 0.00000000000000000e+00,
    0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00,
  };
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c)
      EXPECT_NEAR(sigma[r][c], kS2Ref[r * 3 + c], 1.0e-12);
  }
}

TEST(HbTest, Eg2rnrMatchesReference)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  HbFixture f;
  f.numbers = { 8, 1, 7, 6, 6 };
  f.xyz = { 0.0, 0.0, 0.0, 0.9, 0.3, 0.0, 2.6, 0.2, 0.0, 3.3, 1.3, 0.0,
            3.2, -1.0, 0.0 };
  fillTables(f);
  int n = static_cast<int>(f.numbers.size());
  Environment env;
  double energy = 0.0;
  std::vector<double> gdr(3 * n, 0.0);
  double sigma[3][3] = { { 0 } };
  const std::vector<int> nbrs = { 3, 4 };
  ASSERT_TRUE(abhEg2rnr(0, 2, 1, nbrs, f.numbers, f.xyz, f.qa, f.hbBas,
                        f.hbAci, param.rad, param.repz, testParams(param),
                        0.727406, energy, gdr, sigma, env));
  EXPECT_NEAR(energy, 3.79961504217404814e-26, 1.0e-12);
  const double kG3Ref[15] = {
    -2.03876458480851383e-25, 7.67997591088751203e-26, 0.00000000000000000e+00,
    9.32636247709811651e-26, -1.33044925327796000e-25, 0.00000000000000000e+00,
    1.07128453962931168e-25, 2.05306788917437026e-26, 0.00000000000000000e+00,
    1.66153739887048262e-27, 1.81624272385931883e-26, 0.00000000000000000e+00,
    1.82284234806870267e-27, 1.75520600885839886e-26, 0.00000000000000000e+00,
  };
  for (int i = 0; i < 3 * n; ++i)
    EXPECT_NEAR(gdr[i], kG3Ref[i], 1.0e-12);
  const double kS3Ref[9] = {
    2.71895205869642840e-25, 3.61825816027924459e-26, 0.00000000000000000e+00,
    3.61825816027924344e-26, -2.16390529924172595e-26, 0.00000000000000000e+00,
    0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00,
  };
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c)
      EXPECT_NEAR(sigma[r][c], kS3Ref[r * 3 + c], 1.0e-12);
  }
}

TEST(HbTest, Eg3MatchesReference)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  HbFixture f;
  f.numbers = { 8, 1, 8, 6, 1, 1 };
  f.xyz = { 0.0, 0.0, 0.0, 0.9, 0.3, 0.0, 2.6, 0.1, 0.0, 3.8, -0.1, 0.0,
            4.3, 0.8, 0.0, 4.2, -1.1, 0.0 };
  fillTables(f);
  int n = static_cast<int>(f.numbers.size());
  int npair = n * (n + 1) / 2;
  Environment env;
  double energy = 0.0;
  std::vector<double> gdr(3 * n, 0.0);
  double sigma[3][3] = { { 0 } };
  const std::vector<int> cNbrs = { 2, 4, 5 };
  const std::vector<double> sq(npair, 0.0), sr(npair, 0.0);
  ASSERT_TRUE(abhEg3(0, 2, 1, 3, cNbrs, f.numbers, f.xyz, f.qa, f.hbBas,
                     f.hbAci, sq, sr, param.rad, testParams(param), 0.727406,
                     energy, gdr, sigma, env));
  EXPECT_NEAR(energy, 6.66794807915363877e-10, 1.0e-12);
  const double kG4Ref[18] = {
    -3.68599533867142118e-09, 1.78693100450017754e-09, 0.00000000000000000e+00,
    1.63450029767916486e-09, -3.14926937125757148e-09, 0.00000000000000000e+00,
    2.05054106260055907e-09, 1.37269546818470552e-09, 0.00000000000000000e+00,
    9.53978391697783300e-13, -1.03571014273119132e-11, 0.00000000000000000e+00,
    0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00,
    0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00,
  };
  for (int i = 0; i < 3 * n; ++i)
    EXPECT_NEAR(gdr[i], kG4Ref[i], 1.0e-12);
  const double kS4Ref[9] = {
    4.95078499135627453e-09, 5.05771791317685909e-10, 0.00000000000000000e+00,
    5.05771791317686426e-10, -5.86635157135575489e-10, 0.00000000000000000e+00,
    0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00,
  };
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c)
      EXPECT_NEAR(sigma[r][c], kS4Ref[r * 3 + c], 1.0e-12);
  }
}

TEST(HbTest, RbxMatchesReference)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  HbFixture f;
  f.numbers = { 7, 6, 17, 1, 1, 1 };
  f.xyz = { 0.0, 0.0, 0.0, 5.0, 0.1, 0.0, 3.2, 0.05, 0.0, 5.4, 1.0, 0.0,
            5.4, -0.5, 0.9, 5.4, -0.5, -0.9 };
  fillTables(f);
  Environment env;
  double energy = 0.0;
  double gdr[3][3] = { { 0 } };
  double sigma[3][3] = { { 0 } };
  double cx = param.xbAci[f.numbers[2] - 1];
  ASSERT_TRUE(rbxEg(0, 1, 2, 1.0, cx, f.numbers, f.xyz, f.qa, param.rad,
                    testParams(param), energy, gdr, sigma, env));
  EXPECT_NEAR(energy, -6.13250534080049580e-04, 1.0e-12);
  const double kG5Ref[9] = {
    -2.61516982556254940e-07, 1.87806908103520980e-05, 0.00000000000000000e+00,
    -3.03362493650112964e-03, -5.08722002126678119e-05, 0.00000000000000000e+00,
    3.03388645348368633e-03, 3.20915094023157139e-05, 0.00000000000000000e+00,
  };
  for (int col = 0; col < 3; ++col) {
    for (int c = 0; c < 3; ++c)
      EXPECT_NEAR(gdr[c][col], kG5Ref[col * 3 + c], 1.0e-12);
  }
  const double kS5Ref[9] = {
    -5.45968803135785108e-03, -1.51668170975928767e-04, 0.00000000000000000e+00,
    -1.51668170975928658e-04, -3.48264455115099532e-06, 0.00000000000000000e+00,
    0.00000000000000000e+00, 0.00000000000000000e+00, 0.00000000000000000e+00,
  };
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c)
      EXPECT_NEAR(sigma[r][c], kS5Ref[r * 3 + c], 1.0e-12);
  }
}

} // namespace
} // namespace Xtb
} // namespace Avogadro
