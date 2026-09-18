/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "avogadro/xtb/COPYING.LESSER".
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/xtb/environment.h>
#include <avogadro/xtb/gfnffdata.h>
#include <avogadro/xtb/gfnffdriver.h>
#include <avogadro/xtb/gfnffsetup.h>

#include <cmath>
#include <vector>

namespace Avogadro {
namespace Xtb {

// Full 0d setup (gfnffSetup0d) followed by the driver single point. The
// assembly order mirrors gfnff_ini; every stage kernel is verified
// bitwise in its own batch. Water-dimer list/charge literals below
// reproduce the differentially verified batch-18 (EEQ setup) and
// batch-19 (HB perception) outputs; the end-to-end energies are
// assembly pins (deterministic regression locks, honest until a real
// xtb binary comparison becomes possible).

namespace {

std::vector<int> waterNumbers()
{
  return { 8, 1, 1, 8, 1, 1 };
}

std::vector<double> waterXyz()
{
  return { 0.0, 0.0, 0.0, 1.75, 0.0, 0.0, -0.5, 1.68, 0.0,
           5.0, 0.2, 0.0, 6.2, 1.55, 0.0, 5.8, -1.1, 0.0 };
}

std::vector<int> ethaneNumbers()
{
  return { 6, 6, 1, 1, 1, 1, 1, 1 };
}

std::vector<double> ethaneXyz()
{
  return { 0.0, 0.0, 0.0, 2.9, 0.1, 0.0, -0.7, 1.9, 0.9, -0.7, -0.9, 1.7,
           -0.7, -0.9, -1.7, 3.6, 1.9, -0.9, 3.6, -0.9, 1.7, 3.6, -0.9, -1.7 };
}

} // namespace

TEST(SetupTest, WaterDimerListsMatchVerified)
{
  Environment env;
  DriverInput in;
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(gfnffSetup0d(GffAngewChem2020, waterNumbers(), waterXyz(), 0.0,
                           1.0, in, param, gen, env));
  EXPECT_EQ(in.bonds.size(), 4u);
  EXPECT_EQ(in.bends.size(), 2u);
  EXPECT_EQ(in.torsions.size(), 0u);
  EXPECT_EQ(in.triples.size(), 0u);
  EXPECT_EQ(in.hb1.size(), 0u);
  ASSERT_EQ(in.hb2.size(), 4u);
  EXPECT_EQ(in.hb2[0].a, 0);
  EXPECT_EQ(in.hb2[0].b, 3);
  EXPECT_EQ(in.hb2[0].h, 1);
  EXPECT_EQ(in.hb2[3].a, 3);
  EXPECT_EQ(in.hb2[3].b, 0);
  EXPECT_EQ(in.hb2[3].h, 5);
  EXPECT_EQ(in.xb.size(), 0u);
  ASSERT_EQ(in.nrHb.size(), 4u);
  for (int v : in.nrHb)
    EXPECT_EQ(v, 1);
  // Final EEQ parameters reproduce the batch-18 verified literals.
  ASSERT_EQ(in.chieeq.size(), 6u);
  EXPECT_NEAR(in.chieeq[0], -1.72120099932944770e+00, 0.0);
  EXPECT_NEAR(in.chieeq[1], -1.22705400000000009e+00, 0.0);
}

TEST(SetupTest, WaterDimerEndToEnd)
{
  Environment env;
  DriverInput in;
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(gfnffSetup0d(GffAngewChem2020, waterNumbers(), waterXyz(), 0.0,
                           1.0, in, param, gen, env));
  HbParams hp = makeHbParams(param);
  DriverResult res;
  ASSERT_TRUE(gfnffSinglePoint(in, param, hp, res, env));
  EXPECT_NEAR(res.ees, -2.12024309316180937e-01, 0.0);
  EXPECT_NEAR(res.edisp, -6.71482950307334390e-04, 0.0);
  EXPECT_NEAR(res.erep, 1.35789298101648387e-01, 0.0);
  EXPECT_NEAR(res.ebond, -5.45837078013393762e-01, 0.0);
  EXPECT_NEAR(res.eangl, 3.00835588552948990e-03, 0.0);
  EXPECT_NEAR(res.etors, 0.00000000000000000e+00, 0.0);
  EXPECT_NEAR(res.ehb, -3.84175378489342483e-03, 0.0);
  EXPECT_LT(res.ehb, 0.0);
  EXPECT_NEAR(res.exb, 0.00000000000000000e+00, 0.0);
  EXPECT_NEAR(res.ebatm, 0.00000000000000000e+00, 0.0);
  EXPECT_NEAR(res.etot, -6.23576970077597492e-01, 0.0);
  EXPECT_NEAR(res.gnorm, 4.13434311929300136e-01, 0.0);
  ASSERT_EQ(res.gradient.size(), 18u);
  EXPECT_NEAR(res.gradient[0], 4.38657690312548379e-02, 0.0);
  EXPECT_NEAR(res.gradient[1], 5.03904668589960497e-02, 0.0);
  ASSERT_EQ(res.charges.size(), 6u);
  double qsum = 0.0;
  for (double q : res.charges)
    qsum += q;
  EXPECT_NEAR(qsum, 0.0, 1.0e-12);
  EXPECT_NEAR(res.charges[0], -6.98853824246023936e-01, 0.0);
}

TEST(SetupTest, FormamideEndToEnd)
{
  Environment env;
  std::vector<int> numbers = { 6, 8, 7, 1, 1, 1 };
  std::vector<double> xyz = { 0.0, 0.0, 0.0, 2.3, 0.0, 0.0, -1.2, 1.9, 0.0,
                              -0.5, -1.9, 0.0, -1.0, 2.6, 0.0, -2.2, 1.6,
                              0.0 };
  DriverInput in;
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(gfnffSetup0d(GffAngewChem2020, numbers, xyz, 0.0, 1.0, in,
                           param, gen, env));
  EXPECT_EQ(in.bonds.size(), 6u);
  EXPECT_EQ(in.bends.size(), 8u);
  EXPECT_EQ(in.torsions.size(), 4u);
  EXPECT_EQ(in.triples.size(), 12u);
  HbParams hp = makeHbParams(param);
  DriverResult res;
  ASSERT_TRUE(gfnffSinglePoint(in, param, hp, res, env));
  EXPECT_NEAR(res.ees, 1.47603749621580960e-02, 0.0);
  EXPECT_NEAR(res.edisp, -1.31776089068210778e-03, 0.0);
  EXPECT_NEAR(res.erep, 3.50328828468750997e+00, 0.0);
  EXPECT_NEAR(res.ebond, -7.82368108176616617e-01, 0.0);
  EXPECT_NEAR(res.eangl, 1.88754280357326165e-01, 0.0);
  EXPECT_NEAR(res.etors, 0.00000000000000000e+00, 0.0);
  EXPECT_NEAR(res.ehb, 0.00000000000000000e+00, 0.0);
  EXPECT_NEAR(res.ebatm, 4.71719390465433139e-02, 0.0);
  EXPECT_NEAR(res.etot, 2.97028900998623868e+00, 0.0);
  EXPECT_NEAR(res.gnorm, 1.70721855335344728e+01, 0.0);
  ASSERT_EQ(res.gradient.size(), 18u);
  EXPECT_NEAR(res.gradient[0], -1.97113702567717736e-01, 0.0);
  EXPECT_NEAR(res.gradient[1], 1.91228926084700684e-01, 0.0);
}

TEST(SetupTest, BenzeneEndToEnd)
{
  Environment env;
  std::vector<int> numbers(12, 1);
  for (int i = 0; i < 6; ++i)
    numbers[i] = 6;
  std::vector<double> xyz;
  for (int k = 0; k < 6; ++k) {
    double a = k * 1.04719755119659763 + 1.57079632679489662;
    xyz.push_back(2.64 * std::cos(a));
    xyz.push_back(2.64 * std::sin(a));
    xyz.push_back(0.0);
  }
  for (int k = 0; k < 6; ++k) {
    double a = k * 1.04719755119659763 + 1.57079632679489662;
    xyz.push_back(4.70 * std::cos(a));
    xyz.push_back(4.70 * std::sin(a));
    xyz.push_back(0.0);
  }
  DriverInput in;
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(gfnffSetup0d(GffAngewChem2020, numbers, xyz, 0.0, 1.0, in,
                           param, gen, env));
  EXPECT_EQ(in.bonds.size(), 12u);
  EXPECT_EQ(in.bends.size(), 18u);
  EXPECT_EQ(in.torsions.size(), 24u);
  HbParams hp = makeHbParams(param);
  DriverResult res;
  ASSERT_TRUE(gfnffSinglePoint(in, param, hp, res, env));
  EXPECT_NEAR(res.ees, -4.33668971620166602e-03, 0.0);
  EXPECT_NEAR(res.edisp, -6.25814085022990008e-03, 0.0);
  EXPECT_NEAR(res.erep, 1.42482101897516700e-01, 0.0);
  EXPECT_NEAR(res.ebond, -2.49112906782728372e+00, 0.0);
  EXPECT_NEAR(res.eangl, 1.31654757095640607e-31, 0.0);
  EXPECT_NEAR(res.etors, 0.00000000000000000e+00, 0.0);
  EXPECT_NEAR(res.ehb, 0.00000000000000000e+00, 0.0);
  EXPECT_NEAR(res.ebatm, -3.69267611555710272e-03, 0.0);
  EXPECT_NEAR(res.etot, -2.36293447261175560e+00, 0.0);
  EXPECT_NEAR(res.gnorm, 2.49800102053622693e-02, 0.0);
  ASSERT_EQ(res.gradient.size(), 36u);
  EXPECT_NEAR(res.gradient[0], -5.67698845425250800e-16, 0.0);
  EXPECT_NEAR(res.gradient[1], 5.24995225020835968e-03, 0.0);
}

TEST(SetupTest, EthaneEndToEnd)
{
  Environment env;
  DriverInput in;
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(gfnffSetup0d(GffAngewChem2020, ethaneNumbers(), ethaneXyz(),
                           0.0, 1.0, in, param, gen, env));
  EXPECT_EQ(in.bonds.size(), 7u);
  EXPECT_EQ(in.bends.size(), 12u);
  EXPECT_EQ(in.torsions.size(), 9u);
  EXPECT_EQ(in.triples.size(), 18u);
  EXPECT_EQ(in.hb1.size(), 0u);
  EXPECT_EQ(in.hb2.size(), 0u);
  EXPECT_EQ(in.xb.size(), 0u);
  // Final EEQ parameters reproduce the batch-18 verified literals.
  ASSERT_EQ(in.chieeq.size(), 8u);
  EXPECT_NEAR(in.chieeq[0], -1.31155500000000003e+00, 0.0);
  EXPECT_NEAR(in.chieeq[2], -1.22705400000000009e+00, 0.0);
  HbParams hp = makeHbParams(param);
  DriverResult res;
  ASSERT_TRUE(gfnffSinglePoint(in, param, hp, res, env));
  EXPECT_NEAR(res.ees, -2.36336850883162855e-03, 0.0);
  EXPECT_NEAR(res.edisp, -1.94599033432516248e-03, 0.0);
  EXPECT_NEAR(res.erep, 4.96391868898913916e-02, 0.0);
  EXPECT_NEAR(res.ebond, -1.10918583070999688e+00, 0.0);
  EXPECT_NEAR(res.eangl, 2.30053848513024004e-02, 0.0);
  EXPECT_NEAR(res.etors, 2.50282198278130178e-03, 0.0);
  EXPECT_NEAR(res.ehb, 0.00000000000000000e+00, 0.0);
  EXPECT_NEAR(res.exb, 0.00000000000000000e+00, 0.0);
  EXPECT_NEAR(res.ebatm, -1.32557322099706382e-05, 0.0);
  EXPECT_NEAR(res.etot, -1.03836105156138836e+00, 0.0);
  EXPECT_NEAR(res.gnorm, 9.64202483241317077e-02, 0.0);
  ASSERT_EQ(res.gradient.size(), 24u);
  EXPECT_NEAR(res.gradient[0], -1.84123135164911585e-02, 0.0);
  ASSERT_EQ(res.charges.size(), 8u);
  double qsum = 0.0;
  for (double q : res.charges)
    qsum += q;
  EXPECT_NEAR(qsum, 0.0, 1.0e-12);
  EXPECT_NEAR(res.charges[0], -6.97907019589956296e-02, 0.0);
}

} // namespace Xtb
} // namespace Avogadro
