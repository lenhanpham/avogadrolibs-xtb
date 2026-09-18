/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "avogadro/xtb/COPYING.LESSER".
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/xtb/environment.h>
#include <avogadro/xtb/gfnffdata.h>
#include <avogadro/xtb/gfnffegbond.h>
#include <avogadro/xtb/gfnffegnonbond.h>
#include <avogadro/xtb/gfnffsetup.h>
#include <avogadro/xtb/gfnfftopo.h>

#include <cmath>
#include <vector>

namespace Avogadro {
namespace Xtb {

// Reference values verified bitwise against the unmodified repulsion loop,
// goed_gfnff and ES gradient loop in xtb src/gfnff/gfnff_eg.f90, compiled
// with gfortran -ffp-contract=off (220 energy/charge/gradient outputs over
// an 11-atom water+ethanol input; the harness mirrors solveSymmetric
// op-for-op since production xtb uses LAPACK sytrf here).

namespace {

struct NonbondFixture
{
  std::vector<int> numbers;
  std::vector<double> xyz;
  std::vector<double> alphanb;
  std::vector<int> bpair;
  std::vector<double> chieeq, gameeq, alpeeq, cn;
  std::vector<double> dlogCn;
};

NonbondFixture makeFixture(const GffData& param)
{
  NonbondFixture f;
  f.numbers = { 8, 1, 1, 6, 6, 8, 1, 1, 1, 1, 1 };
  f.xyz = { 0.0, 0.0, 0.0, 0.96, 0.0, 0.0, -0.2485, 0.9269, 0.0, 2.5, 0.1,
            0.0, 3.9, -0.2, 0.1, 5.1, 0.5, -0.1, 2.2, 1.0, 0.2, 4.1, -1.2,
            0.0, 2.3, 0.6, 1.0, 4.5, 0.3, -0.9, 5.6, 1.2, 0.4 };
  int n = static_cast<int>(f.numbers.size());
  f.alphanb.assign(n * n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      double v = 1.5 + 0.13 * ((i * 7 + j * 13) % 11);
      f.alphanb[i * n + j] = v;
      f.alphanb[j * n + i] = v;
    }
  }
  f.bpair.assign(n * n, 0);
  const int bonds[][2] = { { 0, 1 }, { 0, 2 }, { 3, 4 }, { 4, 5 }, { 3, 6 },
                           { 3, 7 }, { 4, 8 }, { 5, 9 }, { 5, 10 } };
  for (const auto& b : bonds)
    f.bpair[b[0] * n + b[1]] = f.bpair[b[1] * n + b[0]] = 1;
  f.chieeq.resize(n);
  f.gameeq.resize(n);
  f.alpeeq.resize(n);
  f.cn.resize(n);
  for (int i = 0; i < n; ++i) {
    int z = f.numbers[i];
    f.chieeq[i] = -param.chi[z - 1] + 0.05 * ((i * 3) % 5);
    f.gameeq[i] = param.gam[z - 1];
    f.alpeeq[i] = param.alp[z - 1] * param.alp[z - 1];
    f.cn[i] = 1.0 + 0.5 * (i % 4);
  }
  f.dlogCn.assign(3 * n * n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int c = 0; c < 3; ++c)
        f.dlogCn[(i * n + j) * 3 + c] =
          0.01 * (((i * 5 + j * 11 + c * 17) % 13) - 6);
    }
  }
  return f;
}

void fillPacked(const NonbondFixture& f, std::vector<double>& sqrab,
                std::vector<double>& srab)
{
  int n = static_cast<int>(f.numbers.size());
  int npair = n * (n + 1) / 2;
  sqrab.assign(npair, 0.0);
  srab.assign(npair, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      double dx = f.xyz[3 * i] - f.xyz[3 * j];
      double dy = f.xyz[3 * i + 1] - f.xyz[3 * j + 1];
      double dz = f.xyz[3 * i + 2] - f.xyz[3 * j + 2];
      double r2 = dx * dx + dy * dy + dz * dz;
      sqrab[packedIndex(i, j)] = r2;
      srab[packedIndex(i, j)] = std::sqrt(r2);
    }
  }
}

TEST(NonbondTest, RepulsionMatchesReference)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  NonbondFixture f = makeFixture(param);
  int n = static_cast<int>(f.numbers.size());
  Environment env;
  double erep = 0.0;
  std::vector<double> grep(3 * n, 0.0);
  double sigma[3][3] = { { 0 } };
  ASSERT_TRUE(repulsionEnergyGradient(n, f.numbers, f.xyz, f.alphanb, f.bpair,
                                      param.repz, param.repScaleN, 100.0, 1.0,
                                      erep, grep, sigma, env));
  EXPECT_NEAR(erep, 6.09848138941394668e-01, 1.0e-12);
  const double kSigmaRef[9] = {
    -5.67254629897657514e-01, 8.72701659154763998e-02, 1.62666051785877425e-01,
    8.72701659154763998e-02, -1.29206282395338290e+00, -2.35264216896597933e-02,
    1.62666051785877425e-01, -2.35264216896597933e-02, -8.74286711417570861e-01,
  };
  const double kGrepRef[33] = {
    5.18687499491521983e-03, 1.30654917311156149e-03, 1.97031617098496854e-03,
    2.03440546385164289e-01, 5.67489855552336109e-02, 1.21359584393775886e-02,
    1.93669572711855330e-02, -1.39753599811784732e-02, 8.99840261063757259e-05,
    -2.40071511057058817e-01, 1.58691782584570007e-01, 3.35863850818813203e-01,
    3.28509597108148588e-01, -8.36014470666675891e-01, -3.10546689909721674e-01,
    -1.33389164635943809e-02, -1.95839595883254870e-02, 1.69308331978246240e-03,
    1.13582023265295321e-02, -1.99725429776465679e-01, 3.33538638488498684e-01,
    -1.74647857430226705e-01, 9.95121019036537380e-01, 8.29226197969096940e-02,
    1.33587517020394442e-02, -5.11616831413192440e-03, -6.87542473827873812e-01,
    -1.33935281689383762e-01, -1.21660817612823297e-01, 2.41675286950227192e-01,
    -1.92273631477189698e-02, -1.57921304098516818e-02, -1.18005742731046859e-02,
  };
  for (int i = 0; i < 3 * n; ++i)
    EXPECT_NEAR(grep[i], kGrepRef[i], 1.0e-12);
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      EXPECT_NEAR(sigma[a][b], kSigmaRef[a * 3 + b], 1.0e-12);
}

TEST(NonbondTest, EeqMatchesReference)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  NonbondFixture f = makeFixture(param);
  int n = static_cast<int>(f.numbers.size());
  std::vector<double> sqrab, srab;
  fillPacked(f, sqrab, srab);
  std::vector<int> frags(n, 1);
  std::vector<double> qfrag(1, 0.0);
  std::vector<double> charges, gamma, erf;
  double ees = 0.0;
  Environment env;
  ASSERT_TRUE(eeqCharges(n, f.numbers, srab, f.chieeq, f.gameeq, f.alpeeq,
                         f.cn, param.cnf, 1, frags, qfrag, 0.0, charges,
                         gamma, erf, ees, env));
  EXPECT_NEAR(ees, 1.79879729694718016e+00, 1.0e-12);
  const double kChargesRef[11] = {
    5.03452084614955275e+00, -4.74078166920610578e+00, -2.03165186727067582e+00,
    4.19958279864708128e+00, -4.52601530081901959e+00, 5.28408034649683511e+00,
    3.74612794488035861e-01, 2.01524195480773116e+00, -7.59398510914621183e-01,
    -1.84315340065403133e+00, -3.00703799172478181e+00,
  };
  for (int i = 0; i < n; ++i)
    EXPECT_NEAR(charges[i], kChargesRef[i], 1.0e-12);
}

TEST(NonbondTest, EsGradientMatchesReference)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  NonbondFixture f = makeFixture(param);
  int n = static_cast<int>(f.numbers.size());
  std::vector<double> sqrab, srab;
  fillPacked(f, sqrab, srab);
  std::vector<int> frags(n, 1);
  std::vector<double> qfrag(1, 0.0);
  std::vector<double> charges, gamma, erf;
  double ees = 0.0;
  Environment env;
  ASSERT_TRUE(eeqCharges(n, f.numbers, srab, f.chieeq, f.gameeq, f.alpeeq,
                         f.cn, param.cnf, 1, frags, qfrag, 0.0, charges,
                         gamma, erf, ees, env));
  std::vector<double> ges(3 * n, 0.0);
  ASSERT_TRUE(esGradient(n, f.xyz, srab, sqrab, gamma, erf, charges, f.cn,
                         param.cnf, f.numbers, f.dlogCn, ges, env));
  const double kGesRef[33] = {
    -5.59172482308776164e+00, -3.59383962418467551e+00, -2.27386280734972662e-01,
    1.24137719840783456e+00, 2.02436664508964226e+00, 5.87675965588500993e-01,
    1.42025632220339415e+00, 1.55307522321763081e+00, 6.15771453689081746e-02,
    1.36601295211810903e+00, -2.11610287908368916e-01, -8.77824722896922149e-01,
    4.74525730958730352e-01, 2.35612379144774131e+00, -1.33447892385409606e-01,
    1.10358934261611075e+00, -2.87856029543832959e+00, -9.29987407353659368e-01,
    3.01965864895814784e-01, -8.18015549226721317e-02, -1.97140385836554799e-01,
    2.92035410237130533e-01, -1.77083429772288170e+00, -5.94714120827587317e-02,
    -2.74473127718348753e-01, 6.54742431230863897e-02, 4.11992828811402978e-01,
    -1.42227073489310740e+00, -1.15472883721782452e-01, -4.72023485741492554e-01,
    1.10533311872433626e+00, 2.67643246602774854e+00, 1.82862035462377692e+00,
  };
  for (int i = 0; i < 3 * n; ++i)
    EXPECT_NEAR(ges[i], kGesRef[i], 1.0e-12);
}

} // namespace
} // namespace Xtb
} // namespace Avogadro
