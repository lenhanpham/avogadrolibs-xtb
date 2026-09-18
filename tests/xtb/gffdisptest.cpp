/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "avogadro/xtb/COPYING.LESSER".
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/xtb/environment.h>
#include <avogadro/xtb/gffgraph.h>
#include <avogadro/xtb/gfnffdata.h>
#include <avogadro/xtb/gfnffdisp.h>
#include <avogadro/xtb/gfnffparams.h>
#include <avogadro/xtb/gfnffsetup.h>
#include <avogadro/xtb/gfnfftopo.h>

#include <cmath>
#include <utility>
#include <vector>

namespace Avogadro {
namespace Xtb {

// Reference values verified bitwise: the model tables against an
// independent replica of newD3Model from xtb include/param_ref.fh
// (33124/33124 C6, 182/182 cnRef) and the loop against the verbatim
// 0d d3_gradient path in src/gfnff/gdisp0.f90 compiled with gfortran
// -ffp-contract=off (46 energy/gradient outputs over a 15-atom
// water+ethanol+methylsilane input).

namespace {

struct DispFixture
{
  std::vector<int> numbers;
  std::vector<double> xyz;
};

DispFixture makeFixture()
{
  DispFixture f;
  f.numbers = { 8, 1, 1, 6, 6, 8, 1, 1, 1, 1, 1, 14, 1, 1, 1 };
  f.xyz = { 0.0, 0.0, 0.0, 0.96, 0.0, 0.0, -0.2485, 0.9269, 0.0, 2.5, 0.1,
            0.0, 3.9, -0.2, 0.1, 5.1, 0.5, -0.1, 2.2, 1.0, 0.2, 4.1, -1.2,
            0.0, 2.3, 0.6, 1.0, 4.5, 0.3, -0.9, 5.6, 1.2, 0.4, 7.5, -0.5,
            0.3, 8.2, 0.4, -0.4, 6.8, -1.3, 1.0, 7.9, 0.2, 1.4 };
  return f;
}

TEST(DispTest, ModelMatchesReference)
{
  DispFixture f = makeFixture();
  Environment env;
  DispModel model;
  ASSERT_TRUE(buildDispersionModel(f.numbers, model, env));
  EXPECT_EQ(model.maxElem, 14);
  // Spot values from the bitwise-verified replica (H-H and C-C blocks).
  EXPECT_NEAR(model.c6[0 + 7 * (0 + 7 * (0 + 14 * 0))], 7.63624129448328937e+00, 0.0);
  EXPECT_NEAR(model.c6[0 + 7 * (1 + 7 * (0 + 14 * 0))], 4.75930570412998843e+00, 0.0);
  EXPECT_NEAR(model.c6[0 + 7 * (0 + 7 * (5 + 14 * 5))], 4.93833448663288266e+01, 0.0);
  EXPECT_EQ(model.nref[0], 2);
  EXPECT_EQ(model.nref[5], 7);
  EXPECT_NEAR(weightCn(4.0, 3.0, 1.98869274482158), 1.67230134465869602e-02, 0.0);
  EXPECT_NEAR(zetaCharge(8, 0.0), 1.00000000000000000e+00, 0.0);
  EXPECT_NEAR(zetaCharge(8, 0.35), 9.28003781176264586e-01, 0.0);
}

TEST(DispTest, GradientMatchesReference)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  DispFixture f = makeFixture();
  int n = static_cast<int>(f.numbers.size());
  int npair = n * (n + 1) / 2;
  std::vector<double> srab(npair, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      double dx = f.xyz[3 * i] - f.xyz[3 * j];
      double dy = f.xyz[3 * i + 1] - f.xyz[3 * j + 1];
      double dz = f.xyz[3 * i + 2] - f.xyz[3 * j + 2];
      srab[packedIndex(i, j)] =
        std::sqrt(dx * dx + dy * dy + dz * dz);
    }
  }
  std::vector<double> rcov(103);
  for (int z = 0; z < 103; ++z)
    rcov[z] = covalentRadiusD3(z + 1);
  std::vector<double> cn, dlogCn;
  gffCoordinationNumber(n, f.numbers, f.xyz, srab, rcov, gen.cnMax, 60.0, cn,
                        dlogCn);
  Environment env;
  DispModel model;
  ASSERT_TRUE(buildDispersionModel(f.numbers, model, env));
  std::vector<double> qa(n);
  for (int i = 0; i < n; ++i)
    qa[i] = -0.3 + 0.07 * ((i * 5) % 11);
  std::vector<double> zetac6;
  ASSERT_TRUE(buildZetaC6(f.numbers, qa, zetac6, env));
  double dispThr, cnThr, repThr, hbThr1, hbThr2;
  gffThresholds(1.0, dispThr, cnThr, repThr, hbThr1, hbThr2);
  std::vector<std::pair<int, int>> pairs;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      double dx = f.xyz[3 * i] - f.xyz[3 * j];
      double dy = f.xyz[3 * i + 1] - f.xyz[3 * j + 1];
      double dz = f.xyz[3 * i + 2] - f.xyz[3 * j + 2];
      if (dx * dx + dy * dy + dz * dz < dispThr)
        pairs.emplace_back(i, j);
    }
  }
  std::vector<double> grad(3 * n, 0.0);
  double edisp = 0.0;
  ASSERT_TRUE(d3Gradient(n, f.numbers, f.xyz, pairs, model, zetac6,
                         param.d3r0, 1.0, cn, dlogCn, edisp, grad, env));
  EXPECT_NEAR(edisp, -7.33795665711284642e-03, 1.0e-12);
  const double kGdispRef[45] = {
    -1.20191153476015535e-04, 1.82767667050689984e-06, -4.56561128486424362e-06,
    -6.19190753826164372e-05, 1.57615737872301167e-06, -3.68484538926353743e-06,
    -5.76147076192247117e-05, 9.40277000132264069e-06, -1.37950698838150875e-06,
    -3.54364957668944469e-05, 8.01833449949049058e-07, -3.40733434607635647e-06,
    -4.54819413198324063e-07, -2.30006827149763398e-06, 2.03619879252132622e-07,
    2.74779732766030470e-05, -6.58657499373800286e-06, 3.74830834670006697e-06,
    -2.91391473949197490e-05, 6.59658773382249096e-06, -1.38498242904600763e-06,
    5.33896746178043293e-07, -7.83716516007203501e-06, -1.33686558772540056e-06,
    -2.86086119407482991e-05, 4.09784402884082336e-06, 3.49134990229002228e-06,
    8.29812519380829011e-06, 6.71073712302053384e-07, -6.88224939393165224e-06,
    2.21523778391305249e-05, 4.96511701042612935e-06, 1.59342442557858538e-06,
    1.28476344683263330e-04, -1.55415911742155944e-05, 4.08006283222131528e-06,
    5.10998634916648894e-05, -3.82282083188137243e-06, 5.16648520979788855e-05,
    3.32645652566490900e-05, 1.59433770281276409e-05, 4.16491996781863261e-06,
    6.20608645063203026e-05, -9.79421658261610083e-06, -4.63051420325509203e-05,
  };
  for (int i = 0; i < 3 * n; ++i)
    EXPECT_NEAR(grad[i], kGdispRef[i], 1.0e-12);
}

} // namespace
} // namespace Xtb
} // namespace Avogadro
