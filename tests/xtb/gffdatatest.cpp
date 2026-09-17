/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "avogadro/xtb/COPYING.LESSER".
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/xtb/gfnffdata.h>
#include <avogadro/xtb/gfnffparams.h>
#include <avogadro/xtb/gfnffsetup.h>

namespace Avogadro {
namespace Xtb {

// Reference values produced by compiling the unmodified xtb sources
// (gfnff_param.f90, data.f90, generator.f90) with gfortran. Tolerances are
// 1e-9 or exact; single-precision literals from the original are mirrored
// with f-suffixed literals here.

TEST(GffDataTest, tables)
{
  EXPECT_NEAR(gffChi[5], 1.311555, 1e-12);
  EXPECT_NEAR(gffGam[7], -0.031236, 1e-12);
  EXPECT_NEAR(gffAlp[0], 0.585069, 1e-12);
  EXPECT_NEAR(gffBond[5], 0.385248, 1e-12);
  EXPECT_NEAR(gffTors[102], 0.286192, 1e-12);
  EXPECT_NEAR(gffEn[0], 2.200000047684, 1e-9); // float(2.2), as in xtb
  EXPECT_NEAR(gffRad[0], 0.32, 1e-12);
  EXPECT_NEAR(gffRepz[25], 8.0, 1e-12);
  EXPECT_EQ(gffMetal[25], 2);
  EXPECT_EQ(gffGroup[1], 8);
  EXPECT_EQ(gffNormCn[7], 2);
  EXPECT_NEAR(sqrtZr4r2Value(1), 2.007348998057, 1e-9);
  EXPECT_NEAR(sqrtZr4r2Value(26), 5.584156058385, 1e-9);
  EXPECT_NEAR(sqrtZr4r2Value(118), 6.859891269321, 1e-9);
  EXPECT_NEAR(sqrtZr4r2Value(0), -1.0, 1e-12);
  EXPECT_NEAR(sqrtZr4r2Value(119), -1.0, 1e-12);
}

TEST(GffDataTest, generatorDefaults)
{
  GffGenerator gen;
  makeDefaultGenerator(gen);
  EXPECT_NEAR(gen.rThr, 1.25, 1e-12);
  EXPECT_NEAR(gen.rShrink, 0.23f, 1e-12);
  EXPECT_NEAR(gen.bstren[1], 1.24, 1e-12);
  EXPECT_NEAR(gen.bsmat[1][0], 1.3234, 1e-12);
  EXPECT_NEAR(gen.split1, 0.33, 1e-12);
  EXPECT_NEAR(gen.hdiag[4], -0.5, 1e-12);
  EXPECT_NEAR(gen.hoffdiag[7], 1.10, 1e-12);
  EXPECT_EQ(gen.torsf[3], 0.0); // never assigned by the reference either
  EXPECT_NEAR(gen.qfacbm[3], 0.70, 1e-12);
  EXPECT_NEAR(gen.d3a2, 4.80, 1e-12);
  EXPECT_NEAR(gen.cnMax, 4.4f, 1e-12);
  EXPECT_NEAR(gen.maxHIter, 5.0, 1e-12);
  EXPECT_NEAR(gen.tdistThr, 12.0, 1e-12);
  EXPECT_NEAR(gen.bstren[8], 3.40, 1e-12); // mean of entries 7 and 8
}

TEST(GffDataTest, thresholds)
{
  double disp, cn, rep, hb1, hb2;
  gffThresholds(1.0, disp, cn, rep, hb1, hb2);
  EXPECT_NEAR(disp, 1500.0, 1e-12);
  EXPECT_NEAR(cn, 100.0, 1e-12);
  EXPECT_NEAR(rep, 400.0, 1e-12);
  EXPECT_NEAR(hb1, 200.0, 1e-12);
  EXPECT_NEAR(hb2, 400.0, 1e-12);
  gffThresholds(100.0, disp, cn, rep, hb1, hb2);
  EXPECT_NEAR(disp, -500.0, 1e-9);
  EXPECT_NEAR(cn, 0.0, 1e-9);
  EXPECT_NEAR(rep, 200.0, 1e-9);
  EXPECT_NEAR(hb1, 100.0, 1e-9);
  EXPECT_NEAR(hb2, 300.0, 1e-9);
}

TEST(GffDataTest, loadParams)
{
  GffData param;
  GffGenerator gen;
  EXPECT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  EXPECT_EQ(param.chi.size(), 103);
  EXPECT_EQ(param.d3r0.size(), 5356);
  EXPECT_NEAR(param.chi[5], 1.311555, 1e-12);
  EXPECT_NEAR(param.rcov[7], 1.587369797410, 1e-9);
  EXPECT_NEAR(param.cnMax, 4.400000095367, 1e-9);
  EXPECT_NEAR(param.repScaleB, 1.758299946785, 1e-9);
  EXPECT_NEAR(param.xhBas[5], 0.80, 1e-12);
  EXPECT_NEAR(param.xhAci[6], 1.600000001490, 1e-9);
  EXPECT_NEAR(param.zb3atm[0], -0.167358237521, 1e-9);
  EXPECT_NEAR(param.zb3atm[25], -17.405256702136, 1e-9);
  EXPECT_NEAR(param.d3r0[0], 46.465512903178, 1e-9);
  EXPECT_NEAR(param.xbAci[16], 0.5, 1e-12);
  EXPECT_NEAR(param.torsHb, 0.939999997616, 1e-9);
  EXPECT_NEAR(param.hbLongCut, 85.0, 1e-12);

  GffData bad;
  GffGenerator badGen;
  EXPECT_FALSE(loadGffParams(99, bad, badGen));
}

} // namespace Xtb
} // namespace Avogadro
