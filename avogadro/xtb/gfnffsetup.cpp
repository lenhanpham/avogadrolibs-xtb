/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_param.f90 (newGFNFFGenerator, gfnff_set_param,
  gfnff_thresholds, gfnff_load_param, loadGFNFFAngewChem2020),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#include "gfnffsetup.h"

#include "constants.h"
#include "gfnffparams.h"

#include <cmath>

namespace Avogadro {
namespace Xtb {

void GffData::init(int ndim)
{
  en.assign(ndim, 0.0);
  rad.assign(ndim, 0.0);
  rcov.assign(ndim, 0.0);
  metal.assign(ndim, 0);
  group.assign(ndim, 0);
  normCn.assign(ndim, 0);
  repa.assign(ndim, 0.0);
  repan.assign(ndim, 0.0);
  repz.assign(ndim, 0.0);
  zb3atm.assign(ndim, 0.0);
  xhAci.assign(ndim, 0.0);
  xhBas.assign(ndim, 0.0);
  xbAci.assign(ndim, 0.0);
  chi.assign(ndim, 0.0);
  gam.assign(ndim, 0.0);
  cnf.assign(ndim, 0.0);
  alp.assign(ndim, 0.0);
  bond.assign(ndim, 0.0);
  angl.assign(ndim, 0.0);
  angl2.assign(ndim, 0.0);
  tors.assign(ndim, 0.0);
  tors2.assign(ndim, 0.0);
  d3r0.assign(ndim * (1 + ndim) / 2, 0.0);
}

void makeDefaultGenerator(GffGenerator& gen)
{
  gen = GffGenerator();
  gen.cnMax = 4.4f;
  gen.linThr = 160.0f;
  gen.fcThr = 1.0e-3;
  gen.tdistThr = 12.0f;
  gen.rThr = 1.25f;
  gen.rThr2 = 1.00f;
  gen.rShrink = 0.23f;
  gen.hQaThr = 0.01f;
  gen.qaBThr = 0.10f;
  gen.srb1 = 0.3731f;
  gen.srb2 = 0.3171f;
  gen.srb3 = 0.2538f;
  gen.qRepScal = 0.3480f;
  gen.nRepScal = -0.1270f;
  gen.hhFac = 0.6290f;
  gen.hh13Rep = 1.4580f;
  gen.hh14Rep = 0.7080f;
  gen.bstren = { 1.00, 1.24, 1.98, 1.22, 1.00, 0.78, 3.40, 3.40, 0.0 };
  gen.qFacBen = -0.54f;
  gen.qFacTor = 12.0;
  gen.fr3 = 0.3f;
  gen.fr4 = 1.0f;
  gen.fr5 = 1.5f;
  gen.fr6 = 5.7f;
  // torsf(4) is never assigned by the reference; it stays 0.0 here while
  // Fortran leaves it indeterminate.
  gen.torsf = { 1.00f, 1.18f, 1.05f, 0.0f, 0.50f, -0.90f, 0.70f, -2.00f };
  gen.fbs1 = 0.50f;
  gen.batmScal = 0.30;
  gen.mchiShift = -0.09;
  gen.rabShift = -0.110f;
  gen.rabShiftH = -0.050f;
  gen.hyperShift = 0.03f;
  gen.hShift3 = -0.11f;
  gen.hShift4 = -0.11f;
  gen.hShift5 = -0.06f;
  gen.metal1Shift = 0.2f;
  gen.metal2Shift = 0.15f;
  gen.metal3Shift = 0.05f;
  gen.etaShift = 0.040f;
  gen.qfacbm = { 1.0, -0.2, -0.2, 0.70, 0.50 };
  gen.qfacbm0 = 0.047f;
  gen.rfgoed1 = 1.175f;
  gen.hTriple = 1.45;
  gen.hueckelP2 = 1.00;
  gen.hueckelP3 = -0.24;
  // hdiag/hoffdiag use 1-based indices in xtb; unset entries stay 0.0.
  gen.hdiag[4] = -0.5;
  gen.hdiag[5] = 0.00;
  gen.hdiag[6] = 0.14;
  gen.hdiag[7] = -0.38;
  gen.hdiag[8] = -0.29;
  gen.hdiag[15] = -0.30;
  gen.hdiag[16] = -0.30;
  gen.hoffdiag[4] = 0.5;
  gen.hoffdiag[5] = 1.00;
  gen.hoffdiag[6] = 0.66;
  gen.hoffdiag[7] = 1.10;
  gen.hoffdiag[8] = 0.23;
  gen.hoffdiag[15] = 0.60;
  gen.hoffdiag[16] = 1.00;
  gen.hIter = 0.700;
  gen.hueckelP = 0.340;
  gen.bzRef = 0.370;
  gen.bzRef2 = 0.315;
  gen.pilpf = 0.530;
  gen.maxHIter = 5.0;
  gen.d3a1 = 0.58;
  gen.d3a2 = 4.80;
  gen.split0 = 0.670;
  gen.fringbo = 0.020;
  gen.aheavy3 = 89.0f;
  gen.aheavy4 = 100.0f;
  gen.split1 = 1.0 - gen.split0;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j)
      gen.bsmat[i][j] = -999.0f;
  }
  gen.bsmat[0][0] = gen.bstren[0];
  gen.bsmat[3][0] = gen.bstren[0];
  gen.bsmat[3][3] = gen.bstren[0];
  gen.bsmat[2][2] = gen.bstren[1];
  gen.bsmat[1][1] = gen.bstren[2];
  gen.bsmat[1][0] = gen.split0 * gen.bstren[0] + gen.split1 * gen.bstren[2];
  gen.bsmat[3][1] = gen.split0 * gen.bstren[0] + gen.split1 * gen.bstren[2];
  gen.bsmat[2][1] = gen.split0 * gen.bstren[1] + gen.split1 * gen.bstren[2];
  gen.bsmat[2][0] = gen.split0 * gen.bstren[0] + gen.split1 * gen.bstren[1];
  gen.bsmat[3][2] = gen.split0 * gen.bstren[0] + gen.split1 * gen.bstren[1];
  gen.bstren[8] = 0.5 * (gen.bstren[6] + gen.bstren[7]);
}

void gffThresholds(double accuracy, double& dispThr, double& cnThr,
                   double& repThr, double& hbThr1, double& hbThr2)
{
  dispThr = 1500.0 - std::log10(accuracy) * 1000.0;
  cnThr = 100.0 - std::log10(accuracy) * 50.0;
  repThr = 400.0 - std::log10(accuracy) * 100.0;
  hbThr1 = 200.0 - std::log10(accuracy) * 50.0;
  hbThr2 = 400.0 - std::log10(accuracy) * 50.0;
}

static void setGffParams(GffData& param, GffGenerator& gen)
{
  param.cnMax = 4.4f;
  param.angleCutA = 0.595;
  param.angleCutT = 0.505;
  param.angleCutANci = 0.395;
  param.angleCutTNci = 0.305;
  param.repScaleB = 1.7583f;
  param.repScaleN = 0.4270;
  param.hbAngleCut = 49.0f;
  param.hbShortCut = 22.0f;
  param.xbAngleCut = 70.0f;
  param.xbShortCut = 5.0f;
  param.hbSf = 1.0f;
  param.hbSt = 15.0f;
  param.xbSf = 0.03f;
  param.xbSt = 15.0f;
  param.hbAlp = 6.0f;
  param.hbLongCut = 85.0f;
  param.hbLongCutXb = 70.0f;
  param.hbAbMix = 0.80f;
  param.hbNbCut = 11.20f;
  param.torsHb = 0.94f;
  param.bendHb = 0.20f;
  param.vbondScale = 0.9f;
  param.xhAciGlobAbH = 0.268f;
  param.xhAciCoh = 0.350f;
  param.xhAciGlob = 1.50f;
  // Element indices below are 0-based (xtb uses 1-based element numbers).
  param.xhBas[5] = 0.80;
  param.xhBas[6] = 1.68;
  param.xhBas[7] = 0.67;
  param.xhBas[8] = 0.52;
  param.xhBas[13] = 4.0;
  param.xhBas[14] = 3.5;
  param.xhBas[15] = 2.0;
  param.xhBas[16] = 1.5;
  param.xhBas[34] = 1.5;
  param.xhBas[52] = 1.9;
  param.xhBas[32] = param.xhBas[14];
  param.xhBas[33] = param.xhBas[15];
  param.xhBas[50] = param.xhBas[14];
  param.xhBas[51] = param.xhBas[15];
  param.xhAci[5] = 0.75f;
  param.xhAci[6] = param.xhAciGlob + 0.1f;
  param.xhAci[7] = param.xhAciGlob;
  param.xhAci[8] = param.xhAciGlob;
  param.xhAci[14] = param.xhAciGlob;
  param.xhAci[15] = param.xhAciGlob;
  param.xhAci[16] = param.xhAciGlob + 1.0f;
  param.xhAci[34] = param.xhAciGlob + 1.0f;
  param.xhAci[52] = param.xhAciGlob + 1.0f;
  param.xbAci[14] = 1.0;
  param.xbAci[15] = 1.0;
  param.xbAci[16] = 0.5;
  param.xbAci[32] = 1.2;
  param.xbAci[33] = 1.2;
  param.xbAci[34] = 0.9;
  param.xbAci[50] = 1.2;
  param.xbAci[51] = 1.2;
  param.xbAci[52] = 1.2;

  // 3-atom bond prefactors and packed D3 R0^2 table.
  double cubeRoot = std::cbrt(gen.batmScal);
  int k = 0;
  for (int i = 0; i < gffElements; ++i) {
    param.zb3atm[i] = -(i + 1) * cubeRoot;
    for (int j = 0; j <= i; ++j) {
      double dum = sqrtZr4r2Value(i + 1) * sqrtZr4r2Value(j + 1) * 3.0;
      param.d3r0[k++] =
        (gen.d3a1 * std::sqrt(dum) + gen.d3a2) * (gen.d3a1 * std::sqrt(dum) + gen.d3a2);
    }
  }
  param.zb3atm[0] = -0.25 * cubeRoot;
}

bool loadGffParams(int version, GffData& param, GffGenerator& gen)
{
  if (version != GffAngewChem2020 && version != GffAngewChem2020_1 &&
      version != GffAngewChem2020_2 && version != GffHarmonic2020 &&
      version != GffMcGfnFF2023) {
    return false;
  }
  param.init(gffElements);
  for (int i = 0; i < gffElements; ++i) {
    param.en[i] = gffEn[i];
    param.rad[i] = gffRad[i];
    // Mirrors covalentRadD3: ([...] * aatoau) * 4/3 with xtb's own
    // Angstrom definition (0.52917726), not the CODATA one.
    param.rcov[i] =
      ((covalentRadD3Angstrom[i] * xtbAngstromToBohr) * 4.0) / 3.0;
    param.metal[i] = gffMetal[i];
    param.group[i] = gffGroup[i];
    param.normCn[i] = gffNormCn[i];
    param.repz[i] = gffRepz[i];
    param.chi[i] = gffChi[i];
    param.gam[i] = gffGam[i];
    param.cnf[i] = gffCnf[i];
    param.alp[i] = gffAlp[i];
    param.bond[i] = gffBond[i];
    param.repa[i] = gffRepa[i];
    param.repan[i] = gffRepan[i];
    param.angl[i] = gffAngl[i];
    param.angl2[i] = gffAngl2[i];
    param.tors[i] = gffTors[i];
    param.tors2[i] = gffTors2[i];
  }
  makeDefaultGenerator(gen);
  setGffParams(param, gen);
  return true;
}

} // namespace Xtb
} // namespace Avogadro
