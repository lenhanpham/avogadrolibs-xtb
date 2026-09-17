/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_param.f90 (newGFNFFGenerator, gfnff_set_param,
  gfnff_thresholds, gfnff_load_param, loadGFNFFAngewChem2020),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#ifndef AVOGADRO_XTB_GFNFFSETUP_H
#define AVOGADRO_XTB_GFNFFSETUP_H

#include "gfnffdata.h"

namespace Avogadro {
namespace Xtb {

// Default generator settings, mirroring newGFNFFGenerator. Unsuffixed
// literals in the original are single precision; they are written with the
// f suffix here so the same rounding applies.
void makeDefaultGenerator(GffGenerator& gen);

// Accuracy-dependent neighbour-list thresholds, mirroring gfnff_thresholds.
void gffThresholds(double accuracy, double& dispThr, double& cnThr,
                   double& repThr, double& hbThr1, double& hbThr2);

// Fill a GffData with tables and derived parameters, mirroring the
// gfnff_load_param + gfnff_set_param sequence (version selects the
// parametrisation; returns false for unknown versions).
bool loadGffParams(int version, GffData& param, GffGenerator& gen);

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_GFNFFSETUP_H
