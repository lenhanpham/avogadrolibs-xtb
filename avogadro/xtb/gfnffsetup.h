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
#include "gfnffdriver.h"

#include <vector>

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

// Damping/scaling parameters for the HB/XB kernels, resolved from the
// loaded tables exactly like the driver test assembles them.
HbParams makeHbParams(const GffData& param);

// Full 0d setup mirroring the gfnff_ini stage order: neighbour lists,
// hybridization, EEQ xi/initial parameters, coordination numbers,
// topology (goedeckera) charges, gamma/final EEQ parameters, HB
// strengths, bpair flags, bonds/types, HB perception + triplets + maps +
// geometry lists, non-bonded tables, vbond terms, bends, torsions and
// bonded-ATM triples into a complete DriverInput. xyz is flat row-major
// in Bohr; totalCharge goes to the driver's EEQ solve and (for the setup
// charges) onto fragment 0. Known deviations: mcf scalings are 1.0
// (angewChem2020) and the setup charges use the topology (rtmp-based)
// EEQ while the driver re-solves with the final parameters. HB
// perception uses the pre-Hueckel pi flags like the reference; bond
// types, vbond, torsions and angles use the post-Hueckel piadr, pibo
// and pbo; bends, torsions and vbond use ring lists over the icase-3
// neighbour lists with per-bond ring fields.
bool gfnffSetup0d(int version, const std::vector<int>& numbers,
                  const std::vector<double>& xyz, double totalCharge,
                  double accuracy, DriverInput& in, GffData& param,
                  GffGenerator& gen, Environment& env);

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_GFNFFSETUP_H
