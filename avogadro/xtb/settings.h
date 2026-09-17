/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/setparam.f90, Copyright (C) 2017-2020 Stefan Grimme.
******************************************************************************/

#ifndef AVOGADRO_XTB_SETTINGS_H
#define AVOGADRO_XTB_SETTINGS_H

namespace Avogadro {
namespace Xtb {

// Run-type selectors, verbatim values of xtb_setparam's p_run_* constants.
enum RunType
{
  RunPreScc = 1,
  RunScc = 2, // single point
  RunGrad = 3, // single point + gradients
  RunOpt = 4, // geometry optimisation
  RunHess = 5, // Hessian / frequencies
  RunMd = 6,
  RunOHess = 7,
  RunOMd = 8,
  RunPath = 10,
  RunScreen = 11,
  RunModef = 13,
  RunMdOpt = 14,
  RunMetaOpt = 15,
  RunBHess = 71,
  RunVip = 100, // vertical ionisation potential
  RunVea = 101, // vertical electron affinity
  RunViPea = 102,
  RunVFukui = 103,
  RunVOmega = 104
};

// Hamiltonian selectors, verbatim values of xtb_setparam's p_ext_* constants.
enum Hamiltonian
{
  ExtGfn0 = 0,
  ExtGfn1 = 1,
  ExtGfn2 = 2,
  ExtGfnFF = 13,
  ExtOniom = 14,
  ExtIff = 15,
  ExtTblite = 16,
  ExtPtb = 17,
  ExtMcGfnFF = 18
};

// Per-calculation options. Replaces xtb's global `set` (TSet in
// src/setparam.f90, owned by src/set_module.f90): instead of one mutable
// module-wide record, every calculation carries its own Settings value.
// Defaults match the TSet initialisation.
struct Settings
{
  Hamiltonian method = ExtGfn2;
  RunType runType = RunScc;
  int maxSccIterations = 250;
  double accuracy = 1.0; // SCC convergence accelerator threshold
  double electronicTemperature = 300.0; // Fermi smearing temperature (K)
  double broydenDamping = 0.40; // SCC mixing damping
  bool solveScc = true;
  bool periodic = false;
  bool newDispersion = true; // use D4 rather than D3
  int verbosity = 1;
};

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_SETTINGS_H
