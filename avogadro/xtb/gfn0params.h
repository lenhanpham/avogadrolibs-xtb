/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfn0param.f90, Copyright (C) 2017-2020 Stefan Grimme.
******************************************************************************/

#ifndef AVOGADRO_XTB_GFN0PARAMS_H
#define AVOGADRO_XTB_GFN0PARAMS_H

#include <array>

namespace Avogadro {
namespace Xtb {

// Number of elements covered by the GFN0 electrostatic parametrisation.
constexpr int gfn0Elements = 86;

// Per-element tables, index 0 = H. Values are exact doubles (the originals
// carry the _wp suffix).
extern const std::array<double, gfn0Elements> electronegativity; // xi
extern const std::array<double, gfn0Elements> hardness; // gamm
extern const std::array<double, gfn0Elements> cnFactor; // cnfak
extern const std::array<double, gfn0Elements> alpha; // alpg

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_GFN0PARAMS_H
