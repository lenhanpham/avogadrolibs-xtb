/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/elem.f90, Copyright (C) 2017-2020 Stefan Grimme.
******************************************************************************/

#ifndef AVOGADRO_XTB_ELEMENTS_H
#define AVOGADRO_XTB_ELEMENTS_H

#include <string>

namespace Avogadro {
namespace Xtb {

// Number of elements tabulated (H..Og).
constexpr int maxElement = 118;

// Resolve an element symbol to its atomic number, mirroring xtb's elem().
// Matching is case-insensitive on up to two leading characters; returns 0
// when the symbol is not recognised.
int elementNumber(const std::string& symbol);

// Canonical lowercase two-character symbol for atomic numbers 1..118
// (e.g. "h ", "he"). Empty string outside the range.
std::string elementSymbol(int atomicNumber);

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_ELEMENTS_H
