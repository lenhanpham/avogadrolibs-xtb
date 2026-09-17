/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_rab.f (gfnffrab, iTabRow6),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#ifndef AVOGADRO_XTB_GFNFFRAB_H
#define AVOGADRO_XTB_GFNFFRAB_H

#include <vector>

namespace Avogadro {
namespace Xtb {

// Legacy periodic-table row for force-constant scaling, mirroring iTabRow6:
// 1 (H-He), 2 (Li-Ne), 3 (Na-Ar), 4 (K-Kr), 5 (Rb-Xe), 6 (Cs and heavier),
// 0 outside 1..118.
int elementRow6(int atomicNumber);

// Bond-length guesses from "normal" coordination numbers, mirroring
// gfnffrab: rabPacked[packedIndex(j, i)] estimates in Bohr. cn holds the
// per-atom input CNs (production passes normal CNs as doubles).
void gfnffBondGuesses(int n, const std::vector<int>& numbers,
                      const std::vector<double>& cn,
                      std::vector<double>& rabPacked);

// Charge- and element-scaled guesses, mirroring the scaling loop in
// gfnff_neigh: rtmp(k) = (rab(k) - qa(i)*f1 - qa(j)*f2) * fat(ai)*fat(aj)
// with f1/f2 doubled for metals (fq = rqshrink). Operates in place.
void scaleBondGuesses(int n, const std::vector<int>& numbers,
                      const std::vector<double>& qa,
                      const std::vector<int>& metal, double fq,
                      std::vector<double>& rabPacked);

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_GFNFFRAB_H
