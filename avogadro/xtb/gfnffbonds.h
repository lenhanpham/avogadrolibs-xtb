/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_ini.f90 (blist/btyp/imetal sections) and
  src/gfnff/gfnff_ini2.F90 (amide),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#ifndef AVOGADRO_XTB_GFNFFBONDS_H
#define AVOGADRO_XTB_GFNFFBONDS_H

#include <vector>

namespace Avogadro {
namespace Xtb {

// A covalent bond: 0-based endpoint atoms plus cell (always 0 outside PBC).
struct Bond
{
  int first = -1;
  int second = -1;
  int cell = 0;
};

// Unique covalent bond list mirroring the bdum/cdum construction: pairs in
// first-seen order (atom outer, list order inner), both directions marked.
std::vector<Bond> buildBondList(int n, const std::vector<int>& nbLists,
                                const std::vector<int>& nbCounts, int numctr);

// Effective metal classes mirroring the imetal setup: metal(Z) cleared when
// the nb count is <= 4 and the group exceeds 3 (Sn/Pb/Bi-like main groups).
std::vector<int> effectiveMetals(int n, const std::vector<int>& numbers,
                                 const std::vector<int>& nbCounts,
                                 int numctr, const int* metal,
                                 const int* group);

// Amide nitrogen predicate mirroring amide(): N with hyb 3 in the pi system
// bonded to exactly one pi carbon that itself bonds a terminal pi oxygen.
// nbLists/nbCounts use the [(cell * n + atom) * maxNeighbors + slot] layout.
bool isAmideNitrogen(int n, const std::vector<int>& numbers,
                     const std::vector<int>& hyb,
                     const std::vector<int>& nbLists,
                     const std::vector<int>& nbCounts, int numctr,
                     const std::vector<int>& piAtoms, int atom);

// Hückel pi-system setup mirroring the piadr section: candidate pi atoms
// (sp/sp2 B/C/N/O/F/S/Cl plus sp3 N/O/F on sp2), SO3/NR3-X exclusions, the
// pi-only neighbour lists and their fragment assignment (via the verified
// fragment finder). piIndex maps pi slots to atoms (piadr), atomPi maps
// atoms to 1-based pi slots with 0 for non-pi atoms (piadr2).
struct PiSystem
{
  std::vector<int> piIndex; // piadr
  std::vector<int> atomPi; // piadr2, 1-based, 0 = not pi
  std::vector<std::vector<int>> piNeighbours; // per pi slot (all cells flat)
  int fragments = 0; // picount
  std::vector<int> fragOfPi; // pimvec, 1-based
};
PiSystem buildPiSystem(int n, const std::vector<int>& numbers,
                       const std::vector<int>& hyb,
                       const std::vector<int>& nbLists,
                       const std::vector<int>& nbCounts, int numctr);

// Bond types mirroring the btyp assignment (piadr is the post-Hückel pi map
// where available, else the setup map; itag marks carbenes/eta atoms):
// 1 single, 2 pi, 3 triple/linear/no-torsion, 4 hypervalent, 5 metal,
// 6 eta, 7 TM-TM.
std::vector<int> assignBondTypes(
  const std::vector<Bond>& bonds, const std::vector<int>& numbers,
  const std::vector<int>& hyb, const std::vector<int>& itag,
  const std::vector<int>& piAtoms, const std::vector<int>& imetal,
  const int* group);

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_GFNFFBONDS_H
