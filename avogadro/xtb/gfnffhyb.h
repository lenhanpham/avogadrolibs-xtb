/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_ini2.F90 (gfnff_neigh hybridization section, nn_nearest_noM)
  and src/gfnff/neighbor.f90 (nbLoc),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#ifndef AVOGADRO_XTB_GFNFFHYB_H
#define AVOGADRO_XTB_GFNFFHYB_H

#include <vector>

namespace Avogadro {
namespace Xtb {

class Environment;

// Neighbor location mirroring nbLoc: for atom indexi, one column per
// non-empty cell with the member list, the cell index and the member count.
// Only the first maxNeighbors - 2 slots are scanned like the reference.
struct NeighborCell
{
  std::vector<int> members; // 0-based atom indices
  int cell = 0;
  int count = 0;
};
std::vector<NeighborCell> locateNeighbors(const std::vector<int>& lists,
                                          const std::vector<int>& counts,
                                          int n, int numctr, int indexi);

// Coordination number of the closest non-metal neighbour, mirroring
// nn_nearest_noM (the unused packed-distance argument is dropped).
int nearestNonMetalCoord(int ii, const std::vector<int>& numbers,
                         const std::vector<double>& xyz,
                         const std::vector<int>& lists,
                         const std::vector<int>& counts, int numctr,
                         const int* metal);

// Bond angle A-B-C in radians from Cartesian coordinates (Bohr), the 0d
// equivalent of banglPBC used by the hybridization geometry checks.
double bondAngle(const std::vector<double>& xyz, int a, int b, int c);

// Hybridization assignment mirroring the gfnff_neigh section of ini2.
// nbf/nb/nbm are the full/reduced/metal-free neighbor lists with counts
// (layouts [(cell * n + atom) * maxNeighbors + slot] / [(cell * n) + atom],
// as produced by fillNeighborList); qa holds topology charges; metal/group
// are per-element tables; linThr is the linearity threshold in degrees.
// hyb/itag are overwritten (0 = unknown). Returns false with an environment
// error when too many atoms are hypervalent (nbCall set).
bool assignHybridization(int n, const std::vector<int>& numbers,
                         const std::vector<double>& xyz,
                         const std::vector<int>& nbfLists,
                         const std::vector<int>& nbfCounts,
                         const std::vector<int>& nbLists,
                         const std::vector<int>& nbCounts,
                         const std::vector<int>& nbmLists,
                         const std::vector<int>& nbmCounts,
                         const std::vector<double>& qa, const int* metal,
                         const int* group, double linThr, int numctr,
                         bool nbCall, std::vector<int>& hyb,
                         std::vector<int>& itag, Environment& env);

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_GFNFFHYB_H
