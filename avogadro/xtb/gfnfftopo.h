/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/lin.f90 (lin), src/gfnff/neighbor.f90 (fillnb),
  src/gfnff/gfnff_ini.f90 (rabd Floyd-Warshall block, rtmp estimator) and
  src/gfnff/gfnff_ini2.F90 (goedeckera),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#ifndef AVOGADRO_XTB_GFNFFTOPO_H
#define AVOGADRO_XTB_GFNFFTOPO_H

#include <vector>

namespace Avogadro {
namespace Xtb {

class Environment;

// Packed symmetric-pair index, mirroring lin() (0-based here, 1-based in
// xtb: lin(i+1,j+1) - 1).
inline int packedIndex(int i, int j)
{
  int a = i > j ? i : j;
  int b = i > j ? j : i;
  return b + a * (a + 1) / 2;
}

// Maximum neighbours per atom in the xtb lists (numnb).
constexpr int maxNeighbors = 42;

// Neighbor-list fill, mirroring fillnb/getnb for one pass:
//   icase 1: full list (metal-scaled radii), written to full
//   icase 2: no hypervalent heavy atoms, written to plain
//   icase 3: no metals / hypervalent atoms, written to noMetal
// Layouts (xtb uses 1-based (slot, atom, cell) arrays with the count in the
// last slot; here 0-based with a separate counts array):
//   dist[(cell * n + i) * n + j] full distance matrix per cell in Bohr
//   radPacked[packedIndex(i, j)] reference pair lengths in Bohr
//   list[(cell * n + atom) * maxNeighbors + slot] neighbour indices
//   counts[(cell * n) + atom] number of neighbours
// fullCounts (summed nbf counts per atom) is required for icase > 1.
// metal/group/normCn are per-element tables indexed by Z - 1.
void fillNeighborList(int n, const std::vector<int>& numbers,
                      const std::vector<double>& radPacked,
                      const std::vector<double>& dist, const double* mchar,
                      int icase, double f, double f2, const int* metal,
                      const int* group, const int* normCn,
                      const int* fullCounts, int numctr,
                      std::vector<int>& list, std::vector<int>& counts);

// Estimated covalent bond lengths, mirroring the rabd Floyd-Warshall block
// and rtmp estimator in gfnff_ini. radAngstrom holds per-atom covalent radii
// in Angstrom; rabdOut is the n x n topology distance matrix in Angstrom
// (row-major). rabdOut uses float storage exactly like xtb's real(sp) rabd.
// rtmpPacked receives rfgoed1-scaled Bohr lengths for j < i pairs.
void estimateBondLengths(int n, const std::vector<int>& numbers,
                         const std::vector<int>& nbCounts,
                         const std::vector<int>& nbList,
                         const std::vector<double>& radAngstrom,
                         double rfgoed1, double tdistThr,
                         std::vector<float>& rabdOut,
                         std::vector<double>& rtmpPacked);

// Topology EEQ charges, mirroring goedeckera (non-periodic): solves the
// (n + nfrag) Lagrangian system with per-fragment charge constraints.
// pairPacked holds packed pair distances in Bohr; chieeq/gameeq/alpeeqSq are
// per-atom EEQ parameters (alpeeqSq = alpha^2); fragOfAtom holds 1-based
// fragment ids; fragCharges the nfrag constraints. charges and energy are
// overwritten (energy is NOT accumulated, unlike eeqCharges).
bool topologyCharges(int n, const std::vector<double>& pairPacked,
                     const std::vector<double>& chieeq,
                     const std::vector<double>& gameeq,
                     const std::vector<double>& alpeeqSq, int nfrag,
                     const std::vector<int>& fragOfAtom,
                     const std::vector<double>& fragCharges, Environment& env,
                     std::vector<double>& charges, double& energy);

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_GFNFFTOPO_H
