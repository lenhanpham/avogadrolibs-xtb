/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/mrec.f90 (mrecgff, mrecgffPBC), src/gfnff/getring36.f90
  (getring36, chkrng), src/gfnff/gfnff_ini2.F90 (ringsatom, ringsbond) and
  src/gfnff/gfnff_eg.f90 (gfnff_dlogcoord),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#ifndef AVOGADRO_XTB_GFFGRAPH_H
#define AVOGADRO_XTB_GFFGRAPH_H

#include <utility>
#include <vector>

namespace Avogadro {
namespace Xtb {

class Environment;

// Fragment assignment by connected components, mirroring mrecgff.
// adjacency[i] lists the neighbours of atom i (0-based). fragOfAtom receives
// 1-based fragment ids; the return value is the fragment count.
int findFragments(const std::vector<std::vector<int>>& adjacency,
                  std::vector<int>& fragOfAtom);

// Periodic variant mirroring mrecgffPBC: cells mark directed (neighbour,
// cell) edges, traversed j-fastest like maxloc over bond(:, i, :).
// neighbours[i] lists (neighbour, cell) pairs per atom.
int findFragmentsPbc(const std::vector<std::vector<std::pair<int, int>>>& neighbours,
                     int numCells, std::vector<int>& fragOfAtom);

// A ring found through an atom: 0-based member list. The hetero code
// (1000 * std-dev of atomic numbers / size, 0 for homogeneous) is a
// port-added convenience; the reference only returns members and sizes.
struct Ring
{
  std::vector<int> members;
  int hetero = 0;
};

// Ring enumeration through atom a0, mirroring the live getring36 in
// gfnff_ini2.F90 (NOT the dead getring36.f90 file): depth-first search for
// 3- to 6-membered rings over the neighbour lists (leaves pruned, branch
// order from the reference's selection sort, branches to atom 0 skipped),
// duplicates removed by atom-set equality over zero-tailed columns.
// At most 500 candidates are collected and 19 rings returned, exactly like
// the reference (truncation, no error).
bool findRingsThrough(int n, const std::vector<int>& numbers,
                      const std::vector<std::vector<int>>& neighbours, int a0,
                      std::vector<Ring>& rings, Environment& env);

// Smallest ring containing atom i (99 when in no ring), mirroring ringsatom.
int smallestRingThrough(const std::vector<Ring>& rings);
// Smallest ring containing bond i-j (0 when in no common ring), mirroring
// ringsbond over the per-atom tables of both ends.
int smallestRingBond(const std::vector<Ring>& ringsOfI,
                     const std::vector<Ring>& ringsOfJ, int i, int j);

// GFN-FF error-function coordination number with log cutoff, mirroring
// gfnff_dlogcoord (non-periodic): rabPacked holds packed pair distances in
// Bohr, rcov per-element radii in Bohr, thr the squared-distance cutoff.
// Like the reference, only the log-cutoff values leave the routine:
// logCn per atom; dlogCn[(i * n + j) * 3 + c] = d(logCn[i])/d(x[j][c]).
// Note kn = -7.5 and the 11-digit sqrtpi literal below match the original.
void gffCoordinationNumber(int n, const std::vector<int>& numbers,
                           const std::vector<double>& xyz,
                           const std::vector<double>& rabPacked,
                           const std::vector<double>& rcov, double cnmax,
                           double thr, std::vector<double>& logCn,
                           std::vector<double>& dlogCn);

// Metallic-character estimator, mirroring the mchar formula in gfnff_ini:
// mchar(i) = exp(-0.005 * en^8) * |dCN| / (cn + 1), with en from the
// (single-precision) GFN-FF en table and dlogCn/cn from gffCoordinationNumber.
void metallicCharacter(int n, const std::vector<int>& numbers,
                       const std::vector<double>& enTable,
                       const std::vector<double>& cn,
                       const std::vector<double>& dlogCn,
                       std::vector<double>& mchar);

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_GFFGRAPH_H
