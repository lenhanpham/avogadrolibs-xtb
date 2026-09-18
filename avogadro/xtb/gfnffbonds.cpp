/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_ini.f90 (blist/btyp/imetal sections) and
  src/gfnff/gfnff_ini2.F90 (amide, Hückel piadr setup),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#include "gfnffbonds.h"

#include "gffgraph.h"
#include "gfnfftopo.h"

namespace Avogadro {
namespace Xtb {

std::vector<Bond> buildBondList(int n, const std::vector<int>& nbLists,
                                const std::vector<int>& nbCounts, int numctr)
{
  std::vector<Bond> bonds;
  std::vector<char> marked(n * n * numctr, 0);
  auto isMarked = [&](int l, int i, int cell) -> char& {
    return marked[(cell * n + i) * n + l];
  };
  for (int i = 0; i < n; ++i) {
    for (int cell = 0; cell < numctr; ++cell) {
      for (int k = 0; k < nbCounts[cell * n + i]; ++k) {
        int l = nbLists[(cell * n + i) * maxNeighbors + k];
        if (l < 0 || l >= n)
          continue;
        if (!isMarked(l, i, cell)) {
          isMarked(l, i, cell) = 1;
          // PBC pairing uses the opposite cell like iTrNeg; outside PBC
          // there is a single cell mapping to itself.
          isMarked(i, l, cell) = 1;
          Bond bond;
          bond.first = l;
          bond.second = i;
          bond.cell = cell;
          bonds.push_back(bond);
        }
      }
    }
  }
  return bonds;
}

std::vector<int> effectiveMetals(int n, const std::vector<int>& numbers,
                                 const std::vector<int>& nbCounts,
                                 int numctr, const int* metal,
                                 const int* group)
{
  std::vector<int> imetal(n);
  for (int i = 0; i < n; ++i) {
    imetal[i] = metal[numbers[i] - 1];
    // Sn/Pb/Bi-like main groups with small CN are better as non-metals.
    int count = 0;
    for (int cell = 0; cell < numctr; ++cell)
      count += nbCounts[cell * n + i];
    if (count <= 4 && group[numbers[i] - 1] > 3)
      imetal[i] = 0;
  }
  return imetal;
}

static bool isPiElement(int z)
{
  return z == 5 || z == 6 || z == 7 || z == 8 || z == 9 || z == 16 || z == 17;
}

static bool isNofElement(int z)
{
  return z == 7 || z == 8 || z == 9 || z == 16 || z == 17;
}

bool isAmideNitrogen(int n, const std::vector<int>& numbers,
                     const std::vector<int>& hyb,
                     const std::vector<int>& nbLists,
                     const std::vector<int>& nbCounts, int numctr,
                     const std::vector<int>& piAtoms, int atom)
{
  if (piAtoms[atom] == 0 || hyb[atom] != 3 || numbers[atom] != 7)
    return false;
  int nc = 0, ic = -1;
  for (int cell = 0; cell < numctr; ++cell) {
    for (int k = 0; k < nbCounts[cell * n + atom]; ++k) {
      int j = nbLists[(cell * n + atom) * maxNeighbors + k];
      if (numbers[j] == 6 && piAtoms[j] != 0) {
        ++nc;
        ic = j;
      }
    }
  }
  if (nc != 1)
    return false;
  int no = 0;
  for (int cell = 0; cell < numctr; ++cell) {
    for (int k = 0; k < nbCounts[cell * n + ic]; ++k) {
      int j = nbLists[(cell * n + ic) * maxNeighbors + k];
      if (numbers[j] == 8 && piAtoms[j] != 0 && nbCounts[cell * n + j] == 1)
        ++no;
    }
  }
  return no == 1;
}

PiSystem buildPiSystem(int n, const std::vector<int>& numbers,
                       const std::vector<int>& hyb,
                       const std::vector<int>& nbLists,
                       const std::vector<int>& nbCounts, int numctr)
{
  PiSystem pi;
  pi.atomPi.assign(n, 0);
  for (int i = 0; i < n; ++i) {
    bool piat =
      (hyb[i] == 1 || hyb[i] == 2) && isPiElement(numbers[i]);
    int kk = 0;
    for (int cell = 0; cell < numctr; ++cell) {
      for (int k = 0; k < nbCounts[cell * n + i]; ++k) {
        int jj = nbLists[(cell * n + i) * maxNeighbors + k];
        if (numbers[i] == 8 && numbers[jj] == 16 &&
            (jj < static_cast<int>(hyb.size()) ? hyb[jj] : 0) == 5)
          piat = false; // SO3 is not pi; skips counting this neighbour
        else if (hyb[jj] == 1 || hyb[jj] == 2)
          ++kk; // attached to sp2 or sp
      }
    }
    bool picon = kk > 0 && isNofElement(numbers[i]);
    if (numbers[i] == 7) {
      int total = 0;
      for (int cell = 0; cell < numctr; ++cell)
        total += nbCounts[cell * n + i];
      if (total > 3)
        continue; // NR3-X is not pi
    }
    if (numbers[i] == 16 && hyb[i] == 5)
      continue; // SO3 is not pi (bare hyb-5 sulphur exclusion)
    if (picon || piat) {
      pi.piIndex.push_back(i);
      pi.atomPi[i] = static_cast<int>(pi.piIndex.size()); // 1-based
    }
  }
  // Pi-only neighbour lists and their fragments.
  int npiall = static_cast<int>(pi.piIndex.size());
  std::vector<std::vector<int>> piAdj(npiall);
  for (int ii = 0; ii < npiall; ++ii) {
    int i = pi.piIndex[ii];
    for (int cell = 0; cell < numctr; ++cell) {
      for (int k = 0; k < nbCounts[cell * n + i]; ++k) {
        int jj = nbLists[(cell * n + i) * maxNeighbors + k];
        if (pi.atomPi[jj] > 0)
          piAdj[ii].push_back(pi.atomPi[jj] - 1);
      }
    }
  }
  pi.piNeighbours = piAdj;
  pi.fragOfPi.assign(npiall, 0);
  pi.fragments = findFragments(piAdj, pi.fragOfPi);
  return pi;
}

std::vector<int> assignBondTypes(
  const std::vector<Bond>& bonds, const std::vector<int>& numbers,
  const std::vector<int>& hyb, const std::vector<int>& itag,
  const std::vector<int>& piAtoms, const std::vector<int>& imetal,
  const int* group)
{
  std::vector<int> btyp;
  btyp.reserve(bonds.size());
  for (const Bond& bond : bonds) {
    int ii = bond.second, jj = bond.first;
    int ia = numbers[ii], ja = numbers[jj];
    int type = 1; // single
    if (hyb[ii] == 2 && hyb[jj] == 2)
      type = 2; // sp2-sp2 = pi
    if (hyb[ii] == 3 && hyb[jj] == 2 && ia == 7)
      type = 2; // N-sp2
    if (hyb[jj] == 3 && hyb[ii] == 2 && ja == 7)
      type = 2;
    if (hyb[ii] == 1 || hyb[jj] == 1)
      type = 3; // sp-X, no torsion
    if ((group[ia - 1] == 7 || ia == 1) && hyb[ii] == 1)
      type = 3; // linear halogen
    if ((group[ja - 1] == 7 || ja == 1) && hyb[jj] == 1)
      type = 3;
    if (hyb[ii] == 5 || hyb[jj] == 5)
      type = 4; // hypervalent
    if (imetal[ii] > 0 || imetal[jj] > 0)
      type = 5; // metal
    if (imetal[ii] == 2 && imetal[jj] == 2)
      type = 7; // TM-TM
    if (imetal[jj] == 2 && itag[ii] == -1 && piAtoms[ii] > 0)
      type = 6; // eta
    if (imetal[ii] == 2 && itag[jj] == -1 && piAtoms[jj] > 0)
      type = 6;
    btyp.push_back(type);
  }
  return btyp;
}

} // namespace Xtb
} // namespace Avogadro
