/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_ini.f90 (tlist/vtors setup) and
  src/gfnff/gfnff_ini2.F90 (ringstors, ringstorl, alphaCO),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#include "gfnfftorsion.h"

#include "constants.h"
#include "environment.h"
#include "gffgraph.h"
#include "gfnffangle.h"
#include "gfnffbonds.h"
#include "gfnffdata.h"
#include "gfnfftopo.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace Avogadro {
namespace Xtb {

bool dihedralNearLinear(const std::vector<double>& xyz, int ii, int jj,
                        int kk, int ll)
{
  // Exact chktors comparison: (angle * 180) / 3.1415926 > 170.
  if ((bondAnglePbc(xyz, jj, ii, kk) * 180.0f) / 3.1415926 > 170.0)
    return true;
  if ((bondAnglePbc(xyz, ii, jj, ll) * 180.0f) / 3.1415926 > 170.0)
    return true;
  return false;
}

// Zero-padded ring tables like c(10,20,n)/s(20,n): s[19] is the count.
struct RingTables
{
  int s[20] = { 0 };
  int c[10][20] = { { 0 } };
};

static RingTables makeTables(const std::vector<Ring>& rings)
{
  RingTables t;
  int nr = std::min<int>(static_cast<int>(rings.size()), 19);
  t.s[19] = nr;
  for (int m = 0; m < nr; ++m) {
    const std::vector<int>& mem = rings[m].members;
    t.s[m] = static_cast<int>(mem.size());
    for (size_t l = 0; l < mem.size() && l < 10; ++l)
      t.c[l][m] = mem[l];
  }
  return t;
}

int smallestRingTorsion(const std::vector<Ring>& ringsI,
                        const std::vector<Ring>& ringsJ,
                        const std::vector<Ring>& ringsK,
                        const std::vector<Ring>& ringsL, int i, int j, int k,
                        int l)
{
  RingTables ti = makeTables(ringsI), tj = makeTables(ringsJ),
             tk = makeTables(ringsK), tl = makeTables(ringsL);
  if (ti.s[19] == 0 || tj.s[19] == 0 || tk.s[19] == 0 || tl.s[19] == 0)
    return 0;
  int rings1 = 99, rings2 = 99, rings3 = 99, rings4 = 99;
  const RingTables* ts[4] = { &ti, &tj, &tk, &tl };
  const int others[4][3] = { { j, k, l }, { i, k, l },
                             { i, j, l }, { i, k, j } };
  int* outs[4] = { &rings1, &rings2, &rings3, &rings4 };
  for (int t = 0; t < 4; ++t) {
    for (int m = 0; m < ts[t]->s[19]; ++m) {
      int itest = 0;
      for (int a = 0; a < ts[t]->s[m]; ++a) {
        if (ts[t]->c[a][m] == others[t][0] ||
            ts[t]->c[a][m] == others[t][1] || ts[t]->c[a][m] == others[t][2])
          ++itest;
      }
      if (itest == 3 && ts[t]->s[m] < *outs[t])
        *outs[t] = ts[t]->s[m];
    }
  }
  int rings = std::min(std::min(rings1, rings2), std::min(rings3, rings4));
  if (rings == 99)
    rings = 0;
  return rings;
}

int largestRingTorsion(const std::vector<Ring>& ringsI,
                       const std::vector<Ring>& ringsJ,
                       const std::vector<Ring>& ringsK,
                       const std::vector<Ring>& ringsL, int i, int j, int k,
                       int l)
{
  RingTables ti = makeTables(ringsI), tj = makeTables(ringsJ),
             tk = makeTables(ringsK), tl = makeTables(ringsL);
  if (ti.s[19] == 0 || tj.s[19] == 0 || tk.s[19] == 0 || tl.s[19] == 0)
    return 0;
  int rings1 = -99, rings2 = -99, rings3 = -99, rings4 = -99;
  const RingTables* ts[4] = { &ti, &tj, &tk, &tl };
  const int others[4][3] = { { j, k, l }, { i, k, l },
                             { i, j, l }, { i, k, j } };
  int* outs[4] = { &rings1, &rings2, &rings3, &rings4 };
  for (int t = 0; t < 4; ++t) {
    for (int m = 0; m < ts[t]->s[19]; ++m) {
      int itest = 0;
      for (int a = 0; a < ts[t]->s[m]; ++a) {
        if (ts[t]->c[a][m] == others[t][0] ||
            ts[t]->c[a][m] == others[t][1] || ts[t]->c[a][m] == others[t][2])
          ++itest;
      }
      if (itest == 3 && ts[t]->s[m] > *outs[t])
        *outs[t] = ts[t]->s[m];
    }
  }
  int ringl = std::max(std::max(rings1, rings2), std::max(rings3, rings4));
  if (ringl == -99)
    ringl = 0;
  return ringl;
}

bool alphaCarbonyl(const std::vector<std::vector<int>>& neighbours,
                   const std::vector<int>& numbers,
                   const std::vector<int>& hyb,
                   const std::vector<int>& piFlags, int a, int b)
{
  if (piFlags[a] != 0 && hyb[b] == 3 && numbers[a] == 6 && numbers[b] == 6) {
    int no = 0;
    for (int j : neighbours[a]) {
      if (numbers[j] == 8 && piFlags[j] != 0 &&
          static_cast<int>(neighbours[j].size()) == 1)
        ++no;
    }
    if (no == 1)
      return true;
  }
  if (piFlags[b] != 0 && hyb[a] == 3 && numbers[b] == 6 && numbers[a] == 6) {
    int no = 0;
    for (int j : neighbours[b]) {
      if (numbers[j] == 8 && piFlags[j] != 0 &&
          static_cast<int>(neighbours[j].size()) == 1)
        ++no;
    }
    if (no == 1)
      return true;
  }
  return false;
}

// Integer power mirroring libgcc __powidf2 for exponent 14 (verified
// against gfortran by disassembly + hex comparison): x2, x4, x6, x8, x14.
static double intPow14(double x)
{
  double x2 = x * x;
  double r = x2;
  double x4 = x2 * x2;
  r = x4 * r;
  double x8 = x4 * x4;
  r = x8 * r;
  return r;
}

bool buildTorsions(const std::vector<TorsionBond>& bonds,
                   const std::vector<TorsionAtom>& atoms,
                   const std::vector<std::vector<int>>& neighbours,
                   const std::vector<double>& xyz,
                   const std::vector<std::vector<Ring>>& ringsPerAtom,
                   const std::vector<int>& group,
                   const std::vector<int>& metalTable,
                   const std::vector<double>& tors,
                   const std::vector<double>& tors2, const GffData& param,
                   const GffGenerator& gen, std::vector<Torsion>& torsions,
                   std::vector<TorsionTerm>& terms, Environment& env)
{
  (void)param;
  int n = static_cast<int>(atoms.size());
  if (n <= 0) {
    env.error("empty torsion setup", "buildTorsions");
    return false;
  }
  // Layout adapter for isAmideNitrogen (numctr = 1).
  std::vector<int> nbLists(n * maxNeighbors, 0), nbCounts(n, 0),
    numbers(n), hyb(n), piAtoms(n);
  for (int i = 0; i < n; ++i) {
    numbers[i] = atoms[i].element;
    hyb[i] = atoms[i].hyb;
    piAtoms[i] = atoms[i].pi;
    nbCounts[i] = static_cast<int>(neighbours[i].size());
    for (size_t s = 0; s < neighbours[i].size(); ++s)
      nbLists[i * maxNeighbors + s] = neighbours[i][s];
  }
  torsions.clear();
  terms.clear();
  for (const TorsionBond& bond : bonds) {
    int ii = bond.first, jj = bond.second;
    int nni = static_cast<int>(neighbours[ii].size());
    int nnj = static_cast<int>(neighbours[jj].size());
    int bbtyp = bond.btype;
    if (bbtyp == 3 || bbtyp == 6)
      continue;
    double fij = tors[atoms[ii].element - 1] * tors[atoms[jj].element - 1];
    if (fij < gen.fcThr)
      continue;
    if (tors[atoms[ii].element - 1] < 0.0 ||
        tors[atoms[jj].element - 1] < 0.0)
      continue;
    if (metalTable[atoms[ii].element - 1] > 1 && nni > 4)
      continue;
    if (metalTable[atoms[jj].element - 1] > 1 && nnj > 4)
      continue;
    double fqq =
      1.0 + std::abs(atoms[ii].charge * atoms[jj].charge) * gen.qFacTor;
    int rings = smallestRingBond(ringsPerAtom[ii], ringsPerAtom[jj], ii, jj);
    bool lring = rings > 0;
    bool sp3ij = atoms[ii].hyb == 3 && atoms[jj].hyb == 3;
    int nhi = 1, nhj = 1;
    for (int k : neighbours[ii]) {
      if (atoms[k].element == 1)
        ++nhi;
    }
    for (int l : neighbours[jj]) {
      if (atoms[l].element == 1)
        ++nhj;
    }
    fij = fij * std::pow(static_cast<double>(nhi) * static_cast<double>(nhj),
                         static_cast<double>(0.07f));
    if (alphaCarbonyl(neighbours, numbers, hyb, piAtoms, ii, jj))
      fij = fij * 1.3;
    if (isAmideNitrogen(n, numbers, hyb, nbLists, nbCounts, 1, piAtoms, ii) &&
        atoms[jj].hyb == 3 && atoms[jj].element == 6)
      fij = fij * 1.3;
    if (isAmideNitrogen(n, numbers, hyb, nbLists, nbCounts, 1, piAtoms, jj) &&
        atoms[ii].hyb == 3 && atoms[ii].element == 6)
      fij = fij * 1.3;
    if (bbtyp == 4)
      fij = fij * 0.2;
    for (int kk : neighbours[ii]) {
      if (kk == jj)
        continue;
      for (int ll : neighbours[jj]) {
        if (ll == ii)
          continue;
        if (ll == kk)
          continue;
        if (dihedralNearLinear(xyz, ii, jj, kk, ll))
          continue;
        double fkl =
          tors2[atoms[kk].element - 1] * tors2[atoms[ll].element - 1];
        if (atoms[kk].element == 7 && atoms[kk].pi == 0)
          fkl = fkl * 0.5;
        if (atoms[ll].element == 7 && atoms[ll].pi == 0)
          fkl = fkl * 0.5;
        if (fkl < gen.fcThr)
          continue;
        if (tors[atoms[kk].element - 1] < 0.0 ||
            tors[atoms[ll].element - 1] < 0.0)
          continue;
        double f1 = gen.torsf[1 - 1];
        double f2 = 0.0;
        double cnkk = static_cast<double>(neighbours[kk].size());
        double cnll = static_cast<double>(neighbours[ll].size());
        fkl = fkl * std::pow(cnkk * cnll, static_cast<double>(-0.14f));
        int nrot = 1;
        double phi = 0.0;
        if (lring) {
          int rings4 = 0;
          if (rings > 3) {
            rings4 = smallestRingTorsion(ringsPerAtom[ii], ringsPerAtom[jj],
                                         ringsPerAtom[kk], ringsPerAtom[ll],
                                         ii, jj, kk, ll);
          } else {
            rings4 = 3;
          }
          nrot = 1;
          if (bbtyp == 2)
            nrot = 2;
          phi = 0.0;
          if (bbtyp == 1 && rings4 > 0) {
            int ringl = largestRingTorsion(
              ringsPerAtom[ii], ringsPerAtom[jj], ringsPerAtom[kk],
              ringsPerAtom[ll], ii, jj, kk, ll);
            bool notpicon = atoms[kk].pi == 0 && atoms[ll].pi == 0;
            if (rings4 == 3 && notpicon) {
              nrot = 1;
              phi = 0.0;
              f1 = gen.fr3;
            }
            if (rings4 == 4 && ringl == rings4 && notpicon) {
              nrot = 6;
              phi = 30.0;
              f1 = gen.fr4;
            }
            if (rings4 == 5 && ringl == rings4 && notpicon) {
              nrot = 6;
              phi = 30.0;
              f1 = gen.fr5;
            }
            if (rings4 == 6 && ringl == rings4 && notpicon) {
              nrot = 3;
              phi = 60.0;
              f1 = gen.fr6;
            }
          }
          if (rings4 == 0 && bbtyp == 1 &&
              static_cast<int>(neighbours[kk].size()) == 1 &&
              static_cast<int>(neighbours[ll].size()) == 1) {
            nrot = 6;
            phi = 30.0;
            f1 = 0.30f;
          }
          if (bbtyp == 2 && rings == 5 &&
              atoms[ii].element * atoms[jj].element == 42) {
            if (isAmideNitrogen(n, numbers, hyb, nbLists, nbCounts, 1,
                                piAtoms, ii) ||
                isAmideNitrogen(n, numbers, hyb, nbLists, nbCounts, 1,
                                piAtoms, jj))
              f1 = 5.0f;
          }
        } else {
          phi = 180.0;
          nrot = 1;
          if (atoms[ii].hyb == 3 && atoms[jj].hyb == 3)
            nrot = 3;
          if (bbtyp == 2)
            nrot = 2;
          if (atoms[ii].pi > 0 &&
              (atoms[jj].pi == 0 && atoms[jj].hyb == 3)) {
            f1 = 0.5;
            if (atoms[ii].element == 7)
              f1 = 0.2;
            phi = 180.0;
            nrot = 3;
          }
          if (atoms[jj].pi > 0 &&
              (atoms[ii].pi == 0 && atoms[ii].hyb == 3)) {
            f1 = 0.5;
            if (atoms[jj].element == 7)
              f1 = 0.2;
            phi = 180.0;
            nrot = 3;
          }
        }
        if (sp3ij) {
          if (group[atoms[ii].element - 1] == 5 &&
              group[atoms[jj].element - 1] == 5) {
            nrot = 3;
            phi = 60.0;
            f1 = 3.0;
          }
          if ((group[atoms[ii].element - 1] == 5 &&
               group[atoms[jj].element - 1] == 6) ||
              (group[atoms[ii].element - 1] == 6 &&
               group[atoms[jj].element - 1] == 5)) {
            nrot = 2;
            phi = 90.0;
            f1 = 1.0;
            if (atoms[ii].element >= 15 && atoms[jj].element >= 15)
              f1 = 20.0;
          }
          if (group[atoms[ii].element - 1] == 6 &&
              group[atoms[jj].element - 1] == 6) {
            nrot = 2;
            phi = 90.0;
            f1 = 5.0;
            if (atoms[ii].element >= 16 && atoms[jj].element >= 16)
              f1 = 25.0;
          }
        }
        if (bond.pibo > 0) {
          double dum = 1.24 - bond.pibo;
          f2 = bond.pibo * std::exp(-2.5 * intPow14(dum));
          if (atoms[kk].pi == 0 && atoms[kk].element > 10)
            f2 = f2 * 1.3f;
          if (atoms[ll].pi == 0 && atoms[ll].element > 10)
            f2 = f2 * 1.3f;
          f1 = f1 * 0.55f;
        }
        if (atoms[kk].hyb == 5 || atoms[ll].hyb == 5)
          fkl = fkl * 1.5f;
        double fctot =
          (f1 + 10.0 * gen.torsf[2 - 1] * f2) * fqq * fij * fkl;
        if (fctot > gen.fcThr) {
          Torsion t;
          t.outer1 = ll;
          t.center1 = ii;
          t.center2 = jj;
          t.outer2 = kk;
          t.multiplicity = nrot;
          torsions.push_back(t);
          TorsionTerm term;
          term.phase = phi * pi / 180.0;
          term.forceConstant = fctot;
          terms.push_back(term);
        }
        bool sp3kl = atoms[kk].hyb == 3 && atoms[ll].hyb == 3;
        if (sp3kl && sp3ij && !lring && bbtyp < 5) {
          double ff = gen.torsf[6 - 1];
          if (atoms[ii].element == 7 || atoms[jj].element == 7)
            ff = gen.torsf[7 - 1];
          if (atoms[ii].element == 8 || atoms[jj].element == 8)
            ff = gen.torsf[8 - 1];
          Torsion t;
          t.outer1 = ll;
          t.center1 = ii;
          t.center2 = jj;
          t.outer2 = kk;
          t.multiplicity = 1;
          torsions.push_back(t);
          TorsionTerm term;
          term.phase = pi;
          term.forceConstant = ff * fij * fkl * fqq;
          terms.push_back(term);
        }
      }
    }
  }
  return true;
}

} // namespace Xtb
} // namespace Avogadro
