/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_ini2.F90 (gfnff_neigh hybridization section,
  nn_nearest_noM) and src/gfnff/neighbor.f90 (nbLoc),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#include "gfnffhyb.h"

#include "constants.h"
#include "environment.h"
#include "gfnfftopo.h"

#include <cmath>
#include <string>

namespace Avogadro {
namespace Xtb {

std::vector<NeighborCell> locateNeighbors(const std::vector<int>& lists,
                                          const std::vector<int>& counts,
                                          int n, int numctr, int indexi)
{
  std::vector<NeighborCell> cols;
  for (int cell = 0; cell < numctr; ++cell) {
    if (counts[cell * n + indexi] == 0)
      continue;
    NeighborCell col;
    col.cell = cell;
    // Only the first maxNeighbors - 2 slots are scanned like the reference.
    for (int k = 0; k < maxNeighbors - 2; ++k) {
      int m = lists[(cell * n + indexi) * maxNeighbors + k];
      if (m != 0 || k < counts[cell * n + indexi]) {
        // NOTE: lists store 0-based atoms, so slot content cannot double as
        // an emptiness flag; membership is driven by the counts instead.
        // Entries beyond the count are ignored like the reference ignores
        // zero slots.
        if (k >= counts[cell * n + indexi])
          break;
        col.members.push_back(m);
      }
    }
    col.count = static_cast<int>(col.members.size());
    cols.push_back(col);
  }
  return cols;
}

int nearestNonMetalCoord(int ii, const std::vector<int>& numbers,
                         const std::vector<double>& xyz,
                         const std::vector<int>& lists,
                         const std::vector<int>& counts, int numctr,
                         const int* metal)
{
  int best = -1;
  double rmin = 1.0e42;
  for (int cell = 0; cell < numctr; ++cell) {
    int nn = counts[cell * numbers.size() + ii];
    for (int k = 0; k < nn; ++k) {
      int jj = lists[(cell * numbers.size() + ii) * maxNeighbors + k];
      if (metal[numbers[jj] - 1] != 0)
        continue;
      double dx = xyz[3 * ii] - xyz[3 * jj];
      double dy = xyz[3 * ii + 1] - xyz[3 * jj + 1];
      double dz = xyz[3 * ii + 2] - xyz[3 * jj + 2];
      double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (dist < rmin) {
        rmin = dist;
        best = jj;
      }
    }
  }
  if (best < 0)
    return 0;
  int total = 0;
  for (int cell = 0; cell < numctr; ++cell)
    total += counts[cell * numbers.size() + best];
  return total;
}

double bondAngle(const std::vector<double>& xyz, int a, int b, int c)
{
  double bax = xyz[3 * a] - xyz[3 * b];
  double bay = xyz[3 * a + 1] - xyz[3 * b + 1];
  double baz = xyz[3 * a + 2] - xyz[3 * b + 2];
  double bcx = xyz[3 * c] - xyz[3 * b];
  double bcy = xyz[3 * c + 1] - xyz[3 * b + 1];
  double bcz = xyz[3 * c + 2] - xyz[3 * b + 2];
  double na = std::sqrt(bax * bax + bay * bay + baz * baz);
  double nc = std::sqrt(bcx * bcx + bcy * bcy + bcz * bcz);
  double cosang = (bax * bcx + bay * bcy + baz * bcz) / (na * nc);
  if (cosang > 1.0)
    cosang = 1.0;
  if (cosang < -1.0)
    cosang = -1.0;
  return std::acos(cosang);
}

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
                         std::vector<int>& itag, Environment& env)
{
  hyb.assign(n, 0);
  itag.assign(n, 0);
  auto fullCount = [&](int i) {
    int s = 0;
    for (int cell = 0; cell < numctr; ++cell)
      s += nbfCounts[cell * n + i];
    return s;
  };
  auto plainCount = [&](int i) {
    int s = 0;
    for (int cell = 0; cell < numctr; ++cell)
      s += nbCounts[cell * n + i];
    return s;
  };
  auto nometalCount = [&](int i) {
    int s = 0;
    for (int cell = 0; cell < numctr; ++cell)
      s += nbmCounts[cell * n + i];
    return s;
  };
  for (int i = 0; i < n; ++i) {
    int ati = numbers[i];
    int numnbm = nometalCount(i);
    int numnbf = fullCount(i);
    int numnb = plainCount(i);
    int nbdiff = numnbf - numnb;
    int nbmdiff = numnbf - numnbm;
    // Reduced (metal-free) working list for eta-coordinated atoms.
    bool etacoord = false;
    if (ati <= 10) {
      if (ati == 6 && numnbf >= 4 && numnbm == 3)
        etacoord = true; // CP case
      if (ati == 6 && numnbf == 3 && numnbm == 2)
        etacoord = true; // alkyne case
      int nm = 0, im = -1;
      for (int cell = 0; cell < numctr; ++cell) {
        for (int k = 0; k < nbfCounts[cell * n + i]; ++k) {
          int kk = nbfLists[(cell * n + i) * maxNeighbors + k];
          if (metal[numbers[kk] - 1] != 0) {
            ++nm;
            im = kk;
          }
        }
      }
      if (nm == 0) {
        etacoord = false;
      } else if (nm == 1) {
        int ncm = 0;
        for (int cell = 0; cell < numctr; ++cell) {
          for (int k = 0; k < nbfCounts[cell * n + i]; ++k) {
            int kk = nbfLists[(cell * n + i) * maxNeighbors + k];
            if (kk == im)
              continue;
            for (int c2 = 0; c2 < numctr; ++c2) {
              for (int l = 0; l < nbfCounts[c2 * n + kk]; ++l) {
                if (nbfLists[(c2 * n + kk) * maxNeighbors + l] == im)
                  ++ncm;
              }
            }
          }
        }
        if (ncm == 0)
          etacoord = false;
      }
    }
    // nbdum member access: reduced list when eta-coordinated.
    auto dumCount = [&](int a) {
      int s = 0;
      for (int cell = 0; cell < numctr; ++cell) {
        s += etacoord ? nbmCounts[cell * n + a] : nbfCounts[cell * n + a];
      }
      return s;
    };
    int nb20i = dumCount(i);
    int nh = 0, no = 0;
    for (int cell = 0; cell < numctr; ++cell) {
      int cnt = etacoord ? nbmCounts[cell * n + i] : nbfCounts[cell * n + i];
      for (int k = 0; k < cnt; ++k) {
        int m = (etacoord ? nbmLists : nbfLists)[(cell * n + i) * maxNeighbors + k];
        if (numbers[m] == 1)
          ++nh;
        if (numbers[m] == 8)
          ++no;
      }
    }
    if (group[ati - 1] == 1) { // H
      if (nb20i == 2)
        hyb[i] = 1; // bridging H
      if (nb20i > 2)
        hyb[i] = 3; // M+ tetra coord
      if (nb20i > 4)
        hyb[i] = 0; // M+ HC
    }
    if (group[ati - 1] == 2) { // Be
      if (nb20i == 2)
        hyb[i] = 1;
      if (nb20i > 2)
        hyb[i] = 3;
      if (nb20i > 4)
        hyb[i] = 0;
    }
    if (group[ati - 1] == 3) { // B
      if (nb20i > 4)
        hyb[i] = 3;
      if (nb20i > 4 && ati > 10 && nbdiff == 0)
        hyb[i] = 5;
      if (nb20i == 4)
        hyb[i] = 3;
      if (nb20i == 3)
        hyb[i] = 2;
      if (nb20i == 2)
        hyb[i] = 1;
    }
    if (group[ati - 1] == 4) { // C
      if (nb20i >= 4)
        hyb[i] = 3;
      if (nb20i > 4 && ati > 10 && nbdiff == 0)
        hyb[i] = 5;
      if (nb20i == 3)
        hyb[i] = 2;
      if (nb20i == 2) {
        // Locate the two neighbours like nbLoc does.
        std::vector<int> members, cells;
        for (int cell = 0; cell < numctr; ++cell) {
          int cnt = etacoord ? nbmCounts[cell * n + i] : nbfCounts[cell * n + i];
          for (int k = 0; k < cnt && k < maxNeighbors - 2; ++k) {
            members.push_back(
              (etacoord ? nbmLists : nbfLists)[(cell * n + i) * maxNeighbors + k]);
            cells.push_back(cell);
          }
        }
        if (members.size() != 2) {
          env.error("Hybridization failed. Neighbors could not be located.",
                    "assignHybridization");
          return false;
        }
        double phi = bondAngle(xyz, members[0], i, members[1]);
        if (phi * 180.0 / pi < 150.0) {
          hyb[i] = 2; // otherwise carbenes will not be recognized
          itag[i] = 1;
        } else {
          hyb[i] = 1; // linear triple bond etc
        }
        if (qa[i] < -0.4) {
          hyb[i] = 2;
          itag[i] = 0;
        }
      }
      if (nb20i == 1)
        hyb[i] = 1; // CO
    }
    if (group[ati - 1] == 5) { // N
      if (nb20i >= 4)
        hyb[i] = 3;
      if (nb20i > 4 && ati > 10 && nbdiff == 0)
        hyb[i] = 5;
      if (nb20i == 3)
        hyb[i] = 3;
      if (nb20i == 3 && ati == 7) {
        int kk = 0, ll = 0, nn = 0;
        for (int cell = 0; cell < numctr; ++cell) {
          int cnt = etacoord ? nbmCounts[cell * n + i] : nbfCounts[cell * n + i];
          for (int k = 0; k < 3 && k < cnt; ++k) {
            int jj = (etacoord ? nbmLists : nbfLists)[(cell * n + i) * maxNeighbors + k];
            if (jj < 0)
              break;
            if (numbers[jj] == 8 && plainCount(jj) == 1)
              ++kk; // NO2 or R2-N=O
            if (numbers[jj] == 5 && plainCount(jj) == 4)
              ++ll; // B-N
            if (numbers[jj] == 16 && plainCount(jj) == 4)
              ++nn; // N-SO2-
          }
        }
        if (nn == 1 && ll == 0 && kk == 0)
          hyb[i] = 3;
        if (ll == 1 && nn == 0)
          hyb[i] = 2;
        if (kk >= 1) {
          hyb[i] = 2;
          itag[i] = 1;
        }
        if (nbmdiff > 0 && nn == 0)
          hyb[i] = 2; // pyridine N
      }
      if (nb20i == 2) {
        hyb[i] = 2;
        std::vector<int> members, cells;
        for (int cell = 0; cell < numctr; ++cell) {
          int cnt = etacoord ? nbmCounts[cell * n + i] : nbfCounts[cell * n + i];
          for (int k = 0; k < cnt && k < maxNeighbors - 2; ++k) {
            members.push_back(
              (etacoord ? nbmLists : nbfLists)[(cell * n + i) * maxNeighbors + k]);
            cells.push_back(cell);
          }
        }
        if (members.size() != 2) {
          env.error("Hybridization failed. Neighbors could not be located.",
                    "assignHybridization");
          return false;
        }
        int jd = members[0], kd = members[1];
        double phi = bondAngle(xyz, jd, i, kd);
        if (plainCount(jd) == 1 && numbers[jd] == 6)
          hyb[i] = 1; // R-N=C
        if (plainCount(kd) == 1 && numbers[kd] == 6)
          hyb[i] = 1;
        if (plainCount(jd) == 1 && numbers[jd] == 7)
          hyb[i] = 1; // R-N=N
        if (plainCount(kd) == 1 && numbers[kd] == 7)
          hyb[i] = 1;
        if (jd >= 0 && metal[numbers[jd] - 1] > 0)
          hyb[i] = 1; // M-NC-R
        if (kd >= 0 && metal[numbers[kd] - 1] > 0)
          hyb[i] = 1;
        if (numbers[jd] == 7 && numbers[kd] == 7 && plainCount(jd) <= 2 &&
            plainCount(kd) <= 2)
          hyb[i] = 1; // N=N=N
        if (phi * 180.0 / pi > linThr)
          hyb[i] = 1; // geometry dependent
      }
      if (nb20i == 1)
        hyb[i] = 1;
    }
    if (group[ati - 1] == 6) { // O
      if (nb20i >= 3)
        hyb[i] = 3;
      if (nb20i > 3 && ati > 10 && nbdiff == 0)
        hyb[i] = 5;
      if (nb20i == 2)
        hyb[i] = 3;
      if (nb20i == 2 && nbmdiff > 0) {
        int j = nearestNonMetalCoord(i, numbers, xyz, nbLists, nbCounts,
                                     numctr, metal);
        if (j == 3)
          hyb[i] = 2; // M-O-X conjugated
        if (j == 4)
          hyb[i] = 3; // M-O-X not conjugated
      }
      if (nb20i == 1)
        hyb[i] = 2;
      if (nb20i == 1 && nbdiff == 0) {
        auto cols = locateNeighbors(nbLists, nbCounts, n, numctr, i);
        if (!cols.empty()) {
          int first = cols[0].members.empty() ? -1 : cols[0].members[0];
          if (first >= 0 && plainCount(first) == 1)
            hyb[i] = 1; // CO
        }
      }
    }
    if (group[ati - 1] == 7) { // F
      if (nb20i == 2)
        hyb[i] = 1;
      if (nb20i > 2 && ati > 10)
        hyb[i] = 5;
    }
    if (group[ati - 1] == 8) { // Ne
      hyb[i] = 0;
      if (nb20i > 0 && ati > 2)
        hyb[i] = 5;
    }
    if (group[ati - 1] <= 0) { // TMs
      int nni = nb20i;
      if (nh != 0 && nh != nni)
        nni -= nh; // don't count Hs
      if (nni <= 2)
        hyb[i] = 1;
      if (nni <= 2 && group[ati - 1] <= -6)
        hyb[i] = 2;
      if (nni == 3)
        hyb[i] = 2;
      if (nni == 4 && group[ati - 1] > -7)
        hyb[i] = 3;
      if (nni == 4 && group[ati - 1] <= -7)
        hyb[i] = 3;
      if (nni == 5 && group[ati - 1] == -3)
        hyb[i] = 3;
    }
    (void)no;
  }

  // Arine special: two bonded carbene carbons lose their tag.
  for (int i = 0; i < n; ++i) {
    int numnb = plainCount(i);
    if (numnb > 12)
      continue;
    for (int cell = 0; cell < numctr; ++cell) {
      for (int k = 0; k < nbCounts[cell * n + i]; ++k) {
        int kk = nbLists[(cell * n + i) * maxNeighbors + k];
        if (numbers[kk] == 6 && numbers[i] == 6 && itag[i] == 1 &&
            itag[kk] == 1) {
          itag[i] = 0;
          itag[kk] = 0;
        }
      }
    }
  }
  int crowded = 0;
  for (int i = 0; i < n; ++i) {
    if (plainCount(i) > 12)
      ++crowded;
  }
  if (static_cast<double>(crowded) / n > 0.3 && nbCall) {
    env.error(" too many atoms with extreme high CN", "assignHybridization");
    return false;
  }
  return true;
}

} // namespace Xtb
} // namespace Avogadro
