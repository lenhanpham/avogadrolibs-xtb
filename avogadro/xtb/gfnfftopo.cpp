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

#include "gfnfftopo.h"

#include "eeq.h"
#include "environment.h"

#include <cmath>
#include <string>

namespace Avogadro {
namespace Xtb {

void fillNeighborList(int n, const std::vector<int>& numbers,
                      const std::vector<double>& radPacked,
                      const std::vector<double>& dist, const double* mchar,
                      int icase, double f, double f2, const int* metal,
                      const int* group, const int* normCn,
                      const int* fullCounts, int numctr,
                      std::vector<int>& list, std::vector<int>& counts)
{
  list.assign(maxNeighbors * n * numctr, 0);
  counts.assign(n * numctr, 0);
  std::vector<char> tags(n * n * numctr, 0);
  auto tagAt = [&](int j, int i, int cell) -> char& {
    return tags[(cell * n + i) * n + j];
  };

  for (int cell = 0; cell < numctr; ++cell) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        int k = packedIndex(j, i);
        double d = dist[(cell * n + i) * n + j];
        if (d <= 0.0)
          continue;
        double fm = 1.0;
        if (icase == 1) {
          if (metal[numbers[i] - 1] == 2)
            fm *= f2;
          if (metal[numbers[j] - 1] == 2)
            fm *= f2;
          if (metal[numbers[i] - 1] == 1)
            fm *= f2 + 0.025;
          if (metal[numbers[j] - 1] == 1)
            fm *= f2 + 0.025;
        }
        if (icase == 2) {
          int nnfi = fullCounts[i];
          int nnfj = fullCounts[j];
          int crit = 6;
          if (group[numbers[i] - 1] <= 2)
            crit = 4;
          if (nnfi > crit)
            continue;
          crit = 6;
          if (group[numbers[j] - 1] <= 2)
            crit = 4;
          if (nnfj > crit)
            continue;
        }
        if (icase == 3) {
          if (mchar[i] > 0.25 || metal[numbers[i] - 1] > 0)
            continue;
          if (mchar[j] > 0.25 || metal[numbers[j] - 1] > 0)
            continue;
          if (fullCounts[i] > normCn[numbers[i] - 1] && numbers[i] > 10)
            continue;
          if (fullCounts[j] > normCn[numbers[j] - 1] && numbers[j] > 10)
            continue;
        }
        double rco = radPacked[k];
        if (d < fm * f * rco)
          tagAt(j, i, cell) = 1;
      }
    }
  }

  for (int cell = 0; cell < numctr; ++cell) {
    for (int i = 0; i < n; ++i) {
      int nn = 0;
      for (int j = 0; j < n; ++j) {
        if (tagAt(j, i, cell) == 1 && i != j) {
          list[(cell * n + i) * maxNeighbors + nn] = j;
          ++nn;
        }
      }
      counts[cell * n + i] = nn;
    }
  }
}

void estimateBondLengths(int n, const std::vector<int>& numbers,
                         const std::vector<int>& nbCounts,
                         const std::vector<int>& nbList,
                         const std::vector<double>& radAngstrom,
                         double rfgoed1, double tdistThr,
                         std::vector<float>& rabdOut,
                         std::vector<double>& rtmpPacked)
{
  constexpr float cutoff = 13.0f;
  rabdOut.assign(n * n, cutoff);
  auto rabd = [&](int i, int j) -> float& { return rabdOut[i * n + j]; };
  for (int i = 0; i < n; ++i) {
    rabd(i, i) = 0.0f;
    // Single cell (numctr = 1) like the molecular path in gfnff_ini.
    for (int k = 0; k < nbCounts[i]; ++k) {
      int j = nbList[i * maxNeighbors + k];
      float r = static_cast<float>(radAngstrom[numbers[i] - 1] +
                                   radAngstrom[numbers[j] - 1]);
      rabd(j, i) = r;
      rabd(i, j) = r;
    }
  }
  float tthr = static_cast<float>(tdistThr);
  for (int k = 0; k < n; ++k) {
    for (int i = 0; i < n; ++i) {
      if (rabd(i, k) > tthr)
        continue;
      for (int j = 0; j < n; ++j) {
        if (rabd(k, j) > tthr)
          continue;
        float via = rabd(i, k) + rabd(k, j);
        if (rabd(i, j) > via)
          rabd(i, j) = via;
      }
    }
  }
  int npair = n * (n + 1) / 2;
  rtmpPacked.assign(npair, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      float r = rabd(j, i);
      if (r > tthr)
        r = cutoff;
      rtmpPacked[packedIndex(j, i)] =
        rfgoed1 * static_cast<double>(r) / 0.52917726;
    }
  }
}

bool topologyCharges(int n, const std::vector<double>& pairPacked,
                     const std::vector<double>& chieeq,
                     const std::vector<double>& gameeq,
                     const std::vector<double>& alpeeqSq, int nfrag,
                     const std::vector<int>& fragOfAtom,
                     const std::vector<double>& fragCharges, Environment& env,
                     std::vector<double>& charges, double& energy)
{
  const char* source = "topologyCharges";
  // tsqrt2pi literal from goedeckera (0.797884560802866, 15 digits).
  constexpr double tsqrt2pi = 0.797884560802866;
  int m = n + nfrag;
  std::vector<double> amat(m * m, 0.0), xvec(m, 0.0);
  auto amatAt = [&](int i, int j) -> double& { return amat[i * m + j]; };
  for (int i = 0; i < n; ++i)
    xvec[i] = chieeq[i];
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      double rij = pairPacked[packedIndex(j, i)];
      double gammij = 1.0 / std::sqrt(alpeeqSq[i] + alpeeqSq[j]);
      double tmp = std::erf(gammij * rij) / rij;
      amatAt(j, i) = tmp;
      amatAt(i, j) = tmp;
    }
    amatAt(i, i) = gameeq[i] + tsqrt2pi / std::sqrt(alpeeqSq[i]);
  }
  for (int f = 0; f < nfrag; ++f) {
    xvec[n + f] = fragCharges[f];
    for (int j = 0; j < n; ++j) {
      if (fragOfAtom[j] == f + 1) {
        amatAt(n + f, j) = 1.0;
        amatAt(j, n + f) = 1.0;
      }
    }
  }
  if (!solveSymmetric(amat, xvec, m)) {
    env.error("Solving linear equations failed", source);
    return false;
  }
  charges.assign(xvec.begin(), xvec.begin() + n);
  if (n == 1)
    charges[0] = fragCharges[0];

  energy = 0.0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      double rij = pairPacked[packedIndex(j, i)];
      double gammij = 1.0 / std::sqrt(alpeeqSq[i] + alpeeqSq[j]);
      double tmp = std::erf(gammij * rij) / rij;
      energy += charges[i] * charges[j] * tmp / rij;
    }
    energy += -charges[i] * chieeq[i] +
              charges[i] * charges[i] * 0.5 *
                (gameeq[i] + tsqrt2pi / std::sqrt(alpeeqSq[i]));
  }
  return true;
}

} // namespace Xtb
} // namespace Avogadro
