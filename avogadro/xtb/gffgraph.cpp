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

#include "gffgraph.h"

#include "environment.h"
#include "gfnfftopo.h"

#include <cmath>
#include <functional>

namespace Avogadro {
namespace Xtb {

int findFragments(const std::vector<std::vector<int>>& adjacency,
                  std::vector<int>& fragOfAtom)
{
  int n = static_cast<int>(adjacency.size());
  // Symmetric bond marks, mirroring mrecgff.
  std::vector<std::vector<char>> bond(n, std::vector<char>(n, 0));
  for (int i = 0; i < n; ++i) {
    for (int j : adjacency[i]) {
      if (j < 0 || j >= n)
        continue;
      bond[i][j] = 1;
      bond[j][i] = 1;
    }
  }
  fragOfAtom.assign(n, 0);
  std::vector<char> taken(n, 0);
  int count = 0;
  std::function<void(int)> visit = [&](int i) {
    // Take and clear max entries like mrecgff2: scan j ascending and clear
    // bond(j, i) as encountered, skipping self and taken atoms.
    for (int k = 0; k < n; ++k) {
      int j = -1;
      for (int jj = 0; jj < n; ++jj) {
        if (bond[jj][i]) {
          j = jj;
          break;
        }
      }
      if (j < 0)
        break;
      bond[j][i] = 0;
      if (i == j)
        continue;
      if (!taken[j]) {
        fragOfAtom[j] = count;
        taken[j] = 1;
        visit(j);
      }
    }
  };
  for (int i = 0; i < n; ++i) {
    if (!taken[i]) {
      ++count;
      fragOfAtom[i] = count;
      taken[i] = 1;
      visit(i);
    }
  }
  return count;
}

int findFragmentsPbc(
  const std::vector<std::vector<std::pair<int, int>>>& neighbours,
  int numCells, std::vector<int>& fragOfAtom)
{
  int n = static_cast<int>(neighbours.size());
  // Directed (neighbour, cell) marks, mirroring mrecgffPBC.
  std::vector<std::vector<std::vector<char>>> bond(
    n, std::vector<std::vector<char>>(n, std::vector<char>(numCells, 0)));
  for (int i = 0; i < n; ++i) {
    for (const auto& nb : neighbours[i]) {
      int j = nb.first, cell = nb.second;
      if (j < 0 || j >= n || cell < 0 || cell >= numCells)
        continue;
      bond[j][i][cell] = 1;
    }
  }
  fragOfAtom.assign(n, 0);
  std::vector<char> taken(n, 0);
  int count = 0;
  std::function<void(int)> visit = [&](int i) {
    // Count of list entries like sum(nb(numnb,i,:)); take marks j-fastest
    // like maxloc over bond(:, i, :).
    int total = 0;
    for (const auto& nb : neighbours[i])
      (void)nb, ++total;
    for (int k = 0; k < total; ++k) {
      int j = -1, jc = -1;
      for (int c = 0; c < numCells && j < 0; ++c) {
        for (int jj = 0; jj < n; ++jj) {
          if (bond[jj][i][c]) {
            j = jj;
            jc = c;
            break;
          }
        }
      }
      if (j < 0)
        break;
      bond[j][i][jc] = 0;
      if (i == j && jc == 0)
        continue;
      if (!taken[j]) {
        fragOfAtom[j] = count;
        taken[j] = 1;
        visit(j);
      }
    }
  };
  for (int i = 0; i < n; ++i) {
    if (!taken[i]) {
      ++count;
      fragOfAtom[i] = count;
      taken[i] = 1;
      visit(i);
    }
  }
  return count;
}

namespace {

// Selection sort ascending with index tracking, mirroring ssort.
void selectionSort(std::vector<double>& values, std::vector<int>& index)
{
  int n = static_cast<int>(values.size());
  for (int ii = 1; ii < n; ++ii) {
    int i = ii - 1;
    int k = i;
    double pp = values[i];
    for (int j = ii; j < n; ++j) {
      if (values[j] > pp)
        continue;
      k = j;
      pp = values[j];
    }
    if (k == i)
      continue;
    values[k] = values[i];
    values[i] = pp;
    int sc = index[i];
    index[i] = index[k];
    index[k] = sc;
  }
}

// True when all n members are distinct, mirroring chkrng.
bool distinctMembers(const std::vector<int>& members, int n)
{
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (members[i] == members[j])
        return false;
    }
  }
  return true;
}

} // namespace

bool findRingsThrough(int n, const std::vector<int>& numbers,
                      const std::vector<std::vector<int>>& neighbours, int a0,
                      std::vector<Ring>& rings, Environment& /*env*/)
{
  if (n <= 2)
    return true;
  // Flag slot check mirroring the early return (slot numnb-1 holds a marker
  // in xtb's lists; a real 41st neighbour pointing at atom 0 would also
  // match, exactly like the reference).
  if (static_cast<int>(neighbours[a0].size()) > 40 &&
      neighbours[a0][40] == 0) {
    return true;
  }
  // Candidate ring storage mirroring cdum/idum exactly: columns hold
  // 1-based atom numbers with 0 for empty slots (atoms are 1-based in
  // the reference); only slots 1:iring are stored, the rest stay zero,
  // so dedup compares member sets (plus shared zeros).
  struct Candidate
  {
    std::vector<int> column = std::vector<int>(10, 0);
    int size = 0;
  };
  std::vector<Candidate> found;
  int nn = static_cast<int>(neighbours[a0].size());
  // Shared path array like the reference's c, reset to empty at each
  // branch start (c = 0 inside the m loop).
  std::vector<int> path(10, 0);
  // Each branch m reorders a0's neighbours with branch m first.
  for (int m = 0; m < nn; ++m) {
    // Leaves pruned like the reference (count set to 0).
    std::vector<std::vector<int>> work(n);
    for (int i = 0; i < n; ++i) {
      work[i] = neighbours[i];
      if (static_cast<int>(work[i].size()) == 1)
        work[i].clear();
    }
    std::vector<double> w(nn);
    std::vector<int> list(nn);
    for (int mm = 0; mm < nn; ++mm) {
      w[mm] = mm + 1;
      list[mm] = mm;
    }
    w[m] = 0.0;
    selectionSort(w, list);
    std::vector<int> order(nn);
    for (int mm = 0; mm < nn; ++mm)
      order[mm] = neighbours[a0][list[mm]];
    // Depth-first walk mirroring the nested i1..i6 loops (rings of 3 to 6
    // members): path[depth] is the current node (1-based); closure to a0
    // recorded from depth 1 to 4. Branches whose m-th neighbour (in the
    // ORIGINAL order) is atom 0 are skipped like the reference
    // (nb(m,a0) == 1 in 1-based numbering).
    if (neighbours[a0][m] == 0)
      continue;
    std::fill(path.begin(), path.end(), 0);
    std::function<void(int, int)> dfs = [&](int depth, int current) {
      for (int next : work[current]) {
        if (next == current)
          continue;
        if (next == a0) {
          // Closure attempt: distinctness covers the prefix plus the
          // closing node, like chkrng over c(1..iring).
          bool distinct = distinctMembers(path, depth + 1);
          for (int d = 0; d <= depth && distinct; ++d) {
            if (path[d] == a0 + 1)
              distinct = false;
          }
          if (depth >= 1 && depth <= 4 && distinct) {
            Candidate c;
            // Store only slots 1:iring like cdum(1:iring,kk); the rest
            // stays zero, so dedup sees member sets (plus shared zeros).
            for (int d = 0; d <= depth; ++d)
              c.column[d] = path[d];
            // The closing node is stored explicitly like c(iring) = a0.
            c.column[depth + 1] = a0 + 1;
            c.size = depth + 2;
            found.push_back(c);
            if (static_cast<int>(found.size()) >= 500)
              return;
          }
          // The reference descends through a0 freely (only chkrng
          // rejects repeats at record time); replicated for identical
          // control flow.
          if (depth + 1 <= 4) {
            path[depth + 1] = a0 + 1;
            dfs(depth + 1, a0);
            if (static_cast<int>(found.size()) >= 500)
              return;
          }
          continue;
        }
        if (depth + 1 <= 4) {
          path[depth + 1] = next + 1;
          dfs(depth + 1, next);
          if (static_cast<int>(found.size()) >= 500)
            return;
        }
      }
    };
    for (int a1 : order) {
      if (a1 == a0)
        continue;
      path[0] = a1 + 1;
      dfs(0, a1);
    }
  }
  // Deduplicate by atom-set equality over the stored columns,
  // mirroring the reference (1-based values with 0 for empty slots).
  std::vector<char> same(found.size(), 0);
  for (size_t i = 0; i < found.size(); ++i) {
    for (size_t j = i + 1; j < found.size(); ++j) {
      if (found[i].size != found[j].size)
        continue;
      if (same[j])
        continue;
      std::vector<char> inI(n + 1, 0), inJ(n + 1, 0);
      for (int d = 0; d < 10; ++d) {
        int a = found[i].column[d], b = found[j].column[d];
        if (a >= 0 && a <= n)
          inI[a] = 1;
        if (b >= 0 && b <= n)
          inJ[b] = 1;
      }
      if (inI == inJ)
        same[j] = 1;
    }
  }
  rings.clear();
  for (size_t i = 0; i < found.size(); ++i) {
    if (same[i])
      continue;
    // At most 19 rings kept, mirroring the m > 19 truncation.
    if (static_cast<int>(rings.size()) >= 19)
      break;
    Ring ring;
    for (int d = 0; d < found[i].size; ++d)
      ring.members.push_back(found[i].column[d] - 1);
    // Hetero code, mirroring the cout(m,19) computation.
    double sum = 0.0;
    for (int a : ring.members)
      sum += numbers[a];
    double av = sum / ring.members.size();
    double sd = 0.0;
    for (int a : ring.members) {
      double d = av - numbers[a];
      sd += d * d;
    }
    if (sd > 1.0e-6) {
      ring.hetero = static_cast<int>(1000.0 * std::sqrt(sd) / ring.members.size());
    }
    rings.push_back(ring);
  }
  return true;
}

int smallestRingThrough(const std::vector<Ring>& rings)
{
  int best = 99;
  for (const Ring& ring : rings) {
    if (static_cast<int>(ring.members.size()) < best)
      best = static_cast<int>(ring.members.size());
  }
  return best;
}

int smallestRingBond(const std::vector<Ring>& ringsOfI,
                     const std::vector<Ring>& ringsOfJ, int i, int j)
{
  int best = 99;
  auto contains = [](const Ring& ring, int a) {
    for (int m : ring.members) {
      if (m == a)
        return true;
    }
    return false;
  };
  for (const Ring& ring : ringsOfI) {
    if (contains(ring, j))
      best = std::min(best, static_cast<int>(ring.members.size()));
  }
  for (const Ring& ring : ringsOfJ) {
    if (contains(ring, i))
      best = std::min(best, static_cast<int>(ring.members.size()));
  }
  if (best == 99)
    return 0;
  return best;
}

void gffCoordinationNumber(int n, const std::vector<int>& numbers,
                           const std::vector<double>& xyz,
                           const std::vector<double>& rabPacked,
                           const std::vector<double>& rcov, double cnmax,
                           double thr, std::vector<double>& logCn,
                           std::vector<double>& dlogCn)
{
  constexpr double kn = -7.5;
  constexpr double sqrtpi = 1.77245385091; // 11-digit literal from the source
  auto logCut = [&](double c) {
    return std::log(1.0 + std::exp(cnmax)) - std::log(1.0 + std::exp(cnmax - c));
  };
  auto dLogCut = [&](double c) {
    return std::exp(cnmax) / (std::exp(cnmax) + std::exp(c));
  };
  std::vector<double> cn(n, 0.0);
  for (int i = 1; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      double r = rabPacked[packedIndex(j, i)];
      if (r > thr)
        continue;
      double r0 = rcov[numbers[i] - 1] + rcov[numbers[j] - 1];
      double dr = (r - r0) / r0;
      double erfCn = 0.5 * (1.0 + std::erf(kn * dr));
      cn[i] += erfCn;
      cn[j] += erfCn;
    }
  }
  logCn.assign(n, 0.0);
  dlogCn.assign(3 * n * n, 0.0);
  auto dlog = [&](int i, int j, int c) -> double& {
    return dlogCn[(i * n + j) * 3 + c];
  };
  for (int i = 0; i < n; ++i) {
    logCn[i] = logCut(cn[i]);
    double dlogdcni = dLogCut(cn[i]);
    for (int j = 0; j < i; ++j) {
      double dlogdcnj = dLogCut(cn[j]);
      double r = rabPacked[packedIndex(j, i)];
      if (r > thr)
        continue;
      double r0 = rcov[numbers[i] - 1] + rcov[numbers[j] - 1];
      double dr = (r - r0) / r0;
      double derivative = kn / sqrtpi * std::exp(-kn * kn * dr * dr) / r0;
      double dx = xyz[3 * j] - xyz[3 * i];
      double dy = xyz[3 * j + 1] - xyz[3 * i + 1];
      double dz = xyz[3 * j + 2] - xyz[3 * i + 2];
      double rlen = std::sqrt(dx * dx + dy * dy + dz * dz);
      double rx = derivative * dx / rlen;
      double ry = derivative * dy / rlen;
      double rz = derivative * dz / rlen;
      dlog(j, j, 0) += dlogdcnj * rx;
      dlog(j, j, 1) += dlogdcnj * ry;
      dlog(j, j, 2) += dlogdcnj * rz;
      dlog(i, j, 0) = -dlogdcnj * rx;
      dlog(i, j, 1) = -dlogdcnj * ry;
      dlog(i, j, 2) = -dlogdcnj * rz;
      dlog(j, i, 0) = dlogdcni * rx;
      dlog(j, i, 1) = dlogdcni * ry;
      dlog(j, i, 2) = dlogdcni * rz;
      dlog(i, i, 0) += -dlogdcni * rx;
      dlog(i, i, 1) += -dlogdcni * ry;
      dlog(i, i, 2) += -dlogdcni * rz;
    }
  }
}

void metallicCharacter(int n, const std::vector<int>& numbers,
                       const std::vector<double>& enTable,
                       const std::vector<double>& cn,
                       const std::vector<double>& dlogCn,
                       std::vector<double>& mchar)
{
  mchar.assign(n, 0.0);
  for (int i = 0; i < n; ++i) {
    double dum2 = 0.0;
    for (int j = 0; j < n; ++j) {
      double dx = dlogCn[(j * n + i) * 3 + 0];
      double dy = dlogCn[(j * n + i) * 3 + 1];
      double dz = dlogCn[(j * n + i) * 3 + 2];
      dum2 += std::sqrt(dx * dx + dy * dy + dz * dz);
    }
    double en = enTable[numbers[i] - 1];
    double e8 = en * en;
    e8 *= e8;
    e8 *= e8;
    mchar[i] = std::exp(-0.005 * e8) * dum2 / (cn[i] + 1.0);
  }
}

} // namespace Xtb
} // namespace Avogadro
