/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_ini.f90 (out-of-plane/improper setup) and
  src/gfnff/gfnff_ini2.F90 (ssort),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#include "gfnffoop.h"

#include "constants.h"
#include "environment.h"
#include "gfnffdata.h"
#include "gfnfftopo.h"

#include <cmath>
#include <vector>

namespace Avogadro {
namespace Xtb {

void sortTripleByDistance(std::vector<double>& dists, std::vector<int>& ids)
{
  // Verbatim ssort for n = 3 (last-minimum selection sort, in place).
  int n = static_cast<int>(dists.size());
  for (int ii = 1; ii < n; ++ii) {
    int i = ii - 1;
    int k = i;
    double pp = dists[i];
    for (int j = ii; j < n; ++j) {
      if (dists[j] > pp)
        continue;
      k = j;
      pp = dists[j];
    }
    if (k == i)
      continue;
    dists[k] = dists[i];
    dists[i] = pp;
    int sc1 = ids[i];
    ids[i] = ids[k];
    ids[k] = sc1;
  }
}

static double distanceTo(const std::vector<double>& xyz, int center, int other)
{
  // gfortran NORM2, replicated exactly (verified by disassembly of the
  // frontend's inlined expansion plus a 3000-vector battery, bitwise):
  // scale starts at 1.0 (not 0), the ratio is squared before folding in.
  double dx = xyz[3 * center] - xyz[3 * other];
  double dy = xyz[3 * center + 1] - xyz[3 * other + 1];
  double dz = xyz[3 * center + 2] - xyz[3 * other + 2];
  double scale = 1.0, ssq = 0.0;
  const double comps[3] = { dx, dy, dz };
  for (int k = 0; k < 3; ++k) {
    if (comps[k] == 0.0)
      continue;
    double ax = std::fabs(comps[k]);
    if (ax > scale) {
      double r = scale / ax;
      ssq = 1.0 + ssq * (r * r);
      scale = ax;
    } else {
      double q = ax / scale;
      ssq = ssq + q * q;
    }
  }
  return scale * std::sqrt(ssq);
}

bool buildOutOfPlane(const std::vector<OopAtom>& atoms,
                     const std::vector<std::vector<int>>& neighbours,
                     const std::vector<double>& xyz,
                     const std::vector<double>& pboPacked,
                     const std::vector<int>& group,
                     const std::vector<double>& repz, const GffData& param,
                     const GffGenerator& gen, std::vector<Oop>& oops,
                     std::vector<OopTerm>& terms, Environment& env)
{
  (void)param;
  int n = static_cast<int>(atoms.size());
  if (n <= 0) {
    env.error("empty out-of-plane setup", "buildOutOfPlane");
    return false;
  }
  oops.clear();
  terms.clear();
  for (int i = 0; i < n; ++i) {
    if (static_cast<int>(neighbours[i].size()) != 3)
      continue;
    if (atoms[i].pi == 0 && atoms[i].element != 7)
      continue;
    int jj = neighbours[i][0], kk = neighbours[i][1],
        ll = neighbours[i][2];
    // Sort neighbours by central distance (two ssort passes like the
    // reference; cells stay central here).
    std::vector<double> sdum3 = { distanceTo(xyz, i, jj),
                                  distanceTo(xyz, i, kk),
                                  distanceTo(xyz, i, ll) };
    std::vector<int> ind3 = { jj, kk, ll };
    sortTripleByDistance(sdum3, ind3);
    sdum3 = { distanceTo(xyz, i, jj), distanceTo(xyz, i, kk),
              distanceTo(xyz, i, ll) };
    jj = ind3[0];
    kk = ind3[1];
    ll = ind3[2];
    Oop oop;
    oop.center = i;
    oop.first = jj;
    oop.second = kk;
    oop.third = ll;
    OopTerm term;
    if (atoms[i].pi == 0 && atoms[i].element == 7) {
      // Saturated nitrogen double well.
      oop.kind = -1;
      term.phase = (80.0 * pi) / 180.0f;
      term.forceConstant = 0.0;
      for (int idum : neighbours[i]) {
        term.forceConstant +=
          0.60 * std::sqrt(repz[atoms[idum].element - 1]);
      }
    } else {
      int ncarbo = 0, nf = 0;
      for (int idum : neighbours[i]) {
        if (atoms[idum].element == 8 || atoms[idum].element == 16)
          ++ncarbo;
        if (group[atoms[idum].element - 1] == 7)
          ++nf;
      }
      double fqq = 1.0 + atoms[i].charge * 5.0;
      oop.kind = 0;
      term.phase = 0.0;
      double sumppi = pboPacked[packedIndex(i, jj)] +
                      pboPacked[packedIndex(i, kk)] +
                      pboPacked[packedIndex(i, ll)];
      double f2 = 1.0 - sumppi * gen.torsf[5 - 1];
      term.forceConstant = gen.torsf[3 - 1] * f2 * fqq;
      // Unsuffixed 38./10. are single precision, widened on use.
      if (atoms[i].element == 5 && ncarbo > 0)
        term.forceConstant = term.forceConstant * 38.0f;
      if (atoms[i].element == 6 && ncarbo > 0)
        term.forceConstant = term.forceConstant * 38.0f;
      if (atoms[i].element == 6 && nf > 0 && ncarbo == 0)
        term.forceConstant = term.forceConstant * 10.0f;
      if (atoms[i].element == 7 && ncarbo > 0)
        term.forceConstant = term.forceConstant * 10.0f / f2;
    }
    oops.push_back(oop);
    terms.push_back(term);
  }
  return true;
}

} // namespace Xtb
} // namespace Avogadro
