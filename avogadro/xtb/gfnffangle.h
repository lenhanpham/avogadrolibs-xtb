/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_ini.f90 (alist/vangl setup), src/constr.f90 (banglPBC,
  mode 1 only) and src/gfnff/gfnff_ini2.F90 (ringsbend),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#ifndef AVOGADRO_XTB_GFNFFANGLE_H
#define AVOGADRO_XTB_GFNFFANGLE_H

#include <vector>

namespace Avogadro {
namespace Xtb {

class Environment;
struct GffData;
struct GffGenerator;
struct Ring;

// Per-atom inputs for the bend setup (topo% data). ring is the
// smallestRingThrough value (99 when in no ring, like ringsatom).
struct AngleAtom
{
  int element = 0;
  int hyb = 0;
  double charge = 0.0; // topo%qa
  int itag = 0;
  int imetal = 0;
  int pi = 0; // post-Hückel piadr flag
  int ring = 99;
};

// One bend: central atom, two neighbours (alist(1..3), 0-based) and the
// actual angle phi in radians. Cells are always the central cell here
// (non-periodic port of the numctr = 1 case).
struct Angle
{
  int center = -1;
  int first = -1;
  int second = -1;
  double phi = 0.0;
};

// One bend term: vangl(1..2) = equilibrium angle (radians), force constant.
struct AngleTerm
{
  double equilibrium = 0.0;
  double forceConstant = 0.0;
};

// Law-of-cosines angle at central, mirroring banglPBC mode 1 with zero
// shifts: acos(clamp(0.5*(d2ij+d2jk-d2ik)/sqrt(d2ij*d2jk))). xyz is flat
// row-major (x,y,z per atom); units cancel.
double bondAnglePbc(const std::vector<double>& xyz, int end1, int central,
                    int end2);

// Smallest ring containing angle i-j-k, mirroring ringsbend verbatim,
// INCLUDING the third loop comparing s(m,j) instead of s(m,k). Each table
// holds the atom's rings (0-based members); missing entries count as zero
// like the zero-initialized Fortran tables.
int smallestRingAngle(const std::vector<Ring>& ringsI,
                      const std::vector<Ring>& ringsJ,
                      const std::vector<Ring>& ringsK, int i, int j, int k);

// Angle enumeration mirroring the alist loop (non-periodic): central atoms
// with 2..6 neighbours, neighbour-slot pairs (j, k < j), angl-table
// product filter, actual-angle eta skip. neighbours holds 0-based
// adjacency; xyz flat coordinates (any consistent unit).
bool buildAngleList(int n, const std::vector<std::vector<int>>& neighbours,
                    const std::vector<int>& numbers,
                    const std::vector<double>& xyz,
                    const std::vector<double>& angl,
                    const std::vector<double>& angl2, double fcThr,
                    const std::vector<int>& metalTable,
                    std::vector<Angle>& angles, Environment& env);

// angl-table product for one angle (fijk in the reference loop).
double fijkOf(const Angle& angle, const std::vector<AngleAtom>& atoms,
              const std::vector<double>& angl,
              const std::vector<double>& angl2);

// Bend terms mirroring the vangl cascade (hyb defaults, B/C/O/N/NR3, ring,
// R-X, triple, carbene, SO3X, halogen, Pb/Sn, METAL rules; fqq/f2/fn/
// fbsmall/feta prefactors). ringsPerAtom holds each atom's rings for
// ringsbend; pboPacked the full-pair Hückel density (packedIndex layout)
// for the NR3 rule.
bool buildAngleTerms(const std::vector<Angle>& angles,
                     const std::vector<AngleAtom>& atoms,
                     const std::vector<std::vector<int>>& neighbours,
                     const std::vector<std::vector<Ring>>& ringsPerAtom,
                     const std::vector<double>& pboPacked,
                     const std::vector<int>& group,
                     const std::vector<int>& metalTable,
                     const std::vector<double>& angl,
                     const std::vector<double>& angl2, const GffData& param,
                     const GffGenerator& gen, std::vector<AngleTerm>& terms,
                     Environment& env);

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_GFNFFANGLE_H
