/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_ini.f90 (out-of-plane/improper setup) and
  src/gfnff/gfnff_ini2.F90 (ssort),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#ifndef AVOGADRO_XTB_GFNFFOOP_H
#define AVOGADRO_XTB_GFNFFOOP_H

#include <vector>

namespace Avogadro {
namespace Xtb {

class Environment;
struct GffData;
struct GffGenerator;

// Per-atom inputs for the out-of-plane setup (topo% data).
struct OopAtom
{
  int element = 0;
  double charge = 0.0; // topo%qa
  int pi = 0; // post-Hückel piadr flag
};

// One improper: central atom with three neighbours (tlist(1..4),
// 0-based) and the kind flag (tlist(5)): -1 saturated N (double well),
// 0 pi case (harmonic). Cells are always the central cell here
// (non-periodic port of the numctr = 1 case).
struct Oop
{
  int center = -1;
  int first = -1;
  int second = -1;
  int third = -1;
  int kind = 0;
};

// One improper term: vtors(1..2) = phase (radians), force constant.
struct OopTerm
{
  double phase = 0.0;
  double forceConstant = 0.0;
};

// Selection sort of three distances with carried indices, mirroring ssort
// exactly (last-minimum selection, in-place). Called twice by the setup
// like the reference; ties resolve identically.
void sortTripleByDistance(std::vector<double>& dists, std::vector<int>& ids);

// Out-of-plane/improper list + terms mirroring the gfnff_ini block
// (non-periodic): three-fold coordinated pi atoms or nitrogens, neighbours
// sorted by central distance, saturated-N double well vs pi harmonic with
// carbonyl/halogen corrections. neighbours holds 0-based adjacency; xyz
// flat coordinates (any consistent unit); pboPacked the full-pair Hückel
// density (packedIndex layout) for sumppi.
bool buildOutOfPlane(const std::vector<OopAtom>& atoms,
                     const std::vector<std::vector<int>>& neighbours,
                     const std::vector<double>& xyz,
                     const std::vector<double>& pboPacked,
                     const std::vector<int>& group,
                     const std::vector<double>& repz, const GffData& param,
                     const GffGenerator& gen, std::vector<Oop>& oops,
                     std::vector<OopTerm>& terms, Environment& env);

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_GFNFFOOP_H
