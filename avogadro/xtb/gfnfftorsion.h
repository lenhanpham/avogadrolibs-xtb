/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_ini.f90 (tlist/vtors setup) and
  src/gfnff/gfnff_ini2.F90 (ringstors, ringstorl, alphaCO),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#ifndef AVOGADRO_XTB_GFNFFTORSION_H
#define AVOGADRO_XTB_GFNFFTORSION_H

#include <vector>

namespace Avogadro {
namespace Xtb {

class Environment;
struct GffData;
struct GffGenerator;
struct Ring;

// Per-atom inputs for the torsion setup (topo% data).
struct TorsionAtom
{
  int element = 0;
  int hyb = 0;
  double charge = 0.0; // topo%qa
  int imetal = 0;
  int pi = 0; // post-Hückel piadr flag
};

// Central bond of a torsion search (ii/jj = blist(2/1,m), 0-based) with
// the assignBondTypes value and the Hückel pi bond order pibo(m).
struct TorsionBond
{
  int first = -1;
  int second = -1;
  int btype = 1;
  double pibo = 0.0;
};

// One torsion: outer/central atoms (ll,ii,jj,kk = tlist(1..4), 0-based)
// with the multiplicity nrot (tlist(5)). Cells are always the central
// cell here (non-periodic port of the numctr = 1 case).
struct Torsion
{
  int outer1 = -1;
  int center1 = -1;
  int center2 = -1;
  int outer2 = -1;
  int multiplicity = 1;
};

// One torsion term: vtors(1..2) = phase (radians), force constant. The
// sp3-sp3 gauche extra shares the layout (multiplicity 1, phase pi).
struct TorsionTerm
{
  double phase = 0.0;
  double forceConstant = 0.0;
};

// True when the dihedral is unusable, mirroring chktors(): either outer
// angle (jj-ii-kk or ii-jj-ll) exceeds 170 degrees under the exact
// (angle * 180) / 3.1415926 > 170 comparison.
bool dihedralNearLinear(const std::vector<double>& xyz, int ii, int jj,
                        int kk, int ll);

// Smallest ring containing torsion i-j-k-l, mirroring ringstors.
int smallestRingTorsion(const std::vector<Ring>& ringsI,
                        const std::vector<Ring>& ringsJ,
                        const std::vector<Ring>& ringsK,
                        const std::vector<Ring>& ringsL, int i, int j, int k,
                        int l);

// Largest ring containing torsion i-j-k-l, mirroring ringstorl.
int largestRingTorsion(const std::vector<Ring>& ringsI,
                       const std::vector<Ring>& ringsJ,
                       const std::vector<Ring>& ringsK,
                       const std::vector<Ring>& ringsL, int i, int j, int k,
                       int l);

// Alpha-carbonyl test mirroring alphaCO(): pi/sp3 C-C pair where the pi
// carbon carries exactly one terminal pi oxygen. neighbours holds 0-based
// adjacency (non-periodic, numctr = 1).
bool alphaCarbonyl(const std::vector<std::vector<int>>& neighbours,
                   const std::vector<int>& numbers,
                   const std::vector<int>& hyb,
                   const std::vector<int>& piFlags, int a, int b);

// Torsion list + terms mirroring the tlist/vtors loop (non-periodic):
// central-bond filters, H-count/nrot/phi rules, ring vs acyclic cases,
// sp3 specials, pi corrections, hypervalent ends, the fcthr cutoff and
// the sp3-sp3 gauche extra. neighbours holds 0-based adjacency;
// ringsPerAtom each atom's rings.
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
                   std::vector<TorsionTerm>& terms, Environment& env);

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_GFNFFTORSION_H
