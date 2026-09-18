/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_ini.f90 (vbond setup loop), Copyright (C) 2019-2020
  Sebastian Ehlert, and src/gfnff/gfnff_ini2.F90 (ctype).
******************************************************************************/

#ifndef AVOGADRO_XTB_GFNFFVBOND_H
#define AVOGADRO_XTB_GFNFFVBOND_H

#include <vector>

namespace Avogadro {
namespace Xtb {

class Environment;
struct GffData;
struct GffGenerator;

// Per-atom inputs for the vbond setup (topo%/neigh% data).
struct VbondAtom
{
  int element = 0;
  int hyb = 0;
  double charge = 0.0; // topo%qa
  int itag = 0;
  int imetal = 0;
  int pi = 0; // post-Hückel piadr flag
  double mchar = 0.0;
};

// Per-bond inputs. first/second are 0-based atoms (ii/jj = blist(2/1,i)).
// guess is rtmp(ij) in Bohr; pibo the Hückel pi bond order; ring the
// smallestRingBond value; ringFirst/ringSecond the smallestRingThrough
// values (99 when in no ring, like ringsatom); coordFirst/coordSecond the
// neighbour counts nni/nnj. The bond type is computed inside (assignBondTypes
// block at the top of the reference loop); outBtypes carries it plus the
// pi-section update to 2.
struct VbondBond
{
  int first = -1;
  int second = -1;
  double guess = 0.0;
  double pibo = 0.0;
  int ring = 0;
  int ringFirst = 99;
  int ringSecond = 99;
  int coordFirst = 0;
  int coordSecond = 0;
};

// One bond term: vbond(1..3,i) = length shift (Bohr), steepness,
// prefactor (kbond, negative).
struct VbondTerm
{
  double shift = 0.0;
  double steepness = 0.0;
  double prefactor = 0.0;
};

// Carbonyl-carbon test mirroring ctype(): 1 when atom a carries exactly one
// pi-bonded oxygen neighbour and is a pi atom itself, else 0. neighbours
// holds 0-based adjacency (non-periodic, numctr = 1).
int carbonylCarbonType(const std::vector<std::vector<int>>& neighbours,
                       const std::vector<int>& numbers,
                       const std::vector<int>& piFlags, int a);

// Projected pi flags for the carbonyl test.
std::vector<int> piOf(const std::vector<VbondAtom>& atoms);

// Bond terms mirroring the vbond loop in gfnff_ini (bbtyp < 5 normal branch
// with bridge/hypervalent/XH/ring/pi/charge corrections, else the metal
// branch with TM/mtyp/heavy-ligand corrections, plus heavy-heavy shifts).
// group/metal/en/bondTab are element tables indexed element - 1; row6 is
// elementRow6 per atom. terms and outBtypes receive one entry per bond.
bool buildVbondTerms(const std::vector<VbondBond>& bonds,
                     const std::vector<VbondAtom>& atoms,
                     const std::vector<std::vector<int>>& neighbours,
                     const std::vector<int>& numbers,
                     const std::vector<int>& group,
                     const std::vector<int>& metalTable,
                     const std::vector<double>& en,
                     const std::vector<double>& bondTab,
                     const std::vector<int>& row6, const GffData& param,
                     const GffGenerator& gen, std::vector<VbondTerm>& terms,
                     std::vector<int>& outBtypes, Environment& env);

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_GFNFFVBOND_H
