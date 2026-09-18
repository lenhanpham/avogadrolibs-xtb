/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_eg.f90 (gfnff_eg, 0d non-periodic path),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#ifndef AVOGADRO_XTB_GFNFFDRIVER_H
#define AVOGADRO_XTB_GFNFFDRIVER_H

#include "gfnffhb.h"

#include <utility>
#include <vector>

namespace Avogadro {
namespace Xtb {

class Environment;
struct GffData;

// Resolved non-bonded repulsion tables (0d, single cell): per-pair
// exponents and exclusion flags. Only flag exactly 1 skips repulsion;
// 1,3/1,4 pairs (flags 2/3) keep scaled repulsion like the reference.
struct NonbondedTables
{
  std::vector<double> alphanb; // n x n, lower triangle read
  std::vector<int> bpair;      // n x n exclusion flags
};

// One bond: 0-based (first, second) endpoints plus vbond shift,
// steepness and prefactor (vbond(1..3) in the reference).
struct DriverBond
{
  int first = -1;
  int second = -1;
  double shift = 0.0;
  double steepness = 0.0;
  double prefactor = 0.0;
};

// One bend: 0-based (center, first, second) with vangl equilibrium
// (radians) and force constant.
struct DriverBend
{
  int center = -1;
  int first = -1;
  int second = -1;
  double equilibrium = 0.0;
  double forceConstant = 0.0;
};

// One torsion in egtors(i,j,k,l) order with tlist(5) multiplicity,
// vtors phase (radians) and force constant.
struct DriverTorsion
{
  int i = -1, j = -1, k = -1, l = -1;
  int multiplicity = 1;
  double phase = 0.0;
  double forceConstant = 0.0;
};

// One bonded-ATM triple.
struct DriverTriple
{
  int iat = -1, jat = -1, kat = -1;
};

// Explicit HB/XB entries (list builders are a later batch; the kernels
// and dispatch below are verified). Neighbours come from the shared
// adjacency (jth_nb order); hb3 carries C explicitly.
struct DriverHb1
{
  int a = -1, b = -1, h = -1;
};
struct DriverHb2
{
  int a = -1, b = -1, h = -1;
};
struct DriverHb3
{
  int a = -1, b = -1, h = -1, c = -1;
};
struct DriverXb
{
  int a = -1, b = -1, x = -1;
};

// Full 0d single-point input. EEQ parameters (chieeq/gameeq/alpeeq),
// topology charges (qa), EEQ fragments (fragOfAtom/fragCharges, single
// fragment when empty), non-bonded tables and bpair flags come from the
// setup pipeline (later batches); the driver only assembles energies.
struct DriverInput
{
  std::vector<int> numbers;
  std::vector<double> xyz; // flat row-major, Bohr
  double accuracy = 1.0;
  double totalCharge = 0.0;
  // mcf scalings (angewChem2020: all 1.0; mc variant scales ees/ehb/rep).
  double mcfEes = 1.0, mcfEhb = 1.0, mcfNrep = 1.0;
  std::vector<DriverBond> bonds;
  std::vector<DriverBend> bends;
  std::vector<DriverTorsion> torsions;
  std::vector<DriverTriple> triples;
  std::vector<DriverHb1> hb1;
  std::vector<DriverHb2> hb2;
  std::vector<DriverHb3> hb3;
  std::vector<DriverXb> xb;
  // 0-based adjacency in jth_nb order (drives the hb2 carbonyl/nitro/rnr
  // dispatch and neighbour loops).
  std::vector<std::vector<int>> neighbours;
  std::vector<int> hbAtH; // not used yet (hbset lists are explicit)
  NonbondedTables nonbonded;
  std::vector<double> chieeq, gameeq, alpeeq; // per-atom EEQ params
  std::vector<double> qa;                     // topology charges
  std::vector<double> hbBas, hbAci;           // per-atom HB strengths
  std::vector<int> nrHb;      // per-bond HB counts (nr_hb, blist order)
  HbBondMaps bondHb;          // A-H...B groups for the egbond_hb branch
  std::vector<int> fragOfAtom;     // 1-based EEQ fragments (topo%fraglist)
  std::vector<double> fragCharges; // per-fragment charges (topo%qfrag)
};

// 0d single-point result: per-term energies, total, gradient, EEQ
// charges, dipole and gradient norm.
struct DriverResult
{
  double ees = 0.0, edisp = 0.0, erep = 0.0, ebond = 0.0, eangl = 0.0;
  double etors = 0.0, ehb = 0.0, exb = 0.0, ebatm = 0.0, etot = 0.0;
  std::vector<double> gradient;
  std::vector<double> charges;
  std::vector<double> cn, dlogCn; // erf CN + derivatives used inside
  double dipole[3] = { 0.0, 0.0, 0.0 };
  double gnorm = 0.0;
};

// GFN-FF 0d single point mirroring the non-periodic gfnff_eg path:
// thresholds, packed distances + D3 list, erf CN, EEQ, D3, ES gradient,
// HB erf CN (gated on the nrHb counts like dncoord_erf), rab estimates,
// bond loop (egbond, or egbond_hb for bonds with nrHb >= 1), bonded
// repulsion, bends, torsions, bonded ATM, HB/XB dispatch and the etot
// sum. Deviations from the reference, all gradient/energy-neutral in
// 0d: the dEdcn->sigma contraction reads uninitialized dcndL there and
// is skipped, as are the sigma updates inside egbond_hb (sigma is
// discarded by the 0d driver); eesinf/De printout, timers, GBSA, efield
// and sTors terms are out of scope.
bool gfnffSinglePoint(const DriverInput& input, const GffData& param,
                      const HbParams& hbpar, DriverResult& result,
                      Environment& env);

} // namespace Xtb
} // namespace Avogadro

#endif
