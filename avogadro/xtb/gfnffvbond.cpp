/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_ini.f90 (vbond setup loop), Copyright (C) 2019-2020
  Sebastian Ehlert, and src/gfnff/gfnff_ini2.F90 (ctype).
******************************************************************************/

#include "gfnffvbond.h"

#include "environment.h"
#include "gfnffdata.h"
#include "gfnffbonds.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace Avogadro {
namespace Xtb {

int carbonylCarbonType(const std::vector<std::vector<int>>& neighbours,
                       const std::vector<int>& numbers,
                       const std::vector<int>& piFlags, int a)
{
  int no = 0;
  for (int j : neighbours[a]) {
    if (numbers[j] == 8 && piFlags[j] != 0)
      ++no;
  }
  if (no == 1 && piFlags[a] != 0)
    return 1;
  return 0;
}

// Projected pi flags for the carbonyl test.
std::vector<int> piOf(const std::vector<VbondAtom>& atoms)
{
  std::vector<int> pi(atoms.size(), 0);
  for (size_t k = 0; k < atoms.size(); ++k)
    pi[k] = atoms[k].pi;
  return pi;
}

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
                     std::vector<int>& outBtypes, Environment& env)
{
  (void)param;
  int nbond = static_cast<int>(bonds.size());
  int n = static_cast<int>(atoms.size());
  if (n <= 0 || nbond <= 0) {
    env.error("empty bond setup", "buildVbondTerms");
    return false;
  }
  terms.assign(nbond, VbondTerm());
  outBtypes.assign(nbond, 1);
  // Bond types from the assignment block at the top of the reference loop.
  // assignBondTypes mirrors it verbatim (conditions symmetric in ii/jj).
  std::vector<Bond> bondList;
  bondList.reserve(bonds.size());
  std::vector<int> atomNumbers(n), hyb(n), itag(n), piAtoms(n), imetal(n);
  for (int i = 0; i < n; ++i) {
    atomNumbers[i] = atoms[i].element;
    hyb[i] = atoms[i].hyb;
    itag[i] = atoms[i].itag;
    piAtoms[i] = atoms[i].pi;
    imetal[i] = atoms[i].imetal;
  }
  for (const VbondBond& bond : bonds) {
    Bond b;
    b.first = bond.first;
    b.second = bond.second;
    bondList.push_back(b);
  }
  std::vector<int> btypes =
    assignBondTypes(bondList, atomNumbers, hyb, itag, piAtoms, imetal,
                    group.data());
  for (int i = 0; i < nbond; ++i) {
    const VbondBond& bond = bonds[i];
    int ii = bond.first, jj = bond.second;
    int ia = atoms[ii].element, ja = atoms[jj].element;
    int nni = bond.coordFirst, nnj = bond.coordSecond;
    int rings = bond.ring;
    double shift = 0.0, fxh = 1.0, ringf = 1.0, fqq = 1.0, fpi = 1.0;
    double fheavy = 1.0, fcn = 1.0, fsrb2 = gen.srb2, bstrength = 0.0;
    int bbtyp = btypes[i];
    outBtypes[i] = bbtyp;
    // Bridge flags recomputed from the btyp assignment (assignBondTypes
    // keeps only the type, so the linear-halogen condition is repeated).
    bool bridge = false;
    if ((group[ia - 1] == 7 || ia == 1) && atoms[ii].hyb == 1)
      bridge = true;
    if ((group[ja - 1] == 7 || ja == 1) && atoms[jj].hyb == 1)
      bridge = true;
    if (bbtyp < 5) {
      // Normal bond.
      int hybi = std::max(atoms[ii].hyb, atoms[jj].hyb);
      int hybj = std::min(atoms[ii].hyb, atoms[jj].hyb);
      if (hybi == 5 || hybj == 5)
        bstrength = gen.bstren[4 - 1];
      else
        bstrength = gen.bsmat[hybi][hybj];
      if (hybi == 3 && hybj == 2 && (ia == 7 || ja == 7))
        // Unsuffixed 1.04 is single precision, widened on use.
        bstrength = gen.bstren[2 - 1] * 1.04f;
      if (bridge) {
        if (group[ia - 1] == 7)
          bstrength = gen.bstren[1 - 1] * 0.50;
        if (group[ja - 1] == 7)
          bstrength = gen.bstren[1 - 1] * 0.50;
        if (ia == 1 || ia == 9)
          bstrength = gen.bstren[1 - 1] * 0.30;
        if (ja == 1 || ja == 9)
          bstrength = gen.bstren[1 - 1] * 0.30;
      }
      if (bbtyp == 4)
        shift = gen.hyperShift;
      if (ia == 1 || ja == 1)
        shift = gen.rabShiftH;
      if (ia == 9 && ja == 9)
        shift = 0.22f;
      if (atoms[ii].hyb == 3 && atoms[jj].hyb == 0)
        shift = shift - 0.022f;
      if (atoms[ii].hyb == 0 && atoms[jj].hyb == 3)
        shift = shift - 0.022f;
      if (atoms[ii].hyb == 1 && atoms[jj].hyb == 0)
        shift = shift + 0.14f;
      if (atoms[ii].hyb == 0 && atoms[jj].hyb == 1)
        shift = shift + 0.14f;
      if (ia == 1 && ja == 6) {
        if (bond.ringSecond == 3)
          fxh = 1.05f;
        if (carbonylCarbonType(neighbours, numbers, piOf(atoms), jj) == 1)
          fxh = 0.95f;
      }
      if (ia == 6 && ja == 1) {
        if (bond.ringFirst == 3)
          fxh = 1.05f;
        if (carbonylCarbonType(neighbours, numbers, piOf(atoms), ii) == 1)
          fxh = 0.95f;
      }
      if (ia == 1 && ja == 5)
        fxh = 1.10f;
      if (ja == 1 && ia == 5)
        fxh = 1.10f;
      if (ia == 1 && ja == 7)
        fxh = 1.06f;
      if (ja == 1 && ia == 7)
        fxh = 1.06f;
      if (ia == 1 && ja == 8)
        fxh = 0.93f;
      if (ja == 1 && ia == 8)
        fxh = 0.93f;
      if (bbtyp == 3 && ia == 6 && ja == 8)
        bstrength = gen.bstren[3 - 1] * 0.90;
      if (bbtyp == 3 && ia == 8 && ja == 6)
        bstrength = gen.bstren[3 - 1] * 0.90;
      // Triple-bond demotions stay local (bbtyp only, like the reference).
      if (bbtyp == 3 &&
          (atoms[ii].hyb == 0 || atoms[jj].hyb == 0))
        bbtyp = 1;
      if (bbtyp == 3 &&
          (atoms[ii].hyb == 3 || atoms[jj].hyb == 3))
        bbtyp = 1;
      if (bbtyp == 3 &&
          (atoms[ii].hyb == 2 || atoms[jj].hyb == 2))
        bbtyp = 2;
      // Pi corrections from the Hückel bond orders.
      if (bond.pibo > 0) {
        shift = gen.hueckelP * (gen.bzRef - bond.pibo);
        if (bbtyp != 3 && bond.pibo > 0.1) {
          outBtypes[i] = 2;
          bbtyp = 2;
        }
        fpi = 1.0 - gen.hueckelP2 * (gen.bzRef2 - bond.pibo);
      }
      if (ia > 10 && ja > 10) {
        double dni = static_cast<double>(nni);
        double dnj = static_cast<double>(nnj);
        // Unsuffixed 0.007 is single precision; dble(nni)**2 = dni*dni.
        fcn = fcn / (1.0 + 0.007f * (dni * dni));
        fcn = fcn / (1.0 + 0.007f * (dnj * dnj));
      }
      double qafac = atoms[ii].charge * atoms[jj].charge * 70.0;
      fqq = 1.0 + gen.qfacbm0 / (1.0 + std::exp(15.0 * qafac));
    } else {
      // Metal-involving bond.
      shift = 0.0;
      bstrength = gen.bstren[bbtyp - 1];
      if (bbtyp == 7) {
        if (row6[ii] > 4 && row6[jj] > 4)
          bstrength = gen.bstren[8 - 1];
        if (row6[ii] == 4 && row6[jj] > 4)
          bstrength = gen.bstren[9 - 1];
        if (row6[jj] == 4 && row6[ii] > 4)
          bstrength = gen.bstren[9 - 1];
        double dum = 2.0 * atoms[ii].mchar + 2.0 * atoms[jj].mchar;
        dum = std::min(dum, 0.5);
        bstrength = bstrength * (1.0 - dum);
      }
      int mtyp1 = 0, mtyp2 = 0;
      if (group[ia - 1] == 1)
        mtyp1 = 1;
      if (group[ia - 1] == 2)
        mtyp1 = 2;
      if (group[ia - 1] > 2 && atoms[ii].imetal == 1)
        mtyp1 = 3;
      if (atoms[ii].imetal == 2)
        mtyp1 = 4;
      if (group[ja - 1] == 1)
        mtyp2 = 1;
      if (group[ja - 1] == 2)
        mtyp2 = 2;
      if (group[ja - 1] > 2 && atoms[jj].imetal == 1)
        mtyp2 = 3;
      if (atoms[jj].imetal == 2)
        mtyp2 = 4;
      double qafac = atoms[ii].charge * atoms[jj].charge * 25.0;
      double dum = 1.0 / (1.0 + std::exp(15.0 * qafac));
      fqq = 1.0 + dum * (gen.qfacbm[mtyp1] + gen.qfacbm[mtyp2]) * 0.5;
      if (atoms[ii].imetal == 2 && ja > 10)
        fheavy = 0.65;
      if (atoms[jj].imetal == 2 && ia > 10)
        fheavy = 0.65;
      if (atoms[ii].imetal == 2 && ja == 15)
        fheavy = 1.60;
      if (atoms[jj].imetal == 2 && ia == 15)
        fheavy = 1.60;
      if (atoms[ii].imetal == 2 && group[ja - 1] == 6)
        fheavy = 0.85;
      if (atoms[jj].imetal == 2 && group[ia - 1] == 6)
        fheavy = 0.85;
      if (atoms[ii].imetal == 2 && group[ja - 1] == 7)
        fheavy = 1.30;
      if (atoms[jj].imetal == 2 && group[ia - 1] == 7)
        fheavy = 1.30;
      if (atoms[ii].imetal == 2 && ja == 1 && row6[ii] <= 5)
        fxh = 0.80;
      if (atoms[jj].imetal == 2 && ia == 1 && row6[jj] <= 5)
        fxh = 0.80;
      if (atoms[ii].imetal == 2 && ja == 1 && row6[ii] > 5)
        fxh = 1.00;
      if (atoms[jj].imetal == 2 && ia == 1 && row6[jj] > 5)
        fxh = 1.00;
      if (atoms[ii].imetal == 1 && ja == 1)
        fxh = 1.20;
      if (atoms[jj].imetal == 1 && ia == 1)
        fxh = 1.20;
      if (atoms[ii].imetal == 2 && atoms[jj].hyb == 1) {
        if (ja == 6) {
          fpi = 1.5;
          shift = -0.45;
        }
        if (ja == 7 && nnj != 1) {
          fpi = 0.4;
          shift = 0.47;
        }
      }
      if (atoms[jj].imetal == 2 && atoms[ii].hyb == 1) {
        if (ia == 6) {
          fpi = 1.5;
          shift = -0.45;
        }
        if (ia == 7 && nni != 1) {
          fpi = 0.4;
          shift = 0.47;
        }
      }
      if (atoms[ii].imetal == 2)
        shift = shift + gen.metal2Shift;
      if (atoms[jj].imetal == 2)
        shift = shift + gen.metal2Shift;
      if (atoms[ii].imetal == 1 && group[ia - 1] <= 2)
        shift = shift + gen.metal1Shift;
      if (atoms[jj].imetal == 1 && group[ja - 1] <= 2)
        shift = shift + gen.metal1Shift;
      if (mtyp1 == 3)
        shift = shift + gen.metal3Shift;
      if (mtyp2 == 3)
        shift = shift + gen.metal3Shift;
      if (bbtyp == 6 && metalTable[ia - 1] == 2)
        shift = shift + gen.etaShift * static_cast<double>(nni);
      if (bbtyp == 6 && metalTable[ja - 1] == 2)
        shift = shift + gen.etaShift * static_cast<double>(nnj);
      // Unsuffixed 0.100/0.030/0.036 are single precision;
      // dble(nni)**2 = dni*dni.
      if (mtyp1 > 0 && mtyp1 < 3) {
        double dni = static_cast<double>(nni);
        fcn = fcn / (1.0 + 0.100f * (dni * dni));
      }
      if (mtyp2 > 0 && mtyp2 < 3) {
        double dnj = static_cast<double>(nnj);
        fcn = fcn / (1.0 + 0.100f * (dnj * dnj));
      }
      if (mtyp1 == 3) {
        double dni = static_cast<double>(nni);
        fcn = fcn / (1.0 + 0.030f * (dni * dni));
      }
      if (mtyp2 == 3) {
        double dnj = static_cast<double>(nnj);
        fcn = fcn / (1.0 + 0.030f * (dnj * dnj));
      }
      if (mtyp1 == 4) {
        double dni = static_cast<double>(nni);
        fcn = fcn / (1.0 + 0.036f * (dni * dni));
      }
      if (mtyp2 == 4) {
        double dnj = static_cast<double>(nnj);
        fcn = fcn / (1.0 + 0.036f * (dnj * dnj));
      }
      // Unsuffixed 0.22/0.28 are single precision, widened on use.
      if (mtyp1 == 4 || mtyp2 == 4)
        fsrb2 = -gen.srb2 * 0.22f;
      else
        fsrb2 = gen.srb2 * 0.28f;
    }
    if (ia > 10 && ja > 10) {
      shift = shift + gen.hShift3;
      if (ia > 18)
        shift = shift + gen.hShift4;
      if (ja > 18)
        shift = shift + gen.hShift4;
      if (ia > 36)
        shift = shift + gen.hShift5;
      if (ja > 36)
        shift = shift + gen.hShift5;
    }
    VbondTerm term;
    term.shift = gen.rabShift + shift;
    if (rings > 0)
      ringf = 1.0 + gen.fringbo * (6.0 - static_cast<double>(rings)) *
                        (6.0 - static_cast<double>(rings));
    double enDiff = en[ia - 1] - en[ja - 1];
    term.steepness =
      gen.srb1 * (1.0 + fsrb2 * enDiff * enDiff + gen.srb3 * bstrength);
    term.prefactor = -bondTab[ia - 1] * bondTab[ja - 1] * ringf * bstrength *
                     fqq * fheavy * fpi * fxh * fcn;
    terms[i] = term;
  }
  return true;
}

} // namespace Xtb
} // namespace Avogadro
