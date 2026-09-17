/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_param.f90, src/param/sqrtzr4r2.f90 and
  src/param/covalentradd3.f90, Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#ifndef AVOGADRO_XTB_GFNFFPARAMS_H
#define AVOGADRO_XTB_GFNFFPARAMS_H

#include "constants.h"

#include <array>

namespace Avogadro {
namespace Xtb {

// Number of elements covered by the GFN-FF parametrisation (H..Lr).
constexpr int gffElements = 103;

// angewChem2020 element tables, index 0 = H. Exact doubles.
extern const std::array<double, gffElements> gffChi;
extern const std::array<double, gffElements> gffGam;
extern const std::array<double, gffElements> gffCnf;
extern const std::array<double, gffElements> gffAlp;
extern const std::array<double, gffElements> gffBond;
extern const std::array<double, gffElements> gffRepa;
extern const std::array<double, gffElements> gffRepan;
extern const std::array<double, gffElements> gffAngl;
extern const std::array<double, gffElements> gffAngl2;
extern const std::array<double, gffElements> gffTors;
extern const std::array<double, gffElements> gffTors2;

// Electronegativities. Unsuffixed (single precision) in the original, hence
// stored as float like the ncoord tables in coordination.cpp.
extern const std::array<float, gffElements> gffEn;
// Covalent radii in Angstrom (D0 suffix, exact doubles).
extern const std::array<double, gffElements> gffRad;
// Element flags: metal class, group, normal coordination number.
extern const std::array<int, gffElements> gffMetal;
extern const std::array<int, gffElements> gffGroup;
extern const std::array<int, gffElements> gffNormCn;
// Repulsion prefactors (Zval). Unsuffixed but integral, hence exact.
extern const std::array<double, gffElements> gffRepz;

// PBE0 r4/r2 expectation values, H..Og (exact doubles).
extern const std::array<double, 118> r4Overr2;
// D3 covalent radii in Angstrom (exact doubles); multiply by
// xtbAngstromToBohr * 4/3 for Bohr, mirroring covalentRadD3.
extern const std::array<double, 118> covalentRadD3Angstrom;

// sqrt(0.5 * r4Overr2 * sqrt(Z)), mirroring getSqrtZr4r2Number
// (returns -1 outside 1..118).
double sqrtZr4r2Value(int atomicNumber);

// D3 covalent radius in Bohr, mirroring covalentRadD3:
// ((Angstrom * aatoau) * 4) / 3 with xtb's own Angstrom definition.
inline double covalentRadiusD3(int atomicNumber)
{
  if (atomicNumber < 1 || atomicNumber > 118)
    return 0.0;
  return ((covalentRadD3Angstrom[atomicNumber - 1] * xtbAngstromToBohr) * 4.0) /
         3.0;
}

// Bond-length guess tables for gfnffrab (single precision originals).
extern const std::array<float, 103> rabEn;
extern const std::array<float, 103> rabR0;
extern const std::array<float, 103> rabCnfak;

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_GFNFFPARAMS_H
