/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/disp/ncoord.f90, Copyright (C) 2017-2020 Stefan Grimme.
******************************************************************************/

#ifndef AVOGADRO_XTB_COORDINATION_H
#define AVOGADRO_XTB_COORDINATION_H

#include <vector>

namespace Avogadro {
namespace Xtb {

// Fractional coordination numbers for D3/D4/GFN methods.
// Ported from xtb_disp_ncoord (non-periodic variants; PBC variants that need
// lattice summation helpers are not ported yet).
//
// Conventions (identical to the Fortran originals):
// - numbers[i] is the 1-based atomic number of atom i.
// - xyz holds Cartesian coordinates in Bohr, flat [3 * nat], xyz[3*i + c].
// - cn[i] receives the coordination number of atom i.
// - dcn[(i * nat + j) * 3 + c] receives d(cn[i])/d(xyz[j][c]).
// - thr is the squared-distance cutoff (default 1600 Bohr^2).

constexpr double coordinationThreshold = 1600.0;

// Covalent radii in Bohr (Pyykko/Atsumi, metals reduced by 10%).
// rcov[z] = 4/3 * rad[z] / 0.52917726 for 1 <= z <= 118.
double covalentRadius(int atomicNumber);

// Pauling electronegativities (dummies of 1.5 past element 86).
double paulingElectronegativity(int atomicNumber);

// Original D3 coordination number (2010) and its derivatives.
void coordinationD3(const std::vector<int>& numbers,
                    const std::vector<double>& xyz, std::vector<double>& cn,
                    double thr = coordinationThreshold);
void coordinationD3Derivative(const std::vector<int>& numbers,
                              const std::vector<double>& xyz,
                              std::vector<double>& cn,
                              std::vector<double>& dcn,
                              double thr = coordinationThreshold);

// Error-function damped coordination number (2018) and its derivatives.
//
// Port note: dncoord_erf's diagonal derivative terms carry the opposite sign
// from finite differences of ncoord_erf (verified against the reference and
// pinned in coordinationtest.cpp). This reproduces xtb exactly; the
// derivative outputs of this routine are not on xtb's production gradient
// paths, so keep the quirk for parity rather than "fixing" it.
void coordinationErf(const std::vector<int>& numbers,
                     const std::vector<double>& xyz, std::vector<double>& cn,
                     double thr = coordinationThreshold);
void coordinationErfDerivative(const std::vector<int>& numbers,
                               const std::vector<double>& xyz,
                               std::vector<double>& cn,
                               std::vector<double>& dcn,
                               double thr = coordinationThreshold);

// Covalent DFT-D4 coordination number and its derivatives.
//
// Port note: dncoord_d4's off-diagonal derivative terms carry the opposite
// sign from finite differences (verified against the reference and pinned in
// coordinationtest.cpp). xtb's D4 gradients hand-inline their own
// derivatives in dftd4's dispgrad instead of using these, so keep the quirk
// for parity rather than "fixing" it.
void coordinationD4(const std::vector<int>& numbers,
                    const std::vector<double>& xyz, std::vector<double>& cn,
                    double thr = coordinationThreshold);
void coordinationD4Derivative(const std::vector<int>& numbers,
                              const std::vector<double>& xyz,
                              std::vector<double>& cn,
                              std::vector<double>& dcn,
                              double thr = coordinationThreshold);

// Doubly damped GFN2-xTB coordination number and its derivatives.
void coordinationGfn(const std::vector<int>& numbers,
                     const std::vector<double>& xyz, std::vector<double>& cn,
                     double thr = coordinationThreshold);
void coordinationGfnDerivative(const std::vector<int>& numbers,
                               const std::vector<double>& xyz,
                               std::vector<double>& cn,
                               std::vector<double>& dcn,
                               double thr = coordinationThreshold);

// Smooth cutoff for large coordination numbers (log_cn_cut / dlog_cn_cut).
double logCoordinationCutoff(double cn, double max = 4.5);
double dLogCoordinationCutoff(double cn, double max = 4.5);
// Applies the cutoff to cn and, when given, its derivatives in place.
void applyLogCoordinationCutoff(std::vector<double>& cn,
                                std::vector<double>* dcn = nullptr,
                                int nat = 0, double max = 4.5);

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_COORDINATION_H
