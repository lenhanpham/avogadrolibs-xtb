/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_eg.f90 (non-bonded repulsion loop, goed_gfnff,
  ES gradient loop), Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#ifndef AVOGADRO_XTB_GFNFFEGNONBOND_H
#define AVOGADRO_XTB_GFNFFEGNONBOND_H

#include <vector>

namespace Avogadro {
namespace Xtb {

class Environment;
struct GffData;
struct GffGenerator;

// Non-bonded repulsion energy + gradient mirroring the REP loop in
// gfnff_eg (non-periodic, single cell): exponential overlap repulsion
// over all pairs above 1e-8 Bohr^2 and below repThr, skipping directly
// bonded pairs (flag exactly 1; 1,3/1,4 pairs keep scaled repulsion).
// alphanb holds per-pair exponents (n x n, only the lower triangle is
// read like the reference); bpair the exclusion flags. xyz is flat
// row-major, Bohr.
bool repulsionEnergyGradient(
  int n, const std::vector<int>& numbers, const std::vector<double>& xyz,
  const std::vector<double>& alphanb, const std::vector<int>& bondedPair,
  const std::vector<double>& repz, double repScaleN, double repThr,
  double mcfNrep, double& energy, std::vector<double>& gradient,
  double sigma[3][3], Environment& env);

// EEQ charges + electrostatic energy mirroring goed_gfnff
// (non-periodic, no solvation): Coulomb matrix with erf-screened
// off-diagonals from the gamma/erf intermediates (also returned, like
// eeqtmp), fragment charge constraints, LU solve (the reference uses
// LAPACK Bunch-Kaufman; eeq.cpp documents the rounding caveat).
// The reference's sqrab argument is unused in its body and dropped here.
// chieeq/gameeq/alpeeq are per-atom EEQ parameters, cn the D3
// coordination numbers, fragOfAtom 1-based fragment ids (all 1 when
// nfrag == 1), qfrag the fragment charges. gammaOut/erfOut receive the
// packed off-diagonal intermediates for the gradient pass.
bool eeqCharges(int n, const std::vector<int>& numbers,
                const std::vector<double>& rabPacked,
                const std::vector<double>& chieeq,
                const std::vector<double>& gameeq,
                const std::vector<double>& alpeeq,
                const std::vector<double>& cn,
                const std::vector<double>& cnf, int nfrag,
                const std::vector<int>& fragOfAtom,
                const std::vector<double>& qfrag, double totalCharge,
                std::vector<double>& charges, std::vector<double>& gammaOut,
                std::vector<double>& erfOut, double& energy,
                Environment& env);

// ES gradient + CN-response term mirroring the driver loop (non-periodic):
// pairwise erf-Coulomb derivatives from the goed intermediates plus the
// charge-response contraction g -= D^T qtmp with qtmp = q*cnf/(2sqrt(cn)).
// dlogCn uses the gffCoordinationNumber layout ([(i*n+j)*3+c] =
// d(cn[i])/d(x[j][c])), transposed vs Fortran dcn, hence the swapped
// pair index below.
bool esGradient(int n, const std::vector<double>& xyz,
                const std::vector<double>& rabPacked,
                const std::vector<double>& sqrabPacked,
                const std::vector<double>& gamma,
                const std::vector<double>& erf,
                const std::vector<double>& charges,
                const std::vector<double>& cn,
                const std::vector<double>& cnf,
                const std::vector<int>& numbers,
                const std::vector<double>& dlogCn,
                std::vector<double>& gradient, Environment& env);

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_GFNFFEGNONBOND_H
