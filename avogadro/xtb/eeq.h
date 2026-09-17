/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/eeq_model.f90 (goedecker_chrgeq) and src/charge_model.f90
  (new_charge_model_2019), Copyright (C) 2017-2020 Stefan Grimme.
******************************************************************************/

#ifndef AVOGADRO_XTB_EEQ_H
#define AVOGADRO_XTB_EEQ_H

#include <vector>

namespace Avogadro {
namespace Xtb {

class Environment;

// Per-atom charge-model parameters, mirroring xtb's chrg_parameter arrays
// (en, gam, kappa, alpha) as filled by new_charge_model_2019 from the GFN0
// tables: en = xi, gam = gamm, kappa = cnfak, alpha = alpg.
struct ChargeParams
{
  std::vector<double> en;
  std::vector<double> gam;
  std::vector<double> kappa;
  std::vector<double> alpha;
};

// Fill params from the GFN0 tables for the given 1-based atomic numbers.
// Records an environment error and returns false for Z outside 1..86.
bool buildChargeParams(const std::vector<int>& numbers, ChargeParams& params,
                       Environment& env);

// Erf-screened Coulomb matrix, mirroring get_coulomb_matrix_0d:
// amat[i][j] = erf(gamij * r) / r with
// gamij = 1 / sqrt(alpha[i]^2 + alpha[j]^2), amat[i][i] = gam[i].
// Row-major (n, n) output.
void coulombMatrix(const std::vector<int>& numbers,
                   const std::vector<double>& xyz, const ChargeParams& params,
                   std::vector<double>& amat);

// Goedecker electronegativity-equilibration charges, mirroring
// goedecker_chrgeq for the non-periodic case.
//
// Inputs: numbers/xyz (Bohr), totalCharge, params, erf-CN per atom (cn, e.g.
// from coordinationErf) and, when gradient or chargeDerivatives are wanted,
// its derivatives dcn in [(i * n + j) * 3 + c] layout.
//
// charges, energy and gradient follow the Fortran inout semantics: energy is
// incremented by the electrostatic energy, gradient is incremented by the
// IES gradient. chargeDerivatives, when given, receives d(q[i])/d(x[j])
// zeroed first, in [(i * n + j) * 3 + c] layout (the Lagrange-multiplier
// column is dropped exactly like xtb's gemm323 reshape does).
//
// Returns false (with an environment error) when the charge constraint is
// violated or the Coulomb matrix is singular.
bool eeqCharges(const std::vector<int>& numbers,
                const std::vector<double>& xyz, double totalCharge,
                const ChargeParams& params, const std::vector<double>& cn,
                const std::vector<double>* dcn, Environment& env,
                std::vector<double>& charges, double& energy,
                std::vector<double>* gradient = nullptr,
                std::vector<double>* chargeDerivatives = nullptr);

// Solve the symmetric linear system A * x = b (row-major A, overwritten)
// by LU factorisation with partial pivoting. Production xtb uses LAPACK
// dsysv (Bunch-Kaufman); results agree up to linear-solver rounding, and the
// real build should use Eigen instead. Returns false when singular.
bool solveSymmetric(std::vector<double>& mat, std::vector<double>& rhs,
                    int n);

// Invert the symmetric matrix A in place (row-major) by Gauss-Jordan
// elimination with partial pivoting, then copy the lower triangle to the
// upper one exactly like xtb's sytri + symmetrisation sequence does.
// Returns false when singular.
bool invertSymmetric(std::vector<double>& mat, int n);

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_EEQ_H
