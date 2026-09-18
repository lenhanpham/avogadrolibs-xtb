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

// N/O/F/S/Cl predicate mirroring nofs().
inline bool isNofs(int atomicNumber)
{
  int z = atomicNumber;
  return z == 7 || z == 8 || z == 9 || z == 16 || z == 17;
}

// Amide predicate mirroring amide() (ini2, 0d): N with hyb 3 and a pi
// flag, bonded to exactly one pi carbon that itself bonds a terminal pi
// oxygen. neighbours holds 0-based adjacency (single cell); piFlags the
// piadr map (1 = pi, 0 = not).
bool isAmide(int n, const std::vector<int>& numbers,
             const std::vector<int>& hyb,
             const std::vector<std::vector<int>>& neighbours,
             const std::vector<int>& piFlags, int atom);

// Amide-H predicate mirroring amideH() (ini2, 0d): degree-1 atom whose
// sole neighbour is an amide nitrogen with exactly one sp3-carbon
// neighbour.
bool isAmideHydrogen(int n, const std::vector<int>& numbers,
                     const std::vector<int>& hyb,
                     const std::vector<std::vector<int>>& neighbours,
                     const std::vector<int>& piFlags, int atom);

// EEQ xi corrections mirroring the ini dxi loop (0d): B/C/O/group-6/
// group-7 rules from hybridization, tags, pi flags, neighbour element
// counts and first-neighbour identities. itag/imetal parallel hyb;
// group/metal are element tables; counts holds neighbour counts.
// Unsuffixed literals are single-precision widened on use (int * float
// products stay in float), like the reference.
bool eeqXiCorrections(
  int n, const std::vector<int>& numbers,
  const std::vector<int>& itag, const std::vector<int>& imetal,
  const std::vector<int>& piFlags,
  const std::vector<std::vector<int>>& neighbours,
  const std::vector<int>& counts, const int* group, std::vector<double>& dxi,
  Environment& env);

// First EEQ parameters mirroring the ini block (0d): chieeq from -chi +
// dxi + cnf * sqrt(min(count, cnmax)), gameeq = gam, alpeeq = alp^2,
// with the mchishift (double) reduction for imetal == 2.
bool eeqInitialParams(int n, const std::vector<int>& numbers,
                      const std::vector<double>& chi,
                      const std::vector<double>& gam,
                      const std::vector<double>& alp,
                      const std::vector<double>& cnf,
                      const std::vector<int>& imetal,
                      const std::vector<int>& counts, double cnMax,
                      double mchiShift, const std::vector<double>& dxi,
                      std::vector<double>& chieeq,
                      std::vector<double>& gameeq,
                      std::vector<double>& alpeeq, Environment& env);

// Third-order gamma corrections mirroring the ini dgam loop (0d):
// element/hyb/metal/group rules scaled by the topology charges qa,
// using isAmide() for the amide branch. ff literals are
// single-precision widened on use.
bool eeqGammaCorrections(
  int n, const std::vector<int>& numbers, const std::vector<int>& hyb,
  const std::vector<int>& imetal, const std::vector<int>& piFlags,
  const std::vector<std::vector<int>>& neighbours, const int* group,
  const std::vector<double>& qa, std::vector<double>& dgam,
  Environment& env);

// Final atomic EEQ parameters mirroring the ini block (0d): chieeq from
// -chi + dxi with the amide-H -0.02 (single-widened) reduction, gameeq
// from gam + dgam, alpeeq from (alp + ff * qa)^2 with the element/metal
// ff rules (single-widened).
bool eeqFinalParams(int n, const std::vector<int>& numbers,
                    const std::vector<int>& hyb,
                    const std::vector<int>& imetal,
                    const std::vector<int>& piFlags,
                    const std::vector<std::vector<int>>& neighbours,
                    const std::vector<double>& chi,
                    const std::vector<double>& gam,
                    const std::vector<double>& alp,
                    const std::vector<double>& dxi,
                    const std::vector<double>& dgam,
                    const std::vector<double>& qa, const int* group,
                    std::vector<double>& chieeq,
                    std::vector<double>& gameeq, std::vector<double>& alpeeq,
                    Environment& env);

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_EEQ_H
