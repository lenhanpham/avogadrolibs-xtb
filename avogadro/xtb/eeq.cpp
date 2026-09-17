/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/eeq_model.f90 (goedecker_chrgeq, get_coulomb_matrix_0d) and
  src/charge_model.f90 (new_charge_model_2019),
  Copyright (C) 2017-2020 Stefan Grimme.
******************************************************************************/

#include "eeq.h"

#include "constants.h"
#include "coordination.h"
#include "environment.h"
#include "gfn0params.h"

#include <cmath>
#include <string>

namespace Avogadro {
namespace Xtb {

namespace {

// sqrt(2/pi), the sqrt2pi module parameter of eeq_model.
constexpr double sqrt2pi = 0.79788456080286535587989211986876374;

inline double& matAt(std::vector<double>& mat, int n, int i, int j)
{
  return mat[i * n + j];
}

} // namespace

bool buildChargeParams(const std::vector<int>& numbers, ChargeParams& params,
                       Environment& env)
{
  const int n = static_cast<int>(numbers.size());
  params.en.assign(n, 0.0);
  params.gam.assign(n, 0.0);
  params.kappa.assign(n, 0.0);
  params.alpha.assign(n, 0.0);
  for (int i = 0; i < n; ++i) {
    int z = numbers[i];
    if (z < 1 || z > gfn0Elements) {
      env.error("element " + std::to_string(z) +
                  " has no GFN0 charge parameters (Z must be 1..86)",
                "buildChargeParams");
      return false;
    }
    params.en[i] = electronegativity[z - 1];
    params.gam[i] = hardness[z - 1];
    params.kappa[i] = cnFactor[z - 1];
    params.alpha[i] = alpha[z - 1];
  }
  return true;
}

void coulombMatrix(const std::vector<int>& numbers,
                   const std::vector<double>& xyz, const ChargeParams& params,
                   std::vector<double>& amat)
{
  const int n = static_cast<int>(numbers.size());
  amat.assign(n * n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      double dx = xyz[3 * j] - xyz[3 * i];
      double dy = xyz[3 * j + 1] - xyz[3 * i + 1];
      double dz = xyz[3 * j + 2] - xyz[3 * i + 2];
      double r = std::sqrt(dx * dx + dy * dy + dz * dz);
      double gamij = 1.0 / std::sqrt(params.alpha[i] * params.alpha[i] +
                                     params.alpha[j] * params.alpha[j]);
      double v = std::erf(gamij * r) / r;
      matAt(amat, n, j, i) = v;
      matAt(amat, n, i, j) = v;
    }
    // Note: xtb squares the tabulated alpha into a scratch array first
    // (alpha(i) = alp(at(i))**2) and uses sqrt(alpha(i)) == alp below.
    matAt(amat, n, i, i) =
      params.gam[i] + sqrt2pi / std::sqrt(params.alpha[i] * params.alpha[i]);
  }
}

bool solveSymmetric(std::vector<double>& mat, std::vector<double>& rhs, int n)
{
  // Doolittle LU with partial pivoting; mirrors the dsysv reference path
  // used in the parity harness (production xtb calls LAPACK dsysv, which
  // agrees up to linear-solver rounding).
  for (int k = 0; k < n; ++k) {
    int p = k;
    double best = std::abs(matAt(mat, n, k, k));
    for (int i = k + 1; i < n; ++i) {
      double v = std::abs(matAt(mat, n, i, k));
      if (v > best) {
        best = v;
        p = i;
      }
    }
    if (best == 0.0)
      return false;
    if (p != k) {
      for (int j = 0; j < n; ++j)
        std::swap(matAt(mat, n, k, j), matAt(mat, n, p, j));
      std::swap(rhs[k], rhs[p]);
    }
    for (int i = k + 1; i < n; ++i) {
      matAt(mat, n, i, k) /= matAt(mat, n, k, k);
      for (int j = k + 1; j < n; ++j)
        matAt(mat, n, i, j) -= matAt(mat, n, i, k) * matAt(mat, n, k, j);
    }
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j)
      rhs[i] -= matAt(mat, n, i, j) * rhs[j];
  }
  for (int i = n - 1; i >= 0; --i) {
    for (int j = i + 1; j < n; ++j)
      rhs[i] -= matAt(mat, n, i, j) * rhs[j];
    rhs[i] /= matAt(mat, n, i, i);
  }
  return true;
}

bool invertSymmetric(std::vector<double>& mat, int n)
{
  // Gauss-Jordan on [A|I] with partial pivoting, then copy the lower
  // triangle to the upper one exactly like xtb's sytri + symmetrisation.
  std::vector<double> aug(n * 2 * n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j)
      aug[i * 2 * n + j] = matAt(mat, n, i, j);
    aug[i * 2 * n + n + i] = 1.0;
  }
  for (int k = 0; k < n; ++k) {
    int p = k;
    double best = std::abs(aug[k * 2 * n + k]);
    for (int i = k + 1; i < n; ++i) {
      double v = std::abs(aug[i * 2 * n + k]);
      if (v > best) {
        best = v;
        p = i;
      }
    }
    if (best == 0.0)
      return false;
    if (p != k) {
      for (int j = 0; j < 2 * n; ++j)
        std::swap(aug[k * 2 * n + j], aug[p * 2 * n + j]);
    }
    double piv = aug[k * 2 * n + k];
    for (int j = k; j < 2 * n; ++j)
      aug[k * 2 * n + j] /= piv;
    for (int i = 0; i < n; ++i) {
      if (i == k)
        continue;
      double f = aug[i * 2 * n + k];
      for (int j = k; j < 2 * n; ++j)
        aug[i * 2 * n + j] -= f * aug[k * 2 * n + j];
    }
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j)
      matAt(mat, n, i, j) = aug[i * 2 * n + n + j];
  }
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j)
      matAt(mat, n, i, j) = matAt(mat, n, j, i);
  }
  return true;
}

bool eeqCharges(const std::vector<int>& numbers,
                const std::vector<double>& xyz, double totalCharge,
                const ChargeParams& params, const std::vector<double>& cn,
                const std::vector<double>* dcn, Environment& env,
                std::vector<double>& charges, double& energy,
                std::vector<double>* gradient,
                std::vector<double>* chargeDerivatives)
{
  const char* source = "eeqCharges";
  const int n = static_cast<int>(numbers.size());
  const int m = n + 1;
  if (static_cast<int>(cn.size()) != n) {
    env.error("coordination number array has wrong size", source);
    return false;
  }
  if ((gradient || chargeDerivatives) && !dcn) {
    env.error("derivatives require coordination number derivatives", source);
    return false;
  }

  // Xvec(i) = -ENi + ki * CNi / (sqrt(CNi) + eps), alpha = alp^2.
  std::vector<double> xvec(m, 0.0), xfac(n, 0.0), alp2(n, 0.0);
  for (int i = 0; i < n; ++i) {
    double tmp = params.kappa[i] / (std::sqrt(cn[i]) + 1e-14);
    xvec[i] = -params.en[i] + tmp * cn[i];
    xfac[i] = 0.5 * tmp;
    alp2[i] = params.alpha[i] * params.alpha[i];
  }

  // Lagrangian Coulomb matrix with the charge-constraint border.
  std::vector<double> amat(m * m, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      double dx = xyz[3 * j] - xyz[3 * i];
      double dy = xyz[3 * j + 1] - xyz[3 * i + 1];
      double dz = xyz[3 * j + 2] - xyz[3 * i + 2];
      double r = std::sqrt(dx * dx + dy * dy + dz * dz);
      double gamij = 1.0 / std::sqrt(alp2[i] + alp2[j]);
      double v = std::erf(gamij * r) / r;
      matAt(amat, m, j, i) = v;
      matAt(amat, m, i, j) = v;
    }
    matAt(amat, m, i, i) = params.gam[i] + sqrt2pi / std::sqrt(alp2[i]);
    matAt(amat, m, m - 1, i) = 1.0;
    matAt(amat, m, i, m - 1) = 1.0;
  }
  matAt(amat, m, m - 1, m - 1) = 0.0;
  xvec[m - 1] = totalCharge;

  // The reference factorises a copy (Atmp) and keeps Amat intact for the
  // energy and gradient assembly, so do the same here.
  std::vector<double> alu = amat;
  std::vector<double> xtmp = xvec;
  if (!solveSymmetric(alu, xtmp, m)) {
    env.error("Coulomb matrix is singular, cannot solve lin. eq.", source);
    return false;
  }
  charges.assign(xtmp.begin(), xtmp.begin() + n);
  double qsum = 0.0;
  for (double q : charges)
    qsum += q;
  if (std::abs(qsum - totalCharge) > 1.0e-6) {
    env.error("charge constraint is not satisfied", source);
    return false;
  }

  // E = q . (1/2 A q - X), accumulated like the reference.
  double es = 0.0;
  for (int i = 0; i < m; ++i) {
    double ax = 0.0;
    for (int j = 0; j < m; ++j)
      ax += matAt(amat, m, i, j) * xtmp[j];
    es += xtmp[i] * (0.5 * ax - xvec[i]);
  }
  energy += es;

  if (!gradient && !chargeDerivatives)
    return true;

  // dXvec(c, j, k) = d(X_k)/d(x_jc), dAmat(c, j, k) as assembled below,
  // both [(k * n + j) * 3 + c]. Mirrors the reference including the
  // erf-CN derivative sign pattern it consumes.
  std::vector<double> dxvec(m * n * 3, 0.0), damat(m * n * 3, 0.0);
  std::vector<double> afac(3 * n, 0.0);
  auto dx = [&](int k, int j, int c) -> double& {
    return dxvec[(k * n + j) * 3 + c];
  };
  auto da = [&](int k, int j, int c) -> double& {
    return damat[(k * n + j) * 3 + c];
  };
  auto dc = [&](int i, int j, int c) -> double {
    return (*dcn)[(i * n + j) * 3 + c];
  };
  for (int i = 0; i < n; ++i) {
    for (int c = 0; c < 3; ++c)
      dx(i, i, c) = dc(i, i, c) * xfac[i];
    for (int j = 0; j < i; ++j) {
      double rij[3] = { xyz[3 * i] - xyz[3 * j],
                         xyz[3 * i + 1] - xyz[3 * j + 1],
                         xyz[3 * i + 2] - xyz[3 * j + 2] };
      double r2 = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];
      double gamij = 1.0 / std::sqrt(alp2[i] + alp2[j]);
      double arg = gamij * gamij * r2;
      double dtmp = 2.0 * gamij * std::exp(-arg) / (sqrtpi * r2) -
                    matAt(amat, m, j, i) / r2;
      for (int c = 0; c < 3; ++c) {
        // Index order follows the reference: first index = matrix column
        // (Lagrangian slot), second = gradient row (atom).
        dx(i, j, c) = dc(i, j, c) * xfac[i];
        dx(j, i, c) = dc(j, i, c) * xfac[j];
        da(j, i, c) = dtmp * rij[c] * xtmp[i];
        da(i, j, c) = -dtmp * rij[c] * xtmp[j];
        afac[i * 3 + c] += dtmp * rij[c] * xtmp[j];
        afac[j * 3 + c] += -dtmp * rij[c] * xtmp[i];
      }
    }
  }

  if (gradient) {
    if (static_cast<int>(gradient->size()) != 3 * n)
      gradient->assign(3 * n, 0.0);
    for (int i = 0; i < n; ++i) {
      for (int c = 0; c < 3; ++c) {
        double gda = 0.0, gdx = 0.0;
        for (int k = 0; k < m; ++k) {
          gda += da(k, i, c) * xtmp[k];
          gdx += dx(k, i, c) * xtmp[k];
        }
        (*gradient)[i * 3 + c] += gda - gdx;
      }
    }
  }

  if (chargeDerivatives) {
    // Invert the intact Lagrangian matrix (amat still holds it, since the
    // solve worked on a copy), mirroring the sytri + symmetrisation path.
    std::vector<double> ainv = amat;
    if (!invertSymmetric(ainv, m)) {
      env.error("Coulomb Matrix is singular, cannot invert", source);
      return false;
    }
    for (int i = 0; i < n; ++i) {
      for (int c = 0; c < 3; ++c)
        da(i, i, c) += afac[i * 3 + c];
    }
    // dq = -dAmat * Ainv + dXvec * Ainv, keeping the first n columns of
    // Ainv (the Lagrange column is dropped, like xtb's gemm323 reshape).
    chargeDerivatives->assign(3 * n * n, 0.0);
    auto dq = [&](int i, int j, int c) -> double& {
      return (*chargeDerivatives)[(i * n + j) * 3 + c];
    };
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        for (int c = 0; c < 3; ++c) {
          double v = 0.0;
          for (int k = 0; k < m; ++k)
            v += (-da(k, i, c) + dx(k, i, c)) * matAt(ainv, m, k, j);
          dq(i, j, c) = v;
        }
      }
    }
  }
  return true;
}

} // namespace Xtb
} // namespace Avogadro
