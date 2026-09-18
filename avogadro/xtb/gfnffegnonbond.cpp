/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_eg.f90 (non-bonded repulsion loop, goed_gfnff,
  ES gradient loop), Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#include "gfnffegnonbond.h"

#include "constants.h"
#include "eeq.h"
#include "environment.h"
#include "gfnffdata.h"
#include "gfnffegbond.h"

#include <cmath>
#include <vector>

namespace Avogadro {
namespace Xtb {

bool repulsionEnergyGradient(
  int n, const std::vector<int>& numbers, const std::vector<double>& xyz,
  const std::vector<double>& alphanb, const std::vector<int>& bondedPair,
  const std::vector<double>& repz, double repScaleN, double repThr,
  double mcfNrep, double& energy, std::vector<double>& gradient,
  double sigma[3][3], Environment& env)
{
  if (n <= 0) {
    env.error("empty repulsion setup", "repulsionEnergyGradient");
    return false;
  }
  for (int iat = 0; iat < n; ++iat) {
    for (int jat = 0; jat <= iat; ++jat) {
      double dxyz[3] = { xyz[3 * iat] - xyz[3 * jat],
                         xyz[3 * iat + 1] - xyz[3 * jat + 1],
                         xyz[3 * iat + 2] - xyz[3 * jat + 2] };
      double n2 = norm2(dxyz);
      double r2 = n2 * n2;
      if (r2 > repThr || r2 < 1.0e-8)
        continue;
      // Only directly bonded pairs (flag 1) are skipped; 1,3/1,4 pairs
      // (flags 2/3) keep scaled repulsion via alphanb.
      if (bondedPair[iat * n + jat] == 1)
        continue;
      int ati = numbers[iat], atj = numbers[jat];
      double rab = std::sqrt(r2);
      double t16 = std::pow(r2, static_cast<double>(0.75f));
      double t19 = t16 * t16;
      double t8 = t16 * alphanb[iat * n + jat];
      double t26 = std::exp(-t8) * repz[ati - 1] * repz[atj - 1] *
                   repScaleN * mcfNrep;
      energy = energy + (t26 / rab);
      double t27 = t26 * (1.5 * t8 + 1.0) / t19;
      double r3[3] = { dxyz[0] * t27, dxyz[1] * t27, dxyz[2] * t27 };
      for (int c = 0; c < 3; ++c) {
        sigma[c][0] -= r3[0] * dxyz[c];
        sigma[c][1] -= r3[1] * dxyz[c];
        sigma[c][2] -= r3[2] * dxyz[c];
      }
      for (int c = 0; c < 3; ++c) {
        gradient[3 * iat + c] -= r3[c];
        gradient[3 * jat + c] += r3[c];
      }
    }
  }
  return true;
}

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
                Environment& env)
{
  if (n <= 0) {
    env.error("empty EEQ setup", "eeqCharges");
    return false;
  }
  const double tsqrt2pi = 0.797884560802866;
  int m = n + nfrag;
  std::vector<double> amat(m * m, 0.0), xvec(m, 0.0);
  for (int i = 0; i < n; ++i)
    xvec[i] = chieeq[i] + cnf[numbers[i] - 1] * std::sqrt(cn[i]);
  auto at = [&](int r, int c) -> double& { return amat[r * m + c]; };
  int npair = n * (n + 1) / 2;
  gammaOut.assign(npair, 0.0);
  erfOut.assign(npair, 0.0);
  for (int i = 0; i < n; ++i) {
    at(i, i) = tsqrt2pi / std::sqrt(alpeeq[i]) + gameeq[i];
    for (int j = 0; j < i; ++j) {
      int ij = j + i * (i + 1) / 2;
      double gammij = 1.0 / std::sqrt(alpeeq[i] + alpeeq[j]);
      double tmp = std::erf(gammij * rabPacked[ij]);
      gammaOut[ij] = gammij;
      erfOut[ij] = tmp;
      at(j, i) = tmp / rabPacked[ij];
      at(i, j) = at(j, i);
    }
  }
  for (int i = 0; i < nfrag; ++i) {
    xvec[n + i] = qfrag[i];
    for (int j = 0; j < n; ++j) {
      if (fragOfAtom[j] == i + 1) {
        at(n + i, j) = 1.0;
        at(j, n + i) = 1.0;
      }
    }
  }
  // No solvation (GBSA) in this port; the reference adds bornMat here.
  if (!solveSymmetric(amat, xvec, m)) {
    env.error("Coulomb matrix is singular, cannot solve lin. eq.",
              "eeqCharges");
    return false;
  }
  charges.assign(xvec.begin(), xvec.begin() + n);
  if (n == 1)
    charges[0] = totalCharge;
  energy = 0.0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      int ij = j + i * (i + 1) / 2;
      energy += charges[i] * charges[j] * erfOut[ij] / rabPacked[ij];
    }
    energy = energy - charges[i] *
                          (chieeq[i] + cnf[numbers[i] - 1] * std::sqrt(cn[i]));
    energy = energy + charges[i] * charges[i] * 0.5 *
                          (gameeq[i] + tsqrt2pi / std::sqrt(alpeeq[i]));
  }
  return true;
}

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
                std::vector<double>& gradient, Environment& env)
{
  if (n <= 0) {
    env.error("empty ES gradient setup", "esGradient");
    return false;
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      int ij = j + i * (i + 1) / 2;
      double r2 = sqrabPacked[ij];
      double rab = rabPacked[ij];
      double gammij = gamma[ij];
      double erff = erf[ij];
      double gg = gammij * gammij;
      double dd = (2.0 * gammij * std::exp(-(gg * r2)) / (sqrtpi * r2) -
                   erff / (rab * r2)) *
                  charges[i] * charges[j];
      for (int c = 0; c < 3; ++c) {
        double r3 = (xyz[3 * i + c] - xyz[3 * j + c]) * dd;
        gradient[3 * i + c] += r3;
        gradient[3 * j + c] -= r3;
      }
    }
  }
  // Charge-response contraction g -= D^T qtmp with the transposed pair
  // index (see header): D[(i,j)][c] lives at [(j*n+i)*3+c].
  for (int i = 0; i < n; ++i) {
    double qtmp = charges[i] * cnf[numbers[i] - 1] /
                  (2.0 * std::sqrt(cn[i]) + 1.0e-16);
    for (int j = 0; j < n; ++j) {
      for (int c = 0; c < 3; ++c)
        gradient[3 * j + c] -= dlogCn[(j * n + i) * 3 + c] * qtmp;
    }
  }
  // NOTE: mctc_gemv form is g += alpha*D*qtmp with alpha = -1, beta = 1.
  (void)env;
  return true;
}

} // namespace Xtb
} // namespace Avogadro
