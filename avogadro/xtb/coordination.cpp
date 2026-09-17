/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/disp/ncoord.f90, Copyright (C) 2017-2020 Stefan Grimme.
******************************************************************************/

#include "coordination.h"

#include "constants.h"
#include "elements.h"

#include <array>
#include <cmath>

namespace Avogadro {
namespace Xtb {

namespace {

// Counting-function parameters, verbatim from xtb_disp_ncoord.
constexpr double k1 = 16.0;
constexpr double ka = 10.0;
constexpr double kb = 20.0;
constexpr double rShift = 2.0;
constexpr double k4 = 4.10451;
constexpr double k5 = 19.08857;
constexpr double k6 = 2.0 * 11.28174 * 11.28174;
constexpr double kn = 7.50;

// Covalent radii in Angstrom (Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
// 188-197), values for metals decreased by 10%.
//
// Stored as float on purpose: in the Fortran original these initialisers are
// unsuffixed literals, which are single precision and only widened to double
// on use. Keeping the same rounding reproduces xtb's values bitwise.
const std::array<float, maxElement> rad = {
  0.32, 0.46, // H,He
  1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67, // Li-Ne
  1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96, // Na-Ar
  1.76, 1.54, // K,Ca
  1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09, // Sc-Zn
  1.12, 1.09, 1.15, 1.10, 1.14, 1.17, // Ga-Kr
  1.89, 1.67, // Rb,Sr
  1.47, 1.39, 1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23, // Y-Cd
  1.28, 1.26, 1.26, 1.23, 1.32, 1.31, // In-Xe
  2.09, 1.76, // Cs,Ba
  1.62, 1.47, 1.58, 1.57, 1.56, 1.55, 1.51, // La-Eu
  1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53, // Gd-Yb
  1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32, // Lu-Hg
  1.30, 1.30, 1.36, 1.31, 1.38, 1.42, // Tl-Rn
  2.01, 1.81, // Fr,Ra
  1.67, 1.58, 1.52, 1.53, 1.54, 1.55, 1.49, // Ac-Am
  1.49, 1.51, 1.51, 1.48, 1.50, 1.56, 1.58, // Cm-No
  1.45, 1.41, 1.34, 1.29, 1.27, 1.21, 1.16, 1.15, 1.09, 1.22, // Lr-Cn
  1.22, 1.29, 1.46, 1.58, 1.48, 1.41 // Nh-Og
};

// Pauling electronegativities (1.5 dummies past element 86).
// float storage for the same single-precision-literal reason as rad above.
const std::array<float, maxElement> en = {
  2.20, 3.00, // H,He
  0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98, 4.50, // Li-Ne
  0.93, 1.31, 1.61, 1.90, 2.19, 2.58, 3.16, 3.50, // Na-Ar
  0.82, 1.00, // K,Ca
  1.36, 1.54, 1.63, 1.66, 1.55, 1.83, 1.88, 1.91, 1.90, 1.65, // Sc-Zn
  1.81, 2.01, 2.18, 2.55, 2.96, 3.00, // Ga-Kr
  0.82, 0.95, // Rb,Sr
  1.22, 1.33, 1.60, 2.16, 1.90, 2.20, 2.28, 2.20, 1.93, 1.69, // Y-Cd
  1.78, 1.96, 2.05, 2.10, 2.66, 2.60, // In-Xe
  0.79, 0.89, // Cs,Ba
  1.10, 1.12, 1.13, 1.14, 1.15, 1.17, 1.18, // La-Eu
  1.20, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, // Gd-Yb
  1.27, 1.30, 1.50, 2.36, 1.90, 2.20, 2.20, 2.28, 2.54, 2.00, // Lu-Hg
  1.62, 2.33, 2.02, 2.00, 2.20, 2.20, // Tl-Rn
  1.50, 1.50, // Fr,Ra
  1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, // Ac-Am
  1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, // Cm-No
  1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, // Rf-Cn
  1.50, 1.50, 1.50, 1.50, 1.50, 1.50 // Nh-Og
};

// Counting functions and their radial derivatives (erf_count, derf_count,
// exp_count, dexp_count in the original).
inline double expCount(double k, double r, double r0)
{
  return 1.0 / (1.0 + std::exp(-k * (r0 / r - 1.0)));
}

inline double dExpCount(double k, double r, double r0)
{
  double expterm = std::exp(-k * (r0 / r - 1.0));
  return (-k * r0 * expterm) / (r * r * (expterm + 1.0) * (expterm + 1.0));
}

inline double erfCount(double k, double r, double r0)
{
  return 0.5 * (1.0 + std::erf(-k * (r - r0) / r0));
}

inline double dErfCount(double k, double r, double r0)
{
  double x = k * (r - r0) / r0;
  return -k / sqrtpi / r0 * std::exp(-x * x);
}

inline double d4Denominator(int ia, int ja)
{
  // Widen before subtracting: the original keeps en in a double array, so
  // the difference is a double operation on single-rounded inputs.
  double den =
    std::abs(static_cast<double>(en[ia - 1]) - static_cast<double>(en[ja - 1])) +
    k5;
  return k4 * std::exp(-den * den / k6);
}

} // namespace

double covalentRadius(int atomicNumber)
{
  if (atomicNumber < 1 || atomicNumber > maxElement)
    return 0.0;
  return 4.0 / 3.0 * static_cast<double>(rad[atomicNumber - 1]) / 0.52917726;
}

double paulingElectronegativity(int atomicNumber)
{
  if (atomicNumber < 1 || atomicNumber > maxElement)
    return 0.0;
  return static_cast<double>(en[atomicNumber - 1]);
}

void coordinationD3(const std::vector<int>& numbers,
                    const std::vector<double>& xyz, std::vector<double>& cn,
                    double thr)
{
  const int nat = static_cast<int>(numbers.size());
  cn.assign(nat, 0.0);
  for (int i = 0; i < nat; ++i) {
    for (int j = 0; j < i; ++j) {
      double dx = xyz[3 * j] - xyz[3 * i];
      double dy = xyz[3 * j + 1] - xyz[3 * i + 1];
      double dz = xyz[3 * j + 2] - xyz[3 * i + 2];
      double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 > thr)
        continue;
      double r = std::sqrt(r2);
      double rco = covalentRadius(numbers[j]) + covalentRadius(numbers[i]);
      double tmp = expCount(k1, r, rco);
      cn[i] += tmp;
      cn[j] += tmp;
    }
  }
}

void coordinationD3Derivative(const std::vector<int>& numbers,
                              const std::vector<double>& xyz,
                              std::vector<double>& cn,
                              std::vector<double>& dcn, double thr)
{
  const int nat = static_cast<int>(numbers.size());
  cn.assign(nat, 0.0);
  dcn.assign(3 * nat * nat, 0.0);
  auto at = [&](int i, int j, int c) -> double& {
    return dcn[(i * nat + j) * 3 + c];
  };
  for (int i = 0; i < nat; ++i) {
    for (int j = 0; j < i; ++j) {
      double rij[3] = { xyz[3 * j] - xyz[3 * i],
                         xyz[3 * j + 1] - xyz[3 * i + 1],
                         xyz[3 * j + 2] - xyz[3 * i + 2] };
      double r2 = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];
      if (r2 > thr)
        continue;
      double r = std::sqrt(r2);
      double rcovij = covalentRadius(numbers[i]) + covalentRadius(numbers[j]);
      double tmp = expCount(k1, r, rcovij);
      double dtmp = dExpCount(k1, r, rcovij);
      cn[i] += tmp;
      cn[j] += tmp;
      for (int c = 0; c < 3; ++c) {
        // Verbatim sign pattern of dncoord_d3 (see the FIXME comments there).
        at(i, i, c) += -dtmp * rij[c] / r;
        at(j, j, c) += dtmp * rij[c] / r;
        at(i, j, c) = dtmp * rij[c] / r;
        at(j, i, c) = -dtmp * rij[c] / r;
      }
    }
  }
}

void coordinationErf(const std::vector<int>& numbers,
                     const std::vector<double>& xyz, std::vector<double>& cn,
                     double thr)
{
  const int nat = static_cast<int>(numbers.size());
  cn.assign(nat, 0.0);
  for (int i = 0; i < nat; ++i) {
    for (int j = 0; j < i; ++j) {
      double dx = xyz[3 * j] - xyz[3 * i];
      double dy = xyz[3 * j + 1] - xyz[3 * i + 1];
      double dz = xyz[3 * j + 2] - xyz[3 * i + 2];
      double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 > thr)
        continue;
      double r = std::sqrt(r2);
      double rco = covalentRadius(numbers[j]) + covalentRadius(numbers[i]);
      double tmp = erfCount(kn, r, rco);
      cn[i] += tmp;
      cn[j] += tmp;
    }
  }
}

void coordinationErfDerivative(const std::vector<int>& numbers,
                               const std::vector<double>& xyz,
                               std::vector<double>& cn,
                               std::vector<double>& dcn, double thr)
{
  const int nat = static_cast<int>(numbers.size());
  cn.assign(nat, 0.0);
  dcn.assign(3 * nat * nat, 0.0);
  auto at = [&](int i, int j, int c) -> double& {
    return dcn[(i * nat + j) * 3 + c];
  };
  for (int i = 0; i < nat; ++i) {
    for (int j = 0; j < i; ++j) {
      double rij[3] = { xyz[3 * j] - xyz[3 * i],
                         xyz[3 * j + 1] - xyz[3 * i + 1],
                         xyz[3 * j + 2] - xyz[3 * i + 2] };
      double r2 = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];
      if (r2 > thr)
        continue;
      double r = std::sqrt(r2);
      double rcovij = covalentRadius(numbers[i]) + covalentRadius(numbers[j]);
      double tmp = erfCount(kn, r, rcovij);
      double dtmp = dErfCount(kn, r, rcovij);
      cn[i] += tmp;
      cn[j] += tmp;
      for (int c = 0; c < 3; ++c) {
        // Verbatim sign pattern of dncoord_erf (differs from dncoord_d3;
        // see port notes — verified against the reference by parity test).
        at(i, i, c) += dtmp * rij[c] / r;
        at(j, j, c) += -dtmp * rij[c] / r;
        at(i, j, c) = dtmp * rij[c] / r;
        at(j, i, c) = -dtmp * rij[c] / r;
      }
    }
  }
}

void coordinationD4(const std::vector<int>& numbers,
                    const std::vector<double>& xyz, std::vector<double>& cn,
                    double thr)
{
  const int nat = static_cast<int>(numbers.size());
  cn.assign(nat, 0.0);
  for (int i = 0; i < nat; ++i) {
    for (int j = 0; j < i; ++j) {
      double dx = xyz[3 * j] - xyz[3 * i];
      double dy = xyz[3 * j + 1] - xyz[3 * i + 1];
      double dz = xyz[3 * j + 2] - xyz[3 * i + 2];
      double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 > thr)
        continue;
      double r = std::sqrt(r2);
      double rco = covalentRadius(numbers[j]) + covalentRadius(numbers[i]);
      double den = d4Denominator(numbers[i], numbers[j]);
      double tmp = den * erfCount(kn, r, rco);
      cn[i] += tmp;
      cn[j] += tmp;
    }
  }
}

void coordinationD4Derivative(const std::vector<int>& numbers,
                              const std::vector<double>& xyz,
                              std::vector<double>& cn,
                              std::vector<double>& dcn, double thr)
{
  const int nat = static_cast<int>(numbers.size());
  cn.assign(nat, 0.0);
  dcn.assign(3 * nat * nat, 0.0);
  auto at = [&](int i, int j, int c) -> double& {
    return dcn[(i * nat + j) * 3 + c];
  };
  for (int i = 0; i < nat; ++i) {
    for (int j = 0; j < i; ++j) {
      double rij[3] = { xyz[3 * j] - xyz[3 * i],
                         xyz[3 * j + 1] - xyz[3 * i + 1],
                         xyz[3 * j + 2] - xyz[3 * i + 2] };
      double r2 = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];
      if (r2 > thr)
        continue;
      double r = std::sqrt(r2);
      double rcovij = covalentRadius(numbers[i]) + covalentRadius(numbers[j]);
      double den = d4Denominator(numbers[i], numbers[j]);
      double tmp = den * erfCount(kn, r, rcovij);
      double dtmp = den * dErfCount(kn, r, rcovij);
      cn[i] += tmp;
      cn[j] += tmp;
      for (int c = 0; c < 3; ++c) {
        // Verbatim sign pattern of dncoord_d4.
        at(i, i, c) += -dtmp * rij[c] / r;
        at(j, j, c) += dtmp * rij[c] / r;
        at(i, j, c) = -dtmp * rij[c] / r;
        at(j, i, c) = dtmp * rij[c] / r;
      }
    }
  }
}

void coordinationGfn(const std::vector<int>& numbers,
                     const std::vector<double>& xyz, std::vector<double>& cn,
                     double thr)
{
  const int nat = static_cast<int>(numbers.size());
  cn.assign(nat, 0.0);
  for (int i = 0; i < nat; ++i) {
    for (int j = 0; j < i; ++j) {
      double dx = xyz[3 * j] - xyz[3 * i];
      double dy = xyz[3 * j + 1] - xyz[3 * i + 1];
      double dz = xyz[3 * j + 2] - xyz[3 * i + 2];
      double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 > thr)
        continue;
      double r = std::sqrt(r2);
      double rcovij = covalentRadius(numbers[i]) + covalentRadius(numbers[j]);
      double tmp = expCount(ka, r, rcovij) * expCount(kb, r, rcovij + rShift);
      cn[i] += tmp;
      cn[j] += tmp;
    }
  }
}

void coordinationGfnDerivative(const std::vector<int>& numbers,
                               const std::vector<double>& xyz,
                               std::vector<double>& cn,
                               std::vector<double>& dcn, double thr)
{
  const int nat = static_cast<int>(numbers.size());
  cn.assign(nat, 0.0);
  dcn.assign(3 * nat * nat, 0.0);
  auto at = [&](int i, int j, int c) -> double& {
    return dcn[(i * nat + j) * 3 + c];
  };
  for (int i = 0; i < nat; ++i) {
    for (int j = 0; j < i; ++j) {
      double rij[3] = { xyz[3 * j] - xyz[3 * i],
                         xyz[3 * j + 1] - xyz[3 * i + 1],
                         xyz[3 * j + 2] - xyz[3 * i + 2] };
      double r2 = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];
      if (r2 > thr)
        continue;
      double r = std::sqrt(r2);
      double rcovij = covalentRadius(numbers[i]) + covalentRadius(numbers[j]);
      // Product rule on 1/(1+e1)/(1+e2), verbatim from dncoord_gfn.
      double expterm1 = std::exp(-ka * (rcovij / r - 1.0));
      double expterm2 = std::exp(-kb * ((rcovij + rShift) / r - 1.0));
      double tmp = 1.0 / (1.0 + expterm1) / (1.0 + expterm2);
      double dtmp =
        (-ka * rcovij * expterm1) /
          (r2 * (expterm1 + 1.0) * (expterm1 + 1.0)) / (1.0 + expterm2) +
        (-kb * (rcovij + rShift) * expterm2) /
          (r2 * (expterm2 + 1.0) * (expterm2 + 1.0)) / (1.0 + expterm1);
      cn[i] += tmp;
      cn[j] += tmp;
      for (int c = 0; c < 3; ++c) {
        at(i, i, c) += -dtmp * rij[c] / r;
        at(j, j, c) += dtmp * rij[c] / r;
        at(i, j, c) = dtmp * rij[c] / r;
        at(j, i, c) = -dtmp * rij[c] / r;
      }
    }
  }
}

double logCoordinationCutoff(double cn, double max)
{
  return std::log(1.0 + std::exp(max)) - std::log(1.0 + std::exp(max - cn));
}

double dLogCoordinationCutoff(double cn, double max)
{
  return std::exp(max) / (std::exp(max) + std::exp(cn));
}

void applyLogCoordinationCutoff(std::vector<double>& cn,
                                std::vector<double>* dcn, int nat, double max)
{
  if (nat <= 0)
    nat = static_cast<int>(cn.size());
  if (dcn) {
    for (int i = 0; i < nat; ++i) {
      double scale = dLogCoordinationCutoff(cn[i], max);
      for (int j = 0; j < nat; ++j) {
        for (int c = 0; c < 3; ++c)
          (*dcn)[(i * nat + j) * 3 + c] *= scale;
      }
    }
  }
  for (int i = 0; i < nat; ++i)
    cn[i] = logCoordinationCutoff(cn[i], max);
}

} // namespace Xtb
} // namespace Avogadro
