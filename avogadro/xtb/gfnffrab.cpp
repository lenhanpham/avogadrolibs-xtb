/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_rab.f (gfnffrab, iTabRow6),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#include "gfnffrab.h"

#include "gfnffparams.h"
#include "gfnfftopo.h"

#include <array>
#include <cmath>

namespace Avogadro {
namespace Xtb {

int elementRow6(int atomicNumber)
{
  if (atomicNumber > 0 && atomicNumber <= 2)
    return 1;
  if (atomicNumber > 2 && atomicNumber <= 10)
    return 2;
  if (atomicNumber > 10 && atomicNumber <= 18)
    return 3;
  if (atomicNumber > 18 && atomicNumber <= 36)
    return 4;
  if (atomicNumber > 36 && atomicNumber <= 54)
    return 5;
  if (atomicNumber > 54)
    return 6;
  return 0;
}

void gfnffBondGuesses(int n, const std::vector<int>& numbers,
                      const std::vector<double>& cn,
                      std::vector<double>& rabPacked)
{
  // p table and scale factors: unsuffixed literals in the original, hence
  // single precision widened on use.
  const float p[6][2] = { { 29.84522887f, -8.87843763f },
                          { -1.70549806f, 2.10878369f },
                          { 6.54013762f, 0.08009374f },
                          { 6.39169003f, -0.85808076f },
                          { 6.00000000f, -1.15000000f },
                          { 5.60000000f, -1.30000000f } };
  float scaleF[103][103];
  for (int i = 0; i < 103; ++i) {
    for (int j = 0; j < 103; ++j)
      scaleF[i][j] = 1.0f;
  }
  for (int z = 57; z <= 71; ++z) {
    scaleF[z - 1][9 - 1] = 0.93533568f;
    scaleF[9 - 1][z - 1] = 0.93533568f;
    scaleF[z - 1][17 - 1] = 1.0190114f;
    scaleF[17 - 1][z - 1] = 1.0190114f;
    scaleF[z - 1][35 - 1] = 1.0425532f;
    scaleF[35 - 1][z - 1] = 1.0425532f;
    scaleF[z - 1][53 - 1] = 1.0551948f;
    scaleF[53 - 1][z - 1] = 1.0551948f;
  }
  for (int z = 89; z <= 103; ++z) {
    scaleF[z - 1][9 - 1] = 0.93533568f;
    scaleF[9 - 1][z - 1] = 0.93533568f;
    scaleF[z - 1][17 - 1] = 1.0190114f;
    scaleF[17 - 1][z - 1] = 1.0190114f;
    scaleF[z - 1][35 - 1] = 1.0425532f;
    scaleF[35 - 1][z - 1] = 1.0425532f;
    scaleF[z - 1][53 - 1] = 1.0551948f;
    scaleF[53 - 1][z - 1] = 1.0551948f;
  }

  rabPacked.assign(n * (n + 1) / 2, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      int ati = numbers[i], atj = numbers[j];
      int ir = elementRow6(ati) - 1, jr = elementRow6(atj) - 1;
      double ra =
        static_cast<double>(rabR0[ati - 1]) + static_cast<double>(rabCnfak[ati - 1]) * cn[i];
      double rb =
        static_cast<double>(rabR0[atj - 1]) + static_cast<double>(rabCnfak[atj - 1]) * cn[j];
      double den = std::abs(static_cast<double>(rabEn[ati - 1]) -
                            static_cast<double>(rabEn[atj - 1]));
      // Widen before adding: the reference keeps p in double arrays, so the
      // sums are double operations on single-rounded inputs.
      double k1 =
        0.005 * (static_cast<double>(p[ir][0]) + static_cast<double>(p[jr][0]));
      double k2 =
        0.005 * (static_cast<double>(p[ir][1]) + static_cast<double>(p[jr][1]));
      double ff = 1.0 - k1 * den - k2 * den * den;
      rabPacked[packedIndex(j, i)] =
        (ra + rb) * ff * scaleF[ati - 1][atj - 1];
    }
  }
}

void scaleBondGuesses(int n, const std::vector<int>& numbers,
                      const std::vector<double>& qa,
                      const std::vector<int>& metal, double fq,
                      std::vector<double>& rabPacked)
{
  // Element specials from gfnff_neigh (fat array, 1-based elements).
  auto fat = [](int z) {
    switch (z) {
      case 1: return 1.02;
      case 4: return 1.03;
      case 5: return 1.02;
      case 8: return 1.02;
      case 9: return 1.05;
      case 10: return 1.10;
      case 11: return 1.01;
      case 12: return 1.02;
      case 15: return 0.97;
      case 18: return 1.10;
      case 19: return 1.02;
      case 20: return 1.02;
      case 38: return 1.02;
      case 34: return 0.99;
      case 50: return 1.01;
      case 51: return 0.99;
      case 52: return 0.95;
      case 53: return 0.98;
      case 56: return 1.02;
      case 76: return 1.02;
      case 82: return 1.06;
      case 83: return 0.95;
      default: return 1.0;
    }
  };
  for (int i = 0; i < n; ++i) {
    int ai = numbers[i];
    double f1 = fq;
    if (metal[ai - 1] > 0)
      f1 *= 2.0;
    for (int j = 0; j <= i; ++j) {
      int aj = numbers[j];
      double f2 = fq;
      if (metal[aj - 1] > 0)
        f2 *= 2.0;
      int k = packedIndex(j, i);
      double rco = rabPacked[k];
      rabPacked[k] =
        (rco - qa[i] * f1 - qa[j] * f2) * fat(ai) * fat(aj);
    }
  }
}

} // namespace Xtb
} // namespace Avogadro
