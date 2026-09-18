/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "avogadro/xtb/COPYING.LESSER".
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/xtb/environment.h>
#include <avogadro/xtb/gfnffhuckel.h>

#include <cmath>
#include <vector>

namespace Avogadro {
namespace Xtb {

// Reference values for Aufbau/Fermi/density verified bitwise against the
// unmodified xtb routines (scc_core.f90 occu, fermismear, dmat) compiled
// with gfortran -ffp-contract=off; hex dumps matched exactly. Eigensolver
// and Hückel values cross-checked against numpy (eigvalsh/eigh) to 1e-14.

TEST(HuckelTest, aufbauEven)
{
  std::vector<double> fa, fb;
  int ha = 0, hb = 0;
  aufbauOccupations(6, 6, 0, fa, fb, ha, hb);
  EXPECT_EQ(ha, 3);
  EXPECT_EQ(hb, 3);
  for (int i = 0; i < 6; ++i) {
    double want = i < 3 ? 1.0 : 0.0;
    EXPECT_DOUBLE_EQ(fa[i], want);
    EXPECT_DOUBLE_EQ(fb[i], want);
  }
}

TEST(HuckelTest, aufbauOdd)
{
  std::vector<double> fa, fb;
  int ha = 0, hb = 0;
  aufbauOccupations(6, 5, 0, fa, fb, ha, hb);
  EXPECT_EQ(ha, 3);
  EXPECT_EQ(hb, 2);
  const double wantA[6] = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };
  const double wantB[6] = { 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
  for (int i = 0; i < 6; ++i) {
    EXPECT_DOUBLE_EQ(fa[i], wantA[i]);
    EXPECT_DOUBLE_EQ(fb[i], wantB[i]);
  }
}

TEST(HuckelTest, fermiSmearing)
{
  // Bitwise reference from the Fortran harness (levels in eV).
  const std::vector<double> eig = { -0.45, -0.30, -0.28, 0.15, 0.17, 0.40 };
  std::vector<double> occ(6, 0.0);
  double fod = 0.0, ef = 0.0, s = 0.0;
  fermiSmearing(6, 3, eig, 4000.0, occ, fod, ef, s);
  const double want[6] = {
    0.75939788905164607, 0.67133106362495387, 0.65840373274357644,
    0.35633705512355940, 0.34314261151814357, 0.21138764793811884
  };
  for (int i = 0; i < 6; ++i)
    EXPECT_DOUBLE_EQ(occ[i], want[i]);
  EXPECT_DOUBLE_EQ(fod, 1.8217346291596455);
  EXPECT_DOUBLE_EQ(ef, -0.053816325699250343);
  EXPECT_DOUBLE_EQ(s, -0.046074055805315171);
}

TEST(HuckelTest, fermiBoundaryLeavesOccUntouched)
{
  const std::vector<double> eig = { -0.45, -0.30, -0.28, 0.15, 0.17, 0.40 };
  // nel+1 > norbs: e_fermi = eig(nel), occ untouched.
  std::vector<double> occ(6, 0.5);
  double fod = 0.0, ef = 0.0, s = 0.0;
  fermiSmearing(6, 6, eig, 4000.0, occ, fod, ef, s);
  for (double v : occ)
    EXPECT_DOUBLE_EQ(v, 0.5);
  EXPECT_DOUBLE_EQ(ef, 0.40);
  // nel = 0: e_fermi = eig(1), occ untouched.
  fermiSmearing(6, 0, eig, 4000.0, occ, fod, ef, s);
  for (double v : occ)
    EXPECT_DOUBLE_EQ(v, 0.5);
  EXPECT_DOUBLE_EQ(ef, -0.45);
}

TEST(HuckelTest, densityMatrix)
{
  std::vector<double> coeffs(9, 0.0);
  coeffs[0] = coeffs[4] = coeffs[8] = 1.0;
  std::vector<double> focc = { 2.0, 1.0, 0.0 };
  std::vector<double> density;
  densityMatrix(coeffs, focc, density, 3);
  const double want[9] = { 2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
  for (int i = 0; i < 9; ++i)
    EXPECT_DOUBLE_EQ(density[i], want[i]);
}

TEST(HuckelTest, eigensystemRing)
{
  // Six-ring with hopping -1: eigenvalues -2,-1,-1,+1,+1,+2 exactly.
  const int n = 6;
  std::vector<double> mat(n * n, 0.0);
  for (int i = 0; i < n; ++i) {
    int j = (i + 1) % n;
    mat[i * n + j] = mat[j * n + i] = -1.0;
  }
  std::vector<double> evals, evecs;
  std::vector<double> work = mat; // eigensolver overwrites its input
  ASSERT_TRUE(symmetricEigensystem(work, n, evals, evecs));
  const double want[6] = { -2.0, -1.0, -1.0, 1.0, 1.0, 2.0 };
  for (int i = 0; i < n; ++i)
    EXPECT_NEAR(evals[i], want[i], 1e-12);
  // Residual ||A V - V L|| and orthonormality.
  for (int k = 0; k < n; ++k) {
    for (int i = 0; i < n; ++i) {
      double av = 0.0;
      for (int j = 0; j < n; ++j)
        av += mat[i * n + j] * evecs[j * n + k];
      EXPECT_NEAR(av, evals[k] * evecs[i * n + k], 1e-11);
    }
  }
  for (int k = 0; k < n; ++k) {
    for (int l = 0; l < n; ++l) {
      double dot = 0.0;
      for (int i = 0; i < n; ++i)
        dot += evecs[i * n + k] * evecs[i * n + l];
      EXPECT_NEAR(dot, k == l ? 1.0 : 0.0, 1e-12);
    }
  }
}

// Benzene pi system from the verified oracle run (numpy cross-check:
// density to 6.7e-15, HOMO/LUMO to 1e-12).
TEST(HuckelTest, benzene)
{
  const double angToBohr = 1.0 / 0.52917726;
  std::vector<HuckelAtom> atoms(6);
  for (int i = 0; i < 6; ++i) {
    atoms[i].element = 6;
    atoms[i].hyb = 2;
    atoms[i].charge = 0.0;
    atoms[i].piElectrons = 1;
  }
  std::vector<HuckelBond> bonds;
  for (int i = 0; i < 6; ++i) {
    HuckelBond b;
    b.first = i;
    b.second = (i + 1) % 6;
    b.distance = 1.39 * angToBohr;
    bonds.push_back(b);
  }
  std::vector<double> hdiag(86, 0.0), hoffdiag(86, 0.0);
  hoffdiag[5] = 1.0;
  Environment env;
  HuckelResult result;
  ASSERT_TRUE(huckelSolve(atoms, bonds, hdiag, -0.24, 0.53, hoffdiag, 0.70,
                          1.45, 5, 6, result, env))
    << env.errorMessage();
  EXPECT_EQ(result.electrons, 6);
  EXPECT_NEAR(result.homo, -2.7206654673464250, 1e-12);
  EXPECT_NEAR(result.lumo, 2.7206654673464232, 1e-12);
  ASSERT_EQ(result.density.size(), 36u);
  // Alternant-hydrocarbon pattern: P(ii) = 1, bonded 2/3, para -1/3.
  for (int i = 0; i < 6; ++i)
    EXPECT_NEAR(result.density[i * 6 + i], 1.0, 1e-12);
  EXPECT_NEAR(result.density[0 * 6 + 1], 0.66641777073276376, 1e-12);
  EXPECT_NEAR(result.density[0 * 6 + 3], -0.33283582023675029, 1e-12);
  EXPECT_NEAR(result.density[0 * 6 + 2], 0.0, 1e-12);
}

// High-lying HOMO triggers the second attempt with one electron less.
TEST(HuckelTest, secondAttempt)
{
  std::vector<HuckelAtom> atoms(2);
  for (int i = 0; i < 2; ++i) {
    atoms[i].element = 6;
    atoms[i].hyb = 2;
    atoms[i].charge = 0.0;
    atoms[i].piElectrons = 1;
  }
  HuckelBond b;
  b.first = 0;
  b.second = 1;
  b.distance = 2.0;
  std::vector<double> hdiag(86, 0.0), hoffdiag(86, 0.0);
  hdiag[5] = 5.0;
  hoffdiag[5] = 1.0;
  Environment env;
  HuckelResult result;
  ASSERT_TRUE(huckelSolve(atoms, std::vector<HuckelBond>{ b }, hdiag, -0.24,
                          0.53, hoffdiag, 0.70, 1.45, 5, 2, result, env))
    << env.errorMessage();
  EXPECT_EQ(result.electrons, 1);
  EXPECT_NEAR(result.homo, 10.2498509725650671, 1e-9);
  EXPECT_NEAR(result.lumo, 16.9615447274349300, 1e-9);
  ASSERT_EQ(result.density.size(), 4u);
  EXPECT_NEAR(result.density[0], 0.5, 1e-12);
  EXPECT_NEAR(result.density[1], 0.49994087254614761, 1e-12);
}

// Subsystem assembly + bond-order save, verified against an independent
// Python replica of the ini Hueckel block (numpy eigh; 43 pibo + 213
// pbo values over formaldehyde, formamide, benzene, nitromethane,
// pyridine and formate agree at 1e-9, covering C=O, amide-N, aromatic,
// nitro-N and odd-electron anion cases).

namespace {

void huckelTables(std::vector<double>& hdiag, std::vector<double>& hoffdiag)
{
  hdiag.assign(17, 0.0);
  hoffdiag.assign(17, 0.0);
  hdiag[4] = -0.5;
  hdiag[5] = 0.0;
  hdiag[6] = 0.14;
  hdiag[7] = -0.38;
  hdiag[8] = -0.29;
  hdiag[15] = -0.30;
  hdiag[16] = -0.30;
  hoffdiag[4] = 0.5;
  hoffdiag[5] = 1.00;
  hoffdiag[6] = 0.66;
  hoffdiag[7] = 1.10;
  hoffdiag[8] = 0.23;
  hoffdiag[15] = 0.60;
  hoffdiag[16] = 1.00;
}

} // namespace

TEST(HuckelPiTest, FormaldehydeBondOrder)
{
  Environment env;
  std::vector<int> numbers = { 6, 8, 1, 1 };
  std::vector<int> hyb = { 2, 2, 0, 0 };
  std::vector<int> itag = { 0, 0, 0, 0 };
  std::vector<double> qa = { 9.36651996035209045e-02, -3.26946970992318675e-01,
                             1.16640885694398830e-01, 1.16640885694398927e-01 };
  std::vector<double> xyz = { 0.0, 0.0, 0.0, 2.28, 0.0, 0.0, -1.03, 1.78, 0.0,
                              -1.03, -1.78, 0.0 };
  std::vector<std::pair<int, int>> blist = { { 1, 0 }, { 2, 0 }, { 3, 0 } };
  std::vector<int> piIndex = { 0, 1 };
  std::vector<int> fragOfPi = { 1, 1 };
  std::vector<double> hdiag, hoffdiag;
  huckelTables(hdiag, hoffdiag);
  std::vector<double> pibo(3, 0.0), pbo(10, 0.0);
  std::vector<int> piadrOut(4, 0);
  ASSERT_TRUE(huckelPiBondOrders(4, numbers, hyb, itag, qa, xyz, blist,
                                 piIndex, fragOfPi, 1, hdiag, -0.24, 0.53,
                                 hoffdiag, 0.70, 1.45, 5.0, pibo, pbo,
                                 piadrOut, env));
  EXPECT_NEAR(pibo[0], 9.94123309805507782e-01, 1e-9);
  EXPECT_NEAR(pibo[1], 0.0, 0.0);
  EXPECT_NEAR(pibo[2], 0.0, 0.0);
  EXPECT_NEAR(pbo[1], 9.94123309805507782e-01, 1e-9);
  EXPECT_EQ(piadrOut[0], 1);
  EXPECT_EQ(piadrOut[1], 1);
  EXPECT_EQ(piadrOut[2], 0);
  EXPECT_EQ(piadrOut[3], 0);
}

TEST(HuckelPiTest, BenzeneAromaticBondOrders)
{
  Environment env;
  std::vector<int> numbers(12, 1);
  for (int i = 0; i < 6; ++i)
    numbers[i] = 6;
  std::vector<int> hyb(12, 0);
  for (int i = 0; i < 6; ++i)
    hyb[i] = 2;
  std::vector<int> itag(12, 0);
  std::vector<double> qa = { -3.43066958971704680e-02,
                             -3.43066958971700517e-02,
                             -3.43066958971705582e-02,
                             -3.43066958971711758e-02,
                             -3.43066958971712868e-02,
                             -3.43066958971708427e-02,
                             3.43066958971710440e-02,
                             3.43066958971691496e-02,
                             3.43066958971710231e-02,
                             3.43066958971712452e-02,
                             3.43066958971732436e-02,
                             3.43066958971686500e-02 };
  std::vector<double> xyz;
  for (int k = 0; k < 6; ++k) {
    double a = k * 1.04719755119659763;
    xyz.push_back(2.64 * std::cos(a));
    xyz.push_back(2.64 * std::sin(a));
    xyz.push_back(0.0);
  }
  for (int k = 0; k < 6; ++k) {
    double a = k * 1.04719755119659763;
    xyz.push_back(4.70 * std::cos(a));
    xyz.push_back(4.70 * std::sin(a));
    xyz.push_back(0.0);
  }
  std::vector<std::pair<int, int>> blist = { { 1, 0 }, { 5, 0 }, { 6, 0 },
                                             { 2, 1 }, { 7, 1 }, { 3, 2 },
                                             { 8, 2 }, { 4, 3 }, { 9, 3 },
                                             { 5, 4 }, { 10, 4 }, { 11, 5 } };
  std::vector<int> piIndex = { 0, 1, 2, 3, 4, 5 };
  std::vector<int> fragOfPi = { 1, 1, 1, 1, 1, 1 };
  std::vector<double> hdiag, hoffdiag;
  huckelTables(hdiag, hoffdiag);
  std::vector<double> pibo(12, 0.0), pbo(78, 0.0);
  std::vector<int> piadrOut(12, 0);
  ASSERT_TRUE(huckelPiBondOrders(12, numbers, hyb, itag, qa, xyz, blist,
                                 piIndex, fragOfPi, 1, hdiag, -0.24, 0.53,
                                 hoffdiag, 0.70, 1.45, 5.0, pibo, pbo,
                                 piadrOut, env));
  // Aromatic C-C orders are 2/3; C-H bonds carry none.
  EXPECT_NEAR(pibo[0], 6.66417770732740555e-01, 1e-9);
  EXPECT_NEAR(pibo[1], 6.66417770732740000e-01, 1e-9);
  EXPECT_NEAR(pibo[2], 0.0, 0.0);
  for (int i = 0; i < 6; ++i)
    EXPECT_EQ(piadrOut[i], 1);
  for (int i = 6; i < 12; ++i)
    EXPECT_EQ(piadrOut[i], 0);
}

} // namespace Xtb
} // namespace Avogadro
