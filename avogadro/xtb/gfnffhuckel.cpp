/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/scc_core.f90 (occu, fermismear, dmat), src/gfnff/gfnff_qm.f90
  (gfnffqmsolve, ZDO path only) and src/gfnff/gfnff_ini.f90 (Hückel loop),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#include "gfnffhuckel.h"

#include "constants.h"
#include "environment.h"
#include "gfnffegbond.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace Avogadro {
namespace Xtb {

bool symmetricEigensystem(std::vector<double>& mat, int n,
                          std::vector<double>& eigenvalues,
                          std::vector<double>& eigenvectors)
{
  if (n <= 0)
    return false;
  std::vector<double> a = mat;
  eigenvectors.assign(n * n, 0.0);
  for (int i = 0; i < n; ++i)
    eigenvectors[i * n + i] = 1.0;
  // Cyclic Jacobi rotations.
  for (int sweep = 0; sweep < 100; ++sweep) {
    double off = 0.0;
    for (int p = 0; p < n; ++p) {
      for (int q = p + 1; q < n; ++q)
        off += a[p * n + q] * a[p * n + q];
    }
    if (off <= 1.0e-24)
      break;
    for (int p = 0; p < n - 1; ++p) {
      for (int q = p + 1; q < n; ++q) {
        double apq = a[p * n + q];
        if (std::abs(apq) < 1.0e-18)
          continue;
        double app = a[p * n + p], aqq = a[q * n + q];
        double theta = (aqq - app) / (2.0 * apq);
        double t = (theta >= 0.0 ? 1.0 : -1.0) /
                   (std::abs(theta) + std::sqrt(theta * theta + 1.0));
        double c = 1.0 / std::sqrt(t * t + 1.0);
        double s = t * c;
        for (int k = 0; k < n; ++k) {
          double akp = a[k * n + p], akq = a[k * n + q];
          a[k * n + p] = c * akp - s * akq;
          a[k * n + q] = s * akp + c * akq;
        }
        for (int k = 0; k < n; ++k) {
          double apk = a[p * n + k], aqk = a[q * n + k];
          a[p * n + k] = c * apk - s * aqk;
          a[q * n + k] = s * apk + c * aqk;
        }
        for (int k = 0; k < n; ++k) {
          double vkp = eigenvectors[k * n + p];
          double vkq = eigenvectors[k * n + q];
          eigenvectors[k * n + p] = c * vkp - s * vkq;
          eigenvectors[k * n + q] = s * vkp + c * vkq;
        }
      }
    }
  }
  eigenvalues.resize(n);
  for (int i = 0; i < n; ++i)
    eigenvalues[i] = a[i * n + i];
  // Ascending order like dsyev, permuting the eigenvectors along.
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(),
            [&](int x, int y) { return eigenvalues[x] < eigenvalues[y]; });
  std::vector<double> sorted(n), sortedVec(n * n);
  for (int i = 0; i < n; ++i) {
    sorted[i] = eigenvalues[idx[i]];
    for (int k = 0; k < n; ++k)
      sortedVec[k * n + i] = eigenvectors[k * n + idx[i]];
  }
  eigenvalues = sorted;
  eigenvectors = sortedVec;
  mat = a;
  return true;
}

void aufbauOccupations(int ndim, int nel, int nopen,
                       std::vector<double>& foccA,
                       std::vector<double>& foccB, int& ihomoa, int& ihomob)
{
  std::vector<double> foc(ndim, 0.0);
  if (nel % 2 == 0) {
    int ihomo = nel / 2;
    for (int i = 0; i < ihomo && i < ndim; ++i)
      foc[i] = 2.0;
    if (2 * ihomo != nel) {
      // Unreachable for even nel (kept for parity with the reference).
      ihomo = ihomo + 1;
      if (ihomo - 1 < ndim)
        foc[ihomo - 1] = 1.0;
    }
    if (nopen > 1) {
      for (int i = 0; i < nopen / 2; ++i) {
        foc[ihomo - i - 1] -= 1.0;
        foc[ihomo + i] += 1.0;
      }
    }
  } else {
    int na = nel / 2 + (nopen - 1) / 2 + 1;
    int nb = nel / 2 - (nopen - 1) / 2;
    for (int i = 0; i < na && i < ndim; ++i)
      foc[i] += 1.0;
    for (int i = 0; i < nb && i < ndim; ++i)
      foc[i] += 1.0;
  }
  foccA.assign(ndim, 0.0);
  foccB.assign(ndim, 0.0);
  for (int i = 0; i < ndim; ++i) {
    if (foc[i] == 2.0) {
      foccA[i] = 1.0;
      foccB[i] = 1.0;
    }
    if (foc[i] == 1.0)
      foccA[i] = 1.0;
  }
  ihomoa = 0;
  ihomob = 0;
  for (int i = 0; i < ndim; ++i) {
    if (foccA[i] > 0.99)
      ihomoa = i + 1;
    if (foccB[i] > 0.99)
      ihomob = i + 1;
  }
}

void fermiSmearing(int norbs, int nel, const std::vector<double>& eig,
                   double t, std::vector<double>& occ, double& fod,
                   double& fermiEnergy, double& entropy)
{
  // boltz = kB * autoev with autoev = 27.21138505 (mctc_convert).
  const double boltz = kB * 27.21138505;
  const double thr = 1e-9;
  const double sqrttiny = std::sqrt(std::numeric_limits<double>::min());
  double bkt = boltz * t;
  double efermi = 0.0;
  if (nel + 1 > norbs) {
    efermi = eig[nel - 1];
  } else if (nel == 0) {
    efermi = eig[0];
  } else {
    efermi = 0.5 * (eig[nel - 1] + eig[nel]);
    double occt = nel;
    for (int cycle = 0; cycle < 200; ++cycle) {
      double total = 0.0, dtotal = 0.0;
      for (int i = 0; i < norbs; ++i) {
        double fermifunct = 0.0, dfermi = 0.0;
        if ((eig[i] - efermi) / bkt < 50) {
          double e = std::exp((eig[i] - efermi) / bkt);
          fermifunct = 1.0 / (e + 1.0);
          dfermi = e / (bkt * (e + 1.0) * (e + 1.0));
        }
        occ[i] = fermifunct;
        total += fermifunct;
        dtotal += dfermi;
      }
      double change = 0.0;
      if (dtotal > sqrttiny)
        change = (occt - total) / dtotal;
      // Redundant re-evaluation kept verbatim from the reference.
      change = (occt - total) / dtotal;
      efermi += change;
      if (std::abs(occt - total) <= thr)
        break;
    }
  }
  // Boundary branches leave occ untouched, like the reference.
  fod = 0.0;
  double s = 0.0;
  for (int i = 0; i < norbs; ++i) {
    if (occ[i] > thr && 1.0 - occ[i] > thr) {
      // Left-associated like the reference s=s+A+B.
      s = s + occ[i] * std::log(occ[i]);
      s = s + (1.0 - occ[i]) * std::log(1.0 - occ[i]);
    }
    if (eig[i] < efermi)
      fod = fod + 1.0 - occ[i]; // (fod+1)-occ like the reference
    else
      fod = fod + occ[i];
  }
  s *= kB * t;
  fermiEnergy = efermi;
  entropy = s;
}

void densityMatrix(const std::vector<double>& coeffs,
                   const std::vector<double>& focc,
                   std::vector<double>& density, int n)
{
  density.assign(n * n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int m = 0; m < n; ++m) {
      double sum = 0.0;
      for (int k = 0; k < n; ++k)
        sum += coeffs[i * n + k] * focc[k] * coeffs[m * n + k];
      density[i * n + m] = sum;
    }
  }
}

// Single ZDO solve mirroring gfnffqmsolve() with ovlp=.false.:
// eigensystem of the Hamiltonian in api, eigenvalue scale, Fermi
// occupations (occu + fermismear, et > 1e-3 always here), electronic
// energy eel, and the density matrix written back into api.
static bool zdoSolve(std::vector<double>& api, int npi, int nel, double et,
                     double& eel, std::vector<double>& focc,
                     std::vector<double>& eig, Environment& env)
{
  std::vector<double> coeffs;
  if (!symmetricEigensystem(api, npi, eig, coeffs)) {
    env.error("diagonalization failed", "huckelSolve");
    return false;
  }
  for (int i = 0; i < npi; ++i)
    eig[i] *= 0.1 * 27.2113957;
  int ihomoa = 0, ihomob = 0;
  std::vector<double> fa(npi, 0.0), fb(npi, 0.0);
  aufbauOccupations(npi, nel, 0, fa, fb, ihomoa, ihomob);
  double fod = 0.0, efa = 0.0, ga = 0.0;
  double fodb = 0.0, efb = 0.0, gb = 0.0;
  if (ihomoa >= 1 && ihomoa <= npi)
    fermiSmearing(npi, ihomoa, eig, et, fa, fod, efa, ga);
  else
    std::fill(fa.begin(), fa.end(), 0.0);
  if (ihomob >= 1 && ihomob <= npi)
    fermiSmearing(npi, ihomob, eig, et, fb, fodb, efb, gb);
  else
    std::fill(fb.begin(), fb.end(), 0.0);
  focc.assign(npi, 0.0);
  for (int i = 0; i < npi; ++i)
    focc[i] = fa[i] + fb[i];
  if (ihomoa + 1 <= npi) {
    if (std::abs(focc[ihomoa - 1] - focc[ihomoa]) < 1.0e-4) {
      // Perfect biradical: break symmetry like the reference
      // (reference also prints a note; kept silent in the library).
      std::fill(focc.begin(), focc.end(), 0.0);
      for (int i = 0; i < nel / 2 && i < npi; ++i)
        focc[i] = 2.0;
    }
  }
  eel = 0.0;
  for (int i = 0; i < npi; ++i)
    eel += focc[i] * eig[i];
  std::vector<double> density;
  densityMatrix(coeffs, focc, density, npi);
  api = density;
  return true;
}

// HOMO/LUMO scan mirroring the pisip/pisea loop: last level with
// occ > 0.5 is the HOMO, the next level the LUMO.
static void homoLumo(const std::vector<double>& focc,
                     const std::vector<double>& eig, double& homo,
                     double& lumo)
{
  int npi = static_cast<int>(eig.size());
  homo = -1e300;
  lumo = 1e300;
  for (int i = 0; i < npi; ++i) {
    if (focc[i] > 0.5) {
      homo = eig[i];
      if (i + 1 < npi)
        lumo = eig[i + 1];
    }
  }
}

bool huckelSolve(const std::vector<HuckelAtom>& atoms,
                 const std::vector<HuckelBond>& bonds,
                 const std::vector<double>& hdiag, double hueckelP3,
                 double pilpf, const std::vector<double>& hoffdiag,
                 double hiter, double htriple, int maxIter, int nelElectrons,
                 HuckelResult& result, Environment& env)
{
  int npi = static_cast<int>(atoms.size());
  if (npi < 2 || nelElectrons < 1) {
    env.error("Hückel system too small", "huckelSolve");
    return false;
  }
  const double et = 4000.0; // electronic temperature for the smearing
  std::vector<double> api(npi * npi, 0.0), apisave(npi * npi, 0.0),
    pold(npi * npi, 2.0 / 3.0);
  std::vector<double> eig(npi), focc(npi, 0.0);
  double eel = 0.0, eold = 0.0, homo = 0.0, lumo = 0.0;
  int nel = nelElectrons;
  for (int nn = 0; nn < maxIter; ++nn) {
    std::fill(api.begin(), api.end(), 0.0);
    for (int i = 0; i < npi; ++i) {
      api[i * npi + i] =
        hdiag[atoms[i].element - 1] + atoms[i].charge * hueckelP3 -
        (atoms[i].piElectrons - 1) * pilpf;
    }
    for (const HuckelBond& bond : bonds) {
      int ia = bond.first, ja = bond.second;
      double dum = std::sqrt(hoffdiag[atoms[ia].element - 1] *
                             hoffdiag[atoms[ja].element - 1]) -
                   1.0e-9 * bond.distance;
      double dum2 = hiter;
      if (atoms[ia].hyb == 1)
        dum2 *= htriple;
      if (atoms[ja].hyb == 1)
        dum2 *= htriple;
      double v = -dum * (1.0 - dum2 * (2.0 / 3.0 - pold[ja * npi + ia]));
      api[ja * npi + ia] = v;
      api[ia * npi + ja] = v;
    }
    apisave = api;
    if (!zdoSolve(api, npi, nel, et, eel, focc, eig, env))
      return false;
    homoLumo(focc, eig, homo, lumo);
    if (std::abs(eel - eold) < 1.0e-4)
      break;
    pold = api; // api now holds the density matrix
    eold = eel;
  }
  if (homo > 0.40) {
    // Second attempt with one electron less on the saved Hamiltonian.
    nel = nelElectrons - 1;
    api = apisave;
    if (!zdoSolve(api, npi, nel, et, eel, focc, eig, env))
      return false;
    homoLumo(focc, eig, homo, lumo);
  }
  result.density = api;
  result.homo = homo;
  result.lumo = lumo;
  result.electrons = nel;
  return true;
}

bool huckelPiBondOrders(int n, const std::vector<int>& numbers,
                        const std::vector<int>& hyb,
                        const std::vector<int>& itag,
                        const std::vector<double>& charges,
                        const std::vector<double>& xyz,
                        const std::vector<std::pair<int, int>>& blist,
                        const std::vector<int>& piIndex,
                        const std::vector<int>& fragOfPi, int piFragments,
                        const std::vector<double>& hdiag, double hueckelP3,
                        double pilpf, const std::vector<double>& hoffdiag,
                        double hiter, double htriple, double maxHIter,
                        std::vector<double>& pibo, std::vector<double>& pbo,
                        std::vector<int>& piadrOut, Environment& env)
{
  int nbond = static_cast<int>(blist.size());
  if (n <= 0 || static_cast<int>(numbers.size()) != n || nbond <= 0) {
    env.error("empty Hueckel setup", "huckelPiBondOrders");
    return false;
  }
  int maxIter = static_cast<int>(std::lround(maxHIter));
  for (int pis = 1; pis <= piFragments; ++pis) {
    // Electron counting over this subsystem's pi atoms (cumulative
    // element/hyb/itag rules like the reference).
    int npi = 0, nelpi = 0;
    std::vector<int> piadr3, piels;
    for (size_t k = 0; k < piIndex.size(); ++k) {
      if (fragOfPi[k] != pis)
        continue;
      int atom = piIndex[k];
      int ati = numbers[atom];
      int hybi = hyb[atom];
      int before = nelpi;
      if (ati == 5 && hybi == 1)
        ++nelpi;
      if (ati == 6 && itag[atom] != 1)
        ++nelpi;
      if (ati == 7 && hybi == 2 && itag[atom] == 1)
        ++nelpi;
      if (ati == 7 && hybi <= 2)
        ++nelpi;
      if (ati == 7 && hybi == 3)
        nelpi += 2;
      if (ati == 8 && hybi == 1)
        ++nelpi;
      if (ati == 8 && hybi == 2)
        ++nelpi;
      if (ati == 8 && hybi == 3)
        nelpi += 2;
      if (ati == 9 && hybi != 1)
        nelpi += 2;
      if (ati == 9 && hybi == 1)
        nelpi += 3;
      if (ati == 16 && hybi == 1)
        ++nelpi;
      if (ati == 16 && hybi == 2)
        ++nelpi;
      if (ati == 16 && hybi == 3)
        nelpi += 2;
      if (ati == 17 && hybi == 0)
        nelpi += 2;
      if (ati == 17 && hybi == 1)
        nelpi += 3;
      int piel = nelpi - before;
      if (piel > 2)
        piel = 2;
      piadr3.push_back(atom);
      piels.push_back(piel);
      ++npi;
    }
    if (npi < 2 || nelpi < 1)
      continue;
    // Slot map for this subsystem (piadr4 equivalent, 0-based, -1 out).
    std::vector<int> slotMap(n, -1);
    for (int s = 0; s < npi; ++s)
      slotMap[piadr3[s]] = s;
    std::vector<HuckelAtom> atoms(npi);
    for (int s = 0; s < npi; ++s) {
      int atom = piadr3[s];
      atoms[s].element = numbers[atom];
      atoms[s].hyb = hyb[atom];
      atoms[s].charge = charges[atom];
      atoms[s].piElectrons = piels[s];
      atoms[s].tagged = itag[atom] == 1;
    }
    std::vector<HuckelBond> hbonds;
    for (const auto& b : blist) {
      int ja = slotMap[b.first], ia = slotMap[b.second];
      if (ia < 0 || ja < 0)
        continue;
      double dx = xyz[3 * b.second] - xyz[3 * b.first];
      double dy = xyz[3 * b.second + 1] - xyz[3 * b.first + 1];
      double dz = xyz[3 * b.second + 2] - xyz[3 * b.first + 2];
      HuckelBond hb;
      hb.first = ja;
      hb.second = ia;
      double d[3] = { dx, dy, dz };
      hb.distance = norm2(d);
      hbonds.push_back(hb);
    }
    HuckelResult result;
    if (!huckelSolve(atoms, hbonds, hdiag, hueckelP3, pilpf, hoffdiag,
                     hiter, htriple, maxIter, nelpi, result, env))
      return false;
    // Bond-order save over all blist bonds (reference order).
    for (int i = 0; i < nbond; ++i) {
      int ja = slotMap[blist[i].first], ia = slotMap[blist[i].second];
      if (ia < 0 || ja < 0)
        continue;
      double bo = result.density[ja * npi + ia];
      pibo[i] = bo;
      int a = blist[i].first, b = blist[i].second;
      int lo = a < b ? a : b, hi = a < b ? b : a;
      pbo[hi * (hi + 1) / 2 + lo] = bo;
      piadrOut[a] = 1;
      piadrOut[b] = 1;
    }
  }
  return true;
}

} // namespace Xtb
} // namespace Avogadro
