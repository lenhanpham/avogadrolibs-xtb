/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/scc_core.f90 (occu, fermismear, dmat), src/gfnff/gfnff_qm.f90
  (gfnffqmsolve, ZDO path only) and src/gfnff/gfnff_ini.f90 (Hückel loop),
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#ifndef AVOGADRO_XTB_GFNFFHUCKEL_H
#define AVOGADRO_XTB_GFNFFHUCKEL_H

#include <vector>

namespace Avogadro {
namespace Xtb {

class Environment;

// Symmetric eigensystem by cyclic Jacobi rotations, replacing LAPACK dsyev
// (unavailable in standalone builds; the real build should use Eigen).
// Eigenvalues ascending with matching orthonormal eigenvectors in columns.
// Row-major matrices throughout.
bool symmetricEigensystem(std::vector<double>& mat, int n,
                          std::vector<double>& eigenvalues,
                          std::vector<double>& eigenvectors);

// Aufbau occupations mirroring occu(): integer alpha/beta weight vectors
// with ihomoa/ihomob = last index above 0.99 (1-based, 0 when empty).
void aufbauOccupations(int ndim, int nel, int nopen,
                       std::vector<double>& foccA,
                       std::vector<double>& foccB, int& ihomoa, int& ihomob);

// Fermi smearing mirroring fermismear() at temperature t (Kelvin):
// Newton search for the Fermi level (max 200 steps), occupations, number
// of fractional_deps electrons (fod) and electronic entropy (s). Keeps the
// reference's redundant Newton-step evaluation verbatim.
void fermiSmearing(int norbs, int nel, const std::vector<double>& eig,
                   double t, std::vector<double>& occ, double& fod,
                   double& fermiEnergy, double& entropy);

// Density matrix P = C diag(focc) C^T mirroring dmat().
void densityMatrix(const std::vector<double>& coeffs,
                   const std::vector<double>& focc,
                   std::vector<double>& density, int n);

// One Hückel pi system: iterative solve with P-dependent off-diagonal
// damping, mirroring the gfnff_ini loop (maxIter = nint(maxhiter)).
// atoms holds per-pi-atom data; bonds holds bonded pi index pairs with
// distances in Bohr; gen bundles hdiag/hueckelp3/pilpf/hoffdiag/hiter/
// htriple. nelElectrons is the (ipis-adjusted) electron count. Outputs the
// density matrix (row-major, for pibo extraction), HOMO/LUMO in eV and the
// electron count used (after a possible second attempt with one less).
struct HuckelAtom
{
  int element = 0;
  int hyb = 0;
  double charge = 0.0;
  int piElectrons = 0; // piel, capped at 2
  bool tagged = false; // itag == 1
};
struct HuckelBond
{
  int first = -1; // pi-slot indices
  int second = -1;
  double distance = 0.0; // Bohr
};
struct HuckelResult
{
  std::vector<double> density;
  double homo = 0.0;
  double lumo = 0.0;
  int electrons = 0;
};
bool huckelSolve(const std::vector<HuckelAtom>& atoms,
                 const std::vector<HuckelBond>& bonds,
                 const std::vector<double>& hdiag, double hueckelP3,
                 double pilpf, const std::vector<double>& hoffdiag,
                 double hiter, double htriple, int maxIter, int nelElectrons,
                 HuckelResult& result, Environment& env);

// Pi bond orders from the gfnff_ini iterative-Hueckel section (0d):
// per-subsystem electron counting (B/C/N/O/F/S/Cl hyb+itag rules with
// piel clamped to 2; ipis is 0 in single-shot setup since the fragment
// charge estimation only runs on qloop repeats), one huckelSolve per
// subsystem with npi >= 2 and nelpi >= 1, then bond-order save into
// pibo (per blist bond) and packed pbo plus the post-Hueckel pi
// participation flags (itmp: atoms in at least one pi-mapped bond).
// blist holds (jj, ii) endpoint pairs in reference order; piIndex and
// fragOfPi (1-based pimvec) come from buildPiSystem; charges are the
// setup (goedeckera) charges; hdiag/
// hoffdiag are 0-based element tables (Z - 1); maxHIter is rounded like
// nint(). pibo/pbo/piadrOut are zero-initialized by the caller and only
// filled for solved subsystems.
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
                        std::vector<int>& piadrOut, Environment& env);

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_GFNFFHUCKEL_H
