/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/type/molecule.f90 (TMolecule), Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#ifndef AVOGADRO_XTB_MOLECULE_H
#define AVOGADRO_XTB_MOLECULE_H

#include <array>
#include <vector>

namespace Avogadro {
namespace Xtb {

// Minimal molecular structure container, mirroring xtb's TMolecule fields
// used by the ported calculators. Coordinates are in Bohr, stored flat
// [3 * n], matching the Fortran xyz(3, nat) layout.
struct Molecule
{
  std::vector<int> numbers; // 1-based atomic numbers
  std::vector<double> positions; // Cartesian coordinates in Bohr
  int charge = 0;
  int unpairedElectrons = 0;
  bool periodic = false;
  std::array<double, 9> lattice = { 0.0, 0.0, 0.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0, 0.0 }; // Bohr, row-major

  int atomCount() const { return static_cast<int>(numbers.size()); }
  bool valid() const
  {
    return positions.size() == 3 * numbers.size();
  }
};

// Calculation results, mirroring xtb's scc_results fields consumed by the
// ported code. Energies in Hartree, gradients in Hartree/Bohr.
struct Results
{
  double energy = 0.0;
  std::vector<double> gradient; // flat [3 * n]
  std::vector<double> charges; // partial charges per atom
  std::array<double, 3> dipole = { 0.0, 0.0, 0.0 }; // e*Bohr
  double homoLumoGap = 0.0; // Hartree
  bool converged = false;
};

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_MOLECULE_H
