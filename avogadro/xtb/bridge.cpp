/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.
******************************************************************************/

#include "bridge.h"

#include "constants.h"

#include <avogadro/core/molecule.h>
#include <avogadro/core/vector.h>

namespace Avogadro {
namespace Xtb {

Molecule fromCoreMolecule(const Core::Molecule& mol, int charge,
                          int unpairedElectrons)
{
  Molecule out;
  out.charge = charge;
  out.unpairedElectrons = unpairedElectrons;
  const auto& numbers = mol.atomicNumbers();
  const auto& positions = mol.atomPositions3d();
  out.numbers.assign(numbers.begin(), numbers.end());
  out.positions.resize(3 * numbers.size());
  for (size_t i = 0; i < numbers.size(); ++i) {
    out.positions[3 * i] = positions[i][0] * angstromToBohr;
    out.positions[3 * i + 1] = positions[i][1] * angstromToBohr;
    out.positions[3 * i + 2] = positions[i][2] * angstromToBohr;
  }
  return out;
}

void toCoreMolecule(const Molecule& mol, Core::Molecule& out)
{
  const int nat = mol.atomCount();
  Core::Array<unsigned char> numbers(nat);
  Core::Array<Vector3> positions(nat);
  for (int i = 0; i < nat; ++i) {
    numbers[i] = static_cast<unsigned char>(mol.numbers[i]);
    positions[i][0] = mol.positions[3 * i] / angstromToBohr;
    positions[i][1] = mol.positions[3 * i + 1] / angstromToBohr;
    positions[i][2] = mol.positions[3 * i + 2] / angstromToBohr;
  }
  out.setAtomicNumbers(numbers);
  out.setAtomPositions3d(positions);
}

} // namespace Xtb
} // namespace Avogadro
