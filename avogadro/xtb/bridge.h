/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.
******************************************************************************/

#ifndef AVOGADRO_XTB_BRIDGE_H
#define AVOGADRO_XTB_BRIDGE_H

#include "molecule.h"

namespace Avogadro {
namespace Core {
class Molecule;
}

namespace Xtb {

// Convert between Avogadro's Core::Molecule (positions in Angstrom) and the
// port's Molecule (positions in Bohr, mirroring xtb's TMolecule).
// The molecular charge is carried separately: Core::Molecule exposes only a
// computed totalCharge(), so callers pass it explicitly.
Molecule fromCoreMolecule(const Core::Molecule& mol, int charge = 0,
                          int unpairedElectrons = 0);
void toCoreMolecule(const Molecule& mol, Core::Molecule& out);

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_BRIDGE_H
