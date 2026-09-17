/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/mctc/mctc_constants.f90, Copyright (C) 2017-2020 Stefan Grimme.
******************************************************************************/

#ifndef AVOGADRO_XTB_CONSTANTS_H
#define AVOGADRO_XTB_CONSTANTS_H

namespace Avogadro {
namespace Xtb {

// Physical constants, matching xtb_mctc_constants (all double precision).
// CODATA 2018 values, atomic units unless suffixed with _SI.
constexpr double pi = 3.1415926535897932384626433832795029;
constexpr double sqrtpi = 1.7724538509055160272981674833411452; // sqrt(pi)
constexpr double twopi = 2.0 * pi;
constexpr double fourpi = 4.0 * pi;
constexpr double pihalf = 0.5 * pi;
constexpr double twothirdpi = 2.0 * pi / 3.0;
// Boltzmann constant in Eh/K
constexpr double kB = 3.166808578545117e-06;
// Speed of light in vacuum in a.u.
constexpr double lightspeed = 137.0359990740;
constexpr double kB_SI = 1.380649e-23;
constexpr double lightspeed_SI = 299792458.0;
constexpr double h_SI = 6.62607015e-34;
constexpr double standardAtmosphere = 101325e0; // Pa
constexpr double bohrRadius = 5.29177210903e-11; // m
constexpr double molarGasConstant = 8.314462618; // J mol^-1 K^-1
// Angstrom -> Bohr conversion (CODATA 2018)
constexpr double angstromToBohr = 1.0e-10 / bohrRadius;
// xtb's internal Angstrom definition (mctc_convert's aatoau = 1/autoaa with
// autoaa = 0.52917726). Differs from the CODATA value above at ~1e-9
// relative; always use this one when reproducing xtb table conversions
// such as covalentRadD3.
constexpr double xtbAngstromToBohr = 1.0 / 0.52917726;

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_CONSTANTS_H
