/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gdisp0.f90 (d3_gradient 0d path, weight_references_d4,
  get_atomic_c6_d4), src/disp/dftd4.F90 (newD3Model, trapzd),
  src/disp/dftd3.f90 (weight_cn) and src/gfnff/gfnff_ini.f90
  (elemental zeta, zetac6 setup), Copyright (C) 2019-2020 Stefan Grimme.
******************************************************************************/

#ifndef AVOGADRO_XTB_GFNFFDISP_H
#define AVOGADRO_XTB_GFNFFDISP_H

#include <vector>

namespace Avogadro {
namespace Xtb {

class Environment;

// Dispersion reference model mirroring TDispersionModel as filled by
// newD3Model: per-element reference counts, reference coordination
// numbers and the integrated C6 table. c6 is [maxElem][maxElem][7][7]
// in the reference; here flat with c6Index().
struct DispModel
{
  int maxElem = 0;
  std::vector<int> nref;    // [maxElem], references per element (Z - 1)
  std::vector<double> cnRef; // [maxElem][7] reference CNs
  std::vector<double> c6;   // [maxElem][maxElem][7][7] C6 table
};

// Gaussian CN weight mirroring weight_cn (dftd3.f90): exp(-wf*(cn-ref)^2)
// with the -200 cutoff (exp(-200) underflows to ~1.4e-87).
double weightCn(double wf, double cn, double ref);

// D4 charge-scaling factor mirroring the elemental zeta in gfnff_ini.f90
// (NOT the dftd4 zeta(a,c,qref,qmod)): exp(3) for qmod < 0, else
// exp(3*(1-exp(c*(1-zeff/qmod)))) with qmod = zeff + q.
double zetaCharge(int atomicNumber, double charge);

// Pair zeta products mirroring the zetac6 loop in gfnff_ini
// (zetac6(ij) = zeta(ati,qa(i))*zeta(atj,qa(j)), packed incl. diagonal).
// qa holds the topology EEQ charges.
bool buildZetaC6(const std::vector<int>& numbers,
                 const std::vector<double>& qa, std::vector<double>& zetac6,
                 Environment& env);

// Reference-model builder mirroring newD3Model (dftd4.F90): per-element
// reference polarizabilities from the param_ref tables, C6 integrated
// with the 23-point trapzd rule (thopi = 3/pi). Elements without REF
// blocks in xtb (Z > 86) keep 0 references like the init-zeroed model.
// maxElem becomes maxval(numbers).
bool buildDispersionModel(const std::vector<int>& numbers, DispModel& model,
                          Environment& env);

// Non-bonded exponents mirroring the ini alphanb loop (0d, single cell):
// sqrt(dum1*dum2) with charge/CN scaling plus the H...H / M...H / C...H /
// O...H specials keyed by the bpair flags (0.85/0.91/1.04 are
// single-precision literals widened on use, like the reference).
// nbCounts holds neighbour counts; pairFlags the bondPairFlags n x n
// table. zetac6 is built alongside via buildZetaC6.
bool buildNonbondedTables(int n, const std::vector<int>& numbers,
                          const std::vector<double>& qa,
                          const std::vector<int>& nbCounts,
                          const std::vector<double>& repan,
                          const std::vector<int>& metal, double nRepScal,
                          double qRepScal, double hhFac, double hh13Rep,
                          double hh14Rep,
                           const std::vector<int>& pairFlags,
                           std::vector<double>& alphanb,
                           std::vector<double>& zetac6, Environment& env);

// D3(0) dispersion energy + gradient mirroring the 0d d3_gradient path
// (zero damping, single pair list, serial pair order): per-pair C6 from
// the CN-weighted references, t6/t8 damping with the packed d3r0 (R0^2)
// radii, zeta_scale/dispscale prefactors, CN-response via the dlogCn
// contraction (direct [(i*n+j)*3+c] layout, matching the dgemv('n')
// call over dcndr memory) and energy = sum(energies) in atom order.
// Like the reference, the virial accumulator is local and discarded,
// so sigma is untouched here. pairs holds 0-based (i,j) with i >= j;
// cn the D3 coordination numbers, dlogCn the gffCoordinationNumber
// derivative layout. xyz is flat row-major, Bohr.
bool d3Gradient(int n, const std::vector<int>& numbers,
                const std::vector<double>& xyz,
                const std::vector<std::pair<int, int>>& pairs,
                const DispModel& model, const std::vector<double>& zetac6,
                const std::vector<double>& d3r0, double dispScale,
                const std::vector<double>& cn,
                const std::vector<double>& dlogCn, double& energy,
                std::vector<double>& gradient, Environment& env);

} // namespace Xtb
} // namespace Avogadro

#endif
