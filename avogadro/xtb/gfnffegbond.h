/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_eg.f90 (egbond, egbend, egtors, gfnffdampa/dampt),
  src/gfnff/math.f (vsub, vlen, valijklffPBC), src/basic_geo.f90
  (crossprod, vecnorm, impsc) and src/constr.f90 (dphidrPBC, omegaPBC,
  domegadrPBC), Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#ifndef AVOGADRO_XTB_GFNFFEGBOND_H
#define AVOGADRO_XTB_GFNFFEGBOND_H

#include <vector>

namespace Avogadro {
namespace Xtb {

class Environment;
struct GffData;
struct GffGenerator;

// Flat 3-vector helpers mirroring vsub/crprod/vlen (math.f, plain
// left-associated sums, NOT the scaled NORM2) and impsc/vecnorm/crossprod
// (basic_geo.f90). norm2 replicates gfortran's NORM2 exactly (verified by
// disassembly plus a 3000-vector battery, bitwise).
double norm2(const double a[3]);
void vecSub(const double a[3], const double b[3], double c[3]);
void vecCross(const double a[3], const double b[3], double c[3]);
double vecLen(const double a[3]);
double vecCos(const double a[3], const double b[3]);
// vecnorm(r, inorm): plain length, normalizing in place for inorm != 0
// when |length| > 1e-14. r has exactly 3 entries.
double vecNorm(double r[3], int inorm);

// Distance damping mirroring gfnffdampa/dampt (angle/torsion variants
// share one formula): rcut = atcuta(t) * (rcov(ia)+rcov(ja))^2.
void bondDamp(double r2, double atcut, const std::vector<double>& rcov,
              int ia, int ja, double& damp, double& ddamp);

// Dihedral angle mirroring valijklffPBC mode 1 with zero shifts.
double dihedralAngle(const std::vector<double>& xyz, int i, int j, int k,
                     int l);
// Dihedral derivatives mirroring dphidrPBC mode 2 with zero shifts:
// dphi/dri..dphi/drl (3-vectors in da..dd).
void dihedralDerivatives(const std::vector<double>& xyz, int i, int j, int k,
                         int l, double phi, double da[3], double db[3],
                         double dc[3], double dd[3]);
// Inversion angle mirroring omegaPBC with zero shifts (no clamping).
double inversionAngle(const std::vector<double>& xyz, int i, int j, int k,
                      int l);
// Inversion derivatives mirroring domegadrPBC with zero shifts.
void inversionDerivatives(const std::vector<double>& xyz, int i, int j, int k,
                          int l, double omega, double da[3], double db[3],
                          double dc[3], double dd[3]);

// Bond energy + gradient mirroring egbond (non-periodic): rab0/rij/drij
// come from the coordination machinery (later batch); dEdcn accumulates
// the 3-body CN derivative part; sigma accumulates stress like the
// reference (unused without PBC). xyz is flat row-major, Bohr.
void bondEnergyGradient(int bond, int iat, int jat, double rab, double rij,
                        const std::vector<double>& drij,
                        const double drijdcn[2], int n,
                        const std::vector<double>& xyz,
                        const std::vector<double>& vbondShift,
                        const std::vector<double>& vbondSteep,
                        const std::vector<double>& vbondPref, double& energy,
                        std::vector<double>& gradient,
                        std::vector<double>& dEdcn, double sigma[3][3],
                        Environment& env);

// Bend energy + gradient mirroring egbend (non-periodic): alist entry
// (center, first, second) with vangl equilibrium/force constant. Returns
// the 3-atom gradient rows and stress in the reference layout.
void bendEnergyGradient(int center, int first, int second, double equilibrium,
                        double forceConstant,
                        const std::vector<int>& numbers,
                        const std::vector<double>& xyz,
                        const std::vector<double>& rcov, double angleCutA,
                        double& energy, double atomGrad[3][3],
                        double stress[3][3]);

// Torsion/improper energy + gradient mirroring egtors (non-periodic).
// Arguments are in egtors(i,j,k,l) order: proper torsions pass
// (ll,ii,jj,kk) = tlist(1..4); impropers pass (center,nb1,nb2,nb3).
// kind > 0 proper torsion (multiplicity, phase, force constant), kind <= 0
// out-of-plane/improper. Returns the 4-atom gradient rows; stress only
// when periodic (reference: nTrans != 1).
void torsionEnergyGradient(int outer1, int center1, int center2, int outer2,
                           int kind, double phase, double forceConstant,
                           const std::vector<int>& numbers,
                           const std::vector<double>& xyz,
                           const std::vector<double>& rcov, double angleCutT,
                           bool periodic, double& energy,
                           double atomGrad[4][3], double stress[3][3]);

// Bond-length estimates mirroring gfnffdrab (non-periodic): reference
// lengths from the r0/cnfak tables with CN dependence, EN-mismatch scaling
// (k1/k2 from the p table via elementRow6) and F/HLn specials, plus full
// gradients over the dlogCn layout ([(i*n+m)*3+c] = d(cn[i])/d(x[m][c])).
// blist holds 0-based (j, i) pairs; rabIn carries the vbond shifts;
// rabOut/grabOut/rabdcnOut receive lengths, gradients ([k][m][c]) and the
// two CN derivatives per bond.
bool bondLengthEstimates(int n, const std::vector<int>& numbers,
                         const std::vector<double>& cn,
                         const std::vector<double>& dlogCn,
                         const std::vector<std::pair<int, int>>& blist,
                         const std::vector<double>& rabIn,
                         std::vector<double>& rabOut,
                         std::vector<double>& grabOut,
                         std::vector<double>& rabdcnOut, Environment& env);

// Bonded repulsion mirroring the nbond loop (non-periodic): exponential
// repulsion over bonded pairs with sqrt(repa*repa) exponents (0.75d0
// power via pow) and repz*repz*repscalb strengths. blist holds 0-based
// (j, i) pairs; xyz is flat row-major, Bohr.
bool bondedRepulsion(int n, const std::vector<int>& numbers,
                     const std::vector<double>& xyz,
                     const std::vector<std::pair<int, int>>& blist,
                     const std::vector<double>& repa,
                     const std::vector<double>& repz, double repScaleB,
                     double& energy, std::vector<double>& gradient,
                     double sigma[3][3], Environment& env);

// Bonded ATM three-body term mirroring batmgfnff_eg (non-periodic,
// 0d image choice): charge-damped C9 Axilrod-Teller-Muto energy over
// (iat, jat, kat) with zb3atm strengths. atomGrad rows hold
// (iat, jat, kat) gradients (row-first like the bend/torsion ports;
// the reference uses columns); stress only when periodic
// (reference: nTrans != 1, skipped here).
void batmEnergyGradient(int iat, int jat, int kat,
                        const std::vector<int>& numbers,
                        const std::vector<double>& xyz,
                        const std::vector<double>& charges,
                        const std::vector<double>& zb3atm, double& energy,
                        double atomGrad[3][3], Environment& env);

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_GFNFFEGBOND_H
