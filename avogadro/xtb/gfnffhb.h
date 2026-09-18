/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_eg.f90 (abhgfnff_eg1/eg2new/eg2_rnr/eg3, rbxgfnff_eg,
  egbend_nci_mul, egtors_nci_mul), src/gfnff/gfnff_ini2.F90 (hbonds),
  src/constr.f90 (dphidrPBC, crprod, impsc) and src/basic_geo.f90
  (crossprod, vecnorm), Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#ifndef AVOGADRO_XTB_GFNFFHB_H
#define AVOGADRO_XTB_GFNFFHB_H

#include <utility>
#include <vector>

namespace Avogadro {
namespace Xtb {

class Environment;

// HB/XB damping and scaling parameters, resolved to double exactly like
// the reference widens its single-precision setup tables on use.
struct HbParams
{
  double hbacut = 0.0;      // HB angle cutoff
  double hblongcut = 0.0;   // HB long-range cutoff (squared distance)
  double hbscut = 0.0;      // HB short-range cutoff factor
  double hbalp = 0.0;       // HB damping exponent
  double hbst = 0.0;        // HB charge-scaling steepness
  double hbsf = 0.0;        // HB charge-scaling shift
  double hbabmix = 0.0;     // HB distance mixing
  double hbnbcut = 0.0;     // HB neighbour-angle parameter
  double xbacut = 0.0;      // XB angle cutoff
  double xbscut = 0.0;      // XB short-range cutoff factor
  double xbst = 0.0;        // XB charge-scaling steepness
  double xbsf = 0.0;        // XB charge-scaling shift
  double hblongcutXb = 0.0; // XB long-range cutoff (squared distance)
  double xhaciGlobAbh = 0.0; // A-H...B global scaling
  double xhaciCoh = 0.0;    // A-H...O=C global scaling
  double torsHb = 0.0;      // HB torsion shift
  double bendHb = 0.0;      // HB bending shift
};

// HB strengths mirroring hbonds (ini2): ca/cb = (hbbas, hbaci) of A/B.
void hbStrengths(int a, int b, const std::vector<double>& hbBas,
                 const std::vector<double>& hbAci, double* ca, double* cb);

// Vector helpers mirroring constr.f90/basic_geo.f90: cross product
// (crprod/crossprod), cosine via norm2 (impsc) and the vecnorm norm with
// its in-place normalisation side effect for flag > 0.
void crossProd(const double* a, const double* b, double* c);
double cosAngle(const double* a, const double* b);
double vecNorm10(double* v, int normalise);

// 0d dihedral value/derivatives mirroring valijklffPBC (mode 2, H in the
// central cell) and dphidrPBC (mode 1). All shift vectors are zero here;
// xyz is flat row-major, Bohr; indices 0-based.
double dihedralValue0d(int r, int b, int c, int h,
                       const std::vector<double>& xyz);
void dihedralDerivs0d(int r, int b, int c, int h, double phi,
                      const std::vector<double>& xyz, double* dda,
                      double* ddb, double* ddc, double* ddd);

// NCI bending factor mirroring egbend_nci_mul (0d): e = 1-ea with the
// linear/nonlinear branches; g columns hold (B, C, H) = (j, k, l).
void bendNciMul(int b, int c, int h, double c0, double fc,
                const std::vector<double>& xyz, double& e, double g[3][3]);

// NCI torsion factor mirroring egtors_nci_mul (0d): e = (1+cos)*fc+tshift
// with phi0/tshift; g columns hold (R, B, C, H).
void torsNciMul(int r, int b, int c, int h, int rn, double phi, double phi0,
                double tshift, const std::vector<double>& xyz, double& e,
                double g[3][4]);

// A-H...B hydrogen bond without neighbour orientation (abhgfnff_eg1,
// 0d): qa holds the topology charges, hbBas/hbAci the per-atom HB
// strengths; gdr columns hold (A, B, H); sigma accumulates the
// mcf-scaled virial (the caller scales energy and gdr by mcf like the
// driver). Returns false on empty input; the expo > 15 early return
// yields zero energy like the reference.
bool abhEg1(int a, int b, int h, const std::vector<int>& numbers,
            const std::vector<double>& xyz, const std::vector<double>& qa,
            const std::vector<double>& hbBas, const std::vector<double>& hbAci,
            const std::vector<double>& rad, const HbParams& par, double mcf,
            double& energy, double gdr[3][3], double sigma[3][3],
            Environment& env);

// A-H...B with neighbour orientation at B (abhgfnff_eg2new, 0d):
// nbrsB lists B's neighbour atom ids (0-based, jth_nb order); gdr is the
// full (3, n) gradient accumulator.
bool abhEg2new(int a, int b, int h, const std::vector<int>& nbrsB,
               const std::vector<int>& numbers,
               const std::vector<double>& xyz,
               const std::vector<double>& qa,
               const std::vector<double>& hbBas,
               const std::vector<double>& hbAci,
               const std::vector<double>& rad, const HbParams& par,
               double mcf, double& energy, std::vector<double>& gdr,
               double sigma[3][3], Environment& env);

// A-H...B with lone-pair orientation at N (abhgfnff_eg2_rnr, 0d):
// nbrsB holds exactly the two neighbours of B.
bool abhEg2rnr(int a, int b, int h, const std::vector<int>& nbrsB,
               const std::vector<int>& numbers,
               const std::vector<double>& xyz,
               const std::vector<double>& qa,
               const std::vector<double>& hbBas,
               const std::vector<double>& hbAci,
               const std::vector<double>& rad,
               const std::vector<double>& repz, const HbParams& par,
               double mcf, double& energy, std::vector<double>& gdr,
               double sigma[3][3], Environment& env);

// A-H...O=C hydrogen bond (abhgfnff_eg3, 0d): cNbrs lists C's neighbour
// atom ids including B (skipped inside like the reference); sqrab/srab
// are packed (unused in the 0d body but kept for signature parity).
bool abhEg3(int a, int b, int h, int c,
            const std::vector<int>& cNbrs, const std::vector<int>& numbers,
            const std::vector<double>& xyz, const std::vector<double>& qa,
            const std::vector<double>& hbBas,
            const std::vector<double>& hbAci,
            const std::vector<double>& sqrab, const std::vector<double>& srab,
            const std::vector<double>& rad, const HbParams& par, double mcf,
            double& energy, std::vector<double>& gdr, double sigma[3][3],
            Environment& env);

// A...X-B halogen bond (rbxgfnff_eg, 0d): cx/cb acidities (cb is 1.0 in
// the reference, kept as argument for parity); gdr columns hold (A,B,X);
// sigma is unscaled by mcf unlike the HB kernels (verbatim).
bool rbxEg(int a, int b, int x, double cb, double cx,
           const std::vector<int>& numbers, const std::vector<double>& xyz,
           const std::vector<double>& qa, const std::vector<double>& rad,
           const HbParams& par, double& energy, double gdr[3][3],
           double sigma[3][3], Environment& env);

// Halogen-bond heavy atom predicate mirroring xatom() (ini2): P, S, Cl,
// As, Se, Br, Sb, Te, I.
bool isXAtom(int atomicNumber);

// HB donor/acceptor perception mirroring the ini HAB-list block (0d,
// single cell): H atoms with positive topology charges (hbatHl),
// AB pairs with negative charges and HB strength (hbatABl), and A-X...B
// halogen triples (xbatABl). neighbours holds 0-based adjacency in
// jth_nb order; firstNb holds the first neighbour per atom (-1 if none);
// piFlags is the piadr map; group/xhbas are element tables (Z - 1);
// bpair holds the n x n integer exclusion flags with pair(a,b) at
// b*bStride + a, i.e. row = second index like pair(j,i). The 0d
// translation sum is 1 (central cell); the getTransVec call is PBC-only
// and skipped. Threshold adjustments (-0.20/+0.05/+0.2, 0.05) are
// single-precision widened on use like the reference.
struct HbDonorLists
{
  std::vector<int> hatH;                    // 0-based H atoms
  std::vector<std::pair<int, int>> hatAB;   // 0-based (i, j), i > j
  std::vector<int> xbatA, xbatB, xbatX;     // 0-based A/B/X triples
};
bool hbDonorLists(int n, const std::vector<int>& numbers,
                  const std::vector<int>& hyb, const std::vector<int>& piFlags,
                  const std::vector<double>& qa,
                  const std::vector<std::vector<int>>& neighbours,
                  const std::vector<int>& firstNb, const int* group,
                  const std::vector<double>& xhbas, double hQaThr,
                  double qaBThr, const std::vector<int>& bpair,
                  const std::vector<double>& hbBas,
                  const std::vector<double>& hbAci, HbDonorLists& lists,
                  Environment& env);

// Bonded A-H...B triplets mirroring bond_hbset0 (count) + bond_hbset
// (fill) in 0d: single central cell, so the two routines coincide and
// one pass fills the triplets (iTrA = iTrB = iTrH = 1 implicitly).
// rab is the squared A...B distance compared against hbThr1 (squared
// threshold); bpair layout as in hbDonorLists.
struct HbBondTriplet
{
  int a = -1, b = -1, h = -1; // 0-based A, B, H
};
bool hbBondTriplets(int n, const std::vector<std::pair<int, int>>& hatAB,
                     const std::vector<int>& hatH,
                     const std::vector<double>& xyz, double hbThr1,
                     const std::vector<int>& bpair,
                     std::vector<HbBondTriplet>& triplets, Environment& env);

// A-H...B group maps mirroring bond_hb_AHB_set0/1/set in 0d (all shifts
// are the central cell): distinct (A,H) groups with their B lists
// (bond_hb_AH/B/Bn), the atom participation flags (isABH), the per-bond
// HB counts in blist order (nr_hb, with the reference quirk that bonds
// without a match keep the stale running B count) and the compact
// atom->HB-index maps (hb_mapABH, 0-based within the AB and H classes).
// bonds holds (blist1, blist2) endpoint pairs in reference blist order;
// only N/O-H bonds participate.
struct HbBondMaps
{
  std::vector<int> ahA, ahH;                 // 0-based A/H per group
  std::vector<std::vector<int>> bAtoms;      // 0-based B per group
  std::vector<char> isAbh;                   // per-atom participation
  std::vector<int> nrHb;                     // per-bond HB count
  std::vector<int> mapAbh;                   // per-atom compact HB index
  int mapNab = 0, mapNh = 0, bMax = 1;
};
bool hbAhbMaps(int n, const std::vector<int>& numbers,
               const std::vector<std::pair<int, int>>& bonds,
               const std::vector<HbBondTriplet>& triplets, HbBondMaps& maps,
               Environment& env);

// Per-atom HB basicity/acidity mirroring the ini hbbas/hbaci blocks
// (0d): element-table values with the carbene, carbonyl/nitro oxygen
// (last neighbour in nb order, single cell) and amide-H neighbour
// (x0.80f) corrections. neighbours/counts parallel hbDonorLists;
// hyb/piFlags feed isAmideHydrogen.
bool hbBasicity(int n, const std::vector<int>& numbers,
                const std::vector<int>& itag,
                const std::vector<std::vector<int>>& neighbours,
                const std::vector<int>& counts,
                const std::vector<double>& xhbas, std::vector<double>& hbBas,
                Environment& env);
bool hbAcidity(int n, const std::vector<int>& numbers,
               const std::vector<std::vector<int>>& neighbours,
               const std::vector<int>& hyb, const std::vector<int>& piFlags,
               const std::vector<double>& xhaci, std::vector<double>& hbAci,
               Environment& env);

// HB erf coordination mirroring dncoord_erf (0d, thr = 900): coordination
// numbers over the A-H...B groups with the 1.78-scaled covalent radii,
// kn = 27.5 and the double-precision 1/sqrt(pi) literal. hbDcn is flat
// [(a * n + b) * 3 + c] = d(cn[b])/d(x[a]) like dlogCn. The strain
// derivative feeds only the discarded 0d sigma and is skipped.
bool hbErfCoordination(int n, const std::vector<int>& numbers,
                       const std::vector<double>& xyz,
                       const std::vector<double>& rcov, const HbBondMaps& maps,
                       std::vector<double>& hbCn, std::vector<double>& hbDcn,
                       Environment& env);

// HB/XB geometry lists mirroring gfnff_hbset0 (counts) + gfnff_hbset
// (fill) in 0d without the OpenMP batching (serial ix-major order) and
// without the rmsd neighbour-list cache (single-point rebuild):
// non-covalent A-H...B (hb1, rab + rAH + rBH below hbThr2), bonded-A
// A-H...B (hb2) and A-X...B halogen entries (xb). All distances squared;
// bpair layout as in hbDonorLists. Triples carry 0-based indices.
struct HbTriple
{
  int a = -1, b = -1, h = -1;
};
struct XbTriple
{
  int a = -1, b = -1, x = -1;
};
bool hbTripletLists(int n, const std::vector<double>& xyz,
                    const std::vector<std::pair<int, int>>& hatAB,
                    const std::vector<int>& hatH,
                    const HbDonorLists& donors,
                    const std::vector<int>& bpair, double hbThr1,
                    double hbThr2, std::vector<HbTriple>& hb1,
                    std::vector<HbTriple>& hb2, std::vector<XbTriple>& xb,
                    Environment& env);

// Bonded HB correction for A-H...B bonds mirroring egbond_hb (0d):
// CN-weakened bond energy/gradient plus the A-H...B group CN-gradient
// part with the once-per-ABH-triple dedup (consideredAbh is flat
// (mapA * mapNab + mapB) * mapNh + mapH over the 0-based compact maps).
// drij/drijdcn come from the rab estimates, steepness/prefactor are
// vbond(2/3) of bond bondIdx, dEdcn accumulates like the reference
// (its sigma contraction is skipped as 0d-neutral). Sigma updates are
// skipped: sigma is discarded by the 0d driver. Bonds without an H
// endpoint are a no-op like the reference early return.
bool egbondHb(int bondIdx, int iat, int jat, double rab, double rij,
              const std::vector<double>& drij, const double drijdcn[2],
              const std::vector<int>& numbers,
              const std::vector<double>& hbCn, const std::vector<double>& hbDcn,
              const std::vector<double>& xyz, double vbondScale,
              double steepness, double prefactor, const HbBondMaps& maps,
              std::vector<char>& consideredAbh, std::vector<double>& dEdcn,
              double& e, std::vector<double>& g, Environment& env);

} // namespace Xtb
} // namespace Avogadro

#endif
