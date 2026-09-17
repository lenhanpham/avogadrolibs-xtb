/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/data.f90 (TGFFData), src/gfnff/generator.f90 (TGFFGenerator),
  src/gfnff/topology.f90 (TGFFTopology) and the version enum in
  src/gfnff/gfnff_param.f90, Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#ifndef AVOGADRO_XTB_GFNFFDATA_H
#define AVOGADRO_XTB_GFNFFDATA_H

#include <array>
#include <vector>

namespace Avogadro {
namespace Xtb {

// Force field versions, mirroring TGFFVersionEnum.
enum GffVersion
{
  GffAngewChem2020 = 1,
  GffAngewChem2020_1 = 2,
  GffAngewChem2020_2 = 3,
  GffHarmonic2020 = -1,
  GffMcGfnFF2023 = 4
};

// Parametrisation data of the force field, mirroring TGFFData.
// Vectors replace the allocatable arrays (sized by initGffData);
// dispm (TDispersionModel) is not ported yet and follows with the D3/D4 port.
struct GffData
{
  double repScaleN = 0.0; // repscaln, non-bonded repulsion scaling
  double repScaleB = 0.0; // repscalb, bonded repulsion scaling
  double angleCutA = 0.0; // atcuta, bend angle damping
  double angleCutT = 0.0; // atcutt, torsion angle damping
  double angleCutANci = 0.0; // atcuta_nci
  double angleCutTNci = 0.0; // atcutt_nci
  double hbAngleCut = 0.0; // hbacut
  double hbShortCut = 0.0; // hbscut
  double xbAngleCut = 0.0; // xbacut
  double xbShortCut = 0.0; // xbscut
  double hbAlp = 0.0; // hbalp, HB/XB damping
  double hbLongCut = 0.0; // hblongcut
  double hbLongCutXb = 0.0; // hblongcut_xb
  double hbSt = 0.0; // hbst, HB charge scaling steepness
  double hbSf = 0.0; // hbsf, HB charge scaling factor
  double xbSt = 0.0; // xbst
  double xbSf = 0.0; // xbsf
  double xhAciGlobAbH = 0.0; // xhaci_globabh
  double xhAciCoh = 0.0; // xhaci_coh
  double xhAciGlob = 0.0; // xhaci_glob, acidity
  double hbAbMix = 0.0; // hbabmix
  double hbNbCut = 0.0; // hbnbcut, neighbour angle parameter
  double torsHb = 0.0; // tors_hb, HB torsion term
  double bendHb = 0.0; // bend_hb, HB bending term
  double vbondScale = 0.0; // vbond_scale
  double cnMax = 0.0; // cnmax, max CN cut-off
  double dispScale = 0.0; // dispscale, D3 scaling

  std::vector<double> en; // electronegativities
  std::vector<double> rad; // radii in Angstrom
  std::vector<double> rcov; // D3 covalent radii in Bohr
  std::vector<int> metal; // metal class
  std::vector<int> group; // periodic group
  std::vector<int> normCn; // normal coordination numbers
  std::vector<double> repa; // rep alpha bond
  std::vector<double> repan;
  std::vector<double> repz; // 3-atom bond prefactor (Zval)
  std::vector<double> zb3atm;
  std::vector<double> xhAci; // HB acidities
  std::vector<double> xhBas; // HB basicities
  std::vector<double> xbAci; // XB acidities
  std::vector<double> chi; // EEQ electronegativities
  std::vector<double> gam; // EEQ hardness
  std::vector<double> cnf; // EEQ CN factors
  std::vector<double> alp; // EEQ exponents
  std::vector<double> bond; // element bond parameters
  std::vector<double> angl; // element angle parameters
  std::vector<double> angl2;
  std::vector<double> tors; // element torsion parameters
  std::vector<double> tors2;
  std::vector<double> d3r0; // Bjorkman radii, packed pairs

  // Allocate all per-element arrays for ndim elements (initGFFData).
  void init(int ndim);
};

// Topology generator settings, mirroring TGFFGenerator. Zero-initialised;
// fields the reference leaves unset (e.g. torsf(4)) stay 0.0 by choice,
// where Fortran leaves them indeterminate.
struct GffGenerator
{
  double linThr = 0.0; // linthr, linear-angle threshold
  double fcThr = 0.0; // fcthr, skip small torsion/bending
  double tdistThr = 0.0; // tdist_thr, Angstrom (sp in the original)
  double rThr = 0.0; // rthr, bond determination threshold
  double rThr2 = 0.0; // rthr2
  double rShrink = 0.0; // rqshrink
  double hQaThr = 0.0; // hqabthr
  double qaBThr = 0.0; // qabthr
  double srb1 = 0.0, srb2 = 0.0, srb3 = 0.0;
  double qRepScal = 0.0; // qrepscal
  double nRepScal = 0.0; // nrepscal
  double hhFac = 0.0; // hhfac, HH repulsion
  double hh13Rep = 0.0, hh14Rep = 0.0;
  std::array<double, 9> bstren = {}; // bond strengths by bond type
  double qFacBen = 0.0; // qfacBEN
  double qFacTor = 0.0; // qfacTOR
  double fr3 = 0.0, fr4 = 0.0, fr5 = 0.0, fr6 = 0.0; // ring torsion factors
  std::array<double, 8> torsf = {}; // torsion factors (index 3 unset = 0)
  double fbs1 = 0.0;
  double batmScal = 0.0; // batmscal
  double mchiShift = 0.0; // mchishift
  double rabShift = 0.0; // rabshift
  double rabShiftH = 0.0; // rabshifth
  double hyperShift = 0.0; // hyper_shift
  double hShift3 = 0.0, hShift4 = 0.0, hShift5 = 0.0;
  double metal1Shift = 0.0, metal2Shift = 0.0, metal3Shift = 0.0;
  double etaShift = 0.0; // eta_shift
  std::array<double, 5> qfacbm = {}; // bond charge dependence (0..4)
  double qfacbm0 = 0.0;
  double rfgoed1 = 0.0;
  double hTriple = 0.0; // htriple
  double hueckelP2 = 0.0; // hueckelp2
  double hueckelP3 = 0.0; // hueckelp3
  std::array<double, 17> hdiag = {}; // Huckel diagonal (1-based in xtb)
  std::array<double, 17> hoffdiag = {}; // Huckel off-diagonal (1-based)
  double hIter = 0.0; // hiter
  double hueckelP = 0.0; // hueckelp
  double bzRef = 0.0; // bzref
  double bzRef2 = 0.0; // bzref2
  double pilpf = 0.0;
  double maxHIter = 0.0; // maxhiter
  double d3a1 = 0.0, d3a2 = 0.0; // D3 parameters
  double split0 = 0.0, split1 = 0.0;
  double fringbo = 0.0;
  double aheavy3 = 0.0, aheavy4 = 0.0;
  std::array<std::array<double, 4>, 4> bsmat = {}; // bond strength matrix
  double cnMax = 0.0; // cnmax
};

// Topology of a system, mirroring TGFFTopology (layout for the upcoming
// energy port; dispm follows with the dispersion port).
struct GffTopology
{
  int nBond = 0, nAngl = 0, nTors = 0, nStors = 0;
  int nAtHbH = 0, nAtHbAB = 0, nAtXbAB = 0;
  int nBatM = 0, nFrag = 0, maxSystem = 0;
  int bondHbNr = 0, bMax = 0;
  int nBondBlist = 0, nBondVbond = 0, nAnglAlloc = 0, nTorsAlloc = 0;
  int readFileType = 0;

  std::vector<int> hyb; // hybridization per atom
  std::vector<int> blist; // bonded atoms, pairs
  std::vector<int> alist; // angles, triples
  std::vector<int> tlist; // torsions, quadruples
  std::vector<int> b3list; // 3-atom terms
  std::vector<int> sTorsl; // special torsions
  std::vector<int> nrHb;
  std::vector<int> bondHbAH;
  std::vector<int> bondHbB;
  std::vector<int> bondHbBn;
  std::vector<int> hbAtABl, xbAtABl, hbAtHl;
  std::vector<int> fraglist;
  std::vector<int> qpdb;
  std::vector<double> vbond, vangl, vtors;
  std::vector<double> chiEeq, gamEeq, alpEeq;
  std::vector<double> alphanb;
  std::vector<double> qa; // topology EEQ charges
  std::vector<double> xyze0;
  std::vector<double> zetac6;
  std::vector<double> qfrag;
  std::vector<double> hbBas, hbAci;
  std::vector<int> hbMapABH;
  std::vector<bool> isABH;
  int hbMapNAB = 0, hbMapNH = 0;
};

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_GFNFFDATA_H
