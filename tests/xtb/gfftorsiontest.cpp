/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "avogadro/xtb/COPYING.LESSER".
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/xtb/environment.h>
#include <avogadro/xtb/gffgraph.h>
#include <avogadro/xtb/gfnffbonds.h>
#include <avogadro/xtb/gfnffdata.h>
#include <avogadro/xtb/gfnffsetup.h>
#include <avogadro/xtb/gfnfftorsion.h>

#include <cmath>
#include <vector>

namespace Avogadro {
namespace Xtb {

// Reference values verified bitwise against the unmodified tlist/vtors
// loop in xtb src/gfnff/gfnff_ini.f90 (1855-2048), compiled with gfortran
// -ffp-contract=off over 32 synthetic fragments (115 torsions).

namespace {

struct Mol
{
  std::vector<TorsionAtom> atoms;
  std::vector<TorsionBond> bonds; // central bonds only
  std::vector<std::vector<int>> neighbours;
  std::vector<std::vector<Ring>> rings;
  std::vector<double> xyz;
};

int addAtom(Mol& m, int z, int hyb, double qa, int imetal, int pi, double x,
            double y, double zc)
{
  TorsionAtom a;
  a.element = z;
  a.hyb = hyb;
  a.charge = qa;
  a.imetal = imetal;
  a.pi = pi;
  m.atoms.push_back(a);
  m.xyz.push_back(x);
  m.xyz.push_back(y);
  m.xyz.push_back(zc);
  m.neighbours.emplace_back();
  m.rings.emplace_back();
  return static_cast<int>(m.atoms.size()) - 1;
}

void link(Mol& m, int a, int b, double pibo = 0.0, int btype = -1)
{
  m.neighbours[a].push_back(b);
  m.neighbours[b].push_back(a);
  TorsionBond bond;
  bond.first = a;
  bond.second = b;
  bond.pibo = pibo;
  bond.btype = btype;
  m.bonds.push_back(bond);
}

void finish(Mol& m, const GffData& param)
{
  // Central-bond types via the verified assignBondTypes (itag zero here).
  int n = static_cast<int>(m.atoms.size());
  std::vector<Bond> blist;
  std::vector<int> numbers(n), hyb(n), itag(n, 0), piF(n), imet(n);
  for (int i = 0; i < n; ++i) {
    numbers[i] = m.atoms[i].element;
    hyb[i] = m.atoms[i].hyb;
    piF[i] = m.atoms[i].pi;
    imet[i] = m.atoms[i].imetal;
  }
  for (const TorsionBond& b : m.bonds) {
    Bond bnd;
    bnd.first = b.first;
    bnd.second = b.second;
    blist.push_back(bnd);
  }
  std::vector<int> btypes = assignBondTypes(blist, numbers, hyb, itag, piF,
                                            imet, param.group.data());
  for (size_t i = 0; i < m.bonds.size(); ++i) {
    if (m.bonds[i].btype < 0)
      m.bonds[i].btype = btypes[i];
  }
}

Mol makeEthane()
{
  Mol m;
  int c1 = addAtom(m, 6, 3, 0.0, 0, 0, -0.77, 0, 0);
  int c2 = addAtom(m, 6, 3, 0.0, 0, 0, 0.77, 0, 0);
  link(m, c1, c2);
  for (int rep = 0; rep < 2; ++rep) {
    int c = rep == 0 ? c1 : c2;
    double s = rep == 0 ? -1.0 : 1.0;
    for (int k = 0; k < 3; ++k) {
      double t = (120 * k + (rep == 0 ? 0 : 30)) * 3.141592653589793 / 180;
      int h = addAtom(m, 1, 0, 0.0, 0, 0, s * 0.77 - 0.36,
                      1.03 * std::cos(t), 1.03 * std::sin(t));
      link(m, c, h);
    }
  }
  return m;
}

Mol makeEthylene()
{
  Mol m;
  int c1 = addAtom(m, 6, 2, 0.0, 0, 1, -0.67, 0, 0);
  int c2 = addAtom(m, 6, 2, 0.0, 0, 1, 0.67, 0, 0);
  link(m, c1, c2, 0.65);
  int h1 = addAtom(m, 1, 0, 0.0, 0, 0, -1.22, 0.95, 0);
  int h2 = addAtom(m, 1, 0, 0.0, 0, 0, -1.22, -0.95, 0);
  int h3 = addAtom(m, 1, 0, 0.0, 0, 0, 1.22, 0.95, 0);
  int h4 = addAtom(m, 1, 0, 0.0, 0, 0, 1.22, -0.95, 0);
  link(m, c1, h1);
  link(m, c1, h2);
  link(m, c2, h3);
  link(m, c2, h4);
  return m;
}

Mol makeBenzene()
{
  Mol m;
  std::vector<int> ids;
  for (int k = 0; k < 6; ++k) {
    double t = k * 3.141592653589793 / 3.0;
    ids.push_back(addAtom(m, 6, 2, 0.0, 0, 1, 1.40 * std::cos(t),
                          1.40 * std::sin(t), 0));
  }
  for (int k = 0; k < 6; ++k) {
    link(m, ids[k], ids[(k + 1) % 6], 0.6);
    double t = k * 3.141592653589793 / 3.0;
    int h = addAtom(m, 1, 0, 0.0, 0, 0, 2.48 * std::cos(t),
                    2.48 * std::sin(t), 0);
    link(m, ids[k], h);
  }
  for (int a : ids) {
    Ring r;
    for (int b : ids)
      r.members.push_back(b);
    m.rings[a].push_back(r);
  }
  return m;
}

Mol makeCyclopropane()
{
  Mol m;
  std::vector<int> ids;
  for (int k = 0; k < 3; ++k) {
    double t = (120 * k + 90) * 3.141592653589793 / 180;
    ids.push_back(addAtom(m, 6, 3, 0.0, 0, 0, 1.51 * std::cos(t),
                          1.51 * std::sin(t), 0));
  }
  for (int k = 0; k < 3; ++k)
    link(m, ids[k], ids[(k + 1) % 3]);
  for (int a : ids) {
    Ring r;
    for (int b : ids)
      r.members.push_back(b);
    m.rings[a].push_back(r);
  }
  for (int k = 0; k < 3; ++k) {
    double t = (120 * k + 90) * 3.141592653589793 / 180;
    int h = addAtom(m, 1, 0, 0.0, 0, 0, 2.6 * std::cos(t),
                    2.6 * std::sin(t), 0);
    link(m, ids[k], h);
  }
  return m;
}

Mol makeLactam()
{
  Mol m;
  int n = addAtom(m, 7, 3, 0.0, 0, 1, 0.0, 1.0, 0);
  int c1 = addAtom(m, 6, 2, 0.0, 0, 1, 1.30, 0.40, 0);
  int o = addAtom(m, 8, 2, 0.0, 0, 1, 2.30, 1.05, 0);
  int c2 = addAtom(m, 6, 3, 0.0, 0, 0, 1.10, -1.05, 0);
  int c3 = addAtom(m, 6, 3, 0.0, 0, 0, -0.15, -1.35, 0);
  int c4 = addAtom(m, 6, 3, 0.0, 0, 0, -1.05, -0.35, 0);
  int h = addAtom(m, 1, 0, 0.0, 0, 0, 0.05, 2.05, 0);
  std::vector<int> ids = { n, c1, c2, c3, c4 };
  for (int k = 0; k < 5; ++k)
    link(m, ids[k], ids[(k + 1) % 5]);
  link(m, c1, o);
  link(m, n, h);
  for (int a : ids) {
    Ring r;
    for (int b : ids)
      r.members.push_back(b);
    m.rings[a].push_back(r);
  }
  return m;
}

Mol makeMchHcch()
{
  // H-C(ring)-C(ring)-H with a 6-ring on the carbons: terminal rule.
  Mol m;
  int ca = addAtom(m, 6, 3, 0.0, 0, 0, -0.77, 0, 0);
  int cb = addAtom(m, 6, 3, 0.0, 0, 0, 0.77, 0, 0);
  int ha = addAtom(m, 1, 0, 0.0, 0, 0, -0.77, 1.05, 0);
  int hb = addAtom(m, 1, 0, 0.0, 0, 0, 0.77, -1.05, 0);
  int x1 = addAtom(m, 6, 3, 0.0, 0, 0, -1.9, -0.6, 0);
  int x2 = addAtom(m, 6, 3, 0.0, 0, 0, 1.9, 0.6, 0);
  int x3 = addAtom(m, 6, 3, 0.0, 0, 0, -1.9, 0.6, 0);
  link(m, ca, cb);
  link(m, ca, ha);
  link(m, cb, hb);
  link(m, ca, x1);
  link(m, cb, x2);
  link(m, ca, x3);
  Ring r;
  r.members = { ca, cb, x1, x2, 100, 101 };
  m.rings[ca].push_back(r);
  m.rings[cb].push_back(r);
  return m;
}

Mol makeAcetone()
{
  Mol m;
  int ck = addAtom(m, 6, 2, 0.0, 0, 1, 0, 0, 0);
  int o = addAtom(m, 8, 2, 0.0, 0, 1, 0, 1.22, 0);
  int m1 = addAtom(m, 6, 3, 0.0, 0, 0, -1.28, -0.5, 0);
  int m2 = addAtom(m, 6, 3, 0.0, 0, 0, 1.28, -0.5, 0);
  link(m, ck, o);
  link(m, ck, m1);
  link(m, ck, m2);
  int h1 = addAtom(m, 1, 0, 0.0, 0, 0, -1.7, -1.3, 0);
  int h2 = addAtom(m, 1, 0, 0.0, 0, 0, 1.7, -1.3, 0);
  link(m, m1, h1);
  link(m, m2, h2);
  return m;
}

Mol makeNMethylacetamide()
{
  Mol m;
  int n = addAtom(m, 7, 3, -0.1, 0, 1, -0.7, 0, 0);
  int c = addAtom(m, 6, 2, 0.0, 0, 1, 0.75, 0, 0);
  link(m, n, c);
  int o = addAtom(m, 8, 2, 0.0, 0, 1, 1.5, 1.0, 0);
  link(m, c, o);
  int me = addAtom(m, 6, 3, 0.0, 0, 0, -1.9, -0.5, 0);
  link(m, n, me);
  int mh1 = addAtom(m, 1, 0, 0.0, 0, 0, -2.5, 0.3, 0);
  int mh2 = addAtom(m, 1, 0, 0.0, 0, 0, -2.0, -1.55, 0);
  link(m, me, mh1);
  link(m, me, mh2);
  int h = addAtom(m, 1, 0, 0.0, 0, 0, -0.7, 1.05, 0);
  link(m, n, h);
  return m;
}

Mol makeEthylamine()
{
  Mol m;
  int c1 = addAtom(m, 6, 3, 0.0, 0, 0, -0.75, 0, 0);
  int c2 = addAtom(m, 6, 3, 0.0, 0, 0, 0.75, 0, 0);
  link(m, c1, c2);
  int n = addAtom(m, 7, 3, -0.2, 0, 0, 1.9, 0.5, 0);
  link(m, c2, n);
  int h = addAtom(m, 1, 0, 0.0, 0, 0, -0.75, 1.05, 0);
  link(m, c1, h);
  return m;
}

Mol makeVinylSilane()
{
  Mol m;
  int c1 = addAtom(m, 6, 2, 0.0, 0, 1, -0.67, 0, 0);
  int c2 = addAtom(m, 6, 2, 0.0, 0, 1, 0.67, 0, 0);
  link(m, c1, c2, 0.7);
  int si = addAtom(m, 14, 3, 0.0, 0, 0, -1.9, -0.6, 0);
  int h1 = addAtom(m, 1, 0, 0.0, 0, 0, -1.2, 0.9, 0);
  link(m, c1, si);
  link(m, c1, h1);
  int h2 = addAtom(m, 1, 0, 0.0, 0, 0, 1.2, 0.9, 0);
  link(m, c2, h2);
  return m;
}

Mol makeHydrazine()
{
  Mol m;
  int n1 = addAtom(m, 7, 3, 0.0, 0, 0, -0.72, 0, 0);
  int n2 = addAtom(m, 7, 3, 0.0, 0, 0, 0.72, 0, 0);
  link(m, n1, n2);
  int h1 = addAtom(m, 1, 0, 0.0, 0, 0, -1.07, 0.9, 0.2);
  int h2 = addAtom(m, 1, 0, 0.0, 0, 0, -1.07, -0.9, -0.2);
  int h3 = addAtom(m, 1, 0, 0.0, 0, 0, 1.07, 0.9, -0.2);
  int h4 = addAtom(m, 1, 0, 0.0, 0, 0, 1.07, -0.9, 0.2);
  link(m, n1, h1);
  link(m, n1, h2);
  link(m, n2, h3);
  link(m, n2, h4);
  return m;
}

Mol makePeroxide()
{
  Mol m;
  int o1 = addAtom(m, 8, 3, 0.0, 0, 0, -0.74, 0, 0);
  int o2 = addAtom(m, 8, 3, 0.0, 0, 0, 0.74, 0, 0);
  link(m, o1, o2);
  int h1 = addAtom(m, 1, 0, 0.0, 0, 0, -1.1, 0.9, 0.3);
  int h2 = addAtom(m, 1, 0, 0.0, 0, 0, 1.1, -0.9, 0.3);
  link(m, o1, h1);
  link(m, o2, h2);
  return m;
}

Mol makeDisulfane()
{
  Mol m;
  int s1 = addAtom(m, 16, 3, 0.0, 0, 0, -1.03, 0, 0);
  int s2 = addAtom(m, 16, 3, 0.0, 0, 0, 1.03, 0, 0);
  link(m, s1, s2);
  int h1 = addAtom(m, 1, 0, 0.0, 0, 0, -1.4, 0.9, 0.3);
  int h2 = addAtom(m, 1, 0, 0.0, 0, 0, 1.4, -0.9, 0.3);
  link(m, s1, h1);
  link(m, s2, h2);
  return m;
}

Mol makeLinearSkip()
{
  Mol m;
  int c1 = addAtom(m, 6, 3, 0.0, 0, 0, -0.77, 0, 0);
  int c2 = addAtom(m, 6, 3, 0.0, 0, 0, 0.77, 0, 0);
  link(m, c1, c2);
  int k = addAtom(m, 1, 0, 0.0, 0, 0, -1.4, 0.9, 0);
  link(m, c1, k);
  int ll = addAtom(m, 1, 0, 0.0, 0, 0, 1.8598339877204666,
                   0.01902312301663903, 0);
  link(m, c2, ll);
  int q = addAtom(m, 1, 0, 0.0, 0, 0, 0.77, -1.05, 0);
  link(m, c2, q);
  return m;
}

Mol makeFused()
{
  Mol m;
  int a = addAtom(m, 6, 3, 0.0, 0, 0, -0.77, 0, 0);
  int b = addAtom(m, 6, 3, 0.0, 0, 0, 0.77, 0, 0);
  link(m, a, b);
  int ha = addAtom(m, 1, 0, 0.0, 0, 0, -0.77, 1.05, 0);
  int hb = addAtom(m, 1, 0, 0.0, 0, 0, 0.77, -1.05, 0);
  link(m, a, ha);
  link(m, b, hb);
  Ring r5, r6;
  r5.members = { a, b, 100, 101, 102 };
  r6.members = { a, b, 200, 201, 202, 203 };
  m.rings[a].push_back(r5);
  m.rings[b].push_back(r5);
  m.rings[a].push_back(r6);
  m.rings[b].push_back(r6);
  return m;
}

Mol makeFused2()
{
  // All four atoms ring-bound: smallest shared ring 5, largest 6, so no
  // size rule fires (defaults nrot=1, phi=0, f1=torsf(1)).
  Mol m;
  int a = addAtom(m, 6, 3, 0.0, 0, 0, -0.77, 0, 0);
  int b = addAtom(m, 6, 3, 0.0, 0, 0, 0.77, 0, 0);
  int k = addAtom(m, 6, 3, 0.0, 0, 0, -1.9, 0.9, 0);
  int l = addAtom(m, 6, 3, 0.0, 0, 0, 1.9, -0.9, 0);
  link(m, a, b);
  link(m, a, k);
  link(m, b, l);
  Ring r5, r6;
  r5.members = { a, b, k, l, 100 };
  r6.members = { a, b, k, l, 200, 201 };
  m.rings[a].push_back(r5);
  m.rings[b].push_back(r5);
  m.rings[k].push_back(r5);
  m.rings[l].push_back(r5);
  m.rings[a].push_back(r6);
  m.rings[b].push_back(r6);
  m.rings[k].push_back(r6);
  m.rings[l].push_back(r6);
  return m;
}

struct Run
{
  std::vector<Torsion> torsions;
  std::vector<TorsionTerm> terms;
};

Run runMol(Mol m, const GffData& param, const GffGenerator& gen)
{
  finish(m, param);
  Run r;
  Environment env;
  int n = static_cast<int>(m.atoms.size());
  std::vector<int> numbers(n);
  for (int i = 0; i < n; ++i)
    numbers[i] = m.atoms[i].element;
  if (!buildTorsions(m.bonds, m.atoms, m.neighbours, m.xyz, m.rings,
                     param.group, param.metal, param.tors, param.tors2,
                     param, gen, r.torsions, r.terms, env))
    ADD_FAILURE() << env.errorMessage();
  return r;
}

} // namespace

TEST(TorsionTest, dihedralNearLinear)
{
  // ii=(0,0,0) jj=(1,0,0): kk=(-1,0,0) gives angle(jj,ii,kk)=180;
  // kk=(0,1,0), ll=(1,1,0) gives 90/90.
  std::vector<double> xyz = { 0, 0, 0, 1, 0, 0, -1, 0, 0, 1, 1, 0,
                              0, 1, 0 };
  EXPECT_TRUE(dihedralNearLinear(xyz, 0, 1, 2, 3));
  EXPECT_FALSE(dihedralNearLinear(xyz, 0, 1, 4, 3));
}

TEST(TorsionTest, ringHelpers)
{
  Ring r5, r6;
  r5.members = { 0, 1, 2, 3, 4 };
  r6.members = { 0, 1, 5, 6, 7, 8 };
  std::vector<Ring> a = { r5, r6 }, b = { r5, r6 }, c = { r5 }, d = { r5 };
  EXPECT_EQ(smallestRingTorsion(a, b, c, d, 0, 1, 2, 3), 5);
  EXPECT_EQ(largestRingTorsion(a, b, c, d, 0, 1, 2, 3), 5);
  std::vector<Ring> none;
  EXPECT_EQ(smallestRingTorsion(a, b, none, d, 0, 1, 2, 3), 0);
  EXPECT_EQ(largestRingTorsion(a, b, none, d, 0, 1, 2, 3), 0);
}

TEST(TorsionTest, alphaCarbonyl)
{
  // Acetone-like: keto C(0) with terminal O(1), methyl C(2).
  std::vector<std::vector<int>> nb(4);
  nb[0] = { 1, 2 };
  nb[1] = { 0 };
  nb[2] = { 0, 3 };
  nb[3] = { 2 };
  std::vector<int> numbers = { 6, 8, 6, 1 };
  std::vector<int> hyb = { 2, 2, 3, 0 };
  std::vector<int> pi = { 1, 1, 0, 0 };
  EXPECT_TRUE(alphaCarbonyl(nb, numbers, hyb, pi, 0, 2));
  EXPECT_TRUE(alphaCarbonyl(nb, numbers, hyb, pi, 2, 0));
  pi[0] = 0;
  EXPECT_FALSE(alphaCarbonyl(nb, numbers, hyb, pi, 0, 2));
}

TEST(TorsionTest, ethane)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeEthane(), param, gen);
  ASSERT_EQ(r.torsions.size(), 9u); // outers are H: no sp3 extra
  { // torsion 0: (5,0,1,2)
    EXPECT_EQ(r.torsions[0].outer1, 5);
    EXPECT_EQ(r.torsions[0].center1, 0);
    EXPECT_EQ(r.torsions[0].center2, 1);
    EXPECT_EQ(r.torsions[0].outer2, 2);
    EXPECT_EQ(r.torsions[0].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[0].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 2.15104531043385117e-01);
  }
  { // torsion 1: (6,0,1,2)
    EXPECT_EQ(r.torsions[1].outer1, 6);
    EXPECT_EQ(r.torsions[1].center1, 0);
    EXPECT_EQ(r.torsions[1].center2, 1);
    EXPECT_EQ(r.torsions[1].outer2, 2);
    EXPECT_EQ(r.torsions[1].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[1].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[1].forceConstant, 2.15104531043385117e-01);
  }
  { // torsion 2: (7,0,1,2)
    EXPECT_EQ(r.torsions[2].outer1, 7);
    EXPECT_EQ(r.torsions[2].center1, 0);
    EXPECT_EQ(r.torsions[2].center2, 1);
    EXPECT_EQ(r.torsions[2].outer2, 2);
    EXPECT_EQ(r.torsions[2].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[2].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[2].forceConstant, 2.15104531043385117e-01);
  }
  { // torsion 3: (5,0,1,3)
    EXPECT_EQ(r.torsions[3].outer1, 5);
    EXPECT_EQ(r.torsions[3].center1, 0);
    EXPECT_EQ(r.torsions[3].center2, 1);
    EXPECT_EQ(r.torsions[3].outer2, 3);
    EXPECT_EQ(r.torsions[3].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[3].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[3].forceConstant, 2.15104531043385117e-01);
  }
  { // torsion 4: (6,0,1,3)
    EXPECT_EQ(r.torsions[4].outer1, 6);
    EXPECT_EQ(r.torsions[4].center1, 0);
    EXPECT_EQ(r.torsions[4].center2, 1);
    EXPECT_EQ(r.torsions[4].outer2, 3);
    EXPECT_EQ(r.torsions[4].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[4].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[4].forceConstant, 2.15104531043385117e-01);
  }
  { // torsion 5: (7,0,1,3)
    EXPECT_EQ(r.torsions[5].outer1, 7);
    EXPECT_EQ(r.torsions[5].center1, 0);
    EXPECT_EQ(r.torsions[5].center2, 1);
    EXPECT_EQ(r.torsions[5].outer2, 3);
    EXPECT_EQ(r.torsions[5].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[5].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[5].forceConstant, 2.15104531043385117e-01);
  }
  { // torsion 6: (5,0,1,4)
    EXPECT_EQ(r.torsions[6].outer1, 5);
    EXPECT_EQ(r.torsions[6].center1, 0);
    EXPECT_EQ(r.torsions[6].center2, 1);
    EXPECT_EQ(r.torsions[6].outer2, 4);
    EXPECT_EQ(r.torsions[6].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[6].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[6].forceConstant, 2.15104531043385117e-01);
  }
  { // torsion 7: (6,0,1,4)
    EXPECT_EQ(r.torsions[7].outer1, 6);
    EXPECT_EQ(r.torsions[7].center1, 0);
    EXPECT_EQ(r.torsions[7].center2, 1);
    EXPECT_EQ(r.torsions[7].outer2, 4);
    EXPECT_EQ(r.torsions[7].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[7].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[7].forceConstant, 2.15104531043385117e-01);
  }
  { // torsion 8: (7,0,1,4)
    EXPECT_EQ(r.torsions[8].outer1, 7);
    EXPECT_EQ(r.torsions[8].center1, 0);
    EXPECT_EQ(r.torsions[8].center2, 1);
    EXPECT_EQ(r.torsions[8].outer2, 4);
    EXPECT_EQ(r.torsions[8].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[8].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[8].forceConstant, 2.15104531043385117e-01);
  }
}

TEST(TorsionTest, ethylene)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeEthylene(), param, gen);
  ASSERT_EQ(r.torsions.size(), 4u);
  { // torsion 0: (4,0,1,2)
    EXPECT_EQ(r.torsions[0].outer1, 4);
    EXPECT_EQ(r.torsions[0].center1, 0);
    EXPECT_EQ(r.torsions[0].center2, 1);
    EXPECT_EQ(r.torsions[0].outer2, 2);
    EXPECT_EQ(r.torsions[0].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[0].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 1.69590891266774912e+00);
  }
  { // torsion 1: (5,0,1,2)
    EXPECT_EQ(r.torsions[1].outer1, 5);
    EXPECT_EQ(r.torsions[1].center1, 0);
    EXPECT_EQ(r.torsions[1].center2, 1);
    EXPECT_EQ(r.torsions[1].outer2, 2);
    EXPECT_EQ(r.torsions[1].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[1].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[1].forceConstant, 1.69590891266774912e+00);
  }
  { // torsion 2: (4,0,1,3)
    EXPECT_EQ(r.torsions[2].outer1, 4);
    EXPECT_EQ(r.torsions[2].center1, 0);
    EXPECT_EQ(r.torsions[2].center2, 1);
    EXPECT_EQ(r.torsions[2].outer2, 3);
    EXPECT_EQ(r.torsions[2].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[2].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[2].forceConstant, 1.69590891266774912e+00);
  }
  { // torsion 3: (5,0,1,3)
    EXPECT_EQ(r.torsions[3].outer1, 5);
    EXPECT_EQ(r.torsions[3].center1, 0);
    EXPECT_EQ(r.torsions[3].center2, 1);
    EXPECT_EQ(r.torsions[3].outer2, 3);
    EXPECT_EQ(r.torsions[3].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[3].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[3].forceConstant, 1.69590891266774912e+00);
  }
}

TEST(TorsionTest, benzene)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeBenzene(), param, gen);
  // 6 ring bonds x (2 H-ends... check count in run
  ASSERT_EQ(r.torsions.size(), 24u);
  { // torsion 0: (2,0,1,6)
    EXPECT_EQ(r.torsions[0].outer1, 2);
    EXPECT_EQ(r.torsions[0].center1, 0);
    EXPECT_EQ(r.torsions[0].center2, 1);
    EXPECT_EQ(r.torsions[0].outer2, 6);
    EXPECT_EQ(r.torsions[0].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[0].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 1.04050586425878988e+00);
  }
  { // torsion 1: (7,0,1,6)
    EXPECT_EQ(r.torsions[1].outer1, 7);
    EXPECT_EQ(r.torsions[1].center1, 0);
    EXPECT_EQ(r.torsions[1].center2, 1);
    EXPECT_EQ(r.torsions[1].outer2, 6);
    EXPECT_EQ(r.torsions[1].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[1].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[1].forceConstant, 1.48279627172098905e+00);
  }
  { // torsion 2: (2,0,1,5)
    EXPECT_EQ(r.torsions[2].outer1, 2);
    EXPECT_EQ(r.torsions[2].center1, 0);
    EXPECT_EQ(r.torsions[2].center2, 1);
    EXPECT_EQ(r.torsions[2].outer2, 5);
    EXPECT_EQ(r.torsions[2].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[2].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[2].forceConstant, 7.30142416867803457e-01);
  }
  { // torsion 3: (7,0,1,5)
    EXPECT_EQ(r.torsions[3].outer1, 7);
    EXPECT_EQ(r.torsions[3].center1, 0);
    EXPECT_EQ(r.torsions[3].center2, 1);
    EXPECT_EQ(r.torsions[3].outer2, 5);
    EXPECT_EQ(r.torsions[3].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[3].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[3].forceConstant, 1.04050586425878988e+00);
  }
  { // torsion 4: (3,1,2,0)
    EXPECT_EQ(r.torsions[4].outer1, 3);
    EXPECT_EQ(r.torsions[4].center1, 1);
    EXPECT_EQ(r.torsions[4].center2, 2);
    EXPECT_EQ(r.torsions[4].outer2, 0);
    EXPECT_EQ(r.torsions[4].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[4].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[4].forceConstant, 7.30142416867803457e-01);
  }
  { // torsion 5: (8,1,2,0)
    EXPECT_EQ(r.torsions[5].outer1, 8);
    EXPECT_EQ(r.torsions[5].center1, 1);
    EXPECT_EQ(r.torsions[5].center2, 2);
    EXPECT_EQ(r.torsions[5].outer2, 0);
    EXPECT_EQ(r.torsions[5].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[5].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[5].forceConstant, 1.04050586425878988e+00);
  }
  { // torsion 6: (3,1,2,7)
    EXPECT_EQ(r.torsions[6].outer1, 3);
    EXPECT_EQ(r.torsions[6].center1, 1);
    EXPECT_EQ(r.torsions[6].center2, 2);
    EXPECT_EQ(r.torsions[6].outer2, 7);
    EXPECT_EQ(r.torsions[6].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[6].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[6].forceConstant, 1.04050586425878988e+00);
  }
  { // torsion 7: (8,1,2,7)
    EXPECT_EQ(r.torsions[7].outer1, 8);
    EXPECT_EQ(r.torsions[7].center1, 1);
    EXPECT_EQ(r.torsions[7].center2, 2);
    EXPECT_EQ(r.torsions[7].outer2, 7);
    EXPECT_EQ(r.torsions[7].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[7].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[7].forceConstant, 1.48279627172098905e+00);
  }
  { // torsion 8: (4,2,3,1)
    EXPECT_EQ(r.torsions[8].outer1, 4);
    EXPECT_EQ(r.torsions[8].center1, 2);
    EXPECT_EQ(r.torsions[8].center2, 3);
    EXPECT_EQ(r.torsions[8].outer2, 1);
    EXPECT_EQ(r.torsions[8].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[8].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[8].forceConstant, 7.30142416867803457e-01);
  }
  { // torsion 9: (9,2,3,1)
    EXPECT_EQ(r.torsions[9].outer1, 9);
    EXPECT_EQ(r.torsions[9].center1, 2);
    EXPECT_EQ(r.torsions[9].center2, 3);
    EXPECT_EQ(r.torsions[9].outer2, 1);
    EXPECT_EQ(r.torsions[9].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[9].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[9].forceConstant, 1.04050586425878988e+00);
  }
  { // torsion 10: (4,2,3,8)
    EXPECT_EQ(r.torsions[10].outer1, 4);
    EXPECT_EQ(r.torsions[10].center1, 2);
    EXPECT_EQ(r.torsions[10].center2, 3);
    EXPECT_EQ(r.torsions[10].outer2, 8);
    EXPECT_EQ(r.torsions[10].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[10].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[10].forceConstant, 1.04050586425878988e+00);
  }
  { // torsion 11: (9,2,3,8)
    EXPECT_EQ(r.torsions[11].outer1, 9);
    EXPECT_EQ(r.torsions[11].center1, 2);
    EXPECT_EQ(r.torsions[11].center2, 3);
    EXPECT_EQ(r.torsions[11].outer2, 8);
    EXPECT_EQ(r.torsions[11].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[11].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[11].forceConstant, 1.48279627172098905e+00);
  }
  { // torsion 12: (5,3,4,2)
    EXPECT_EQ(r.torsions[12].outer1, 5);
    EXPECT_EQ(r.torsions[12].center1, 3);
    EXPECT_EQ(r.torsions[12].center2, 4);
    EXPECT_EQ(r.torsions[12].outer2, 2);
    EXPECT_EQ(r.torsions[12].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[12].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[12].forceConstant, 7.30142416867803457e-01);
  }
  { // torsion 13: (10,3,4,2)
    EXPECT_EQ(r.torsions[13].outer1, 10);
    EXPECT_EQ(r.torsions[13].center1, 3);
    EXPECT_EQ(r.torsions[13].center2, 4);
    EXPECT_EQ(r.torsions[13].outer2, 2);
    EXPECT_EQ(r.torsions[13].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[13].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[13].forceConstant, 1.04050586425878988e+00);
  }
  { // torsion 14: (5,3,4,9)
    EXPECT_EQ(r.torsions[14].outer1, 5);
    EXPECT_EQ(r.torsions[14].center1, 3);
    EXPECT_EQ(r.torsions[14].center2, 4);
    EXPECT_EQ(r.torsions[14].outer2, 9);
    EXPECT_EQ(r.torsions[14].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[14].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[14].forceConstant, 1.04050586425878988e+00);
  }
  { // torsion 15: (10,3,4,9)
    EXPECT_EQ(r.torsions[15].outer1, 10);
    EXPECT_EQ(r.torsions[15].center1, 3);
    EXPECT_EQ(r.torsions[15].center2, 4);
    EXPECT_EQ(r.torsions[15].outer2, 9);
    EXPECT_EQ(r.torsions[15].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[15].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[15].forceConstant, 1.48279627172098905e+00);
  }
  { // torsion 16: (0,4,5,3)
    EXPECT_EQ(r.torsions[16].outer1, 0);
    EXPECT_EQ(r.torsions[16].center1, 4);
    EXPECT_EQ(r.torsions[16].center2, 5);
    EXPECT_EQ(r.torsions[16].outer2, 3);
    EXPECT_EQ(r.torsions[16].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[16].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[16].forceConstant, 7.30142416867803457e-01);
  }
  { // torsion 17: (11,4,5,3)
    EXPECT_EQ(r.torsions[17].outer1, 11);
    EXPECT_EQ(r.torsions[17].center1, 4);
    EXPECT_EQ(r.torsions[17].center2, 5);
    EXPECT_EQ(r.torsions[17].outer2, 3);
    EXPECT_EQ(r.torsions[17].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[17].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[17].forceConstant, 1.04050586425878988e+00);
  }
  { // torsion 18: (0,4,5,10)
    EXPECT_EQ(r.torsions[18].outer1, 0);
    EXPECT_EQ(r.torsions[18].center1, 4);
    EXPECT_EQ(r.torsions[18].center2, 5);
    EXPECT_EQ(r.torsions[18].outer2, 10);
    EXPECT_EQ(r.torsions[18].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[18].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[18].forceConstant, 1.04050586425878988e+00);
  }
  { // torsion 19: (11,4,5,10)
    EXPECT_EQ(r.torsions[19].outer1, 11);
    EXPECT_EQ(r.torsions[19].center1, 4);
    EXPECT_EQ(r.torsions[19].center2, 5);
    EXPECT_EQ(r.torsions[19].outer2, 10);
    EXPECT_EQ(r.torsions[19].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[19].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[19].forceConstant, 1.48279627172098905e+00);
  }
  { // torsion 20: (1,5,0,4)
    EXPECT_EQ(r.torsions[20].outer1, 1);
    EXPECT_EQ(r.torsions[20].center1, 5);
    EXPECT_EQ(r.torsions[20].center2, 0);
    EXPECT_EQ(r.torsions[20].outer2, 4);
    EXPECT_EQ(r.torsions[20].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[20].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[20].forceConstant, 7.30142416867803457e-01);
  }
  { // torsion 21: (6,5,0,4)
    EXPECT_EQ(r.torsions[21].outer1, 6);
    EXPECT_EQ(r.torsions[21].center1, 5);
    EXPECT_EQ(r.torsions[21].center2, 0);
    EXPECT_EQ(r.torsions[21].outer2, 4);
    EXPECT_EQ(r.torsions[21].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[21].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[21].forceConstant, 1.04050586425878988e+00);
  }
  { // torsion 22: (1,5,0,11)
    EXPECT_EQ(r.torsions[22].outer1, 1);
    EXPECT_EQ(r.torsions[22].center1, 5);
    EXPECT_EQ(r.torsions[22].center2, 0);
    EXPECT_EQ(r.torsions[22].outer2, 11);
    EXPECT_EQ(r.torsions[22].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[22].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[22].forceConstant, 1.04050586425878988e+00);
  }
  { // torsion 23: (6,5,0,11)
    EXPECT_EQ(r.torsions[23].outer1, 6);
    EXPECT_EQ(r.torsions[23].center1, 5);
    EXPECT_EQ(r.torsions[23].center2, 0);
    EXPECT_EQ(r.torsions[23].outer2, 11);
    EXPECT_EQ(r.torsions[23].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[23].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[23].forceConstant, 1.48279627172098905e+00);
  }
}

TEST(TorsionTest, cyclopropane)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeCyclopropane(), param, gen);
  ASSERT_EQ(r.torsions.size(), 9u);
  { // torsion 0: (4,0,1,2)
    EXPECT_EQ(r.torsions[0].outer1, 4);
    EXPECT_EQ(r.torsions[0].center1, 0);
    EXPECT_EQ(r.torsions[0].center2, 1);
    EXPECT_EQ(r.torsions[0].outer2, 2);
    EXPECT_EQ(r.torsions[0].multiplicity, 1);
    EXPECT_DOUBLE_EQ(r.terms[0].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 4.10950655831036904e-02);
  }
  { // torsion 1: (2,0,1,3)
    EXPECT_EQ(r.torsions[1].outer1, 2);
    EXPECT_EQ(r.torsions[1].center1, 0);
    EXPECT_EQ(r.torsions[1].center2, 1);
    EXPECT_EQ(r.torsions[1].outer2, 3);
    EXPECT_EQ(r.torsions[1].multiplicity, 1);
    EXPECT_DOUBLE_EQ(r.terms[1].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[1].forceConstant, 4.10950655831036904e-02);
  }
  { // torsion 2: (4,0,1,3)
    EXPECT_EQ(r.torsions[2].outer1, 4);
    EXPECT_EQ(r.torsions[2].center1, 0);
    EXPECT_EQ(r.torsions[2].center2, 1);
    EXPECT_EQ(r.torsions[2].outer2, 3);
    EXPECT_EQ(r.torsions[2].multiplicity, 1);
    EXPECT_DOUBLE_EQ(r.terms[2].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[2].forceConstant, 5.85634469981228950e-02);
  }
  { // torsion 3: (5,1,2,0)
    EXPECT_EQ(r.torsions[3].outer1, 5);
    EXPECT_EQ(r.torsions[3].center1, 1);
    EXPECT_EQ(r.torsions[3].center2, 2);
    EXPECT_EQ(r.torsions[3].outer2, 0);
    EXPECT_EQ(r.torsions[3].multiplicity, 1);
    EXPECT_DOUBLE_EQ(r.terms[3].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[3].forceConstant, 4.10950655831036904e-02);
  }
  { // torsion 4: (0,1,2,4)
    EXPECT_EQ(r.torsions[4].outer1, 0);
    EXPECT_EQ(r.torsions[4].center1, 1);
    EXPECT_EQ(r.torsions[4].center2, 2);
    EXPECT_EQ(r.torsions[4].outer2, 4);
    EXPECT_EQ(r.torsions[4].multiplicity, 1);
    EXPECT_DOUBLE_EQ(r.terms[4].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[4].forceConstant, 4.10950655831036904e-02);
  }
  { // torsion 5: (5,1,2,4)
    EXPECT_EQ(r.torsions[5].outer1, 5);
    EXPECT_EQ(r.torsions[5].center1, 1);
    EXPECT_EQ(r.torsions[5].center2, 2);
    EXPECT_EQ(r.torsions[5].outer2, 4);
    EXPECT_EQ(r.torsions[5].multiplicity, 1);
    EXPECT_DOUBLE_EQ(r.terms[5].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[5].forceConstant, 5.85634469981228950e-02);
  }
  { // torsion 6: (3,2,0,1)
    EXPECT_EQ(r.torsions[6].outer1, 3);
    EXPECT_EQ(r.torsions[6].center1, 2);
    EXPECT_EQ(r.torsions[6].center2, 0);
    EXPECT_EQ(r.torsions[6].outer2, 1);
    EXPECT_EQ(r.torsions[6].multiplicity, 1);
    EXPECT_DOUBLE_EQ(r.terms[6].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[6].forceConstant, 4.10950655831036904e-02);
  }
  { // torsion 7: (1,2,0,5)
    EXPECT_EQ(r.torsions[7].outer1, 1);
    EXPECT_EQ(r.torsions[7].center1, 2);
    EXPECT_EQ(r.torsions[7].center2, 0);
    EXPECT_EQ(r.torsions[7].outer2, 5);
    EXPECT_EQ(r.torsions[7].multiplicity, 1);
    EXPECT_DOUBLE_EQ(r.terms[7].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[7].forceConstant, 4.10950655831036904e-02);
  }
  { // torsion 8: (3,2,0,5)
    EXPECT_EQ(r.torsions[8].outer1, 3);
    EXPECT_EQ(r.torsions[8].center1, 2);
    EXPECT_EQ(r.torsions[8].center2, 0);
    EXPECT_EQ(r.torsions[8].outer2, 5);
    EXPECT_EQ(r.torsions[8].multiplicity, 1);
    EXPECT_DOUBLE_EQ(r.terms[8].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[8].forceConstant, 5.85634469981228950e-02);
  }
}

TEST(TorsionTest, lactamCb7)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeLactam(), param, gen);
  ASSERT_EQ(r.torsions.size(), 10u);
  { // torsion 0: (3,0,1,5)
    EXPECT_EQ(r.torsions[0].outer1, 3);
    EXPECT_EQ(r.torsions[0].center1, 0);
    EXPECT_EQ(r.torsions[0].center2, 1);
    EXPECT_EQ(r.torsions[0].outer2, 5);
    EXPECT_EQ(r.torsions[0].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[0].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 4.38969743638367771e-01);
  }
  { // torsion 1: (2,0,1,5)
    EXPECT_EQ(r.torsions[1].outer1, 2);
    EXPECT_EQ(r.torsions[1].center1, 0);
    EXPECT_EQ(r.torsions[1].center2, 1);
    EXPECT_EQ(r.torsions[1].outer2, 5);
    EXPECT_EQ(r.torsions[1].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[1].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[1].forceConstant, 5.39894257232641728e-01);
  }
  { // torsion 2: (3,0,1,6)
    EXPECT_EQ(r.torsions[2].outer1, 3);
    EXPECT_EQ(r.torsions[2].center1, 0);
    EXPECT_EQ(r.torsions[2].center2, 1);
    EXPECT_EQ(r.torsions[2].outer2, 6);
    EXPECT_EQ(r.torsions[2].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[2].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[2].forceConstant, 5.91042572400507504e-01);
  }
  { // torsion 3: (2,0,1,6)
    EXPECT_EQ(r.torsions[3].outer1, 2);
    EXPECT_EQ(r.torsions[3].center1, 0);
    EXPECT_EQ(r.torsions[3].center2, 1);
    EXPECT_EQ(r.torsions[3].outer2, 6);
    EXPECT_EQ(r.torsions[3].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[3].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[3].forceConstant, 7.26930489500714483e-01);
  }
  { // torsion 4: (4,1,3,0)
    EXPECT_EQ(r.torsions[4].outer1, 4);
    EXPECT_EQ(r.torsions[4].center1, 1);
    EXPECT_EQ(r.torsions[4].center2, 3);
    EXPECT_EQ(r.torsions[4].outer2, 0);
    EXPECT_EQ(r.torsions[4].multiplicity, 1);
    EXPECT_DOUBLE_EQ(r.terms[4].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[4].forceConstant, 9.77954733009831073e-02);
  }
  { // torsion 5: (4,1,3,2)
    EXPECT_EQ(r.torsions[5].outer1, 4);
    EXPECT_EQ(r.torsions[5].center1, 1);
    EXPECT_EQ(r.torsions[5].center2, 3);
    EXPECT_EQ(r.torsions[5].outer2, 2);
    EXPECT_EQ(r.torsions[5].multiplicity, 1);
    EXPECT_DOUBLE_EQ(r.terms[5].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[5].forceConstant, 1.56246493278137971e-01);
  }
  { // torsion 6: (5,3,4,1)
    EXPECT_EQ(r.torsions[6].outer1, 5);
    EXPECT_EQ(r.torsions[6].center1, 3);
    EXPECT_EQ(r.torsions[6].center2, 4);
    EXPECT_EQ(r.torsions[6].outer2, 1);
    EXPECT_EQ(r.torsions[6].multiplicity, 1);
    EXPECT_DOUBLE_EQ(r.terms[6].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[6].forceConstant, 9.23294034462790597e-02);
  }
  { // torsion 7: (0,4,5,3)
    EXPECT_EQ(r.torsions[7].outer1, 0);
    EXPECT_EQ(r.torsions[7].center1, 4);
    EXPECT_EQ(r.torsions[7].center2, 5);
    EXPECT_EQ(r.torsions[7].outer2, 3);
    EXPECT_EQ(r.torsions[7].multiplicity, 1);
    EXPECT_DOUBLE_EQ(r.terms[7].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[7].forceConstant, 7.52272871546023902e-02);
  }
  { // torsion 8: (1,5,0,4)
    EXPECT_EQ(r.torsions[8].outer1, 1);
    EXPECT_EQ(r.torsions[8].center1, 5);
    EXPECT_EQ(r.torsions[8].center2, 0);
    EXPECT_EQ(r.torsions[8].outer2, 4);
    EXPECT_EQ(r.torsions[8].multiplicity, 1);
    EXPECT_DOUBLE_EQ(r.terms[8].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[8].forceConstant, 1.07833861824351951e-01);
  }
  { // torsion 9: (6,5,0,4)
    EXPECT_EQ(r.torsions[9].outer1, 6);
    EXPECT_EQ(r.torsions[9].center1, 5);
    EXPECT_EQ(r.torsions[9].center2, 0);
    EXPECT_EQ(r.torsions[9].outer2, 4);
    EXPECT_EQ(r.torsions[9].multiplicity, 1);
    EXPECT_DOUBLE_EQ(r.terms[9].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[9].forceConstant, 1.53671068824131951e-01);
  }
}

TEST(TorsionTest, terminalMethylRule)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeMchHcch(), param, gen);
  ASSERT_EQ(r.torsions.size(), 6u);
  { // torsion 0: (3,0,1,2)
    EXPECT_EQ(r.torsions[0].outer1, 3);
    EXPECT_EQ(r.torsions[0].center1, 0);
    EXPECT_EQ(r.torsions[0].center2, 1);
    EXPECT_EQ(r.torsions[0].outer2, 2);
    EXPECT_EQ(r.torsions[0].multiplicity, 6);
    EXPECT_DOUBLE_EQ(r.terms[0].phase, 5.23598775598298816e-01);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 5.85634469981228950e-02);
  }
  { // torsion 1: (5,0,1,2)
    EXPECT_EQ(r.torsions[1].outer1, 5);
    EXPECT_EQ(r.torsions[1].center1, 0);
    EXPECT_EQ(r.torsions[1].center2, 1);
    EXPECT_EQ(r.torsions[1].outer2, 2);
    EXPECT_EQ(r.torsions[1].multiplicity, 6);
    EXPECT_DOUBLE_EQ(r.terms[1].phase, 5.23598775598298816e-01);
    EXPECT_DOUBLE_EQ(r.terms[1].forceConstant, 4.79277072459354872e-02);
  }
  { // torsion 2: (3,0,1,4)
    EXPECT_EQ(r.torsions[2].outer1, 3);
    EXPECT_EQ(r.torsions[2].center1, 0);
    EXPECT_EQ(r.torsions[2].center2, 1);
    EXPECT_EQ(r.torsions[2].outer2, 4);
    EXPECT_EQ(r.torsions[2].multiplicity, 6);
    EXPECT_DOUBLE_EQ(r.terms[2].phase, 5.23598775598298816e-01);
    EXPECT_DOUBLE_EQ(r.terms[2].forceConstant, 4.79277072459354872e-02);
  }
  { // torsion 3: (5,0,1,4)
    EXPECT_EQ(r.torsions[3].outer1, 5);
    EXPECT_EQ(r.torsions[3].center1, 0);
    EXPECT_EQ(r.torsions[3].center2, 1);
    EXPECT_EQ(r.torsions[3].outer2, 4);
    EXPECT_EQ(r.torsions[3].multiplicity, 6);
    EXPECT_DOUBLE_EQ(r.terms[3].phase, 5.23598775598298816e-01);
    EXPECT_DOUBLE_EQ(r.terms[3].forceConstant, 3.92235300276249851e-02);
  }
  { // torsion 4: (3,0,1,6)
    EXPECT_EQ(r.torsions[4].outer1, 3);
    EXPECT_EQ(r.torsions[4].center1, 0);
    EXPECT_EQ(r.torsions[4].center2, 1);
    EXPECT_EQ(r.torsions[4].outer2, 6);
    EXPECT_EQ(r.torsions[4].multiplicity, 6);
    EXPECT_DOUBLE_EQ(r.terms[4].phase, 5.23598775598298816e-01);
    EXPECT_DOUBLE_EQ(r.terms[4].forceConstant, 4.79277072459354872e-02);
  }
  { // torsion 5: (5,0,1,6)
    EXPECT_EQ(r.torsions[5].outer1, 5);
    EXPECT_EQ(r.torsions[5].center1, 0);
    EXPECT_EQ(r.torsions[5].center2, 1);
    EXPECT_EQ(r.torsions[5].outer2, 6);
    EXPECT_EQ(r.torsions[5].multiplicity, 6);
    EXPECT_DOUBLE_EQ(r.terms[5].phase, 5.23598775598298816e-01);
    EXPECT_DOUBLE_EQ(r.terms[5].forceConstant, 3.92235300276249851e-02);
  }
}

TEST(TorsionTest, acetoneAlpha)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeAcetone(), param, gen);
  ASSERT_EQ(r.torsions.size(), 4u);
  { // torsion 0: (4,0,2,1)
    EXPECT_EQ(r.torsions[0].outer1, 4);
    EXPECT_EQ(r.torsions[0].center1, 0);
    EXPECT_EQ(r.torsions[0].center2, 2);
    EXPECT_EQ(r.torsions[0].outer2, 1);
    EXPECT_EQ(r.torsions[0].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[0].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 1.10417156777052522e-01);
  }
  { // torsion 1: (4,0,2,3)
    EXPECT_EQ(r.torsions[1].outer1, 4);
    EXPECT_EQ(r.torsions[1].center1, 0);
    EXPECT_EQ(r.torsions[1].center2, 2);
    EXPECT_EQ(r.torsions[1].outer2, 3);
    EXPECT_EQ(r.torsions[1].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[1].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[1].forceConstant, 8.97764522485269972e-02);
  }
  { // torsion 2: (5,0,3,1)
    EXPECT_EQ(r.torsions[2].outer1, 5);
    EXPECT_EQ(r.torsions[2].center1, 0);
    EXPECT_EQ(r.torsions[2].center2, 3);
    EXPECT_EQ(r.torsions[2].outer2, 1);
    EXPECT_EQ(r.torsions[2].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[2].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[2].forceConstant, 1.10417156777052522e-01);
  }
  { // torsion 3: (5,0,3,2)
    EXPECT_EQ(r.torsions[3].outer1, 5);
    EXPECT_EQ(r.torsions[3].center1, 0);
    EXPECT_EQ(r.torsions[3].center2, 3);
    EXPECT_EQ(r.torsions[3].outer2, 2);
    EXPECT_EQ(r.torsions[3].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[3].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[3].forceConstant, 8.97764522485269972e-02);
  }
}

TEST(TorsionTest, amideAdjacenct)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeNMethylacetamide(), param, gen);
  ASSERT_EQ(r.torsions.size(), 6u);
  { // torsion 0: (2,0,1,3)
    EXPECT_EQ(r.torsions[0].outer1, 2);
    EXPECT_EQ(r.torsions[0].center1, 0);
    EXPECT_EQ(r.torsions[0].center2, 1);
    EXPECT_EQ(r.torsions[0].outer2, 3);
    EXPECT_EQ(r.torsions[0].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[0].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 1.02020142842162431e-01);
  }
  { // torsion 1: (2,0,1,6)
    EXPECT_EQ(r.torsions[1].outer1, 2);
    EXPECT_EQ(r.torsions[1].center1, 0);
    EXPECT_EQ(r.torsions[1].center2, 1);
    EXPECT_EQ(r.torsions[1].outer2, 6);
    EXPECT_EQ(r.torsions[1].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[1].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[1].forceConstant, 1.45386097900142897e-01);
  }
  { // torsion 2: (4,0,3,1)
    EXPECT_EQ(r.torsions[2].outer1, 4);
    EXPECT_EQ(r.torsions[2].center1, 0);
    EXPECT_EQ(r.torsions[2].center2, 3);
    EXPECT_EQ(r.torsions[2].outer2, 1);
    EXPECT_EQ(r.torsions[2].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[2].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[2].forceConstant, 3.31910197363283477e-02);
  }
  { // torsion 3: (5,0,3,1)
    EXPECT_EQ(r.torsions[3].outer1, 5);
    EXPECT_EQ(r.torsions[3].center1, 0);
    EXPECT_EQ(r.torsions[3].center2, 3);
    EXPECT_EQ(r.torsions[3].outer2, 1);
    EXPECT_EQ(r.torsions[3].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[3].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[3].forceConstant, 3.31910197363283477e-02);
  }
  { // torsion 4: (4,0,3,6)
    EXPECT_EQ(r.torsions[4].outer1, 4);
    EXPECT_EQ(r.torsions[4].center1, 0);
    EXPECT_EQ(r.torsions[4].center2, 3);
    EXPECT_EQ(r.torsions[4].outer2, 6);
    EXPECT_EQ(r.torsions[4].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[4].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[4].forceConstant, 4.46894255694229720e-02);
  }
  { // torsion 5: (5,0,3,6)
    EXPECT_EQ(r.torsions[5].outer1, 5);
    EXPECT_EQ(r.torsions[5].center1, 0);
    EXPECT_EQ(r.torsions[5].center2, 3);
    EXPECT_EQ(r.torsions[5].outer2, 6);
    EXPECT_EQ(r.torsions[5].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[5].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[5].forceConstant, 4.46894255694229720e-02);
  }
}

TEST(TorsionTest, ethylamineHalving)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeEthylamine(), param, gen);
  ASSERT_EQ(r.torsions.size(), 1u);
  { // torsion 0: (2,0,1,3)
    EXPECT_EQ(r.torsions[0].outer1, 2);
    EXPECT_EQ(r.torsions[0].center1, 0);
    EXPECT_EQ(r.torsions[0].center2, 1);
    EXPECT_EQ(r.torsions[0].outer2, 3);
    EXPECT_EQ(r.torsions[0].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[0].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 6.20009947528828639e-02);
  }
}

TEST(TorsionTest, vinylSilaneHeavy)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeVinylSilane(), param, gen);
  ASSERT_EQ(r.torsions.size(), 2u);
  { // torsion 0: (4,0,1,2)
    EXPECT_EQ(r.torsions[0].outer1, 4);
    EXPECT_EQ(r.torsions[0].center1, 0);
    EXPECT_EQ(r.torsions[0].center2, 1);
    EXPECT_EQ(r.torsions[0].outer2, 2);
    EXPECT_EQ(r.torsions[0].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[0].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 4.54862992019890389e-01);
  }
  { // torsion 1: (4,0,1,3)
    EXPECT_EQ(r.torsions[1].outer1, 4);
    EXPECT_EQ(r.torsions[1].center1, 0);
    EXPECT_EQ(r.torsions[1].center2, 1);
    EXPECT_EQ(r.torsions[1].outer2, 3);
    EXPECT_EQ(r.torsions[1].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[1].phase, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[1].forceConstant, 1.71909058439274398e+00);
  }
}

TEST(TorsionTest, hydrazine)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeHydrazine(), param, gen);
  ASSERT_EQ(r.torsions.size(), 4u);
  { // torsion 0: (4,0,1,2)
    EXPECT_EQ(r.torsions[0].outer1, 4);
    EXPECT_EQ(r.torsions[0].center1, 0);
    EXPECT_EQ(r.torsions[0].center2, 1);
    EXPECT_EQ(r.torsions[0].outer2, 2);
    EXPECT_EQ(r.torsions[0].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[0].phase, 1.04719755119659763e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 4.54023930861501723e-01);
  }
  { // torsion 1: (5,0,1,2)
    EXPECT_EQ(r.torsions[1].outer1, 5);
    EXPECT_EQ(r.torsions[1].center1, 0);
    EXPECT_EQ(r.torsions[1].center2, 1);
    EXPECT_EQ(r.torsions[1].outer2, 2);
    EXPECT_EQ(r.torsions[1].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[1].phase, 1.04719755119659763e+00);
    EXPECT_DOUBLE_EQ(r.terms[1].forceConstant, 4.54023930861501723e-01);
  }
  { // torsion 2: (4,0,1,3)
    EXPECT_EQ(r.torsions[2].outer1, 4);
    EXPECT_EQ(r.torsions[2].center1, 0);
    EXPECT_EQ(r.torsions[2].center2, 1);
    EXPECT_EQ(r.torsions[2].outer2, 3);
    EXPECT_EQ(r.torsions[2].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[2].phase, 1.04719755119659763e+00);
    EXPECT_DOUBLE_EQ(r.terms[2].forceConstant, 4.54023930861501723e-01);
  }
  { // torsion 3: (5,0,1,3)
    EXPECT_EQ(r.torsions[3].outer1, 5);
    EXPECT_EQ(r.torsions[3].center1, 0);
    EXPECT_EQ(r.torsions[3].center2, 1);
    EXPECT_EQ(r.torsions[3].outer2, 3);
    EXPECT_EQ(r.torsions[3].multiplicity, 3);
    EXPECT_DOUBLE_EQ(r.terms[3].phase, 1.04719755119659763e+00);
    EXPECT_DOUBLE_EQ(r.terms[3].forceConstant, 4.54023930861501723e-01);
  }
}

TEST(TorsionTest, peroxide)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makePeroxide(), param, gen);
  ASSERT_EQ(r.torsions.size(), 1u);
  { // torsion 0: (3,0,1,2)
    EXPECT_EQ(r.torsions[0].outer1, 3);
    EXPECT_EQ(r.torsions[0].center1, 0);
    EXPECT_EQ(r.torsions[0].center2, 1);
    EXPECT_EQ(r.torsions[0].outer2, 2);
    EXPECT_EQ(r.torsions[0].multiplicity, 2);
    EXPECT_DOUBLE_EQ(r.terms[0].phase, 1.57079632679489656e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 9.06706195855715658e-01);
  }
}

TEST(TorsionTest, linearSkip)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeLinearSkip(), param, gen);
  ASSERT_EQ(r.torsions.size(), 1u); // (kk,ll)=(k,normal-H); linear pair out
  ASSERT_EQ(r.torsions[0].outer1, 4);
  ASSERT_EQ(r.torsions[0].center1, 0);
  ASSERT_EQ(r.torsions[0].center2, 1);
  ASSERT_EQ(r.torsions[0].outer2, 2);
  ASSERT_EQ(r.torsions[0].multiplicity, 3);
  EXPECT_DOUBLE_EQ(r.terms[0].phase, 3.14159265358979312e+00);
  EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 2.00831460807167905e-01);
}

TEST(TorsionTest, fusedDefaults)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeFused(), param, gen);
  ASSERT_EQ(r.torsions.size(), 1u);
  ASSERT_EQ(r.torsions[0].outer1, 3);
  ASSERT_EQ(r.torsions[0].center1, 0);
  ASSERT_EQ(r.torsions[0].center2, 1);
  ASSERT_EQ(r.torsions[0].outer2, 2);
  ASSERT_EQ(r.torsions[0].multiplicity, 6);
  EXPECT_DOUBLE_EQ(r.terms[0].phase, 5.23598775598298816e-01);
  EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 5.85634469981228950e-02);
}

TEST(TorsionTest, fusedRingMismatch)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeFused2(), param, gen);
  ASSERT_EQ(r.torsions.size(), 1u);
  ASSERT_EQ(r.torsions[0].outer1, 3);
  ASSERT_EQ(r.torsions[0].center1, 0);
  ASSERT_EQ(r.torsions[0].center2, 1);
  ASSERT_EQ(r.torsions[0].outer2, 2);
  ASSERT_EQ(r.torsions[0].multiplicity, 1);
  EXPECT_DOUBLE_EQ(r.terms[0].phase, 0.0);
  EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 1.18653678033530011e-01);
}

TEST(TorsionTest, skips)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  // Acetylene C#C, Sr-C (tors<0), Mg-O (tors=0) give no torsions.
  Mol m;
  int c1 = addAtom(m, 6, 1, 0.0, 0, 1, -0.60, 0, 0);
  int c2 = addAtom(m, 6, 1, 0.0, 0, 1, 0.60, 0, 0);
  link(m, c1, c2);
  int sr = addAtom(m, 38, 0, 0.4, 1, 0, 3.0, 0, 0);
  int cs = addAtom(m, 6, 3, 0.0, 0, 0, 4.5, 0, 0);
  link(m, sr, cs);
  int mg = addAtom(m, 12, 0, 0.4, 1, 0, 6.0, 0, 0);
  int ox = addAtom(m, 8, 3, -0.3, 0, 0, 7.5, 0, 0);
  link(m, mg, ox);
  Run r = runMol(m, param, gen);
  EXPECT_EQ(r.torsions.size(), 0u);
}

} // namespace Xtb
} // namespace Avogadro
