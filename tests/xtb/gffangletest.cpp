/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "avogadro/xtb/COPYING.LESSER".
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/xtb/environment.h>
#include <avogadro/xtb/gffgraph.h>
#include <avogadro/xtb/gfnffangle.h>
#include <avogadro/xtb/gfnffdata.h>
#include <avogadro/xtb/gfnffsetup.h>
#include <avogadro/xtb/gfnfftopo.h>

#include <cmath>
#include <vector>

namespace Avogadro {
namespace Xtb {

// Reference values verified bitwise against the unmodified alist/vangl
// loops in xtb src/gfnff/gfnff_ini.f90 (1499-1540, 1552-1825), compiled
// with gfortran -ffp-contract=off over 40 synthetic molecules (163
// angles). Includes the ringsbend s(m,j) quirk verbatim.

namespace {

struct Mol
{
  std::vector<AngleAtom> atoms;
  std::vector<std::vector<int>> neighbours;
  std::vector<std::vector<Ring>> rings;
  std::vector<double> xyz;
  std::vector<double> pbo;
};

int addAtom(Mol& m, int z, int hyb, double qa, int itag, int imetal, int pi,
            double x, double y, double zc)
{
  AngleAtom a;
  a.element = z;
  a.hyb = hyb;
  a.charge = qa;
  a.itag = itag;
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

void link(Mol& m, int a, int b)
{
  m.neighbours[a].push_back(b);
  m.neighbours[b].push_back(a);
}

void finish(Mol& m)
{
  int n = static_cast<int>(m.atoms.size());
  m.pbo.assign(n * (n + 1) / 2, 0.0);
  for (int i = 0; i < n; ++i) {
    int best = 99;
    for (const Ring& r : m.rings[i])
      best = std::min(best, static_cast<int>(r.members.size()));
    if (!m.rings[i].empty())
      m.atoms[i].ring = best;
  }
}

Mol makeWater()
{
  Mol m;
  int o = addAtom(m, 8, 3, -0.6, 0, 0, 0, 0, 0, 0);
  int h1 = addAtom(m, 1, 0, 0.3, 0, 0, 0, 0.96, 0, 0);
  int h2 = addAtom(m, 1, 0, 0.3, 0, 0, 0, 0.96 * -0.25881904510252074,
                   0.96 * 0.96592582628906831, 0);
  link(m, o, h1);
  link(m, o, h2);
  finish(m);
  return m;
}

Mol makeFormaldehyde()
{
  Mol m;
  int c = addAtom(m, 6, 2, 0.0, 0, 0, 1, 0, 0, 0);
  int o = addAtom(m, 8, 2, 0.0, 0, 0, 1, 1.21, 0, 0);
  int h1 = addAtom(m, 1, 0, 0.0, 0, 0, 0, -0.48, 0.88, 0);
  int h2 = addAtom(m, 1, 0, 0.0, 0, 0, 0, -0.48, -0.88, 0);
  link(m, c, o);
  link(m, c, h1);
  link(m, c, h2);
  finish(m);
  return m;
}

Mol makeMethylcyclopropane()
{
  Mol m;
  // Ring triangle in xy-plane, methyl on atom 0.
  int r0 = addAtom(m, 6, 3, 0.0, 0, 0, 0, 0.0, 1.0, 0);
  int r1 = addAtom(m, 6, 3, 0.0, 0, 0, 0, 0.86602540378443871, -0.5, 0);
  int r2 = addAtom(m, 6, 3, 0.0, 0, 0, 0, -0.86602540378443871, -0.5, 0);
  int me = addAtom(m, 6, 3, 0.0, 0, 0, 0, 0.0, 2.5, 0);
  link(m, r0, r1);
  link(m, r1, r2);
  link(m, r2, r0);
  link(m, r0, me);
  Ring ring;
  ring.members = { r0, r1, r2 };
  m.rings[r0].push_back(ring);
  m.rings[r1].push_back(ring);
  m.rings[r2].push_back(ring);
  finish(m);
  return m;
}

Mol makeFormamide()
{
  Mol m;
  int c = addAtom(m, 6, 2, 0.0, 0, 0, 1, 0, 0, 0);
  int o = addAtom(m, 8, 2, 0.0, 0, 0, 1, 1.23, 0, 0);
  int n = addAtom(m, 7, 3, 0.0, 0, 0, 1, -1.35, 0.2, 0);
  int h1 = addAtom(m, 1, 0, 0.0, 0, 0, 0, -1.7, 1.05, 0);
  int h2 = addAtom(m, 1, 0, 0.0, 0, 0, 0, -1.9, -0.55, 0);
  int hc = addAtom(m, 1, 0, 0.0, 0, 0, 0, 0.0, -1.1, 0);
  link(m, c, o);
  link(m, c, n);
  link(m, n, h1);
  link(m, n, h2);
  link(m, c, hc);
  finish(m);
  return m;
}

Mol makeAnilineN()
{
  Mol m;
  int c = addAtom(m, 6, 2, 0.0, 0, 0, 1, 0, 0, 0);
  int cc = addAtom(m, 6, 2, 0.0, 0, 0, 1, 1.40, 0, 0);
  int n = addAtom(m, 7, 3, -0.2, 0, 0, 1, -1.40, 0.2, 0);
  int h1 = addAtom(m, 1, 0, 0.0, 0, 0, 0, -1.75, 1.0, 0);
  int h2 = addAtom(m, 1, 0, 0.0, 0, 0, 0, -1.9, -0.6, 0);
  link(m, c, cc);
  link(m, c, n);
  link(m, n, h1);
  link(m, n, h2);
  finish(m);
  m.pbo[packedIndex(c, n)] = 0.3;
  return m;
}

Mol makeSf6()
{
  Mol m;
  int s = addAtom(m, 16, 5, 0.0, 0, 0, 0, 0, 0, 0);
  const double pts[6][3] = { { 1.56, 0, 0 }, { -1.56, 0, 0 },
                             { 0, 1.56, 0 }, { 0, -1.56, 0 },
                             { 0, 0, 1.56 }, { 0, 0, -1.56 } };
  for (const auto& p : pts) {
    int f = addAtom(m, 9, 1, 0.0, 0, 0, 0, p[0], p[1], p[2]);
    link(m, s, f);
  }
  finish(m);
  return m;
}

Mol makeFeCo()
{
  Mol m;
  int fe = addAtom(m, 26, 0, 0.2, 0, 2, 0, -1.8, 0, 0);
  int c = addAtom(m, 6, 1, 0.0, 0, 0, 1, 0, 0, 0);
  int o = addAtom(m, 8, 2, 0.0, 0, 0, 1, 1.15, 0, 0);
  link(m, fe, c);
  link(m, c, o);
  finish(m);
  return m;
}

Mol makeMom()
{
  Mol m;
  int m1 = addAtom(m, 26, 0, 0.2, 0, 2, 0, -1.80, 0, 0);
  int o = addAtom(m, 8, 3, -0.3, 0, 0, 0, 0, 0, 0);
  int m2 = addAtom(m, 26, 0, 0.2, 0, 2, 0, 1.80 * 0.99619469809174555,
                   1.80 * 0.087155742747658179, 0);
  link(m, m1, o);
  link(m, o, m2);
  finish(m);
  return m;
}

Mol makeFeCp()
{
  Mol m;
  int fe = addAtom(m, 26, 0, 0.2, 0, 2, 0, 0, 0, 0);
  int c1 = addAtom(m, 6, 2, 0.0, -1, 0, 1, 2.00, 0, 0);
  int c2 = addAtom(m, 6, 2, 0.0, -1, 0, 1, 0, 2.00, 0);
  link(m, fe, c1);
  link(m, fe, c2);
  finish(m);
  return m;
}

Mol makeCaSkip()
{
  Mol m;
  int ca = addAtom(m, 20, 0, 0.5, 0, 1, 0, 0, 0, 0);
  int c1 = addAtom(m, 6, 2, 0.0, 0, 0, 1, 2.40, 0, 0);
  int c2 = addAtom(m, 6, 2, 0.0, 0, 0, 1, 2.40 * 0.64278760968653936,
                   2.40 * 0.76604444311897801, 0);
  link(m, ca, c1);
  link(m, ca, c2);
  finish(m);
  return m;
}

struct Run
{
  std::vector<Angle> angles;
  std::vector<AngleTerm> terms;
};

Run runMol(const Mol& m, const GffData& param, const GffGenerator& gen)
{
  Run r;
  Environment env;
  int n = static_cast<int>(m.atoms.size());
  std::vector<int> numbers(n);
  for (int i = 0; i < n; ++i)
    numbers[i] = m.atoms[i].element;
  if (!buildAngleList(n, m.neighbours, numbers, m.xyz, param.angl,
                      param.angl2, gen.fcThr, param.metal, r.angles, env))
    ADD_FAILURE() << env.errorMessage();
  if (!buildAngleTerms(r.angles, m.atoms, m.neighbours, m.rings, m.pbo,
                       param.group, param.metal, param.angl, param.angl2,
                       param, gen, r.terms, env))
    ADD_FAILURE() << env.errorMessage();
  return r;
}

} // namespace

TEST(AngleTest, bondAnglePbc)
{
  // (1,0,0)-(0,0,0)-(0,1,0): right angle; linear case gives pi.
  std::vector<double> xyz = { 1, 0, 0, 0, 0, 0, 0, 1, 0 };
  EXPECT_DOUBLE_EQ(bondAnglePbc(xyz, 0, 1, 2), 1.5707963267948966);
  xyz = { -1, 0, 0, 0, 0, 0, 1, 0, 0 };
  EXPECT_DOUBLE_EQ(bondAnglePbc(xyz, 0, 1, 2), 3.1415926535897931);
}

TEST(AngleTest, smallestRingAngle)
{
  // i in 6-ring {i,j,...}, j ringless, k with 5-ring {i,j,...} at slot 0:
  // loop 3 sees itest == 2 with sj[0] == 0 < 99, so rings3 = 5 verbatim.
  Ring ri, rk;
  ri.members = { 0, 1, 2, 3, 4, 5 };
  rk.members = { 0, 1, 6, 7, 8 };
  std::vector<Ring> vi = { ri }, vj, vk = { rk };
  EXPECT_EQ(smallestRingAngle(vi, vj, vk, 0, 1, 6), 0);
  std::vector<Ring> vj2 = { ri };
  EXPECT_EQ(smallestRingAngle(vi, vj2, vk, 0, 1, 6), 5);
  // Shared 3-ring on all ends.
  Ring r3;
  r3.members = { 0, 1, 6 };
  std::vector<Ring> v3 = { r3 };
  EXPECT_EQ(smallestRingAngle(v3, v3, v3, 0, 1, 6), 3);
}

TEST(AngleTest, water)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeWater(), param, gen);
  ASSERT_EQ(r.angles.size(), 1u);
  EXPECT_EQ(r.angles[0].center, 0);
  { // angle 0: (0,2,1)
    EXPECT_EQ(r.angles[0].center, 0);
    EXPECT_EQ(r.angles[0].first, 2);
    EXPECT_EQ(r.angles[0].second, 1);
    EXPECT_DOUBLE_EQ(r.terms[0].equilibrium, 1.74532925199432953e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 2.85657714295727283e-01);
  }
}

TEST(AngleTest, formaldehyde)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeFormaldehyde(), param, gen);
  ASSERT_EQ(r.angles.size(), 3u);
  { // angle 0: (0,2,1)
    EXPECT_EQ(r.angles[0].center, 0);
    EXPECT_EQ(r.angles[0].first, 2);
    EXPECT_EQ(r.angles[0].second, 1);
    EXPECT_DOUBLE_EQ(r.terms[0].equilibrium, 2.09439510239319526e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 1.52118514500338642e-01);
  }
  { // angle 1: (0,3,1)
    EXPECT_EQ(r.angles[1].center, 0);
    EXPECT_EQ(r.angles[1].first, 3);
    EXPECT_EQ(r.angles[1].second, 1);
    EXPECT_DOUBLE_EQ(r.terms[1].equilibrium, 2.09439510239319526e+00);
    EXPECT_DOUBLE_EQ(r.terms[1].forceConstant, 1.52118514500338642e-01);
  }
  { // angle 2: (0,3,2)
    EXPECT_EQ(r.angles[2].center, 0);
    EXPECT_EQ(r.angles[2].first, 3);
    EXPECT_EQ(r.angles[2].second, 2);
    EXPECT_DOUBLE_EQ(r.terms[2].equilibrium, 2.09439510239319526e+00);
    EXPECT_DOUBLE_EQ(r.terms[2].forceConstant, 1.66711857358269266e-01);
  }
}

TEST(AngleTest, methylcyclopropane)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeMethylcyclopropane(), param, gen);
  ASSERT_EQ(r.angles.size(), 5u);
  { // angle 0: (0,2,1)
    EXPECT_EQ(r.angles[0].center, 0);
    EXPECT_EQ(r.angles[0].first, 2);
    EXPECT_EQ(r.angles[0].second, 1);
    EXPECT_DOUBLE_EQ(r.terms[0].equilibrium, 1.43116998663535000e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 3.00013637909886999e-01);
  }
  { // angle 1: (0,3,1)
    EXPECT_EQ(r.angles[1].center, 0);
    EXPECT_EQ(r.angles[1].first, 3);
    EXPECT_EQ(r.angles[1].second, 1);
    EXPECT_DOUBLE_EQ(r.terms[1].equilibrium, 1.98094870101356380e+00);
    EXPECT_DOUBLE_EQ(r.terms[1].forceConstant, 2.56383156322592343e-01);
  }
  { // angle 2: (0,3,2)
    EXPECT_EQ(r.angles[2].center, 0);
    EXPECT_EQ(r.angles[2].first, 3);
    EXPECT_EQ(r.angles[2].second, 2);
    EXPECT_DOUBLE_EQ(r.terms[2].equilibrium, 1.98094870101356380e+00);
    EXPECT_DOUBLE_EQ(r.terms[2].forceConstant, 2.56383156322592343e-01);
  }
  { // angle 3: (1,2,0)
    EXPECT_EQ(r.angles[3].center, 1);
    EXPECT_EQ(r.angles[3].first, 2);
    EXPECT_EQ(r.angles[3].second, 0);
    EXPECT_DOUBLE_EQ(r.terms[3].equilibrium, 1.43116998663535000e+00);
    EXPECT_DOUBLE_EQ(r.terms[3].forceConstant, 1.66724446368596846e-01);
  }
  { // angle 4: (2,0,1)
    EXPECT_EQ(r.angles[4].center, 2);
    EXPECT_EQ(r.angles[4].first, 0);
    EXPECT_EQ(r.angles[4].second, 1);
    EXPECT_DOUBLE_EQ(r.terms[4].equilibrium, 1.43116998663535000e+00);
    EXPECT_DOUBLE_EQ(r.terms[4].forceConstant, 1.66724446368596846e-01);
  }
}

TEST(AngleTest, formamideAmide)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeFormamide(), param, gen);
  ASSERT_EQ(r.angles.size(), 6u);
  { // angle 0: (0,2,1)
    EXPECT_EQ(r.angles[0].center, 0);
    EXPECT_EQ(r.angles[0].first, 2);
    EXPECT_EQ(r.angles[0].second, 1);
    EXPECT_DOUBLE_EQ(r.terms[0].equilibrium, 2.09439510239319526e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 1.85592059616190741e-01);
  }
  { // angle 1: (0,5,1)
    EXPECT_EQ(r.angles[1].center, 0);
    EXPECT_EQ(r.angles[1].first, 5);
    EXPECT_EQ(r.angles[1].second, 1);
    EXPECT_DOUBLE_EQ(r.terms[1].equilibrium, 2.09439510239319526e+00);
    EXPECT_DOUBLE_EQ(r.terms[1].forceConstant, 1.52118514500338642e-01);
  }
  { // angle 2: (0,5,2)
    EXPECT_EQ(r.angles[2].center, 0);
    EXPECT_EQ(r.angles[2].first, 5);
    EXPECT_EQ(r.angles[2].second, 2);
    EXPECT_DOUBLE_EQ(r.terms[2].equilibrium, 2.09439510239319526e+00);
    EXPECT_DOUBLE_EQ(r.terms[2].forceConstant, 2.03396654714853325e-01);
  }
  { // angle 3: (2,3,0)
    EXPECT_EQ(r.angles[3].center, 2);
    EXPECT_EQ(r.angles[3].first, 3);
    EXPECT_EQ(r.angles[3].second, 0);
    EXPECT_DOUBLE_EQ(r.terms[3].equilibrium, 2.00712863979347889e+00);
    EXPECT_DOUBLE_EQ(r.terms[3].forceConstant, 3.43523404649629660e-01);
  }
  { // angle 4: (2,4,0)
    EXPECT_EQ(r.angles[4].center, 2);
    EXPECT_EQ(r.angles[4].first, 4);
    EXPECT_EQ(r.angles[4].second, 0);
    EXPECT_DOUBLE_EQ(r.terms[4].equilibrium, 2.00712863979347889e+00);
    EXPECT_DOUBLE_EQ(r.terms[4].forceConstant, 3.43523404649629660e-01);
  }
  { // angle 5: (2,4,3)
    EXPECT_EQ(r.angles[5].center, 2);
    EXPECT_EQ(r.angles[5].first, 4);
    EXPECT_EQ(r.angles[5].second, 3);
    EXPECT_DOUBLE_EQ(r.terms[5].equilibrium, 1.81514242207410281e+00);
    EXPECT_DOUBLE_EQ(r.terms[5].forceConstant, 1.97921945712523389e-01);
  }
}

TEST(AngleTest, anilineNonAmide)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeAnilineN(), param, gen);
  ASSERT_EQ(r.angles.size(), 4u);
  { // angle 0: (0,2,1)
    EXPECT_EQ(r.angles[0].center, 0);
    EXPECT_EQ(r.angles[0].first, 2);
    EXPECT_EQ(r.angles[0].second, 1);
    EXPECT_DOUBLE_EQ(r.terms[0].equilibrium, 2.09439510239319526e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 1.36872628493815263e-01);
  }
  { // angle 1: (2,3,0)
    EXPECT_EQ(r.angles[1].center, 2);
    EXPECT_EQ(r.angles[1].first, 3);
    EXPECT_EQ(r.angles[1].second, 0);
    EXPECT_DOUBLE_EQ(r.terms[1].equilibrium, 1.97222205475359225e+00);
    EXPECT_DOUBLE_EQ(r.terms[1].forceConstant, 2.29341765811006709e-01);
  }
  { // angle 2: (2,4,0)
    EXPECT_EQ(r.angles[2].center, 2);
    EXPECT_EQ(r.angles[2].first, 4);
    EXPECT_EQ(r.angles[2].second, 0);
    EXPECT_DOUBLE_EQ(r.terms[2].equilibrium, 1.97222205475359225e+00);
    EXPECT_DOUBLE_EQ(r.terms[2].forceConstant, 2.29341765811006709e-01);
  }
  { // angle 3: (2,4,3)
    EXPECT_EQ(r.angles[3].center, 2);
    EXPECT_EQ(r.angles[3].first, 4);
    EXPECT_EQ(r.angles[3].second, 3);
    EXPECT_DOUBLE_EQ(r.terms[3].equilibrium, 1.81514242207410281e+00);
    EXPECT_DOUBLE_EQ(r.terms[3].forceConstant, 1.97921945712523389e-01);
  }
}

TEST(AngleTest, sf6)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeSf6(), param, gen);
  ASSERT_EQ(r.angles.size(), 15u);
  { // angle 0: (0,2,1)
    EXPECT_EQ(r.angles[0].center, 0);
    EXPECT_EQ(r.angles[0].first, 2);
    EXPECT_EQ(r.angles[0].second, 1);
    EXPECT_DOUBLE_EQ(r.terms[0].equilibrium, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 1.51820557444507315e-01);
  }
  { // angle 1: (0,3,1)
    EXPECT_EQ(r.angles[1].center, 0);
    EXPECT_EQ(r.angles[1].first, 3);
    EXPECT_EQ(r.angles[1].second, 1);
    EXPECT_DOUBLE_EQ(r.terms[1].equilibrium, 1.57079632679489656e+00);
    EXPECT_DOUBLE_EQ(r.terms[1].forceConstant, 2.72342851555634591e-01);
  }
  { // angle 2: (0,3,2)
    EXPECT_EQ(r.angles[2].center, 0);
    EXPECT_EQ(r.angles[2].first, 3);
    EXPECT_EQ(r.angles[2].second, 2);
    EXPECT_DOUBLE_EQ(r.terms[2].equilibrium, 1.57079632679489656e+00);
    EXPECT_DOUBLE_EQ(r.terms[2].forceConstant, 2.72342851555634591e-01);
  }
  { // angle 3: (0,4,1)
    EXPECT_EQ(r.angles[3].center, 0);
    EXPECT_EQ(r.angles[3].first, 4);
    EXPECT_EQ(r.angles[3].second, 1);
    EXPECT_DOUBLE_EQ(r.terms[3].equilibrium, 1.57079632679489656e+00);
    EXPECT_DOUBLE_EQ(r.terms[3].forceConstant, 2.72342851555634591e-01);
  }
  { // angle 4: (0,4,2)
    EXPECT_EQ(r.angles[4].center, 0);
    EXPECT_EQ(r.angles[4].first, 4);
    EXPECT_EQ(r.angles[4].second, 2);
    EXPECT_DOUBLE_EQ(r.terms[4].equilibrium, 1.57079632679489656e+00);
    EXPECT_DOUBLE_EQ(r.terms[4].forceConstant, 2.72342851555634591e-01);
  }
  { // angle 5: (0,4,3)
    EXPECT_EQ(r.angles[5].center, 0);
    EXPECT_EQ(r.angles[5].first, 4);
    EXPECT_EQ(r.angles[5].second, 3);
    EXPECT_DOUBLE_EQ(r.terms[5].equilibrium, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[5].forceConstant, 1.51820557444507315e-01);
  }
  { // angle 6: (0,5,1)
    EXPECT_EQ(r.angles[6].center, 0);
    EXPECT_EQ(r.angles[6].first, 5);
    EXPECT_EQ(r.angles[6].second, 1);
    EXPECT_DOUBLE_EQ(r.terms[6].equilibrium, 1.57079632679489656e+00);
    EXPECT_DOUBLE_EQ(r.terms[6].forceConstant, 2.72342851555634591e-01);
  }
  { // angle 7: (0,5,2)
    EXPECT_EQ(r.angles[7].center, 0);
    EXPECT_EQ(r.angles[7].first, 5);
    EXPECT_EQ(r.angles[7].second, 2);
    EXPECT_DOUBLE_EQ(r.terms[7].equilibrium, 1.57079632679489656e+00);
    EXPECT_DOUBLE_EQ(r.terms[7].forceConstant, 2.72342851555634591e-01);
  }
  { // angle 8: (0,5,3)
    EXPECT_EQ(r.angles[8].center, 0);
    EXPECT_EQ(r.angles[8].first, 5);
    EXPECT_EQ(r.angles[8].second, 3);
    EXPECT_DOUBLE_EQ(r.terms[8].equilibrium, 1.57079632679489656e+00);
    EXPECT_DOUBLE_EQ(r.terms[8].forceConstant, 2.72342851555634591e-01);
  }
  { // angle 9: (0,5,4)
    EXPECT_EQ(r.angles[9].center, 0);
    EXPECT_EQ(r.angles[9].first, 5);
    EXPECT_EQ(r.angles[9].second, 4);
    EXPECT_DOUBLE_EQ(r.terms[9].equilibrium, 1.57079632679489656e+00);
    EXPECT_DOUBLE_EQ(r.terms[9].forceConstant, 2.72342851555634591e-01);
  }
  { // angle 10: (0,6,1)
    EXPECT_EQ(r.angles[10].center, 0);
    EXPECT_EQ(r.angles[10].first, 6);
    EXPECT_EQ(r.angles[10].second, 1);
    EXPECT_DOUBLE_EQ(r.terms[10].equilibrium, 1.57079632679489656e+00);
    EXPECT_DOUBLE_EQ(r.terms[10].forceConstant, 2.72342851555634591e-01);
  }
  { // angle 11: (0,6,2)
    EXPECT_EQ(r.angles[11].center, 0);
    EXPECT_EQ(r.angles[11].first, 6);
    EXPECT_EQ(r.angles[11].second, 2);
    EXPECT_DOUBLE_EQ(r.terms[11].equilibrium, 1.57079632679489656e+00);
    EXPECT_DOUBLE_EQ(r.terms[11].forceConstant, 2.72342851555634591e-01);
  }
  { // angle 12: (0,6,3)
    EXPECT_EQ(r.angles[12].center, 0);
    EXPECT_EQ(r.angles[12].first, 6);
    EXPECT_EQ(r.angles[12].second, 3);
    EXPECT_DOUBLE_EQ(r.terms[12].equilibrium, 1.57079632679489656e+00);
    EXPECT_DOUBLE_EQ(r.terms[12].forceConstant, 2.72342851555634591e-01);
  }
  { // angle 13: (0,6,4)
    EXPECT_EQ(r.angles[13].center, 0);
    EXPECT_EQ(r.angles[13].first, 6);
    EXPECT_EQ(r.angles[13].second, 4);
    EXPECT_DOUBLE_EQ(r.terms[13].equilibrium, 1.57079632679489656e+00);
    EXPECT_DOUBLE_EQ(r.terms[13].forceConstant, 2.72342851555634591e-01);
  }
  { // angle 14: (0,6,5)
    EXPECT_EQ(r.angles[14].center, 0);
    EXPECT_EQ(r.angles[14].first, 6);
    EXPECT_EQ(r.angles[14].second, 5);
    EXPECT_DOUBLE_EQ(r.terms[14].equilibrium, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[14].forceConstant, 1.51820557444507315e-01);
  }
}

TEST(AngleTest, feCo)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeFeCo(), param, gen);
  ASSERT_EQ(r.angles.size(), 1u);
  { // angle 0: (1,2,0)
    EXPECT_EQ(r.angles[0].center, 1);
    EXPECT_EQ(r.angles[0].first, 2);
    EXPECT_EQ(r.angles[0].second, 0);
    EXPECT_DOUBLE_EQ(r.terms[0].equilibrium, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 2.23324684824972253e-01);
  }
}

TEST(AngleTest, momLinear)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeMom(), param, gen);
  ASSERT_EQ(r.angles.size(), 1u);
  { // angle 0: (1,2,0)
    EXPECT_EQ(r.angles[0].center, 1);
    EXPECT_EQ(r.angles[0].first, 2);
    EXPECT_EQ(r.angles[0].second, 0);
    EXPECT_DOUBLE_EQ(r.terms[0].equilibrium, 3.14159265358979312e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 1.71242307529582863e-03);
  }
}

TEST(AngleTest, feCpFeta)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeFeCp(), param, gen);
  ASSERT_EQ(r.angles.size(), 1u);
  { // angle 0: (0,2,1)
    EXPECT_EQ(r.angles[0].center, 0);
    EXPECT_EQ(r.angles[0].first, 2);
    EXPECT_EQ(r.angles[0].second, 1);
    EXPECT_DOUBLE_EQ(r.terms[0].equilibrium, 1.57079632679489656e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 5.58165811230479630e-03);
  }
}

TEST(AngleTest, caEtaSkip)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeCaSkip(), param, gen);
  EXPECT_EQ(r.angles.size(), 0u);
}

} // namespace Xtb
} // namespace Avogadro
