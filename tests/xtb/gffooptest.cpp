/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "avogadro/xtb/COPYING.LESSER".
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/xtb/environment.h>
#include <avogadro/xtb/gfnffdata.h>
#include <avogadro/xtb/gfnffoop.h>
#include <avogadro/xtb/gfnffsetup.h>
#include <avogadro/xtb/gfnfftopo.h>

#include <cmath>
#include <vector>

namespace Avogadro {
namespace Xtb {

// Reference values verified bitwise against the unmodified out-of-plane
// block in xtb src/gfnff/gfnff_ini.f90 (2061-2171), compiled with gfortran
// -ffp-contract=off over 12 synthetic fragments (15 impropers).

namespace {

struct Mol
{
  std::vector<OopAtom> atoms;
  std::vector<std::vector<int>> neighbours;
  std::vector<double> xyz;
  std::vector<double> pbo;
};

int addAtom(Mol& m, int z, double qa, int pi, double x, double y, double zc)
{
  OopAtom a;
  a.element = z;
  a.charge = qa;
  a.pi = pi;
  m.atoms.push_back(a);
  m.xyz.push_back(x);
  m.xyz.push_back(y);
  m.xyz.push_back(zc);
  m.neighbours.emplace_back();
  return static_cast<int>(m.atoms.size()) - 1;
}

struct BondP
{
  int a, b;
  double p;
};

struct MolEx
{
  Mol m;
  std::vector<BondP> bonds;
};

void linkEx(MolEx& e, int a, int b, double p = 0.0)
{
  e.m.neighbours[a].push_back(b);
  e.m.neighbours[b].push_back(a);
  e.bonds.push_back({ a, b, p });
}

void finishEx(MolEx& e)
{
  int n = static_cast<int>(e.m.atoms.size());
  e.m.pbo.assign(n * (n + 1) / 2, 0.0);
  for (const BondP& b : e.bonds)
    e.m.pbo[packedIndex(b.a, b.b)] = b.p;
}

MolEx makeBenzene()
{
  MolEx e;
  std::vector<int> ids;
  for (int k = 0; k < 6; ++k) {
    double t = k * 3.141592653589793 / 3.0;
    ids.push_back(addAtom(e.m, 6, 0.0, 1, 1.40 * std::cos(t),
                          1.40 * std::sin(t), 0));
  }
  for (int k = 0; k < 6; ++k) {
    linkEx(e, ids[k], ids[(k + 1) % 6], 0.65);
    double t = k * 3.141592653589793 / 3.0;
    int h = addAtom(e.m, 1, 0.0, 0, 2.48 * std::cos(t), 2.48 * std::sin(t),
                    0);
    linkEx(e, ids[k], h, 0.02);
  }
  finishEx(e);
  return e;
}

MolEx makeFormaldehyde()
{
  MolEx e;
  int c = addAtom(e.m, 6, 0.0, 1, 0, 0, 0);
  int o = addAtom(e.m, 8, 0.0, 1, 1.21, 0, 0);
  int h1 = addAtom(e.m, 1, 0.0, 0, -0.48, 0.88, 0);
  int h2 = addAtom(e.m, 1, 0.0, 0, -0.48, -0.88, 0);
  linkEx(e, c, o, 0.8);
  linkEx(e, c, h1, 0.05);
  linkEx(e, c, h2, 0.05);
  finishEx(e);
  return e;
}

MolEx makeAmmonia()
{
  MolEx e;
  int n = addAtom(e.m, 7, 0.0, 0, 0, 0, 0);
  for (int k = 0; k < 3; ++k) {
    double t = k * 2.0 * 3.141592653589793 / 3.0;
    int h = addAtom(e.m, 1, 0.0, 0, 1.01 * std::cos(t) * 0.83,
                    1.01 * std::sin(t) * 0.83, -0.38);
    linkEx(e, n, h);
  }
  finishEx(e);
  return e;
}

MolEx makeMethylCation()
{
  MolEx e;
  int c = addAtom(e.m, 6, 1.0, 1, 0, 0, 0);
  for (int k = 0; k < 3; ++k) {
    double t = k * 2.0 * 3.141592653589793 / 3.0;
    int h = addAtom(e.m, 1, 0.0, 0, 1.09 * std::cos(t),
                    1.09 * std::sin(t), 0);
    linkEx(e, c, h, 0.03);
  }
  finishEx(e);
  return e;
}

MolEx makeNOxide()
{
  MolEx e;
  int n = addAtom(e.m, 7, 0.0, 1, 0, 0, 0);
  int o = addAtom(e.m, 8, 0.0, 1, 1.25, 0, 0);
  int c1 = addAtom(e.m, 6, 0.0, 0, -1.2, 0.7, 0);
  int c2 = addAtom(e.m, 6, 0.0, 0, -1.2, -0.7, 0);
  linkEx(e, n, o, 0.5);
  linkEx(e, n, c1, 0.1);
  linkEx(e, n, c2, 0.1);
  finishEx(e);
  return e;
}

MolEx makeIsobutane()
{
  MolEx e;
  int c0 = addAtom(e.m, 6, 0.0, 0, 0, 0, 0);
  for (int k = 0; k < 3; ++k) {
    double t = k * 2.0 * 3.141592653589793 / 3.0;
    int c = addAtom(e.m, 6, 0.0, 0, 1.54 * std::cos(t),
                    1.54 * std::sin(t), 0.3);
    linkEx(e, c0, c);
  }
  finishEx(e);
  return e;
}

struct Run
{
  std::vector<Oop> oops;
  std::vector<OopTerm> terms;
};

Run runMol(const MolEx& e, const GffData& param, const GffGenerator& gen)
{
  Run r;
  Environment env;
  if (!buildOutOfPlane(e.m.atoms, e.m.neighbours, e.m.xyz, e.m.pbo,
                       param.group, param.repz, param, gen, r.oops, r.terms,
                       env))
    ADD_FAILURE() << env.errorMessage();
  return r;
}

} // namespace

TEST(OopTest, sortTripleByDistance)
{
  // Distinct values sort ascending with carried ids.
  std::vector<double> d = { 2.0, 0.5, 1.0 };
  std::vector<int> ids = { 10, 20, 30 };
  sortTripleByDistance(d, ids);
  EXPECT_DOUBLE_EQ(d[0], 0.5);
  EXPECT_DOUBLE_EQ(d[1], 1.0);
  EXPECT_DOUBLE_EQ(d[2], 2.0);
  EXPECT_EQ(ids[0], 20);
  EXPECT_EQ(ids[1], 30);
  EXPECT_EQ(ids[2], 10);
  // All-equal tie: last-minimum selection gives [ll,jj,kk] order.
  d = { 1.0, 1.0, 1.0 };
  ids = { 10, 20, 30 };
  sortTripleByDistance(d, ids);
  EXPECT_EQ(ids[0], 30);
  EXPECT_EQ(ids[1], 10);
  EXPECT_EQ(ids[2], 20);
}

TEST(OopTest, benzene)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeBenzene(), param, gen);
  ASSERT_EQ(r.oops.size(), 6u);
  { // oop 0: (0,6,5,1) kind 0
    EXPECT_EQ(r.oops[0].center, 0);
    EXPECT_EQ(r.oops[0].first, 6);
    EXPECT_EQ(r.oops[0].second, 5);
    EXPECT_EQ(r.oops[0].third, 1);
    EXPECT_EQ(r.oops[0].kind, 0);
    EXPECT_DOUBLE_EQ(r.terms[0].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 3.56999983787536612e-01);
  }
  { // oop 1: (1,7,0,2) kind 0
    EXPECT_EQ(r.oops[1].center, 1);
    EXPECT_EQ(r.oops[1].first, 7);
    EXPECT_EQ(r.oops[1].second, 0);
    EXPECT_EQ(r.oops[1].third, 2);
    EXPECT_EQ(r.oops[1].kind, 0);
    EXPECT_DOUBLE_EQ(r.terms[1].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[1].forceConstant, 3.56999983787536612e-01);
  }
  { // oop 2: (2,8,1,3) kind 0
    EXPECT_EQ(r.oops[2].center, 2);
    EXPECT_EQ(r.oops[2].first, 8);
    EXPECT_EQ(r.oops[2].second, 1);
    EXPECT_EQ(r.oops[2].third, 3);
    EXPECT_EQ(r.oops[2].kind, 0);
    EXPECT_DOUBLE_EQ(r.terms[2].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[2].forceConstant, 3.56999983787536612e-01);
  }
  { // oop 3: (3,9,4,2) kind 0
    EXPECT_EQ(r.oops[3].center, 3);
    EXPECT_EQ(r.oops[3].first, 9);
    EXPECT_EQ(r.oops[3].second, 4);
    EXPECT_EQ(r.oops[3].third, 2);
    EXPECT_EQ(r.oops[3].kind, 0);
    EXPECT_DOUBLE_EQ(r.terms[3].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[3].forceConstant, 3.56999983787536612e-01);
  }
  { // oop 4: (4,10,3,5) kind 0
    EXPECT_EQ(r.oops[4].center, 4);
    EXPECT_EQ(r.oops[4].first, 10);
    EXPECT_EQ(r.oops[4].second, 3);
    EXPECT_EQ(r.oops[4].third, 5);
    EXPECT_EQ(r.oops[4].kind, 0);
    EXPECT_DOUBLE_EQ(r.terms[4].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[4].forceConstant, 3.56999983787536612e-01);
  }
  { // oop 5: (5,11,0,4) kind 0
    EXPECT_EQ(r.oops[5].center, 5);
    EXPECT_EQ(r.oops[5].first, 11);
    EXPECT_EQ(r.oops[5].second, 0);
    EXPECT_EQ(r.oops[5].third, 4);
    EXPECT_EQ(r.oops[5].kind, 0);
    EXPECT_DOUBLE_EQ(r.terms[5].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[5].forceConstant, 3.56999983787536612e-01);
  }
}

TEST(OopTest, formaldehyde)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeFormaldehyde(), param, gen);
  ASSERT_EQ(r.oops.size(), 1u);
  { // oop 0: (0,3,2,1) kind 0
    EXPECT_EQ(r.oops[0].center, 0);
    EXPECT_EQ(r.oops[0].first, 3);
    EXPECT_EQ(r.oops[0].second, 2);
    EXPECT_EQ(r.oops[0].third, 1);
    EXPECT_EQ(r.oops[0].kind, 0);
    EXPECT_DOUBLE_EQ(r.terms[0].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 2.19449990034103415e+01);
  }
}

TEST(OopTest, ammonia)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeAmmonia(), param, gen);
  ASSERT_EQ(r.oops.size(), 1u);
  { // oop 0: (0,3,1,2) kind -1
    EXPECT_EQ(r.oops[0].center, 0);
    EXPECT_EQ(r.oops[0].first, 3);
    EXPECT_EQ(r.oops[0].second, 1);
    EXPECT_EQ(r.oops[0].third, 2);
    EXPECT_EQ(r.oops[0].kind, -1);
    EXPECT_DOUBLE_EQ(r.terms[0].phase, 1.39626340159546358e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 1.79999999999999982e+00);
  }
}

TEST(OopTest, methylCation)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeMethylCation(), param, gen);
  ASSERT_EQ(r.oops.size(), 1u);
  { // oop 0: (0,3,1,2) kind 0
    EXPECT_EQ(r.oops[0].center, 0);
    EXPECT_EQ(r.oops[0].first, 3);
    EXPECT_EQ(r.oops[0].second, 1);
    EXPECT_EQ(r.oops[0].third, 2);
    EXPECT_EQ(r.oops[0].kind, 0);
    EXPECT_DOUBLE_EQ(r.terms[0].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 6.01649972677230771e+00);
  }
}

TEST(OopTest, nOxide)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeNOxide(), param, gen);
  ASSERT_EQ(r.oops.size(), 1u);
  { // oop 0: (0,1,3,2) kind 0
    EXPECT_EQ(r.oops[0].center, 0);
    EXPECT_EQ(r.oops[0].first, 1);
    EXPECT_EQ(r.oops[0].second, 3);
    EXPECT_EQ(r.oops[0].third, 2);
    EXPECT_EQ(r.oops[0].kind, 0);
    EXPECT_DOUBLE_EQ(r.terms[0].phase, 0.00000000000000000e+00);
    EXPECT_DOUBLE_EQ(r.terms[0].forceConstant, 1.04999995231628418e+01);
  }
}

TEST(OopTest, isobutaneSkipped)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Run r = runMol(makeIsobutane(), param, gen);
  EXPECT_EQ(r.oops.size(), 0u); // pi == 0 and not N
}

} // namespace Xtb
} // namespace Avogadro
