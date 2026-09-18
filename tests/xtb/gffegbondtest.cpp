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
#include <avogadro/xtb/gfnffegbond.h>
#include <avogadro/xtb/gfnffoop.h>
#include <avogadro/xtb/gfnffsetup.h>
#include <avogadro/xtb/gfnfftopo.h>
#include <avogadro/xtb/gfnffvbond.h>

#include <cmath>
#include <vector>

namespace Avogadro {
namespace Xtb {

// Reference values verified bitwise against the unmodified egbond/egbend/
// egtors in xtb src/gfnff/gfnff_eg.f90, compiled with gfortran
// -ffp-contract=off (17 energy/gradient/stress outputs over ethane,
// water, SF6, ethylene and NH3 inputs).

namespace {

struct Mol
{
  std::vector<int> numbers;
  std::vector<double> xyz;
  std::vector<std::vector<int>> neighbours;
};

int addAtom(Mol& m, int z, double x, double y, double zc)
{
  m.numbers.push_back(z);
  m.xyz.push_back(x);
  m.xyz.push_back(y);
  m.xyz.push_back(zc);
  m.neighbours.emplace_back();
  return static_cast<int>(m.numbers.size()) - 1;
}

void link(Mol& m, int a, int b)
{
  m.neighbours[a].push_back(b);
  m.neighbours[b].push_back(a);
}






Mol makeEthane()
{
  Mol m;
  const double xyz[8][3] = {
    { -0.77, 0, 0 }, { 0.77, 0, 0 }, { -1.13, 1.03, 0 },
    { -1.13, -0.515, 0.892 }, { -1.13, -0.515, -0.892 },
    { 1.13, 1.03, 0 }, { 1.13, -0.515, 0.892 }, { 1.13, -0.515, -0.892 },
  };
  const int z[8] = { 6, 6, 1, 1, 1, 1, 1, 1 };
  for (int i = 0; i < 8; ++i)
    addAtom(m, z[i], xyz[i][0], xyz[i][1], xyz[i][2]);
  link(m, 0, 1);
  link(m, 0, 2);
  link(m, 0, 3);
  link(m, 0, 4);
  link(m, 1, 5);
  link(m, 1, 6);
  link(m, 1, 7);
  return m;
}

Mol makeWater()
{
  Mol m;
  int o = addAtom(m, 8, 0, 0, 0);
  int h1 = addAtom(m, 1, 0.96, 0, 0);
  int h2 = addAtom(m, 1, -0.2485, 0.9269, 0);
  link(m, o, h1);
  link(m, o, h2);
  return m;
}

Mol makeSf6Trans()
{
  Mol m;
  int s = addAtom(m, 16, 0, 0, 0);
  int f1 = addAtom(m, 9, 1.56, 0, 0);
  int f2 = addAtom(m, 9, -1.56, 0, 0);
  link(m, s, f1);
  link(m, s, f2);
  return m;
}

Mol makeEthylene()
{
  Mol m;
  int c1 = addAtom(m, 6, -0.67, 0, 0);
  int c2 = addAtom(m, 6, 0.67, 0, 0);
  link(m, c1, c2);
  int h1 = addAtom(m, 1, -1.22, 0.95, 0);
  int h2 = addAtom(m, 1, -1.22, -0.95, 0);
  int h3 = addAtom(m, 1, 1.22, 0.95, 0);
  int h4 = addAtom(m, 1, 1.22, -0.95, 0);
  link(m, c1, h1);
  link(m, c1, h2);
  link(m, c2, h3);
  link(m, c2, h4);
  return m;
}

Mol makeAmmonia()
{
  Mol m;
  int n = addAtom(m, 7, 0, 0, 0);
  int h1 = addAtom(m, 1, 0.84, 0.5, -0.38);
  int h2 = addAtom(m, 1, -0.42, 0.85, -0.38);
  int h3 = addAtom(m, 1, -0.42, -0.35, 0.76);
  link(m, n, h1);
  link(m, n, h2);
  link(m, n, h3);
  return m;
}

Mol makeFormaldehydeC()
{
  Mol m;
  int c = addAtom(m, 6, 0, 0, 0);
  int o = addAtom(m, 8, 1.21, 0, 0);
  int h1 = addAtom(m, 1, -0.48, 0.88, 0);
  int h2 = addAtom(m, 1, -0.48, -0.88, 0);
  link(m, c, o);
  link(m, c, h1);
  link(m, c, h2);
  return m;
}

} // namespace

TEST(EgBondTest, vecHelpers)
{
  double a[3] = { 1.0, 0.0, 0.0 }, b[3] = { 0.0, 1.0, 0.0 }, c[3];
  vecSub(a, b, c);
  EXPECT_DOUBLE_EQ(c[0], 1.0);
  EXPECT_DOUBLE_EQ(c[1], -1.0);
  EXPECT_DOUBLE_EQ(c[2], 0.0);
  vecCross(a, b, c);
  EXPECT_DOUBLE_EQ(c[0], 0.0);
  EXPECT_DOUBLE_EQ(c[1], 0.0);
  EXPECT_DOUBLE_EQ(c[2], 1.0);
  EXPECT_DOUBLE_EQ(vecLen(a), 1.0);
  EXPECT_DOUBLE_EQ(vecCos(a, a), 1.0);
  EXPECT_DOUBLE_EQ(vecCos(a, b), 0.0);
  double r[3] = { 3.0, 4.0, 0.0 };
  EXPECT_DOUBLE_EQ(vecNorm(r, 0), 5.0);
  EXPECT_DOUBLE_EQ(vecNorm(r, 1), 5.0);
  EXPECT_DOUBLE_EQ(r[0], 6.00000000000000089e-01);
  EXPECT_DOUBLE_EQ(r[1], 0.8);
}

TEST(EgBondTest, bondDamp)
{
  std::vector<double> rcov(103, 0.0);
  rcov[0] = 0.8062830717001458;
  rcov[7] = 1.5873697974096619;
  double damp = 0.0, ddamp = 0.0;
  bondDamp(0.9216, 0.595, rcov, 1, 8, damp, ddamp);
  EXPECT_DOUBLE_EQ(damp, 9.31895831454505941e-01);
  EXPECT_DOUBLE_EQ(ddamp, -2.75460029393320704e-01);
}

TEST(EgBondTest, ethaneBonds)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Mol m = makeEthane();
  int n = static_cast<int>(m.numbers.size());
  // vbond via the verified builder (C-C then C-H)
  std::vector<VbondAtom> va(n);
  for (int i = 0; i < n; ++i) {
    va[i].element = m.numbers[i];
    va[i].hyb = m.numbers[i] == 6 ? 3 : 0;
  }
  std::vector<VbondBond> vb;
  for (const auto& pr : std::vector<std::pair<int, int>>{ { 0, 1 },
                                                          { 0, 2 } }) {
    VbondBond b;
    b.first = pr.first;
    b.second = pr.second;
    b.guess = 2.9;
    b.coordFirst = static_cast<int>(m.neighbours[pr.first].size());
    b.coordSecond = static_cast<int>(m.neighbours[pr.second].size());
    vb.push_back(b);
  }
  std::vector<int> row6(n, 2);
  row6[0] = row6[1] = 2;
  std::vector<VbondTerm> vt;
  std::vector<int> vbt;
  Environment env;
  ASSERT_TRUE(buildVbondTerms(vb, va, m.neighbours, m.numbers, param.group,
                              param.metal, param.en, param.bond, row6, param,
                              gen, vt, vbt, env))
    << env.errorMessage();
  double ebond = 0.0;
  std::vector<double> gbond(3 * n, 0.0), dEdcn(n, 0.0);
  double sigma[3][3] = { { 0 } };
  // synthetic rij/drij (seeded): recomputed here identically to the dump
  unsigned long long seed = 12345;
  auto rnd = [&seed]() {
    seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(seed >> 11) * (1.0 / 9007199254740992.0) - 0.5;
  };
  for (int bi = 0; bi < 2; ++bi) {
    int ia = vb[bi].first, ja = vb[bi].second;
    double dx = m.xyz[3 * ia] - m.xyz[3 * ja];
    double dy = m.xyz[3 * ia + 1] - m.xyz[3 * ja + 1];
    double dz = m.xyz[3 * ia + 2] - m.xyz[3 * ja + 2];
    double rab = std::sqrt(dx * dx + dy * dy + dz * dz);
    std::vector<double> drij(3 * n);
    for (int k = 0; k < 3 * n; ++k)
      drij[k] = rnd() * 0.01;
    double dcn[2] = { rnd() * 0.01, rnd() * 0.01 };
    std::vector<double> vsh(2, 0.0), vst(2), vpr(2);
    vst[bi] = vt[bi].steepness;
    vpr[bi] = vt[bi].prefactor;
    bondEnergyGradient(bi, ia, ja, rab, rab * 0.98, drij, dcn, n, m.xyz,
                       vsh, vst, vpr, ebond, gbond, dEdcn, sigma, env);
  }
  EXPECT_DOUBLE_EQ(ebond, -3.16615326816609199e-01);
  EXPECT_DOUBLE_EQ(gbond[0], -3.22615395274875255e-03);
  EXPECT_DOUBLE_EQ(gbond[1], -3.25548324468713426e-03);
  EXPECT_DOUBLE_EQ(gbond[2], -9.69905044684064383e-06);
  EXPECT_DOUBLE_EQ(gbond[3], 4.36225796169284908e-03);
  EXPECT_DOUBLE_EQ(gbond[4], -7.97723018998806792e-06);
  EXPECT_DOUBLE_EQ(gbond[5], 1.35490324915156817e-05);
  EXPECT_DOUBLE_EQ(gbond[6], -1.14159827351677733e-03);
  EXPECT_DOUBLE_EQ(gbond[7], 3.26490743031206641e-03);
  EXPECT_DOUBLE_EQ(gbond[8], -2.29187925472272615e-05);
  EXPECT_DOUBLE_EQ(gbond[9], 2.11454864270007245e-05);
  EXPECT_DOUBLE_EQ(gbond[10], -2.09451802433224886e-05);
  EXPECT_DOUBLE_EQ(gbond[11], -1.12019466239608738e-05);
  EXPECT_DOUBLE_EQ(gbond[12], 2.57515713206376780e-05);
  EXPECT_DOUBLE_EQ(gbond[13], 6.80024796337323252e-07);
  EXPECT_DOUBLE_EQ(gbond[14], -1.08730021810582675e-05);
  EXPECT_DOUBLE_EQ(gbond[15], 1.45583816551213301e-05);
  EXPECT_DOUBLE_EQ(gbond[16], 4.61202919377626713e-06);
  EXPECT_DOUBLE_EQ(gbond[17], 2.09016488272843990e-05);
  EXPECT_DOUBLE_EQ(gbond[18], -5.16869439360704442e-06);
  EXPECT_DOUBLE_EQ(gbond[19], 1.21292188745447704e-05);
  EXPECT_DOUBLE_EQ(gbond[20], 8.38130714249828702e-06);
  EXPECT_DOUBLE_EQ(gbond[21], 3.49125790563342480e-05);
  EXPECT_DOUBLE_EQ(gbond[22], 1.15490116957362737e-05);
  EXPECT_DOUBLE_EQ(gbond[23], -9.73381109047320351e-06);
  EXPECT_DOUBLE_EQ(dEdcn[0], -2.13821118489582734e-05);
  EXPECT_DOUBLE_EQ(dEdcn[1], -5.93075014655581776e-07);
  EXPECT_DOUBLE_EQ(dEdcn[2], -8.96031434738523657e-06);
  EXPECT_DOUBLE_EQ(dEdcn[3], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(dEdcn[4], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(dEdcn[5], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(dEdcn[6], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(dEdcn[7], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(sigma[0][0], 7.08873796106385817e-03);
  EXPECT_DOUBLE_EQ(sigma[0][1], -1.17165559846630751e-03);
  EXPECT_DOUBLE_EQ(sigma[0][2], 4.77920866967083307e-05);
  EXPECT_DOUBLE_EQ(sigma[1][0], -1.21706118808986657e-03);
  EXPECT_DOUBLE_EQ(sigma[1][1], 3.36068117055887293e-03);
  EXPECT_DOUBLE_EQ(sigma[1][2], -7.88024836394766952e-06);
  EXPECT_DOUBLE_EQ(sigma[2][0], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(sigma[2][1], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(sigma[2][2], 0.00000000000000000e+00);
}

TEST(EgBondTest, waterBend)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Mol m = makeWater();
  std::vector<Angle> angles;
  Environment env;
  ASSERT_TRUE(buildAngleList(3, m.neighbours, m.numbers, m.xyz, param.angl,
                             param.angl2, gen.fcThr, param.metal, angles,
                             env));
  ASSERT_EQ(angles.size(), 1u);
  std::vector<AngleAtom> aa(3);
  aa[0].element = 8;
  aa[0].hyb = 3;
  aa[0].charge = -0.6;
  aa[1].element = aa[2].element = 1;
  aa[1].charge = aa[2].charge = 0.3;
  std::vector<std::vector<Ring>> rings(3);
  std::vector<AngleTerm> at;
  ASSERT_TRUE(buildAngleTerms(angles, aa, m.neighbours, rings,
                              std::vector<double>(), param.group,
                              param.metal, param.angl, param.angl2, param,
                              gen, at, env));
  double e = 0.0, ag[3][3] = { { 0 } }, ds[3][3] = { { 0 } };
  const Angle& a = angles[0];
  bendEnergyGradient(a.center, a.first, a.second, at[0].equilibrium,
                     at[0].forceConstant, m.numbers, m.xyz, param.rcov,
                     param.angleCutA, e, ag, ds);
  EXPECT_DOUBLE_EQ(e, 1.80540304959001864e-03);
  EXPECT_DOUBLE_EQ(ag[0][0], 4.15307503114280205e-02);
  EXPECT_DOUBLE_EQ(ag[0][1], 5.41147310134281259e-02);
  EXPECT_DOUBLE_EQ(ag[0][2], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[1][0], -4.10184358381704070e-02);
  EXPECT_DOUBLE_EQ(ag[1][1], -1.15268130695000900e-02);
  EXPECT_DOUBLE_EQ(ag[1][2], -0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[2][0], -5.12314473257615755e-04);
  EXPECT_DOUBLE_EQ(ag[2][1], -4.25879179439280359e-02);
  EXPECT_DOUBLE_EQ(ag[2][2], -0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ds[0][0], 9.70125941145803458e-03);
  EXPECT_DOUBLE_EQ(ds[0][1], -3.80199881784001464e-02);
  EXPECT_DOUBLE_EQ(ds[0][2], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ds[1][0], -3.80199881784001464e-02);
  EXPECT_DOUBLE_EQ(ds[1][1], -1.06842030341196322e-02);
  EXPECT_DOUBLE_EQ(ds[1][2], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ds[2][0], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ds[2][1], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ds[2][2], 0.00000000000000000e+00);
}

TEST(EgBondTest, sf6LinearBend)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  std::vector<int> numbers = { 16, 9, 9 };
  std::vector<double> xyz = { 0, 0, 0, 1.56, 0, 0, -1.56, 0, 0 };
  double e = 0.0, ag[3][3] = { { 0 } }, ds[3][3] = { { 0 } };
  bendEnergyGradient(0, 1, 2, 3.141592653589793, 0.27234285155563459,
                     numbers, xyz, param.rcov, param.angleCutA, e, ag, ds);
  EXPECT_DOUBLE_EQ(e, 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[0][0], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[0][1], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[0][2], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[1][0], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[1][1], -0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[1][2], -0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[2][0], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[2][1], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[2][2], -0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ds[0][0], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ds[0][1], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ds[0][2], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ds[1][0], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ds[1][1], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ds[1][2], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ds[2][0], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ds[2][1], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ds[2][2], 0.00000000000000000e+00);
}

TEST(EgBondTest, ethyleneTorsion)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Mol m = makeEthylene();
  double e = 0.0, ag[4][3] = { { 0 } }, ds[3][3] = { { 0 } };
  // trans H-C-C-H (atoms 2,0,1,4 in builder order)
  torsionEnergyGradient(2, 0, 1, 4, 2, 3.141592653589793,
                        1.69590891266774912, m.numbers, m.xyz, param.rcov,
                        param.angleCutT, false, e, ag, ds);
  EXPECT_DOUBLE_EQ(e, 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[0][0], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[0][1], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[0][2], -0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[1][0], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[1][1], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[1][2], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[2][0], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[2][1], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[2][2], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[3][0], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[3][1], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[3][2], 0.00000000000000000e+00);
}

TEST(EgBondTest, crossTorsion)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Mol m = makeEthane();
  double e = 0.0, ag[4][3] = { { 0 } }, ds[3][3] = { { 0 } };
  // H-C-C-H (atoms 2,0,1,5 in builder order), sp3 main only
  torsionEnergyGradient(2, 0, 1, 5, 3, 3.141592653589793,
                        0.21510453104338512, m.numbers, m.xyz, param.rcov,
                        param.angleCutT, false, e, ag, ds);
  EXPECT_DOUBLE_EQ(e, 3.17853942992015415e-01);
  EXPECT_DOUBLE_EQ(ag[0][0], 3.65942527896341158e-02);
  EXPECT_DOUBLE_EQ(ag[0][1], -1.04700223259230979e-01);
  EXPECT_DOUBLE_EQ(ag[0][2], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[1][0], 4.39392783391122782e-02);
  EXPECT_DOUBLE_EQ(ag[1][1], 1.04700223259230979e-01);
  EXPECT_DOUBLE_EQ(ag[1][2], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[2][0], -4.39392783391122782e-02);
  EXPECT_DOUBLE_EQ(ag[2][1], 1.04700223259230993e-01);
  EXPECT_DOUBLE_EQ(ag[2][2], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[3][0], -3.65942527896341158e-02);
  EXPECT_DOUBLE_EQ(ag[3][1], -1.04700223259230993e-01);
  EXPECT_DOUBLE_EQ(ag[3][2], 0.00000000000000000e+00);
}

TEST(EgBondTest, ammoniaOop)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Mol m = makeAmmonia();
  std::vector<OopAtom> oa(4);
  oa[0].element = 7;
  oa[1].element = oa[2].element = oa[3].element = 1;
  std::vector<Oop> oops;
  std::vector<OopTerm> ot;
  Environment env;
  std::vector<double> pbo(10, 0.0);
  ASSERT_TRUE(buildOutOfPlane(oa, m.neighbours, m.xyz, pbo, param.group,
                              param.repz, param, gen, oops, ot, env));
  ASSERT_EQ(oops.size(), 1u);
  double e = 0.0, ag[4][3] = { { 0 } }, ds[3][3] = { { 0 } };
  const Oop& o = oops[0];
  torsionEnergyGradient(o.center, o.first, o.second, o.third, o.kind,
                        ot[0].phase, ot[0].forceConstant, m.numbers, m.xyz,
                        param.rcov, param.angleCutT, false, e, ag, ds);
  EXPECT_DOUBLE_EQ(e, 1.68109890230921030e-02);
  EXPECT_DOUBLE_EQ(ag[0][0], 2.43727779489243124e-02);
  EXPECT_DOUBLE_EQ(ag[0][1], 2.71756510552775242e-02);
  EXPECT_DOUBLE_EQ(ag[0][2], 4.74866188111812346e-02);
  EXPECT_DOUBLE_EQ(ag[1][0], 4.79999416660118611e-03);
  EXPECT_DOUBLE_EQ(ag[1][1], 2.28673224800665653e-02);
  EXPECT_DOUBLE_EQ(ag[1][2], -6.27597464661684978e-02);
  EXPECT_DOUBLE_EQ(ag[2][0], -3.40468666446705568e-03);
  EXPECT_DOUBLE_EQ(ag[2][1], -2.71288483037145679e-02);
  EXPECT_DOUBLE_EQ(ag[2][2], 1.94086283932249776e-02);
  EXPECT_DOUBLE_EQ(ag[3][0], -2.57680854510584506e-02);
  EXPECT_DOUBLE_EQ(ag[3][1], -2.29141252316295216e-02);
  EXPECT_DOUBLE_EQ(ag[3][2], -4.13550073823771430e-03);
}

TEST(EgBondTest, formaldehydeOop)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));
  Mol m = makeFormaldehydeC();
  std::vector<OopAtom> oa(4);
  oa[0].element = 6;
  oa[0].pi = 1;
  oa[1].element = 8;
  oa[1].pi = 1;
  oa[2].element = oa[3].element = 1;
  std::vector<Oop> oops;
  std::vector<OopTerm> ot;
  Environment env;
  std::vector<double> pbo(10, 0.0);
  pbo[packedIndex(0, 1)] = 0.8;
  pbo[packedIndex(0, 2)] = 0.05;
  pbo[packedIndex(0, 3)] = 0.05;
  ASSERT_TRUE(buildOutOfPlane(oa, m.neighbours, m.xyz, pbo, param.group,
                              param.repz, param, gen, oops, ot, env));
  ASSERT_EQ(oops.size(), 1u);
  double e = 0.0, ag[4][3] = { { 0 } }, ds[3][3] = { { 0 } };
  const Oop& o = oops[0];
  torsionEnergyGradient(o.center, o.first, o.second, o.third, o.kind,
                        ot[0].phase, ot[0].forceConstant, m.numbers, m.xyz,
                        param.rcov, param.angleCutT, false, e, ag, ds);
  EXPECT_DOUBLE_EQ(e, 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[0][0], -0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[0][1], -0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[0][2], 4.30502948189318843e-16);
  EXPECT_DOUBLE_EQ(ag[1][0], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[1][1], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[1][2], -1.54114960742330115e-16);
  EXPECT_DOUBLE_EQ(ag[2][0], 0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[2][1], -0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[2][2], -1.54114960742330115e-16);
  EXPECT_DOUBLE_EQ(ag[3][0], -0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[3][1], -0.00000000000000000e+00);
  EXPECT_DOUBLE_EQ(ag[3][2], -1.22273026704658612e-16);
}

} // namespace Xtb
} // namespace Avogadro
