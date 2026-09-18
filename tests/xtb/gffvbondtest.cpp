/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "avogadro/xtb/COPYING.LESSER".
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/xtb/environment.h>
#include <avogadro/xtb/gfnffdata.h>
#include <avogadro/xtb/gfnffsetup.h>
#include <avogadro/xtb/gfnffvbond.h>

#include <vector>

namespace Avogadro {
namespace Xtb {

// Reference values verified bitwise against the unmodified vbond loop in
// xtb src/gfnff/gfnff_ini.f90 (lines 1238-1439, print dropped), compiled
// with gfortran -ffp-contract=off and driven with identical inputs
// (tables, generator scalars, per-atom/per-bond data). All 38 synthetic
// cases cover bond types 1-7, bridge/hypervalent/XH/ring/pi/charge
// corrections and every TM-ligand/mtyp branch.

// (element, hyb, imetal, pi, itag, charge, mchar)
static const int kElements[] = {
  6, 6, 6, 1, 8, 8, 7, 7, 5, 9, 17, 11, 12, 13, 26, 15, 16, 6, 1, 1,
  6, 46, 78, 1, 6, 7, 17, 16, 6, 1, 7,
};
static const int kHyb[] = {
  3, 2, 1, 0, 2, 3, 3, 1, 3, 1, 3, 0, 0, 0, 0, 3, 2, 2, 0, 0,
  3, 0, 0, 0, 1, 1, 0, 5, 3, 0, 1,
};
static const int kImetal[] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 0, 0, 0, 0, 0,
  0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0,
};
static const int kPi[] = {
  0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
  0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1,
};
static const int kItag[] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
};
static const double kCharges[] = {
  0.0, 0.0, 0.0, 0.0, 0.0, -0.3, 0.1, 0.0, 0.0, -0.2, 0.0, 0.5, 0.4, 0.3,
  0.2, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.1, 0.1, 0.0, 0.0, 0.0, -0.1, 0.0,
  0.0, 0.0, 0.0,
};
static const double kMchar[] = {
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1,
  0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.05, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0,
};

// (first, second, ring, ringFirst, ringSecond, pibo)
static const int kBondAtoms[][2] = {
  { 0, 0 }, { 1, 1 }, { 2, 2 }, { 2, 0 }, { 2, 1 }, { 0, 3 }, { 19, 20 },
  { 18, 17 }, { 8, 3 }, { 6, 3 }, { 5, 3 }, { 9, 9 }, { 9, 0 }, { 2, 4 },
  { 1, 6 }, { 10, 10 }, { 27, 9 }, { 11, 10 }, { 12, 5 }, { 13, 10 },
  { 13, 3 }, { 14, 2 }, { 14, 7 }, { 14, 6 }, { 14, 15 }, { 14, 16 },
  { 14, 26 }, { 14, 23 }, { 22, 23 }, { 14, 14 }, { 14, 21 }, { 21, 21 },
  { 14, 24 }, { 1, 4 }, { 28, 29 }, { 2, 7 }, { 14, 30 }, { 30, 14 },
};
static const int kBondRings[][3] = {
  { 0, 99, 99 }, { 6, 5, 5 }, { 0, 99, 99 }, { 0, 99, 99 }, { 0, 99, 99 },
  { 0, 99, 99 }, { 3, 3, 3 }, { 0, 99, 99 }, { 0, 99, 99 }, { 0, 99, 99 },
  { 0, 99, 99 }, { 0, 99, 99 }, { 0, 99, 99 }, { 0, 99, 99 }, { 0, 99, 99 },
  { 0, 99, 99 }, { 0, 99, 99 }, { 0, 99, 99 }, { 0, 99, 99 }, { 0, 99, 99 },
  { 0, 99, 99 }, { 0, 99, 99 }, { 0, 99, 99 }, { 0, 99, 99 }, { 0, 99, 99 },
  { 0, 99, 99 }, { 0, 99, 99 }, { 0, 99, 99 }, { 0, 99, 99 }, { 0, 99, 99 },
  { 0, 99, 99 }, { 0, 99, 99 }, { 0, 99, 99 }, { 0, 99, 99 }, { 5, 5, 99 },
  { 0, 99, 99 }, { 0, 99, 99 }, { 0, 99, 99 },
};
static const double kBondPibo[] = {
  0.0, 0.6, 0.9, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.7, 0.0, 0.85, 0.0, 0.0,
};

// __EXPECTATIONS__

TEST(VbondTest, carbonylCarbonType)
{
  // Formaldehyde-like: C(0) bonded to pi O(1), pi C itself.
  std::vector<std::vector<int>> neighbours(3);
  neighbours[0].push_back(1);
  neighbours[0].push_back(2);
  neighbours[1].push_back(0);
  neighbours[2].push_back(0);
  std::vector<int> numbers = { 6, 8, 1 };
  std::vector<int> pi = { 1, 1, 0 };
  EXPECT_EQ(carbonylCarbonType(neighbours, numbers, pi, 0), 1);
  EXPECT_EQ(carbonylCarbonType(neighbours, numbers, pi, 2), 0);
  pi[0] = 0;
  EXPECT_EQ(carbonylCarbonType(neighbours, numbers, pi, 0), 0);
}

TEST(VbondTest, terms)
{
  GffData param;
  GffGenerator gen;
  ASSERT_TRUE(loadGffParams(GffAngewChem2020, param, gen));

  const int n = sizeof(kElements) / sizeof(kElements[0]);
  std::vector<VbondAtom> atoms(n);
  std::vector<int> numbers(n);
  for (int i = 0; i < n; ++i) {
    atoms[i].element = kElements[i];
    atoms[i].hyb = kHyb[i];
    atoms[i].charge = kCharges[i];
    atoms[i].itag = kItag[i];
    atoms[i].imetal = kImetal[i];
    atoms[i].pi = kPi[i];
    atoms[i].mchar = kMchar[i];
    numbers[i] = kElements[i];
  }
  // Carbonyl-C adjacency for the aldehyde test (atom 17: O 4, H 18).
  std::vector<std::vector<int>> neighbours(n);
  neighbours[17].push_back(4);
  neighbours[17].push_back(18);
  neighbours[4].push_back(17);
  neighbours[18].push_back(17);

  const int nbond = sizeof(kBondAtoms) / sizeof(kBondAtoms[0]);
  // Coordination counts: per-atom maxima like sum(nb(...)).
  const int kCoords[][2] = {
    { 0, 4 }, { 1, 3 }, { 2, 2 }, { 3, 1 }, { 4, 1 }, { 5, 2 }, { 6, 3 },
    { 7, 2 }, { 8, 4 }, { 9, 1 }, { 10, 1 }, { 11, 6 }, { 12, 6 },
    { 13, 4 }, { 14, 6 }, { 15, 4 }, { 16, 2 }, { 17, 3 }, { 18, 1 },
    { 19, 2 }, { 20, 4 }, { 21, 6 }, { 22, 4 }, { 23, 1 }, { 24, 2 },
    { 25, 2 }, { 26, 1 }, { 27, 6 }, { 28, 4 }, { 29, 1 }, { 30, 1 },
  };
  std::vector<int> coordOf(n, 0);
  for (const auto& kc : kCoords)
    coordOf[kc[0]] = kc[1];

  std::vector<VbondBond> bonds;
  for (int i = 0; i < nbond; ++i) {
    VbondBond b;
    b.first = kBondAtoms[i][0];
    b.second = kBondAtoms[i][1];
    b.guess = 2.5;
    b.pibo = kBondPibo[i];
    b.ring = kBondRings[i][0];
    b.ringFirst = kBondRings[i][1];
    b.ringSecond = kBondRings[i][2];
    b.coordFirst = coordOf[b.first];
    b.coordSecond = coordOf[b.second];
    bonds.push_back(b);
  }
  auto row6of = [](int z) {
    if (z <= 2)
      return 1;
    if (z <= 10)
      return 2;
    if (z <= 18)
      return 3;
    if (z <= 36)
      return 4;
    if (z <= 54)
      return 5;
    return 6;
  };
  std::vector<int> row6(n);
  for (int i = 0; i < n; ++i)
    row6[i] = row6of(numbers[i]);

  Environment env;
  std::vector<VbondTerm> terms;
  std::vector<int> outBtypes;
  ASSERT_TRUE(buildVbondTerms(bonds, atoms, neighbours, numbers, param.group,
                              param.metal, param.en, param.bond, row6, param,
                              gen, terms, outBtypes, env))
    << env.errorMessage();
  ASSERT_EQ(terms.size(), static_cast<size_t>(nbond));
  ASSERT_EQ(outBtypes.size(), static_cast<size_t>(nbond));

  static const int kWantBtype[] = { 1, 2, 3, 3, 2, 1, 1, 1, 1, 1, 1, 3, 3, 3, 2, 1, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 7, 7, 7, 6, 2, 1, 3, 5, 5 };
  for (int i = 0; i < nbond; ++i)
    EXPECT_EQ(outBtypes[i], kWantBtype[i]) << "bond " << i;
  static const double kWantShift[] = { -1.09999999403953552e-01, -1.88199999403953544e-01, -2.90199999403953579e-01, -1.09999999403953552e-01, -1.54199999403953569e-01, -1.82000000029802322e-01, -1.82000000029802322e-01, -1.60000000149011612e-01, -1.82000000029802322e-01, -1.82000000029802322e-01, -1.82000000029802322e-01, 1.09999999403953552e-01, -1.09999999403953552e-01, -1.09999999403953552e-01, -1.20199999403953567e-01, -2.19999998807907104e-01, -8.00000000745058060e-02, -1.99999958276748657e-02, 9.00000035762786865e-02, -1.69999998062849045e-01, -5.99999986588954926e-02, -4.09999993443489086e-01, 5.10000006556510899e-01, 4.00000065565109253e-02, -1.79999992251396179e-01, -1.79999992251396179e-01, -1.79999992251396179e-01, 4.00000065565109253e-02, 4.00000065565109253e-02, -1.39999985694885254e-01, -1.99999984353780746e-01, -2.59999983012676239e-01, -1.69999998807907116e-01, -2.22199999403953546e-01, -1.82000000029802322e-01, -2.73199999403953564e-01, 4.00000065565109253e-02, 4.00000065565109253e-02 };
  static const double kWantSteep[] = { 4.67792797484965384e-01, 4.90519065871560966e-01, 5.60591726730230611e-01, 4.98416444135902859e-01, 5.13643043954921885e-01, 4.82285765771657993e-01, 4.82285765771657993e-01, 4.89785434339234493e-01, 4.70821536982033551e-01, 5.51272323249272445e-01, 6.49706271153752790e-01, 4.01507848024061698e-01, 6.43440009389355416e-01, 6.35555936063904481e-01, 5.23622057953419651e-01, 4.67792797484965384e-01, 6.33395498572032367e-01, 6.32529086997132284e-01, 6.18085805480545702e-01, 5.47379948674885375e-01, 4.79324238861461638e-01, 4.54299780970791878e-01, 4.29684912062278435e-01, 4.29684912062278435e-01, 4.64419542239412786e-01, 4.53151938521677344e-01, 4.21751508019528432e-01, 4.64229536541149512e-01, 4.67792797484965384e-01, 6.30664387588900088e-01, 6.27101126645084217e-01, 6.30664387588900088e-01, 4.33467368283079191e-01, 5.84232446624176283e-01, 4.82285765771657993e-01, 5.88997961158041039e-01, 4.29684912062278435e-01, 4.29684912062278435e-01 };
  static const double kWantPref[] = { -1.51903797889919145e-01, -2.42043511557797159e-01, -4.76719688917933260e-01, -2.01029486127519003e-01, -2.67164905941438324e-01, -1.64816772070439127e-01, -2.04207971321564236e-01, -1.68976745277118284e-01, -1.08811166975994175e-01, -1.71988915763077510e-01, -1.34977860366525648e-01, -3.28099375307999985e-02, -3.91193444625152625e-02, -2.38371602051396791e-01, -2.09240905848818298e-01, -4.82948925629607603e-02, -2.95138581538734804e-02, -7.26163497270772718e-03, -1.79289996918290256e-02, -3.87101948206184723e-02, -6.37072107110298186e-02, -7.91871221793751662e-02, -2.07881817401661881e-02, -3.78123639768018632e-02, -5.46033105076288439e-02, -3.38545186985727692e-02, -4.95832353872095247e-02, -3.25686037592071570e-02, -3.69840187237760465e-02, -1.75599286660182777e-02, -1.46169356157753467e-02, -1.23001662217386530e-02, -6.17659552999126268e-02, -2.29730297061205874e-01, -1.68113107511847909e-01, -4.54501598375569082e-01, -5.19704543504154667e-02, -5.19704543504154667e-02 };
  for (int i = 0; i < nbond; ++i) {
    EXPECT_DOUBLE_EQ(terms[i].shift, kWantShift[i]) << "shift " << i;
    EXPECT_DOUBLE_EQ(terms[i].steepness, kWantSteep[i]) << "steep " << i;
    EXPECT_DOUBLE_EQ(terms[i].prefactor, kWantPref[i]) << "pref " << i;
  }
}

} // namespace Xtb
} // namespace Avogadro
