/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "avogadro/xtb/COPYING.LESSER".
******************************************************************************/

#include <gtest/gtest.h>

#include <avogadro/xtb/elements.h>

using Avogadro::Xtb::elementNumber;
using Avogadro::Xtb::elementSymbol;
using Avogadro::Xtb::maxElement;

TEST(ElementsTest, resolveSymbols)
{
  EXPECT_EQ(elementNumber("H"), 1);
  EXPECT_EQ(elementNumber("h"), 1);
  EXPECT_EQ(elementNumber("He"), 2);
  EXPECT_EQ(elementNumber("HE"), 2);
  EXPECT_EQ(elementNumber("C"), 6);
  EXPECT_EQ(elementNumber("Og"), 118);
  EXPECT_EQ(elementNumber("Xx"), 0);
  EXPECT_EQ(elementNumber(""), 0);
  EXPECT_EQ(elementNumber("   "), 0);
}

TEST(ElementsTest, canonicalSymbols)
{
  EXPECT_EQ(elementSymbol(1), "h ");
  EXPECT_EQ(elementSymbol(2), "he");
  EXPECT_EQ(elementSymbol(6), "c ");
  EXPECT_EQ(elementSymbol(118), "og");
  EXPECT_EQ(elementSymbol(0), "");
  EXPECT_EQ(elementSymbol(119), "");
}

TEST(ElementsTest, roundTrip)
{
  for (int z = 1; z <= maxElement; ++z)
    EXPECT_EQ(elementNumber(elementSymbol(z)), z);
}
