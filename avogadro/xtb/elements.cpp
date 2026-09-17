/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/elem.f90, Copyright (C) 2017-2020 Stefan Grimme.
******************************************************************************/

#include "elements.h"

#include <array>
#include <cctype>

namespace Avogadro {
namespace Xtb {

// Lowercase two-character symbols, H..Og, in atomic-number order.
static const std::array<std::string, maxElement> symbols = {
  "h ", "he", "li", "be", "b ", "c ", "n ", "o ", "f ", "ne", "na", "mg",
  "al", "si", "p ", "s ", "cl", "ar", "k ", "ca", "sc", "ti", "v ", "cr",
  "mn", "fe", "co", "ni", "cu", "zn", "ga", "ge", "as", "se", "br", "kr",
  "rb", "sr", "y ", "zr", "nb", "mo", "tc", "ru", "rh", "pd", "ag", "cd",
  "in", "sn", "sb", "te", "i ", "xe", "cs", "ba", "la", "ce", "pr", "nd",
  "pm", "sm", "eu", "gd", "tb", "dy", "ho", "er", "tm", "yb", "lu", "hf",
  "ta", "w ", "re", "os", "ir", "pt", "au", "hg", "tl", "pb", "bi", "po",
  "at", "rn", "fr", "ra", "ac", "th", "pa", "u ", "np", "pu", "am", "cm",
  "bk", "cf", "es", "fm", "md", "no", "lr", "rf", "db", "sg", "bh", "hs",
  "mt", "ds", "rg", "cn", "nh", "fl", "mc", "lv", "ts", "og"
};

int elementNumber(const std::string& symbol)
{
  // Mirror elem(): scan for the last non-blank, then take up to two
  // alphabetic characters, lowercased; stop at a blank/tab after content.
  std::string e;
  std::string::size_type last = symbol.find_last_not_of(" \t");
  if (last == std::string::npos)
    return 0;
  for (std::string::size_type j = 0; j <= last && e.size() < 2; ++j) {
    unsigned char c = static_cast<unsigned char>(symbol[j]);
    if (std::isalpha(c))
      e += static_cast<char>(std::tolower(c));
    else if (!e.empty())
      break;
  }
  while (e.size() < 2)
    e += ' ';
  for (int i = 0; i < maxElement; ++i) {
    if (e == symbols[i])
      return i + 1;
  }
  return 0;
}

std::string elementSymbol(int atomicNumber)
{
  if (atomicNumber < 1 || atomicNumber > maxElement)
    return std::string();
  return symbols[atomicNumber - 1];
}

} // namespace Xtb
} // namespace Avogadro
