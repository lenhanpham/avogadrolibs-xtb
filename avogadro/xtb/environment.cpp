/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/type/environment.f90, Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#include "environment.h"

namespace Avogadro {
namespace Xtb {

Environment::Environment() : m_verbosity(1) {}

void Environment::error(const std::string& message, const std::string& source)
{
  m_errors.push_back("[" + source + "] " + message);
}

std::string Environment::errorMessage() const
{
  std::string out;
  for (size_t i = 0; i < m_errors.size(); ++i) {
    if (i > 0)
      out += "\n";
    out += m_errors[i];
  }
  return out;
}

void Environment::clearErrors()
{
  m_errors.clear();
}

} // namespace Xtb
} // namespace Avogadro
