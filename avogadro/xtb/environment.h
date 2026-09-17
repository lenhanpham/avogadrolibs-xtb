/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/type/environment.f90, Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#ifndef AVOGADRO_XTB_ENVIRONMENT_H
#define AVOGADRO_XTB_ENVIRONMENT_H

#include <string>
#include <vector>

namespace Avogadro {
namespace Xtb {

// Computational environment: verbosity, error collection and output routing.
// Ported from xtb's TEnvironment. Avogadro code must never crash, so errors
// are collected here instead of aborting (xtb's env%error / env%terminate).
class Environment
{
public:
  Environment();
  ~Environment() = default;

  // Record an error with the calling context; further compute calls that
  // check failed() must be skipped by the caller.
  void error(const std::string& message, const std::string& source);
  // True once any error was recorded.
  bool failed() const { return !m_errors.empty(); }
  // Human-readable listing of all recorded errors.
  std::string errorMessage() const;
  void clearErrors();

  int verbosity() const { return m_verbosity; }
  void setVerbosity(int level) { m_verbosity = level; }

private:
  int m_verbosity;
  std::vector<std::string> m_errors;
};

} // namespace Xtb
} // namespace Avogadro

#endif // AVOGADRO_XTB_ENVIRONMENT_H
