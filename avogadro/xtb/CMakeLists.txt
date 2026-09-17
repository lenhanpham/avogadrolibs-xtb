add_library(Xtb)

avogadro_headers(Xtb
  bridge.h
  constants.h
  coordination.h
  eeq.h
  elements.h
  environment.h
  gfn0params.h
  gfnffdata.h
  gfnffhyb.h
  gfnffparams.h
  gfnffsetup.h
  gfnfftopo.h
  gffgraph.h
  gfnffrab.h
  molecule.h
  settings.h
)

target_sources(Xtb PRIVATE
  bridge.cpp
  coordination.cpp
  eeq.cpp
  elements.cpp
  environment.cpp
  gfn0params.cpp
  gfnffparams.cpp
  gfnffsetup.cpp
  gfnfftopo.cpp
  gfnffrab.cpp
  gfnffhyb.cpp
  gffgraph.cpp
)

avogadro_add_library(Xtb)

target_link_libraries(Xtb
  PUBLIC Avogadro::Core)

# The ported xtb module stays LGPL-licensed (see COPYING.LESSER); keep it a
# shared library in binary distributions so end users can relink a modified
# build, as the license requires for static linking.
if(NOT BUILD_SHARED_LIBS)
  message(WARNING
    "Building Avogadro::Xtb statically: binary distributions must ship "
    "object files or build scripts to satisfy the LGPL relinking clause.")
endif()
