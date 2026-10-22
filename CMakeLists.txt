cmake_minimum_required(VERSION 3.14.0 FATAL_ERROR)

project(KFParticleLibrary)

add_compile_options(-O3 -Wall -Wextra -Wsign-promo -Wctor-dtor-privacy -Wreorder -Wno-deprecated -Wno-parentheses)
# add_compile_options(-Weffc++ -Wnon-virtual-dtor -Woverloaded-virtual -Wold-style-cast) # for the future
If(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wno-pmf-conversions -Wstrict-null-sentinel -Wno-non-template-friend -Wmissing-field-initializers)
EndIf()


if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
message(STATUS "===== KFParticle project is built standalone ======")
  add_compile_options(-msse4.2)
else()
  message(STATUS "===== KFParticle project is embedded, use SIMD flags from the host ======")
endif()


macro(add_target_property _target _prop _value)
  get_target_property(_oldprop "${_target}" ${_prop})
  if(NOT _oldprop)
    set_target_properties("${_target}" PROPERTIES ${_prop} "${_value}")
  else(NOT _oldprop)
    set_target_properties("${_target}" PROPERTIES ${_prop} "${_oldprop} ${_value}")
  endif(NOT _oldprop)
endmacro(add_target_property)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS} ${ROOTSYS})

set(FIXTARGET FALSE CACHE BOOL "Compile for fix target geometry.")

find_package(ROOT REQUIRED COMPONENTS Core Hist Matrix Physics EG Gpad Graf Graf3d RIO MathCore)

Set(INCLUDE_DIRECTORIES
  ${PROJECT_SOURCE_DIR}/KFParticle
  ${PROJECT_SOURCE_DIR}/KFParticlePerformance
  ${PROJECT_SOURCE_DIR}/KFParticleTest
)

include(${ROOT_USE_FILE})
set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)
Include_Directories( ${CMAKE_SOURCE_DIR} ${INCLUDE_DIRECTORIES} )
set(CMAKE_INCLUDE_DIRECTORIES_BEFORE OFF)

set(LINK_DIRECTORIES
  ${ROOT_LIBRARY_DIR}
)

link_directories( ${LINK_DIRECTORIES})

set (SOURCES
  KFParticle/KFParticleTopoReconstructor.cxx
  KFParticle/KFVertex.cxx
  KFParticle/KFPTrack.cxx
  KFParticle/KFPTrackVector.cxx
  KFParticle/KFPVertex.cxx
  KFParticle/KFParticlePVReconstructor.cxx
  KFParticle/KFParticleDatabase.cxx
#  KFParticle/KFParticleBase.cxx
#  KFParticle/KFParticleBaseSIMD.cxx
  KFParticle/KFParticle.cxx
  KFParticle/KFParticleSIMD.cxx
  KFParticle/KFParticleFinder.cxx
  KFParticle/KFPEmcCluster.cxx
  KFParticlePerformance/KFMCVertex.cxx
  KFParticlePerformance/KFParticlePerformanceBase.cxx
  KFParticlePerformance/KFTopoPerformance.cxx
  KFParticlePerformance/KFMCParticle.cxx
  KFParticleTest/KFParticleTest.cxx
)

set (HEADERS
#  KFParticle/KFParticleBase.h
  KFParticle/KFParticle.h
  KFParticle/KFVertex.h
  KFParticlePerformance/KFPartEfficiencies.h
  KFParticlePerformance/KFMCParticle.h
  KFParticleTest/KFParticleTest.h
)

Message("CMAKE_BUILD_DIRECTORY: ${CMAKE_BUILD_DIRECTORY}")

if(FIXTARGET)

  find_program(BIN_SED sed)
  if(BIN_SED)
    add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/G__KFParticle_mod.cxx
                       COMMAND sed -e "s#${CMAKE_SOURCE_DIR}/##g" ${CMAKE_BINARY_DIR}/G__KFParticle.cxx > ${CMAKE_BINARY_DIR}/G__KFParticle_mod.cxx
                       DEPENDS ${CMAKE_BINARY_DIR}/G__KFParticle.cxx
                      )
  else()
    Message("Can't correct KFParticle dictionary. Running KFParticle from an installed CbmRoot will not be possible")
  endif()

#  set(SourceFile "${CMAKE_BINARY_DIR}/G__KFParticle.cxx")
#  set(SourceFile_new "${CMAKE_BINARY_DIR}/G__KFParticle_mod.cxx")
#  file(READ ${SourceFile} Contents)
#  STRING(REGEX REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/" "" Contents_Modified "${Contents}" )
#  file(WRITE ${SourceFile_new} ${Contents_Modified})

  add_library(KFParticle SHARED ${SOURCES} G__KFParticle_mod.cxx)
  target_include_directories(KFParticle
    PUBLIC
        ${PROJECT_SOURCE_DIR}
        ${PROJECT_SOURCE_DIR}/KFParticle
        ${PROJECT_SOURCE_DIR}/KFParticlePerformance
        ${PROJECT_SOURCE_DIR}/KFParticleTest
)
  target_link_libraries(KFParticle ${ROOT_LIBRARIES} )
  add_target_property(KFParticle COMPILE_FLAGS "-DDO_TPCCATRACKER_EFF_PERFORMANCE -DNonhomogeneousField -DCBM -DUSE_TIMERS")
#  add_target_property(KFParticle COMPILE_FLAGS "-DDO_TPCCATRACKER_EFF_PERFORMANCE -DHomogeneousField -DCBM -DUSE_TIMERS")


  ROOT_GENERATE_DICTIONARY(G__KFParticle ${HEADERS} LINKDEF KFLinkDef.h OPTIONS "-DDO_TPCCATRACKER_EFF_PERFORMANCE" "-DNonhomogeneousField" "-DCBM" "-DUSE_TIMERS")
#  ROOT_GENERATE_DICTIONARY(G__KFParticle ${HEADERS} LINKDEF KFLinkDef.h OPTIONS "-DDO_TPCCATRACKER_EFF_PERFORMANCE" "-DHomogeneousField" "-DCBM" "-DUSE_TIMERS")

else(FIXTARGET)
  ROOT_GENERATE_DICTIONARY(G__KFParticle ${HEADERS} LINKDEF KFLinkDef.h OPTIONS "-DDO_TPCCATRACKER_EFF_PERFORMANCE" "-DHomogeneousField" "-DUSE_TIMERS")
  add_library(KFParticle SHARED ${SOURCES} G__KFParticle.cxx)
  target_link_libraries(KFParticle ${ROOT_LIBRARIES} )
  add_target_property(KFParticle COMPILE_FLAGS "-DDO_TPCCATRACKER_EFF_PERFORMANCE -DHomogeneousField -DUSE_TIMERS")
endif(FIXTARGET)

if (ROOT_VERSION_MAJOR LESS 6)
    add_custom_target(libKFParticle.rootmap ALL DEPENDS KFParticle COMMAND rlibmap -o libKFParticle.rootmap -l libKFParticle.so -c ${PROJECT_SOURCE_DIR}/KFLinkDef.h)
endif (ROOT_VERSION_MAJOR LESS 6)


Set(NODICT_HEADERS
  KFParticle/KFParticleTopoReconstructor.h
  KFParticle/KFParticlePVReconstructor.h
  KFParticle/KFPVertex.h
  KFParticle/KFPTrack.h
  KFParticle/KFParticleFinder.h
  KFParticle/KFParticleSIMD.h
  KFParticle/KFParticleDef.h
  KFParticle/KFPSimd.h
  KFParticle/KFPSimdPaddedVector.h
  KFParticle/KFPTrackVector.h
  KFParticle/KFPEmcCluster.h
  KFParticle/KFParticleField.h
  KFParticle/KFPTrackVector.h
  KFParticle/KFParticleDatabase.h
  KFParticle/KFParticleMath.h
  KFParticlePerformance/KFTopoPerformance.h
  KFParticlePerformance/KFParticlePerformanceBase.h
  KFParticlePerformance/KFMCCounter.h
  KFParticlePerformance/KFPVEfficiencies.h
  KFParticlePerformance/KFMCVertex.h
  KFParticlePerformance/KFMCTrack.h
  KFParticlePerformance/KFPartMatch.h
)

install(TARGETS KFParticle DESTINATION lib)

install(DIRECTORY KFParticle
        DESTINATION include
        FILES_MATCHING
        PATTERN "*.h")

install(DIRECTORY KFParticlePerformance
        DESTINATION include
        FILES_MATCHING
        PATTERN "*.h")

install(DIRECTORY KFParticleTest
        DESTINATION include
        FILES_MATCHING
        PATTERN "*.h")


#if(EXISTS ${CMAKE_BINARY_DIR}/lib/libKFParticle.rootmap)     
  install(FILES ${CMAKE_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}KFParticle.rootmap DESTINATION lib)
  install(FILES ${CMAKE_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}KFParticle_rdict.pcm DESTINATION lib)
#endif()

#install(FILES ${HEADERS} ${NODICT_HEADERS} DESTINATION include/KFParticle)

#install(FILES ${CMAKE_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}KFParticle_rdict.pcm DESTINATION lib OPTIONAL)
#install(FILES ${CMAKE_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}KFParticle.rootmap DESTINATION lib)
#install(TARGETS KFParticle DESTINATION lib)
