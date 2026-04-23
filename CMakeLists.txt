cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)

project(KFParticleLibrary)

# add at the top of your KFParticle CMakeLists.txt
if(NOT CMAKE_CXX_FLAGS MATCHES "-msse4.1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
endif()

macro(add_target_property _target _prop _value)
  get_target_property(_oldprop "${_target}" ${_prop})
  if(NOT _oldprop)
    set_target_properties("${_target}" PROPERTIES ${_prop} "${_value}")
  else(NOT _oldprop)
    set_target_properties("${_target}" PROPERTIES ${_prop} "${_oldprop} ${_value}")
  endif(NOT _oldprop)
endmacro(add_target_property)

if (CMAKE_VERSION VERSION_LESS "3.1")
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND NOT ${CMAKE_CXX_FLAGS} MATCHES "-std=")
  set (CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
endif ()
else ()
  if (NOT DEFINED CMAKE_CXX_STANDARD AND NOT ${CMAKE_CXX_FLAGS} MATCHES "-std=")
    set (CMAKE_CXX_STANDARD 11)
  endif()
endif ()

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS} ${ROOTSYS})

set(FIXTARGET FALSE CACHE BOOL "Compile for fix target geometry.")

find_package(ROOT REQUIRED COMPONENTS Core Hist Matrix Physics EG Gpad Graf Graf3d RIO MathCore)
If(DEFINED Vc_INCLUDE_DIR AND Vc_LIBRARIES)
  Message("Vc found")
Else()
  find_package(Vc REQUIRED)
EndIf()

Set(INCLUDE_DIRECTORIES
  ${PROJECT_SOURCE_DIR}/KFParticle
  ${PROJECT_SOURCE_DIR}/KFParticlePerformance
  ${PROJECT_SOURCE_DIR}/KFParticleTest
  ${PROJECT_SOURCE_DIR}/KFParticle/KFPSimd
  ${PROJECT_SOURCE_DIR}/KFParticle/KFPSimd/Utils
  ${PROJECT_SOURCE_DIR}/KFParticle/KFPSimd/SSE
)

include(${ROOT_USE_FILE})
set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)
Include_Directories( ${CMAKE_SOURCE_DIR} ${Vc_INCLUDE_DIR} ${INCLUDE_DIRECTORIES} )
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

set (SIMD_HEADERS
  KFParticle/KFPSimd/simd.h
)

set (SIMD_UTILS_HEADERS
  KFParticle/KFPSimd/Utils/macros.h
  KFParticle/KFPSimd/Utils/memory.h
  KFParticle/KFPSimd/Utils/simd_detect.h
  KFParticle/KFPSimd/Utils/tag.h
)

set (SIMD_SSE_HEADERS
  KFParticle/KFPSimd/SSE/constants.h
  KFParticle/KFPSimd/SSE/float32.h
  KFParticle/KFPSimd/SSE/int32.h
  KFParticle/KFPSimd/SSE/mask32.h
  KFParticle/KFPSimd/SSE/types.h
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
  ROOT_GENERATE_DICTIONARY(G__KFParticle ${HEADERS} LINKDEF KFLinkDef.h OPTIONS "-DDO_TPCCATRACKER_EFF_PERFORMANCE" "-DNonhomogeneousField" "-DCBM" "-DUSE_TIMERS")
#  ROOT_GENERATE_DICTIONARY(G__KFParticle ${HEADERS} LINKDEF KFLinkDef.h OPTIONS "-DDO_TPCCATRACKER_EFF_PERFORMANCE" "-DHomogeneousField" "-DCBM" "-DUSE_TIMERS")

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
  target_link_libraries(KFParticle ${ROOT_LIBRARIES} ${Vc_LIBRARIES} )
  add_target_property(KFParticle COMPILE_FLAGS "-DDO_TPCCATRACKER_EFF_PERFORMANCE -DNonhomogeneousField -DCBM -DUSE_TIMERS")
#  add_target_property(KFParticle COMPILE_FLAGS "-DDO_TPCCATRACKER_EFF_PERFORMANCE -DHomogeneousField -DCBM -DUSE_TIMERS")
else(FIXTARGET)
  ROOT_GENERATE_DICTIONARY(G__KFParticle ${HEADERS} LINKDEF KFLinkDef.h OPTIONS "-DDO_TPCCATRACKER_EFF_PERFORMANCE" "-DHomogeneousField" "-DUSE_TIMERS")
  add_library(KFParticle SHARED ${SOURCES} G__KFParticle.cxx)
  target_link_libraries(KFParticle ${ROOT_LIBRARIES} ${Vc_LIBRARIES} )
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
#  KFParticle/KFParticleBaseSIMD.h
  KFParticle/KFParticleDef.h
  KFParticle/KFPSimdAllocator.h
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

install(FILES ${HEADERS} ${NODICT_HEADERS} DESTINATION include/KFParticle)

install(FILES ${SIMD_HEADERS} DESTINATION include/KFParticle)
install(FILES ${SIMD_UTILS_HEADERS} DESTINATION include/KFParticle/Utils) 
install(FILES ${SIMD_SSE_HEADERS} DESTINATION include/KFParticle/SSE) 

install(FILES ${CMAKE_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}KFParticle_rdict.pcm DESTINATION lib OPTIONAL)
install(FILES ${CMAKE_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}KFParticle.rootmap DESTINATION lib)
install(TARGETS KFParticle DESTINATION lib)
