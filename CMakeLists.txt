cmake_minimum_required(VERSION 3.10)

project(KFParticleLibrary)

# Используем C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-msse4.1)
# Подключаем ROOT
find_package(ROOT REQUIRED COMPONENTS Core Hist Matrix Physics EG Gpad Graf Graf3d RIO MathCore)
#find_package(Vc REQUIRED)

# Указываем директории с заголовками
include_directories(${ROOT_INCLUDE_DIRS} ${Vc_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR})
include_directories(
    ${PROJECT_SOURCE_DIR}/KFParticle
    ${PROJECT_SOURCE_DIR}/KFParticlePerformance
    ${PROJECT_SOURCE_DIR}/KFParticleTest
)
# Указываем файлы исходного кода
set(SOURCES
    KFParticle/KFParticleTopoReconstructor.cxx
    KFParticle/KFVertex.cxx
    KFParticle/KFPTrack.cxx
    KFParticle/KFPTrackVector.cxx
    KFParticle/KFPVertex.cxx
    KFParticle/KFParticlePVReconstructor.cxx
    KFParticle/KFParticleDatabase.cxx
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

# Указываем заголовочные файлы (для установки)
set(HEADERS
    KFParticle/KFParticle.h
    KFParticle/KFVertex.h
    KFParticlePerformance/KFPartEfficiencies.h
    KFParticlePerformance/KFMCParticle.h
    KFParticleTest/KFParticleTest.h
)

# Генерация словаря ROOT (если нужно)
ROOT_GENERATE_DICTIONARY(G__KFParticle ${HEADERS} LINKDEF KFLinkDef.h)

# Создаём библиотеку
add_library(KFParticle SHARED ${SOURCES} G__KFParticle.cxx)
target_link_libraries(KFParticle ${ROOT_LIBRARIES} ${Vc_LIBRARIES})

target_compile_definitions(KFParticle PRIVATE -DDO_TPCCATRACKER_EFF_PERFORMANCE -DHomogeneousField -DUSE_TIMERS)

# Установка файлов
install(FILES ${HEADERS} DESTINATION include)
install(FILES ${CMAKE_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}KFParticle_rdict.pcm DESTINATION lib OPTIONAL)
install(FILES ${CMAKE_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}KFParticle.rootmap DESTINATION lib)
install(TARGETS KFParticle DESTINATION lib)
