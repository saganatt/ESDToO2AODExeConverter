cmake_minimum_required(VERSION 3.8)
project(ESDConverter)

if(NOT DEFINED ENV{ENVIRONMENT})
    message(FATAL_ERROR "You must set ENVIRONMENT environment variable")
endif()

if(NOT DEFINED ENV{ARCHITECTURE})
    message(FATAL_ERROR "You must set ARCHITECTURE environment variable")
endif()

if(NOT DEFINED ENV{ALIBUILD_WORK_DIR})
    message(FATAL_ERROR "You must set ALIBUILD_WORK_DIR environment variable")
endif()

if(NOT DEFINED ENV{ROOT_VERSION_UBUNTU})
    message(FATAL_ERROR "You must set ROOT_VERSION_UBUNTU environment variable")
endif()

if(NOT DEFINED ENV{ROOT_VERSION})
    message(FATAL_ERROR "You must set ROOT_VERSION environment variable")
endif()

if(NOT DEFINED ENV{BUILD_LINK})
    message(FATAL_ERROR "You must set BUILD_LINK environment variable")
endif()

set(ALI5_VERSION_UBUNTU $ENV{ENVIRONMENT})
set(BUILD_LINK $ENV{BUILD_LINK})
set(ROOT_VERSION $ENV{ROOT_VERSION})
set(ROOT_VERSION_UBUNTU $ENV{ROOT_VERSION_UBUNTU})
set(BASE_DIR $ENV{ALIBUILD_WORK_DIR}/..)
set(SW_UBUNTU $ENV{ALIBUILD_WORK_DIR}/$ENV{ARCHITECTURE})
set(SW_SOURCES_ROOT $ENV{ALIBUILD_WORK_DIR}/SOURCES/ROOT/${ROOT_VERSION}/${ROOT_VERSION})
set(ALI5_EVE "AliRoot/EVE") # NOTE: Create a symlink to AliRoot or change this to proper path!
set(ALI5_SW_UBUNTU $ENV{ALIBUILD_WORK_DIR}/$ENV{ARCHITECTURE}/AliRoot/${ALI5_VERSION_UBUNTU})

message("ALI5_VERSION_UBUNTU: " ${ALI5_VERSION_UBUNTU})
message("ARCHITECTURE: " $ENV{ARCHITECTURE})
message("BASE_DIR: " ${BASE_DIR})
message("ROOT_VERSION: " ${ROOT_VERSION})
message("ROOT_VERSION_UBUNTU: " ${ROOT_VERSION_UBUNTU})
message("SW_UBUNTU: " ${SW_UBUNTU})
message("ALI5_EVE: " ${ALI5_EVE})
message("ALI5_SW_UBUNTU: " ${ALI5_SW_UBUNTU})

message("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
list(APPEND CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
message("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS "-Wall -Wextra -O0")

include_directories(${ALI5_SW_UBUNTU}/include)
include_directories(${SW_UBUNTU}/ROOT/${ROOT_VERSION_UBUNTU}/include)

link_directories(${ALI5_SW_UBUNTU}/lib)
link_directories(${SW_UBUNTU}/ROOT/${ROOT_VERSION_UBUNTU}/lib)

add_executable(ESDConverter extractO2.cxx )
target_compile_features(ESDConverter PRIVATE cxx_std_11)
target_link_libraries(ESDConverter
        STEERBase STEER ESD CDB AOD RAWDatabase RAWDatarec ANALYSIS ANALYSISalice ITSbase ITSrec MUONrec TPCbase TPCrec TPCsim AliHLTUtil EG Eve Ged Gpad Graf Graf3d Gui Minuit Net RGL Rint HLTbase ASImage Core MathCore Matrix RIO Tree
        )
