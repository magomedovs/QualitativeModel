##
#  CMake script for the TravelingWave program:
##

# Set the name of the project and target:
SET(TARGET "calculate")

SET(TARGET_SRC
  ${TARGET}.cc SolutionStruct.cc WaveConstructor.cc Parameters.cc LimitSolution.cc InitialGuess.cc 
  )

# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 3.13.4)

FIND_PACKAGE(deal.II 9.5.0
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})

DEAL_II_INVOKE_AUTOPILOT()

#  FIND_PACKAGE(Eigen3 3.3 REQUIRED NO_MODULE)
#  INCLUDE_DIRECTORIES(/usr/include/eigen3/)
#
#  ADD_EXECUTABLE(${TARGET} ${TARGET_SRC})
#  DEAL_II_SETUP_TARGET(${TARGET})
#
#  TARGET_LINK_LIBRARIES(${TARGET} Eigen3::Eigen)

