# Load appropriate system files into the configuration file

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  MESSAGE(STATUS "This is a Darwin system")
  # 1. My Mac laptop 
  IF( $ENV{MACHINE_NAME} MATCHES "poincare")
    MESSAGE(STATUS "The machine is my Mac laptop banach")
    INCLUDE(${CMAKE_CURRENT_LIST_DIR}/poincare_PETSc_VTK.cmake)
  ENDIF( $ENV{MACHINE_NAME} MATCHES "poincare")
ELSEIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  MESSAGE(STATUS "This is a Linux system ")
  MESSAGE(STATUS $ENV{MACHINE_NAME})
  IF( $ENV{MACHINE_NAME} MATCHES "linux-hilbert")
    MESSAGE(STATUS "The Machine is my home Linux desktop hilbert")
    INCLUDE(${CMAKE_CURRENT_LIST_DIR}/hilbert_PETSc_VTK.cmake)
  ELSEIF( $ENV{MACHINE_NAME} MATCHES "STAMPEDE")
    MESSAGE(STATUS "The Machine is TACC Stampede")
    INCLUDE(${CMAKE_CURRENT_LIST_DIR}/stampede_PETSc_v36p0_shared_VTK.cmake)
  ELSEIF( $ENV{MACHINE_NAME} MATCHES "sp2")
    MESSAGE(STATUS "The Machine is TACC Stampede2")
    INCLUDE(${CMAKE_CURRENT_LIST_DIR}/stampede2_PETSc_v36p0_shared_VTK.cmake)
  ELSEIF( $ENV{MACHINE_NAME} MATCHES "LS5")
    MESSAGE(STATUS "The Machine is TACC Lonestar 5")
    INCLUDE(${CMAKE_CURRENT_LIST_DIR}/ls5_PETSc_v36p0_shared_VTK.cmake)
  ELSEIF( $ENV{MACHINE_NAME} MATCHES "chern")
    MESSAGE(STATUS "The Machine is Stanford desktop chern")
    INCLUDE(${CMAKE_CURRENT_LIST_DIR}/stanford_bacon_v2018.cmake)
  ELSEIF( $ENV{MACHINE_NAME} MATCHES "lax")
    MESSAGE(STATUS "The Machine is my desktop Lax")
    INCLUDE(${CMAKE_CURRENT_LIST_DIR}/lax_PETSc_VTK.cmake)
  ELSEIF( $ENV{MACHINE_NAME} MATCHES "comet")
    MESSAGE(STATUS "The Machine is COMET@SDSC")
    INCLUDE(${CMAKE_CURRENT_LIST_DIR}/comet_PETSc_v375_shared_VTK.cmake)
  ELSEIF( $ENV{MACHINE_NAME} MATCHES "sherlock")
    MESSAGE(STATUS "The Machine is sherlock@Stanford")
    INCLUDE(${CMAKE_CURRENT_LIST_DIR}/stanford_sherlock.cmake)
  ELSEIF( $ENV{MACHINE_NAME} MATCHES "ingridxlan")
    MESSAGE(STATUS "The Machine is Stanford ingridxlan")
    INCLUDE(${CMAKE_CURRENT_LIST_DIR}/marsden_linux_ingridxlan.cmake)
  ELSEIF( $ENV{MACHINE_NAME} MATCHES "oguz")
    MESSAGE(STATUS "The Machine is Oguz's personal ubuntu laptop")
    INCLUDE(${CMAKE_CURRENT_LIST_DIR}/oguz-ubuntu.cmake)
  ELSE($ENV{MACHINE_NAME} MATCHES "linux-hilbert")
    MESSAGE(STATUS "This Machine cannot be identified.")
  ENDIF( $ENV{MACHINE_NAME} MATCHES "linux-hilbert")
ELSE($ENV{CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  MESSAGE(STATUS "This system cannot be identified.")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# End of the file
