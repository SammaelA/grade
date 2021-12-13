# Install script for directory: /media/sammael/Data/programming/bachelor_work/grade/dependencies/Hydra/HydraCore/hydra_app

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/hydra/hydra" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/hydra/hydra")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/hydra/hydra"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/hydra" TYPE EXECUTABLE FILES "/media/sammael/Data/programming/bachelor_work/grade/dependencies/Hydra/HydraCore/cmake-build-debug/hydra_app/hydra")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/hydra/hydra" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/hydra/hydra")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/hydra/hydra"
         OLD_RPATH "/media/sammael/Data/programming/bachelor_work/grade/dependencies/Hydra/HydraCore/../HydraAPI/bin:/media/sammael/Data/programming/bachelor_work/grade/dependencies/Hydra/HydraCore/LIBRARY/lib_x64_linux:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/hydra/hydra")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/hydra/shaders" TYPE FILE FILES
    "/media/sammael/Data/programming/bachelor_work/grade/dependencies/Hydra/HydraCore/hydra_drv/cfetch.h"
    "/media/sammael/Data/programming/bachelor_work/grade/dependencies/Hydra/HydraCore/hydra_drv/cglobals.h"
    "/media/sammael/Data/programming/bachelor_work/grade/dependencies/Hydra/HydraCore/hydra_drv/shaders/texproc.cl"
    "/media/sammael/Data/programming/bachelor_work/grade/dependencies/Hydra/HydraCore/hydra_drv/shaders/image.xx"
    "/media/sammael/Data/programming/bachelor_work/grade/dependencies/Hydra/HydraCore/hydra_drv/shaders/light.xx"
    "/media/sammael/Data/programming/bachelor_work/grade/dependencies/Hydra/HydraCore/hydra_drv/shaders/material.xx"
    "/media/sammael/Data/programming/bachelor_work/grade/dependencies/Hydra/HydraCore/hydra_drv/shaders/mlt.xx"
    "/media/sammael/Data/programming/bachelor_work/grade/dependencies/Hydra/HydraCore/hydra_drv/shaders/screen.xx"
    "/media/sammael/Data/programming/bachelor_work/grade/dependencies/Hydra/HydraCore/hydra_drv/shaders/sort.xx"
    "/media/sammael/Data/programming/bachelor_work/grade/dependencies/Hydra/HydraCore/hydra_drv/shaders/trace.xx"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS ${ENV}${CMAKE_INSTALL_PREFIX}/hydra/shadercache)
                FILE(REMOVE_RECURSE ${ENV}${CMAKE_INSTALL_PREFIX}/hydra/shadercache)
              ENDIF()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  FILE(MAKE_DIRECTORY ${ENV}${CMAKE_INSTALL_PREFIX}/hydra/shadercache)
endif()

