##
## Author: Hank Anderson <hank@statease.com>
## Description: Ported from portion of OpenBLAS/Makefile.system
##              Sets C related variables.

if (${CMAKE_C_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_C_COMPILER_ID} STREQUAL "LSB" OR ${CMAKE_C_COMPILER_ID} MATCHES "Clang")

  set(CCOMMON_OPT "${CCOMMON_OPT} -Wall")
  set(COMMON_PROF "${COMMON_PROF} -fno-inline")
  set(NO_UNINITIALIZED_WARN "-Wno-uninitialized")

  if (QUIET_MAKE)
    set(CCOMMON_OPT "${CCOMMON_OPT} ${NO_UNINITIALIZED_WARN} -Wno-unused")
  endif ()

  if (NO_BINARY_MODE)

    if (MIPS32)
        set(CCOMMON_OPT "${CCOMMON_OPT} -mabi=32")
      set(BINARY_DEFINED 1)
    endif ()

    if (MIPS64)
      if (BINARY64)
        set(CCOMMON_OPT "${CCOMMON_OPT} -mabi=64")
      else ()
        set(CCOMMON_OPT "${CCOMMON_OPT} -mabi=n32")
      endif ()
      set(BINARY_DEFINED 1)
    endif ()

    if (${CORE} STREQUAL "LOONGSON3A" OR ${CORE} STREQUAL "LOONGSON3B")
      set(CCOMMON_OPT "${CCOMMON_OPT} -march=mips64")
      set(FCOMMON_OPT "${FCOMMON_OPT} -march=mips64")
    endif ()

    if (LOONGARCH64)
      if (BINARY64)
        set(CCOMMON_OPT "${CCOMMON_OPT} -mabi=lp64")
      else ()
        set(CCOMMON_OPT "${CCOMMON_OPT} -mabi=lp32")
      endif ()
      set(BINARY_DEFINED 1)
    endif ()

    if (CMAKE_SYSTEM_NAME STREQUAL "AIX")
      set(BINARY_DEFINED 1)
    endif ()
  endif ()

  if (NOT BINARY_DEFINED)
    if (BINARY64)
      set(CCOMMON_OPT "${CCOMMON_OPT} -m64")
    else ()
      set(CCOMMON_OPT "${CCOMMON_OPT} -m32")
    endif ()
  endif ()
endif ()

if (${CMAKE_C_COMPILER_ID} STREQUAL "PGI")
  if (BINARY64)
    set(CCOMMON_OPT "${CCOMMON_OPT} -tp p7-64")
  else ()
    set(CCOMMON_OPT "${CCOMMON_OPT} -tp p7")
  endif ()
endif ()

if (${CMAKE_C_COMPILER_ID} STREQUAL "PATHSCALE")
  if (BINARY64)
    set(CCOMMON_OPT "${CCOMMON_OPT} -m64")
  else ()
    set(CCOMMON_OPT "${CCOMMON_OPT} -m32")
  endif ()
endif ()

if (${CMAKE_C_COMPILER_ID} STREQUAL "OPEN64")

  if (MIPS64)

    if (NOT BINARY64)
      set(CCOMMON_OPT "${CCOMMON_OPT} -n32")
    else ()
      set(CCOMMON_OPT "${CCOMMON_OPT} -n64")
    endif ()

    if (${CORE} STREQUAL "LOONGSON3A")
      set(CCOMMON_OPT "${CCOMMON_OPT} -loongson3 -static")
    endif ()

    if (${CORE} STREQUAL "LOONGSON3B")
      set(CCOMMON_OPT "${CCOMMON_OPT} -loongson3 -static")
    endif ()

  else ()

    if (BINARY64)
      set(CCOMMON_OPT "${CCOMMON_OPT} -m32")
    else ()
      set(CCOMMON_OPT "${CCOMMON_OPT} -m64")
    endif ()
  endif ()
endif ()

if (${CMAKE_C_COMPILER_ID} STREQUAL "SUN")
  set(CCOMMON_OPT "${CCOMMON_OPT} -w")
  if (X86)
    set(CCOMMON_OPT "${CCOMMON_OPT} -m32")
  else ()
    set(CCOMMON_OPT "${CCOMMON_OPT} -m64")
  endif ()
endif ()

if (${CORE} STREQUAL SKYLAKEX)
  if (NOT DYNAMIC_ARCH)
    if (NOT NO_AVX512)
      set (CCOMMON_OPT "${CCOMMON_OPT} -march=skylake-avx512")
    endif ()
  endif ()
endif ()

if (${CORE} STREQUAL COOPERLAKE)
  if (NOT DYNAMIC_ARCH)
    if (NOT NO_AVX512)
      execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
      if (${GCC_VERSION} VERSION_GREATER 10.1 OR ${GCC_VERSION} VERSION_EQUAL 10.1)
        set (CCOMMON_OPT  "${CCOMMON_OPT} -march=cooperlake")
      else ()
        set (CCOMMON_OPT "${CCOMMON_OPT} -march=skylake-avx512")
      endif()  
    endif ()
  endif ()
endif ()

if (${CORE} STREQUAL SAPPHIRERAPIDS)
  if (NOT DYNAMIC_ARCH)
    if (NOT NO_AVX512)
      execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
      if (${GCC_VERSION} VERSION_GREATER 11.0 OR ${GCC_VERSION} VERSION_EQUAL 11.0)
        set (CCOMMON_OPT  "${CCOMMON_OPT} -march=sapphirerapids")
      else ()
        set (CCOMMON_OPT "${CCOMMON_OPT} -march=skylake-avx512")
      endif()  
    endif ()
  endif ()
endif ()

if (${CORE} STREQUAL A64FX)
  if (NOT DYNAMIC_ARCH)
    execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (${GCC_VERSION} VERSION_GREATER 11.0 OR ${GCC_VERSION} VERSION_EQUAL 11.0)
      set (CCOMMON_OPT  "${CCOMMON_OPT} -march=armv8.2-a+sve -mtune=a64fx")
    else ()
      set (CCOMMON_OPT "${CCOMMON_OPT} -march=armv8.2-a+sve")
    endif()
  endif ()
endif ()

if (${CORE} STREQUAL ARMV8SVE)
  if (NOT DYNAMIC_ARCH)
    set (CCOMMON_OPT "${CCOMMON_OPT} -march=armv8-a+sve")
  endif ()
endif ()

if (${CORE} STREQUAL POWER10)
  if (NOT DYNAMIC_ARCH)
    execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (${GCC_VERSION} VERSION_GREATER 10.2 OR ${GCC_VERSION} VERSION_EQUAL 10.2)
      set (CCOMMON_OPT  "${CCOMMON_OPT} -mcpu=power10 -mtune=power10 -mvsx -fno-fast-math")
    else ()
      message(FATAL_ERROR "Compiler GCC.${GCC_VERSION} does not support Power10." )
    endif()
  endif ()
endif ()

if (${CORE} STREQUAL POWER9)
  if (NOT DYNAMIC_ARCH)
    execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (${GCC_VERSION} VERSION_GREATER 5.0 OR ${GCC_VERSION} VERSION_EQUAL 5.0)
      set (CCOMMON_OPT  "${CCOMMON_OPT} -mcpu=power9 -mtune=power9 -mvsx -fno-fast-math")
    else ()
      set (CCOMMON_OPT  "${CCOMMON_OPT} -mcpu=power8 -mtune=power8 -mvsx -fno-fast-math")
      message(WARNING "Compiler GCC.${GCC_VERSION} does not fully support Power9.")
    endif ()
  endif ()
endif ()

if (${CORE} STREQUAL POWER8)
  if (NOT DYNAMIC_ARCH)
    set (CCOMMON_OPT  "${CCOMMON_OPT} -mcpu=power8 -mtune=power8 -mvsx -fno-fast-math")
  endif ()
endif ()

if (NOT DYNAMIC_ARCH)
	if (HAVE_AVX2)
        set (CCOMMON_OPT  "${CCOMMON_OPT} -mavx2")
	endif ()
	if (HAVE_AVX)
        set (CCOMMON_OPT  "${CCOMMON_OPT} -mavx")
	endif ()
	#	if (HAVE_FMA3)
	#set (CCOMMON_OPT  "${CCOMMON_OPT} -mfma")
	#endif ()
	if (HAVE_SSE)
	set (CCOMMON_OPT  "${CCOMMON_OPT} -msse")
	endif ()
	if (HAVE_SSE2)
        set (CCOMMON_OPT  "${CCOMMON_OPT} -msse2")
	endif ()
	if (HAVE_SSE3)
        set (CCOMMON_OPT  "${CCOMMON_OPT} -msse3")
	endif ()
	if (HAVE_SSSE3)
        set (CCOMMON_OPT  "${CCOMMON_OPT} -mssse3")
	endif ()
	if (HAVE_SSE4_1)
	set (CCOMMON_OPT  "${CCOMMON_OPT} -msse4.1")
	endif ()
endif()