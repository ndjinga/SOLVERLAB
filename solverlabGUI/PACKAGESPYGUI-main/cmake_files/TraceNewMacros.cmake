
# utilities trace and color and log


# function to print all current variables contents
function(printCmakeTrace)
  unset(_variableNames  PARENT_SCOPE)
  get_cmake_property(_variableNames VARIABLES)
  list(REMOVE_DUPLICATES _variableNames)
  message( STATUS "${ColBlue}BEGIN TRACE" )
  foreach(_variableName ${_variableNames})
    # avoid multi-lines codes as too much
    if ( "${${_variableName}}" MATCHES ".#?define ." ) # too long anyway
      continue()
    endif()
    if ( "${${_variableName}}" MATCHES ".#?include ." ) # too long anyway
      continue()
    endif()
    if ( "${${_variableName}}" MATCHES ".DOXYFILE_ENCODING." ) # illisible anyway
      continue()
    endif()
    if ( _variableName MATCHES ".?_REGEX$" )  # illisible anyway
      continue()
    endif()
    if ( _variableName MATCHES "Col." )  # avoid setcolors below
      continue()
    endif()
    message("TRACE:   ${_variableName}=${${_variableName}}")
  endforeach()
  message( STATUS "END TRACE${ColReset}" )
endfunction()


function(printCmakeTraceVariable regExp)
  unset(_variableNames  PARENT_SCOPE)
  get_cmake_property(_variableNames VARIABLES)
  list(REMOVE_DUPLICATES _variableNames)
  message( STATUS "${ColBlue}" )
  foreach(_variableName ${_variableNames})
    if ( _variableName MATCHES "${regExp}" )
      message("TRACE:   ${_variableName}=${${_variableName}}")
    endif()
  endforeach()
  message( STATUS "${ColReset}" )
endfunction()


function(setcolors)
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColReset       "${Esc}[m"     PARENT_SCOPE)
  set(ColBold        "${Esc}[1m"    PARENT_SCOPE)
  set(ColRed         "${Esc}[31m"   PARENT_SCOPE)
  set(ColGreen       "${Esc}[32m"   PARENT_SCOPE)
  set(ColYellow      "${Esc}[33m"   PARENT_SCOPE)
  set(ColBlue        "${Esc}[34m"   PARENT_SCOPE)
  set(ColMagenta     "${Esc}[35m"   PARENT_SCOPE)
  set(ColCyan        "${Esc}[36m"   PARENT_SCOPE)
  set(ColWhite       "${Esc}[37m"   PARENT_SCOPE)
  set(ColBoldRed     "${Esc}[1;31m" PARENT_SCOPE)
  set(ColBoldGreen   "${Esc}[1;32m" PARENT_SCOPE)
  set(ColBoldYellow  "${Esc}[1;33m" PARENT_SCOPE)
  set(ColBoldBlue    "${Esc}[1;34m" PARENT_SCOPE)
  set(ColBoldMagenta "${Esc}[1;35m" PARENT_SCOPE)
  set(ColBoldCyan    "${Esc}[1;36m" PARENT_SCOPE)
  set(ColBoldWhite   "${Esc}[1;37m" PARENT_SCOPE)
endif()
endfunction()

function(unsetcolors)
  unset(ColReset       PARENT_SCOPE)
  unset(ColBold        PARENT_SCOPE)
  unset(ColRed         PARENT_SCOPE)
  unset(ColGreen       PARENT_SCOPE)
  unset(ColYellow      PARENT_SCOPE)
  unset(ColBlue        PARENT_SCOPE)
  unset(ColMagenta     PARENT_SCOPE)
  unset(ColCyan        PARENT_SCOPE)
  unset(ColWhite       PARENT_SCOPE)
  unset(ColBoldRed     PARENT_SCOPE)
  unset(ColBoldGreen   PARENT_SCOPE)
  unset(ColBoldYellow  PARENT_SCOPE)
  unset(ColBoldBlue    PARENT_SCOPE)
  unset(ColBoldMagenta PARENT_SCOPE)
  unset(ColBoldCyan    PARENT_SCOPE)
  unset(ColBoldWhite   PARENT_SCOPE)
endfunction()


function(test_colormessage)
  message("This is normal")
  message("${ColRed}This is Red${ColReset}")
  message("${ColGreen}This is Green${ColReset}")
  message("${ColYellow}This is Yellow${ColReset}")
  message("${ColBlue}This is Blue${ColReset}")
  message("${ColMagenta}This is Magenta${ColReset}")
  message("${ColCyan}This is Cyan${ColReset}")
  message("${ColWhite}This is White${ColReset}")
  message("${ColBoldRed}This is BoldRed${ColReset}")
  message("${ColBoldGreen}This is BoldGreen${ColReset}")
  message("${ColBoldYellow}This is BoldYellow${ColReset}")
  message("${ColBoldBlue}This is BoldBlue${ColReset}")
  message("${ColBoldMagenta}This is BoldMagenta${ColReset}")
  message("${ColBoldCyan}This is BoldCyan${ColReset}")
  message("${ColBoldWhite}This is BoldWhite\n\n${ColReset}")
endfunction()


function(log_fatal msg)
  message( FATAL_ERROR  "${ColBoldRed}FATAL:  ${msg}${ColReset}" )
endfunction()


function(log_warning msg)
  message( "${ColRed}WARNING: ${msg}${ColReset}" )
endfunction()


function(log_info msg)
  message( "${ColGreen}INFO:    ${msg}${ColReset}" )
endfunction()

function(log_title msg)
  set(stars "*****************************************************")
  message( "${ColBoldGreen}${stars}\n  ${msg}\n${stars}${ColReset}" )
endfunction()
