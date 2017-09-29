function(boost_python_module NAME)
  find_package(PythonLibs 2 REQUIRED)
  set(py_version ${PYTHONLIBS_VERSION_STRING})
  STRING( REGEX REPLACE "[^0-9]" "" boost_py_version ${py_version})
  find_package(Boost QUIET COMPONENTS "python-py${boost_py_version}")
  set(Boost_PYTHON_FOUND ${Boost_PYTHON-PY${boost_py_version}_FOUND})
  set(Boost_PYTHON_LIBRARY ${Boost_PYTHON-PY${boost_py_version}_LIBRARY})
  while(NOT "${py_version}" STREQUAL "" AND NOT Boost_PYTHON_FOUND)
    STRING( REGEX REPLACE "([0-9.]+).[0-9]+" "\\1" py_version ${py_version})
    STRING( REGEX REPLACE "[^0-9]" "" boost_py_version ${py_version})
    find_package(Boost QUIET COMPONENTS "python-py${boost_py_version}")
    set(Boost_PYTHON_FOUND ${Boost_PYTHON-PY${boost_py_version}_FOUND})
    set(Boost_PYTHON_LIBRARY ${Boost_PYTHON-PY${boost_py_version}_LIBRARY})
    STRING( REGEX MATCHALL "([0-9.]+).[0-9]+" has_more_versions ${py_version})
    if("${has_more_versions}" STREQUAL "")
      break()
    endif()
  endwhile()
  if(NOT Boost_PYTHON_FOUND)
    find_package(Boost QUIET COMPONENTS python3)
    set(Boost_PYTHON_FOUND ${Boost_PYTHON3_FOUND})
    set(Boost_PYTHON_LIBRARY ${Boost_PYTHON3_LIBRARY})
  endif()
  if(NOT Boost_PYTHON_FOUND)
    find_package(Boost COMPONENTS python REQUIRED)
  endif()
  find_package(Numpy)



  set(DEP_LIBS
    ${Boost_PYTHON_LIBRARY}
    ${PYTHON_LIBRARIES}
    )
  #these are required includes for every ecto module
  include_directories(
    ${PYTHON_INCLUDE_PATH}
    ${Boost_INCLUDE_DIRS}
    )
  add_library(${NAME} SHARED
    ${ARGN}
    )
  set_target_properties(${NAME}
    PROPERTIES
    OUTPUT_NAME ${NAME}
    COMPILE_FLAGS "${FASTIDIOUS_FLAGS}"
    LINK_FLAGS -dynamic
    PREFIX ""
  )
  if( WIN32 )
    set_target_properties(${NAME} PROPERTIES SUFFIX ".pyd")
  elseif( APPLE OR ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # on mac osx, python cannot import libraries with .dylib extension
    set_target_properties(${NAME} PROPERTIES SUFFIX ".so")
  endif()  
  target_link_libraries(${NAME}
    ${DEP_LIBS}
    )
endfunction()
