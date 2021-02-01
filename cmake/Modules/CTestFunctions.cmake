FUNCTION (MANGLE_FORTRAN_NAME CNAME FNAME)
    SET (TMP)
    IF (WIN32)
        STRING (TOUPPER "${FNAME}" TMP)
    ELSE ()
        STRING (TOLOWER "${FNAME}_" TMP)
    ENDIF ()
    SET (${CNAME} ${TMP} PARENT_SCOPE)
ENDFUNCTION ()


FUNCTION (MANGLE_FORTRAN_FILENAME_LIST MANGLED)
    SET (TMP)
    FOREACH (TFILE ${ARGN})
        string (REGEX REPLACE ".f90$" "" TESTNAME ${TFILE})
        MANGLE_FORTRAN_NAME (C_TESTNAME ${TESTNAME})
        list (APPEND TMP ${C_TESTNAME})
    ENDFOREACH ()
    SET (${MANGLED} ${TMP} PARENT_SCOPE)
ENDFUNCTION()


FUNCTION (ADD_FORTRAN_TEST_EXECUTABLE TARGET)
    SET (TEST_FILES ${ARGN})
    MANGLE_FORTRAN_FILENAME_LIST (TEST_FILES_MANGLED ${TEST_FILES})

    create_test_sourcelist (_ main.c ${TEST_FILES_MANGLED})

    add_library (${TARGET}_fortran ${TEST_FILES} $<TARGET_OBJECTS:lib>)
    add_executable (${TARGET} main.c)
    target_link_libraries (${TARGET} ${TARGET}_fortran ${MPI_Fortran_LIBRARIES})

    SET (INDEX 0)
    list (LENGTH TEST_FILES LEN)
    WHILE (${LEN} GREATER ${INDEX})
        list (GET TEST_FILES ${INDEX} TEST_FILE)
        list (GET TEST_FILES_MANGLED ${INDEX} TEST_FILE_MANGLED)
        add_test (
            NAME ${TEST_FILE}
            COMMAND $<TARGET_FILE:${TARGET}>  ${TEST_FILE_MANGLED})
        math (EXPR INDEX "${INDEX} + 1")
    ENDWHILE ()
ENDFUNCTION ()


