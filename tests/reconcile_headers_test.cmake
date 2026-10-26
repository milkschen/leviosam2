if(NOT DEFINED LEVIOSAM2 OR NOT DEFINED TEST_TMP_DIR)
    message(FATAL_ERROR "LEVIOSAM2 and TEST_TMP_DIR are required")
endif()
file(REMOVE_RECURSE "${TEST_TMP_DIR}")
file(MAKE_DIRECTORY "${TEST_TMP_DIR}")
set(HEADER "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200\n")
set(RECORD "read1\t0\tchr2\t10\t60\t1M\t*\t0\t0\tA\tI\tAS:i:100\n")
file(WRITE "${TEST_TMP_DIR}/first.sam" "${HEADER}read1\t0\tchr1\t1\t10\t1M\t*\t0\t0\tA\tI\tAS:i:1\n")
file(WRITE "${TEST_TMP_DIR}/same.sam" "${HEADER}@PG\tID:second\tPN:other-aligner\n${RECORD}")
execute_process(COMMAND "${LEVIOSAM2}" reconcile
    -s "first:${TEST_TMP_DIR}/first.sam" -s "second:${TEST_TMP_DIR}/same.sam"
    -o "${TEST_TMP_DIR}/valid.sam"
    RESULT_VARIABLE RESULT ERROR_VARIABLE STDERR)
if(NOT RESULT EQUAL 0)
    message(FATAL_ERROR "compatible inputs failed (${RESULT}): ${STDERR}")
endif()
file(READ "${TEST_TMP_DIR}/valid.sam" OUTPUT)
string(FIND "${OUTPUT}" "read1\t0\tchr2\t10\t60\t1M" POS)
if(POS EQUAL -1)
    message(FATAL_ERROR "selected record lost its reference: ${OUTPUT}")
endif()
file(WRITE "${TEST_TMP_DIR}/reordered.sam" "@SQ\tSN:chr2\tLN:200\n@SQ\tSN:chr1\tLN:100\n${RECORD}")
file(WRITE "${TEST_TMP_DIR}/renamed.sam" "@SQ\tSN:other\tLN:100\n@SQ\tSN:chr2\tLN:200\n${RECORD}")
file(WRITE "${TEST_TMP_DIR}/length.sam" "@SQ\tSN:chr1\tLN:101\n@SQ\tSN:chr2\tLN:200\n${RECORD}")
file(WRITE "${TEST_TMP_DIR}/count.sam" "@SQ\tSN:chr2\tLN:200\n${RECORD}")
file(WRITE "${TEST_TMP_DIR}/invalid.sam" "invalid header\n")
foreach(CASE reordered renamed length count missing invalid)
    execute_process(COMMAND "${LEVIOSAM2}" reconcile
        -s "first:${TEST_TMP_DIR}/first.sam" -s "second:${TEST_TMP_DIR}/${CASE}.sam"
        -o "${TEST_TMP_DIR}/${CASE}-output.sam"
        RESULT_VARIABLE RESULT ERROR_VARIABLE STDERR)
    if(NOT RESULT EQUAL 1 OR EXISTS "${TEST_TMP_DIR}/${CASE}-output.sam")
        message(FATAL_ERROR "${CASE} input did not fail cleanly before output (${RESULT}): ${STDERR}")
    endif()
    if(CASE MATCHES "^(reordered|renamed|length|count)$" AND
       NOT STDERR MATCHES "Incompatible reference dictionary")
        message(FATAL_ERROR "${CASE} input lacked a useful diagnostic: ${STDERR}")
    endif()
endforeach()
execute_process(COMMAND "${LEVIOSAM2}" reconcile
    -s "first:${TEST_TMP_DIR}/first.sam" -s "second:${TEST_TMP_DIR}/same.sam"
    -o "${TEST_TMP_DIR}/missing-directory/output.sam"
    RESULT_VARIABLE RESULT ERROR_VARIABLE STDERR)
if(NOT RESULT EQUAL 1 OR NOT STDERR MATCHES "Failed to open output file")
    message(FATAL_ERROR "invalid output did not fail cleanly (${RESULT}): ${STDERR}")
endif()
