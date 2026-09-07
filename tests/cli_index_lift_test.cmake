if(NOT DEFINED LEVIOSAM2 OR NOT DEFINED TEST_TMP_DIR)
    message(FATAL_ERROR "LEVIOSAM2 and TEST_TMP_DIR are required")
endif()

file(REMOVE_RECURSE "${TEST_TMP_DIR}")
file(MAKE_DIRECTORY "${TEST_TMP_DIR}")

set(CHAIN "${TEST_TMP_DIR}/tiny.chain")
set(FAI "${TEST_TMP_DIR}/target.fai")
set(SAM "${TEST_TMP_DIR}/input.sam")
set(INDEX_PREFIX "${TEST_TMP_DIR}/tiny")

file(WRITE "${CHAIN}"
    "chain 1 src 1000 + 100 200 dst 1000 + 300 400 1\n"
    "50 0 0\n"
    "50\n\n"
    "chain 1 revsrc 1000 + 100 200 revdst 1000 - 600 700 2\n"
    "50 0 0\n"
    "50\n")
file(WRITE "${FAI}" "dst\t1000\nrevdst\t1000\n")
file(WRITE "${SAM}"
    "@HD\tVN:1.6\tSO:unsorted\n"
    "@SQ\tSN:src\tLN:1000\n"
    "@SQ\tSN:revsrc\tLN:1000\n"
    "forward\t0\tsrc\t111\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n"
    "forward-cross\t0\tsrc\t146\t60\t20M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\n"
    "reverse\t0\trevsrc\t111\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n"
    "reverse-cross\t0\trevsrc\t146\t60\t20M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\n")

execute_process(
    COMMAND "${LEVIOSAM2}" index -c "${CHAIN}" -F "${FAI}" -p "${INDEX_PREFIX}"
    RESULT_VARIABLE INDEX_RESULT
    OUTPUT_VARIABLE INDEX_STDOUT
    ERROR_VARIABLE INDEX_STDERR)
if(NOT INDEX_RESULT EQUAL 0 OR NOT EXISTS "${INDEX_PREFIX}.clft")
    message(FATAL_ERROR "index failed (${INDEX_RESULT}): ${INDEX_STDERR}")
endif()

set(LIFT_PREFIX "${TEST_TMP_DIR}/lifted")
execute_process(
    COMMAND "${LEVIOSAM2}" lift -C "${INDEX_PREFIX}.clft" -a "${SAM}"
            -p "${LIFT_PREFIX}" -O sam
    RESULT_VARIABLE LIFT_RESULT
    OUTPUT_VARIABLE LIFT_STDOUT
    ERROR_VARIABLE LIFT_STDERR)
if(NOT LIFT_RESULT EQUAL 0)
    message(FATAL_ERROR "indexed lift failed (${LIFT_RESULT}): ${LIFT_STDERR}")
endif()
file(READ "${LIFT_PREFIX}.sam" LIFTED)

foreach(EXPECTED
        "forward\t0\tdst\t311\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tLO:Z:L"
        "forward-cross\t0\tdst\t346\t60\t20M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\tLO:Z:L"
        "reverse\t16\trevdst\t381\t60\t10M\t*\t0\t0\tTTTTTTTTTT\tIIIIIIIIII\tLO:Z:L"
        "reverse-cross\t16\trevdst\t336\t60\t20M\t*\t0\t0\tTTTTTTTTTTTTTTTTTTTT\tIIIIIIIIIIIIIIIIIIII\tLO:Z:L")
    string(FIND "${LIFTED}" "${EXPECTED}" MATCH_POS)
    if(MATCH_POS EQUAL -1)
        message(FATAL_ERROR "missing expected lifted record: ${EXPECTED}\n${LIFTED}")
    endif()
endforeach()

# Exercise direct chain construction and guard against -B falling through to -c.
set(DIRECT_PREFIX "${TEST_TMP_DIR}/direct")
execute_process(
    COMMAND "${LEVIOSAM2}" lift -c "${CHAIN}" -B 0.5 -F "${FAI}"
            -a "${SAM}" -p "${DIRECT_PREFIX}" -O sam
    RESULT_VARIABLE DIRECT_RESULT
    OUTPUT_VARIABLE DIRECT_STDOUT
    ERROR_VARIABLE DIRECT_STDERR)
if(NOT DIRECT_RESULT EQUAL 0)
    message(FATAL_ERROR "direct lift failed (${DIRECT_RESULT}): ${DIRECT_STDERR}")
endif()
file(READ "${DIRECT_PREFIX}.sam" DIRECT_LIFTED)
string(FIND "${DIRECT_LIFTED}" "forward\t0\tdst\t311\t60\t10M" DIRECT_POS)
if(DIRECT_POS EQUAL -1)
    message(FATAL_ERROR "direct lift produced an unexpected record: ${DIRECT_LIFTED}")
endif()

set(MISSING_PREFIX "${TEST_TMP_DIR}/missing")
execute_process(
    COMMAND "${LEVIOSAM2}" index -c "${TEST_TMP_DIR}/missing.chain"
            -F "${FAI}" -p "${MISSING_PREFIX}"
    RESULT_VARIABLE MISSING_RESULT
    OUTPUT_QUIET ERROR_QUIET)
if(MISSING_RESULT EQUAL 0 OR EXISTS "${MISSING_PREFIX}.clft")
    message(FATAL_ERROR "missing chain input was accepted")
endif()

execute_process(
    COMMAND "${LEVIOSAM2}" index -c "${CHAIN}"
            -F "${TEST_TMP_DIR}/missing.fai" -p "${TEST_TMP_DIR}/missing-fai"
    RESULT_VARIABLE MISSING_FAI_RESULT
    OUTPUT_QUIET ERROR_QUIET)
if(MISSING_FAI_RESULT EQUAL 0 OR EXISTS "${TEST_TMP_DIR}/missing-fai.clft")
    message(FATAL_ERROR "missing FAI input was accepted")
endif()

file(WRITE "${TEST_TMP_DIR}/malformed.chain" "chain malformed\n")
execute_process(
    COMMAND "${LEVIOSAM2}" index -c "${TEST_TMP_DIR}/malformed.chain"
            -F "${FAI}" -p "${TEST_TMP_DIR}/malformed"
    RESULT_VARIABLE MALFORMED_RESULT
    OUTPUT_QUIET ERROR_QUIET)
if(MALFORMED_RESULT EQUAL 0)
    message(FATAL_ERROR "malformed chain input was accepted")
endif()

execute_process(
    COMMAND "${LEVIOSAM2}" unknown-command
    RESULT_VARIABLE UNKNOWN_RESULT
    OUTPUT_QUIET ERROR_QUIET)
if(UNKNOWN_RESULT EQUAL 0)
    message(FATAL_ERROR "unknown command returned success")
endif()

execute_process(
    COMMAND "${LEVIOSAM2}" --version
    RESULT_VARIABLE VERSION_RESULT
    OUTPUT_VARIABLE VERSION_STDOUT
    ERROR_VARIABLE VERSION_STDERR)
if(NOT VERSION_RESULT EQUAL 0 OR VERSION_STDOUT STREQUAL "")
    message(FATAL_ERROR "global --version failed (${VERSION_RESULT}): ${VERSION_STDERR}")
endif()

file(WRITE "${TEST_TMP_DIR}/empty.clft" "")
foreach(MISSING_INDEX "${TEST_TMP_DIR}/missing.clft" "${TEST_TMP_DIR}/empty.clft")
    execute_process(
        COMMAND "${LEVIOSAM2}" lift -C "${MISSING_INDEX}" -a "${SAM}"
                -p "${TEST_TMP_DIR}/bad-index"
        RESULT_VARIABLE BAD_INDEX_RESULT
        OUTPUT_QUIET ERROR_QUIET)
    if(BAD_INDEX_RESULT EQUAL 0)
        message(FATAL_ERROR "missing/empty ChainMap index was accepted: ${MISSING_INDEX}")
    endif()
endforeach()

execute_process(
    COMMAND "${LEVIOSAM2}" lift -C "${INDEX_PREFIX}.clft"
            -a "${TEST_TMP_DIR}/missing.sam" -p "${TEST_TMP_DIR}/bad-input"
            --hts_threads 1
    RESULT_VARIABLE BAD_INPUT_RESULT
    OUTPUT_QUIET ERROR_QUIET)
if(BAD_INPUT_RESULT EQUAL 0)
    message(FATAL_ERROR "missing alignment input was accepted")
endif()

foreach(INVALID_ARGS "-T;0" "-t;not-a-number" "-G;-1")
    execute_process(
        COMMAND "${LEVIOSAM2}" lift -C "${INDEX_PREFIX}.clft" -a "${SAM}"
                -p "${TEST_TMP_DIR}/invalid" ${INVALID_ARGS}
        RESULT_VARIABLE INVALID_RESULT
        OUTPUT_QUIET ERROR_QUIET)
    if(INVALID_RESULT EQUAL 0)
        message(FATAL_ERROR "invalid arguments were accepted: ${INVALID_ARGS}")
    endif()
endforeach()
