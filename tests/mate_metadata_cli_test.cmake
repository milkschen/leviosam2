if(NOT DEFINED LEVIOSAM2 OR NOT DEFINED TEST_TMP_DIR)
    message(FATAL_ERROR "LEVIOSAM2 and TEST_TMP_DIR are required")
endif()
file(REMOVE_RECURSE "${TEST_TMP_DIR}")
file(MAKE_DIRECTORY "${TEST_TMP_DIR}")
set(CHAIN "${TEST_TMP_DIR}/tiny.chain")
set(FAI "${TEST_TMP_DIR}/target.fai")
file(WRITE "${CHAIN}"
    "chain 1 src 1000 + 100 210 dst 1000 + 300 400 1\n50 10 0\n50\n\n"
    "chain 1 revsrc 1000 + 100 200 revdst 1000 - 600 700 2\n100\n")
file(WRITE "${FAI}" "dst\t1000\nrevdst\t1000\n")
set(HEADER "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:src\tLN:1000\n@SQ\tSN:revsrc\tLN:1000\n@RG\tID:a\tSM:sample\tPL:ILLUMINA\n@RG\tID:b\tSM:sample\tPL:ILLUMINA\n")
set(INPUT "${TEST_TMP_DIR}/input.sam")
file(WRITE "${INPUT}" "${HEADER}"
    "gap\t99\tsrc\t111\t60\t10M\t=\t156\t65\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:a\tMC:Z:20M\n"
    "known\t73\tsrc\t111\t60\t10M\t=\t111\t0\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:a\n"
    "known\t133\tsrc\t111\t0\t*\t=\t111\t0\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:a\n"
    "shared\t99\tsrc\t111\t60\t10M\t=\t181\t80\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:a\n"
    "shared\t99\tsrc\t121\t60\t10M\t=\t191\t80\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:b\n"
    "gap\t355\tsrc\t121\t60\t10M\t=\t156\t0\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:a\n"
    "gap\t2147\tsrc\t131\t60\t10M\t=\t156\t0\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:a\n"
    "gap\t147\tsrc\t156\t60\t20M\t=\t111\t-65\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\tRG:Z:a\tMC:Z:10M\n"
    "shared\t147\tsrc\t181\t60\t10M\t=\t111\t-80\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:a\n"
    "shared\t147\tsrc\t191\t60\t10M\t=\t121\t-80\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:b\n"
    "reverse\t99\trevsrc\t111\t60\t30M\t=\t186\t85\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tRG:Z:a\n"
    "reverse\t147\trevsrc\t186\t60\t10M\t=\t111\t-85\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:a\n")

execute_process(COMMAND "${LEVIOSAM2}" index -c "${CHAIN}" -F "${FAI}"
    -p "${TEST_TMP_DIR}/tiny" RESULT_VARIABLE RESULT ERROR_VARIABLE STDERR)
if(NOT RESULT EQUAL 0)
    message(FATAL_ERROR "index failed: ${STDERR}")
endif()

function(check_record PATH NAME RG FLAG REF POS CIGAR MREF MPOS TLEN)
    file(STRINGS "${PATH}" LINES)
    set(MATCH_COUNT 0)
    foreach(LINE IN LISTS LINES)
        if(LINE MATCHES "^${NAME}\t${FLAG}\t" AND LINE MATCHES "\tRG:Z:${RG}(\t|$)")
            math(EXPR MATCH_COUNT "${MATCH_COUNT} + 1")
            string(REPLACE "\t" ";" FIELDS "${LINE}")
            foreach(PAIR "2;${REF}" "3;${POS}" "5;${CIGAR}" "6;${MREF}" "7;${MPOS}" "8;${TLEN}")
                list(GET PAIR 0 INDEX)
                list(GET PAIR 1 EXPECTED)
                list(GET FIELDS ${INDEX} ACTUAL)
                if(NOT ACTUAL STREQUAL EXPECTED)
                    message(FATAL_ERROR "${PATH}: field ${INDEX}: ${ACTUAL} != ${EXPECTED}: ${LINE}")
                endif()
            endforeach()
            if(LINE MATCHES "\tMC:")
                message(FATAL_ERROR "stale mate CIGAR emitted: ${LINE}")
            endif()
        endif()
    endforeach()
    if(NOT MATCH_COUNT EQUAL 1)
        message(FATAL_ERROR "${PATH}: expected one ${NAME}/${RG}/${FLAG}, got ${MATCH_COUNT}")
    endif()
endfunction()

function(check_complete PATH)
    check_record("${PATH}" gap a 99 dst 311 10M = 351 55)
    check_record("${PATH}" gap a 147 dst 351 5S15M = 311 -55)
    check_record("${PATH}" gap a 355 dst 321 10M = 351 0)
    check_record("${PATH}" gap a 2147 dst 331 10M = 351 0)
    check_record("${PATH}" known a 73 dst 311 10M = 311 0)
    check_record("${PATH}" known a 133 dst 311 * = 311 0)
    check_record("${PATH}" shared a 99 dst 311 10M = 371 70)
    check_record("${PATH}" shared a 147 dst 371 10M = 311 -70)
    check_record("${PATH}" shared b 99 dst 321 10M = 381 70)
    check_record("${PATH}" shared b 147 dst 381 10M = 321 -70)
    check_record("${PATH}" reverse a 83 revdst 361 30M = 306 -85)
    check_record("${PATH}" reverse a 163 revdst 306 10M = 361 85)
endfunction()

foreach(MODE direct indexed)
    if(MODE STREQUAL "direct")
        set(MAPPING -c "${CHAIN}" -F "${FAI}")
    else()
        set(MAPPING -C "${TEST_TMP_DIR}/tiny.clft")
    endif()
    foreach(THREADS 1 4)
        set(PREFIX "${TEST_TMP_DIR}/${MODE}-${THREADS}")
        execute_process(COMMAND "${LEVIOSAM2}" lift ${MAPPING} -a "${INPUT}"
            -p "${PREFIX}" -O sam -G 100 -T 1 -t ${THREADS}
            RESULT_VARIABLE RESULT ERROR_VARIABLE STDERR)
        if(NOT RESULT EQUAL 0)
            message(FATAL_ERROR "${MODE}/${THREADS} failed: ${STDERR}")
        endif()
        check_complete("${PREFIX}.sam")
    endforeach()
endforeach()

# Neither input nor output has to be seekable. Every primary pair crosses
# worker chunks, and secondary/supplementary records use the actual primary mate.
execute_process(COMMAND "${LEVIOSAM2}" lift -C "${TEST_TMP_DIR}/tiny.clft"
    -a - -p - -O sam -G 100 -T 1 -t 4 INPUT_FILE "${INPUT}"
    OUTPUT_FILE "${TEST_TMP_DIR}/stdio.sam" RESULT_VARIABLE RESULT ERROR_VARIABLE STDERR)
if(NOT RESULT EQUAL 0)
    message(FATAL_ERROR "stdin/stdout lift failed: ${STDERR}")
endif()
check_complete("${TEST_TMP_DIR}/stdio.sam")

# The corrected 55-base span must be used for the split decision, not the
# inconsistent 75/-70 estimates that existed before final mate reconciliation.
execute_process(COMMAND "${LEVIOSAM2}" lift -C "${TEST_TMP_DIR}/tiny.clft"
    -a "${INPUT}" -p "${TEST_TMP_DIR}/split" -O sam -G 100 -T 1 -t 4 -S isize:60
    RESULT_VARIABLE RESULT ERROR_VARIABLE STDERR)
if(NOT RESULT EQUAL 0)
    message(FATAL_ERROR "split lift failed: ${STDERR}")
endif()
check_record("${TEST_TMP_DIR}/split-committed.sam" gap a 99 dst 311 10M = 351 55)
check_record("${TEST_TMP_DIR}/split-committed.sam" gap a 147 dst 351 5S15M = 311 -55)
check_record("${TEST_TMP_DIR}/split-deferred.sam" reverse a 83 revdst 361 30M = 306 -85)

# Incomplete region-subset BAMs: use MC exactly when possible; without MC,
# preserve source-known mappedness and leave unknown target coordinates unset.
file(WRITE "${TEST_TMP_DIR}/partial-input.sam" "${HEADER}"
    "mc\t99\tsrc\t111\t60\t10M\t=\t156\t65\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:a\tMC:Z:20M\tMQ:i:42\n"
    "unknown\t99\tsrc\t111\t60\t10M\t=\t156\t65\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:a\tMQ:i:42\n"
    "unmapped\t73\tsrc\t111\t60\t10M\t=\t111\t0\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:a\n"
    "self-unmapped\t133\tsrc\t111\t0\t*\t=\t156\t0\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:a\tMC:Z:20M\n")
execute_process(COMMAND "${LEVIOSAM2}" lift -C "${TEST_TMP_DIR}/tiny.clft"
    -a "${TEST_TMP_DIR}/partial-input.sam" -p "${TEST_TMP_DIR}/partial" -O sam -G 100
    RESULT_VARIABLE RESULT ERROR_VARIABLE STDERR)
if(NOT RESULT EQUAL 0 OR NOT STDERR MATCHES "target mate coordinates and TLEN are unknown")
    message(FATAL_ERROR "incomplete-mate handling failed: ${STDERR}")
endif()
check_record("${TEST_TMP_DIR}/partial.sam" mc a 99 dst 311 10M = 351 55)
check_record("${TEST_TMP_DIR}/partial.sam" unknown a 65 dst 311 10M * 0 0)
check_record("${TEST_TMP_DIR}/partial.sam" unmapped a 73 dst 311 10M = 311 0)
check_record("${TEST_TMP_DIR}/partial.sam" self-unmapped a 133 dst 351 * = 351 0)
file(STRINGS "${TEST_TMP_DIR}/partial.sam" PARTIAL_LINES)
foreach(LINE IN LISTS PARTIAL_LINES)
    if(LINE MATCHES "^mc\t" AND NOT LINE MATCHES "\tMQ:i:42(\t|$)")
        message(FATAL_ERROR "virtual MC mate lost known MQ: ${LINE}")
    elseif(LINE MATCHES "^unknown\t" AND LINE MATCHES "\tMQ:")
        message(FATAL_ERROR "unknown target mate retained stale MQ: ${LINE}")
    endif()
endforeach()

file(WRITE "${TEST_TMP_DIR}/duplicate.sam" "${HEADER}"
    "same\t65\tsrc\t111\t60\t10M\t=\t181\t0\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:a\n"
    "same\t65\tsrc\t121\t60\t10M\t=\t191\t0\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:a\n")
execute_process(COMMAND "${LEVIOSAM2}" lift -C "${TEST_TMP_DIR}/tiny.clft"
    -a "${TEST_TMP_DIR}/duplicate.sam" -p "${TEST_TMP_DIR}/duplicate-output" -O sam -G 100
    RESULT_VARIABLE RESULT ERROR_VARIABLE STDERR)
if(RESULT EQUAL 0 OR NOT STDERR MATCHES "ambiguous duplicate primary")
    message(FATAL_ERROR "ambiguous primary identity was not rejected: ${STDERR}")
endif()

file(WRITE "${TEST_TMP_DIR}/bad-mc.sam" "${HEADER}"
    "bad\t65\tsrc\t111\t60\t10M\t=\t181\t0\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:a\tMC:i:10\n")
execute_process(COMMAND "${LEVIOSAM2}" lift -C "${TEST_TMP_DIR}/tiny.clft"
    -a "${TEST_TMP_DIR}/bad-mc.sam" -p "${TEST_TMP_DIR}/bad-mc-output" -O sam -G 100
    RESULT_VARIABLE RESULT ERROR_VARIABLE STDERR)
if(RESULT EQUAL 0 OR NOT STDERR MATCHES "MC tag must be a CIGAR string")
    message(FATAL_ERROR "malformed MC tag was not rejected clearly: ${STDERR}")
endif()
