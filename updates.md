1. Parallelization

    Now `liftover` can be run in parallel, with little overhead in both CPU time and peak memory usage.
    Two arguments (`--threads`/`-t` and `--chunk_size`/`-T`) are added for `liftover lift`.

2. SAM format fixes

    - Fix invalid CIGAR format such as `1I1D1I`, which results from an indel between the REF and the ALT genomes overlaps with an indel between the read and the reference genome (ALT).

    - We are looking into other SAM formatting issues, if any, and will test the lifted SAM with Picard and GATK.

3. Warning output fix

    Now only shows the `chromosome not found` warning once for each chromosome unmatched with the lift map.
