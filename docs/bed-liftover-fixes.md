# BED interval liftover fixes

BED coordinates are zero-based, half-open boundaries. The previous mapper
mixed base coordinates with boundaries, which shifted forward endpoints and
could emit reversed bounds. It also retained source-coordinate BED6 strand and
BED12 coding/exon fields after lifting the outer interval.

The BED mapper now projects the interval through the actual overlapping chain
blocks, transforms both boundaries, and requires compatible target contig,
strand and order. Reverse mappings flip BED6 strand. BED12 coding bounds and
exon blocks move with the outer interval, are reordered on reversal, and are
validated against it. BED3-6 and BED12 with extra columns are supported;
BED7-11 are explicitly rejected instead of emitting stale coordinate fields.
An empty coding range is anchored at the new start. Zero-length insertion
points map only strictly inside an aligned block.

The inclusive `-G` limit bounds each endpoint clip, each source/target gap, and
the absolute net length change over the retained source span. Thus `-G 0`
accepts gapless mappings, a 500-base gap passes at 500, and a 501-base gap does
not. Invalid records are written to the unmapped BED with a reason; invalid
options, missing inputs/indexes, and I/O failures cause a nonzero exit.
Comments, blank lines, and track/browser directives are skipped.

Shared coordinate helpers now reject negative clipping-error sentinels before
arithmetic. No alignment reverse-position formula or CLFT serialization format
is changed. The convenience CIGAR guard depends on the first-interval lookup
fix in [PR #118](https://github.com/milkschen/leviosam2/pull/118), so this PR is
stacked on that branch.

## Compatibility and mapping policy

Existing CLFT indexes do not retain original chain IDs. Compatible, monotonic
blocks may therefore be joined within `-G` even when they originated in separate
chains. Incompatible mappings are rejected. This policy differs from UCSC's
single-chain selection, and the absolute gap tolerance differs from UCSC
liftOver's mapped-fraction threshold. Cross-gap agreement alone is not a shared
mapping truth. Enforcing original chain identity would require a separate
index-format design.

## Validation

The isolated PR checkout passes all **10 CTest targets**, including the new BED
CLI regression, and all **5 existing pysam/Picard alignment integration tests**.
The combined working-tree validation also passed 12 CTest targets and the same
5 alignment integration tests; its two extra targets concern separate mate
metadata work and are not included in this PR.

The scoped CLI regression covers exact forward/reverse single bases and longer
intervals, both block boundaries, clipped endpoints, gap thresholds of 0/500/501,
adjacent gapless blocks, inconsistent intermediate mappings, BED12 transforms,
invalid records/options, missing files, metadata lines, help and verbosity.
Parsed and serialized indexes are tested for negative clipping sentinels.

Fresh isolated-binary experiment reruns pass all **18 hand-defined nonempty
synthetic cases** in BED3, BED6 and BED12, with separate strand and coding/block
checks. On a deterministic broader real-chain sample, **1,250/1,250 within-block
intervals** (625 forward and 625 reverse) match both the independent block oracle
and UCSC liftOver exactly. Baseline exact counts were 505/625 forward and 0/625
reverse. These are selected unambiguous single-block strata, not genome-wide
accuracy estimates. Raw experiment outputs are intentionally excluded from this
PR.

Reproduce the checked-in regression and existing alignment checks:

```sh
pixi run test
pixi run -e integration test-integration
```
