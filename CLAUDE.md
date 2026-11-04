# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

leeHom trims adapters and merges overlapping paired-end reads using a Bayesian reconstruction (aimed at ancient DNA / short fragments). Input is BAM or FASTQ; output is the same format. See README.md for output files, log fields, priors and `FF` tags.

## Build and test

- `make` (top level) builds htslib (static, from the `htslib/` submodule, only zlib needed), libgab, then `src/leeHom`. Requires `autoconf` for the first htslib build. `make clean` only cleans `src/`.
- `htslib/` and `libgab/` are git submodules (`git clone --recursive`). Override with `HTSLIBINC`/`HTSLIBLIB` (system htslib), `LIBGABINC`, `LIBGABLIB`, `CFLAGS`, `LDFLAGS`, `CXX`. bamtools is no longer used.
- `make test` runs `test/test.sh`: md5 regression checks, a fuzz test of the SIMD kernels against the original scalar code (`test/fuzz_kernels.cpp`, built by `make -C src fuzz_kernels`), and a check that output is identical and in input order for `-t 1` vs `-t 4`. Expected md5s are in `test/*.md5sum_`.
- Debug env vars: `LEEHOM_REFERENCE=1` (original scalar likelihood code), `LEEHOM_NOAVX2=1` (portable SIMD path), `LEEHOM_LEGACY_RNG=1` (global `rand()` stream, only comparable to the pre-optimization binary at `-t 1`), `LEEHOM_STATS=1` (CPU seconds per stage).

## Architecture

- `src/MergeTrimReads.{h,cpp}`: the algorithm (format independent): chimera/adapter detection, likelihood over candidate fragment lengths, merge/trim decision, consensus. `process_PE`/`process_SR` work on strings; `processPair`/`processSingle` on `BamRecord`s. The scalar originals are kept as `*Reference`.
- `src/Kernels.{h,cpp}`: vectorized `adapterScan`, `overlapScan`, `chimeraScan` (AVX2 with runtime dispatch, portable fallback). They must stay **bit-identical** to the scalar code: each SIMD lane runs the same sequential floating-point additions as the original, inactive lanes add +0.0. Do not reassociate sums. Use the fuzz test after any change.
- `src/BamRecord.{h,cpp}`: plain record converted from/to htslib `bam1_t`; tag semantics mirror what bamtools did.
- `src/Pipeline.h`: `OrderedPipeline`, reader -> workers -> writer with condition variables; results are written in the order chunks were read, whatever the thread count.
- `src/leeHom.cpp`: CLI, chunking, BAM (htslib thread pool) and FASTQ (workers gzip their own chunk as an independent gzip member) run functions. Fastq input with `-o`/`-u` goes through the BAM run function: the fastq records are turned into unaligned `bam1_t` (`fastqToBamRecord`) and the random context of a pair is its pair index, as in fastq output. Each worker has its own `MergeTrimReads` copy; tallies are summed at the end. Random tie-breaks are seeded per cluster index (`setRandomContext`), so results do not depend on `-t` or chunk size.
