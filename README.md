# phmm-tree (improved)

Original software: phmm-tree (Yin et al.)
- Citation: Yin Y, Mao X, Yang J, Chen X, Mao F, Xu Y. Phylogenetic tree-based HMM search for homologs of carbohydrate-active enzymes in metagenomes. Bioinformatics. 2016;33(3):453–460. DOI: https://doi.org/10.1093/bioinformatics/btw779
- This repository integrates the original phmm-tree workflow with modern build tooling and concurrency orchestration to improve usability and performance.

Primary documentation of changes
- See `IMPLEMENTATION_SUMMARY.md` for the complete rationale and technical details of the changes in this fork.

Concise changelog:
- Concurrency, control and orchestration
  - Added options to perform only some phylogenetic analyses: `-fitch` `-kitsch`, `-upgma`, `-nj` and to limit Fitch and Kitsch to `-fm` or `-min` only
  - Multi-threaded pairwise PRC analysis; added `-prc_threads` for PRC distance calculations; auto-detect max threads.
  - Added `-phylo_concurrent_threads` to run PHYLIP-based tree builders (Fitch, Kitsch, Neighbor-Joining, UPGMA) in parallel as separate processes; default auto-detect (0) or 1.
  - Moved PHYLIP algorithms to process-isolated workers to avoid global-state conflicts; tried parallelization with in-process OpenMP there but these algorithms are sequential in nature. `phylo_threads` exist exists to control this but it does not seem to offer any speedups.
- PRC backend
  - Integrated prcX as an alternative backend with `-prc_backend auto|legacy|prcx`; auto-selects when available and suitable, and can be forced even for HMMER2.
  - Suppressed legacy HMMER version error when prcX is explicitly chosen; print detected HMM profile format.
- Output and UX polish
  - Standardized output names to `<algo>_<variant>_{report,tree.nwk}`; ensured `.nwk` extension.
  - Prefixed output folder with `prcx_` when using prcX; streamlined logs; full usage/help shown when args are missing or invalid.
  - Help now appears when phmm-tree is run without any parameters or fails.
- Bugfixes:
  - Fixed a bug where tree tips with the same name would be generated from similarly named HMMs, due to an insufficient name collision detection algorithm.
- Minor fixes
  - Safer file operations and error messages; standardized string find checks; removed unused variables; cleaned nested comments that caused warnings.
  - Updated Makefile to use the C++ toolchain for linking, separate `CFLAGS/CXXFLAGS`, and detect OpenMP via C++ compiler; set C++17 (`gnu++17`).

Versioning
- These changes do not alter the core algorithms’ logic or outputs. The upstream PHYLIP and phmm-tree algorithmic behavior is preserved; therefore, this does not constitute a new upstream phmm-tree release.

Authorship and permissions
- Integration and modernization changes (c) 2025, Łukasz F. Sobala (lukasz.sobala@hirszfeld.pl). These changes were made with very little programming knowledge, mostly with the help of coding agents.
- Performed with the permission of the original phmm-tree creator: Yanbin Yin (yanbin.yin@gmail.com).
- Original PHYLIP components retain their licenses and attributions; see source headers.

Build
- Requires a C/C++ toolchain (GCC 13+ recommended) and OpenMP (optional for PRC).
- From the repo root:
  - `make -j`

Usage (example)
- HMMER3 dataset with auto backend and concurrent trees (auto-detect threads):
  - `./phmm-tree -prc -hmms ./path/to/hmms`
