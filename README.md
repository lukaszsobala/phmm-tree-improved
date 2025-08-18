# phmm-tree (improved)

Original software: phmm-tree (Yin et al.)
- Citation: Yin Y, Mao X, Yang J, Chen X, Mao F, Xu Y. Phylogenetic tree-based HMM search for homologs of carbohydrate-active enzymes in metagenomes. Bioinformatics. 2016;33(3):453–460. DOI: https://doi.org/10.1093/bioinformatics/btw779
- This repository integrates the original phmm-tree workflow with modern build tooling and concurrency orchestration to improve usability and performance.

Primary documentation of changes
- See `IMPLEMENTATION_SUMMARY.md` for the complete rationale and technical details of the changes in this fork.

Concise changelog (highlights)
- Concurrency and orchestration
  - Added `-phylo_concurrent_threads` to run PHYLIP-based tree builders in parallel as separate processes; default auto-detect (0) or 1.
  - Moved PHYLIP algorithms (Fitch, Kitsch, Neighbor-Joining, UPGMA) to process-isolated workers to avoid global-state conflicts; tried parallelization with in-process OpenMP there but these algorithms are sequential in nature.
  - Added `-prc_threads` for PRC distance calculations; clarified auto-detect behavior.
- PRC backend
  - Integrated prcX as an alternative backend with `-prc_backend auto|legacy|prcx`; auto-selects when available and suitable, and can be forced even for HMMER2.
  - Suppressed legacy HMMER version error when prcX is explicitly chosen; print detected HMM profile format.
- Output and UX polish
  - Standardized output names to `<algo>_<variant>_{report,tree.nwk}`; ensured `.nwk` extension.
  - Prefixed output folder with `prcx_` when using prcX; streamlined logs; full usage/help shown when args are missing or invalid.
- Build and toolchain modernization
  - Updated Makefile to use the C++ toolchain for linking, separate `CFLAGS/CXXFLAGS`, and detect OpenMP via C++ compiler; set C++17 (`gnu++17`).
  - Replaced deprecated timing APIs (`ftime/timeb`) with `gettimeofday` helper; removed legacy headers; fixed minor GCC 13 warnings and unsafe patterns.
- Minor fixes
  - Safer file operations and error messages; standardized string find checks; removed unused variables; cleaned nested comments that caused warnings.

Versioning
- These changes do not alter the core algorithms’ logic or outputs. The upstream PHYLIP and phmm-tree algorithmic behavior is preserved; therefore, this does not constitute a new upstream phmm-tree release.

Authorship and permissions
- Integration and modernization changes (c) 2025, Łukasz F. Sobala (lukasz.sobala@hirszfeld.pl).
- Performed with the permission of the original phmm-tree creator: Yanbin Yin (yanbin.yin@gmail.com).
- Original PHYLIP components retain their licenses and attributions; see source headers.

Build
- Requires a C/C++ toolchain (GCC 13+ recommended) and OpenMP (optional for PRC).
- From the repo root:
  - `make`

Usage (example)
- HMMER3 dataset with auto backend and concurrent trees (auto-detect threads):
  - `./phmm-tree -prc -hmms ./path/to/hmms`

Notes
- For full details on options and behavior, see `IMPLEMENTATION_SUMMARY.md`.
