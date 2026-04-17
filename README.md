Parallel Global Pair-wise Sequence Alignment (GPSA)
This project implements and analyzes parallel versions of the Global Pair-wise Sequence Alignment algorithm using OpenMP. It was developed as part of the Parallel Architectures and Programming Models course at the University of Vienna.

The core objective was to optimize the computation of large-scale sequence alignment matrices using two different OpenMP tasking approaches and comparing their efficiency against a sequential baseline.

Project Overview
The project parallelizes the dynamic programming approach to sequence alignment (Needleman-Wunsch style). The computation is performed over large matrices, such as 51480 x 53640, requiring careful management of task granularity and dependencies.

Implementation Methods
Sequential: Baseline implementation for performance comparison.

OpenMP Taskloop: Parallelization using the taskloop construct, processing the matrix along anti-diagonals with synchronization barriers.

OpenMP Explicit Tasks: A more advanced version using fine-grained task constructs with depend clauses to allow a wavefront execution pattern.

Key Performance Results
Based on testing on the Alma cluster, the following results were achieved:

Best Speedup (Taskloop): 12.09x

Best Speedup (Explicit Tasks): 13.14x

Optimal Configuration: 256 x 256 block size (42,420 total tasks).

The explicit tasks version consistently outperformed the taskloop version by reducing global synchronization and allowing tasks to start as soon as their specific dependencies were met.

File Structure
main.cpp: Entry point, argument parsing, and timing logic.

implementation.hpp: Core parallel logic for Taskloop and Explicit Task versions.

helpers.hpp: Memory management (contiguous 2D allocation) and sequence loading.

OpenMP_Taskloop_vs_Explicit_Detailed_Report.pdf: A comprehensive technical analysis of speedup, cache effects, and granularity.

a1a.pdf: Project requirements and algorithmic description.

How to Build and Run
Compilation
Use g++ with the OpenMP flag and C++20 standard:

Bash
g++ -O2 -std=c++20 -fopenmp -o gpsa main.cpp
Execution
Run the executable without arguments to use default settings, or specify parameters:

Bash
./gpsa --x X.txt --y Y.txt --exec-mode 0 --block-size-x 256 --block-size-y 256
Arguments:

--x, --y: Input sequence files.

--exec-mode: 0 (All), 1 (Sequential), 2 (Taskloop), 3 (Explicit Tasks).

--block-size-x, --block-size-y: Adjust task granularity.
