#include <unordered_map>
#include <algorithm>
#include <math.h>
#include <omp.h>
#include "helpers.hpp"

// ===============================================================
// Sequential version (reference baseline)
// ===============================================================
unsigned long SequenceInfo::gpsa_sequential(float** S) {
    unsigned long visited = 0;

    // Boundary
    for (unsigned int i = 1; i < rows; i++) {
        S[i][0] = i * gap_penalty;
        visited++;
    }

    for (unsigned int j = 0; j < cols; j++) {
        S[0][j] = j * gap_penalty;
        visited++;
    }

    // Main part
    for (unsigned int i = 1; i < rows; i++) {
        for (unsigned int j = 1; j < cols; j++) {
            float match = S[i - 1][j - 1] + (X[i - 1] == Y[j - 1] ? match_score : mismatch_score);
            float del = S[i - 1][j] + gap_penalty;
            float insert = S[i][j - 1] + gap_penalty;
            S[i][j] = std::max({match, del, insert});
            visited++;
        }
    }

    return visited;
}

unsigned long SequenceInfo::gpsa_taskloop(float** S, long grain_size, int block_size_x, int block_size_y) {
    unsigned long visited = 0;

    // === 1) Initialize matrix borders (sequential part) ===
    for (unsigned int i = 1; i < rows; ++i) {
        S[i][0] = i * gap_penalty;
        ++visited;
    }
    for (unsigned int j = 0; j < cols; ++j) {
        S[0][j] = j * gap_penalty;
        ++visited;
    }

    // === 2) Determine block sizes (from arguments or grain_size) ===
    int bx = block_size_x;
    int by = block_size_y;

    // If user didn’t specify block size, derive it from grain_size
    if (bx == 1 && by == 1 && grain_size > 1) {
        int side = (int)std::sqrt((double)grain_size);
        bx = std::max(8, std::min(side, (int)rows));
        by = std::max(8, std::min(side, (int)cols));
    }

    // Safety limits (avoid crash for edge cases)
    bx = std::max(8, std::min(bx, (int)rows));
    by = std::max(8, std::min(by, (int)cols));

    // === 3) Compute number of blocks ===
    const int NBX = (rows + bx - 1) / bx;
    const int NBY = (cols + by - 1) / by;

    // === 4) Anti-diagonal parallel traversal using taskloop ===
    #pragma omp parallel
    {
        #pragma omp single
        {
            // Iterate over all anti-diagonals of blocks
            for (int diag = 0; diag < NBX + NBY - 1; ++diag) {
                const int start_bi = std::max(0, diag - (NBY - 1));
                const int end_bi   = std::min(diag, NBX - 1);

                // One taskloop per anti-diagonal
                #pragma omp taskloop firstprivate(diag, start_bi, end_bi) \
                                     shared(S, X, Y, rows, cols, match_score, mismatch_score, gap_penalty, bx, by, NBY) \
                                     reduction(+:visited) grainsize(1)
                for (int bi = start_bi; bi <= end_bi; ++bi) {
                    const int bj = diag - bi;
                    if (bj < 0 || bj >= NBY) continue;

                    const int i0 = bi * bx + 1;
                    const int i1 = std::min((int)rows, i0 + bx);
                    const int j0 = bj * by + 1;
                    const int j1 = std::min((int)cols, j0 + by);

                    // Compute block (sequential within task)
                    for (int i = i0; i < i1; ++i) {
                        for (int j = j0; j < j1; ++j) {
                            float match_val = S[i-1][j-1] + (X[i-1]==Y[j-1]?match_score:mismatch_score);
                            float del = S[i-1][j] + gap_penalty;
                            float ins = S[i][j-1] + gap_penalty;
                            S[i][j] = std::max({match_val, del, ins});
                            ++visited;
                        }
                    }
                }
                // Wait for all tasks in this diagonal to complete
                #pragma omp taskwait
            } // end for diag
        } // end single
    } // end parallel

    return visited;
}

unsigned long SequenceInfo::gpsa_tasks(float** S, long grain_size, int block_size_x, int block_size_y) {
    unsigned long visited = 0;

    // === 1) Inicijalizacija granica (uvek sekvencijalno, definisano ponašanje) ===
    for (unsigned int i = 1; i < rows; ++i) { 
        S[i][0] = i * gap_penalty; 
        ++visited; 
    }
    for (unsigned int j = 0; j < cols; ++j) { 
        S[0][j] = j * gap_penalty; 
        ++visited; 
    }

    // === 2) Odredi veličinu blokova na osnovu argumenata ===
    int bx = block_size_x;
    int by = block_size_y;

    // Ako korisnik nije uneo block size, koristi grain_size kao kvadratni blok
    if (bx == 1 && by == 1 && grain_size > 1) {
        int side = (int)std::sqrt((double)grain_size);
        bx = std::max(8, std::min(side, (int)rows));
        by = std::max(8, std::min(side, (int)cols));
    }

    // Sigurnosna ograničenja (da ne crashuje)
    bx = std::max(8, std::min(bx, (int)rows));
    by = std::max(8, std::min(by, (int)cols));

    // === 3) Broj blokova po dimenziji ===
    const int NBX = (rows + bx - 1) / bx;
    const int NBY = (cols + by - 1) / by;

    // === 4) Tokeni za depend klauze ===
    std::unique_ptr<unsigned char[]> tokens(new unsigned char[NBX * NBY]());
    auto token = [&](int bi, int bj) -> unsigned char& {
        return tokens[(size_t)bi * NBY + bj];
    };
    std::cout <<  "Nbx: " << NBX << " Nby: " << NBY << std::endl;
     std::cout <<  "bx: " << bx << " by: " << by << std::endl;

    // === 4) Paralelizacija po anti-dijagonalama ===
    #pragma omp parallel 
    #pragma omp single
    {
        #pragma omp taskgroup task_reduction(+:visited)
        {
            for (int diag = 0; diag < NBX + NBY - 1; ++diag) {
                const int start_bi = std::max(0, diag - (NBY - 1));
                const int end_bi   = std::min(diag, NBX - 1);

                for (int bi = start_bi; bi <= end_bi; ++bi) {
                    const int bj = diag - bi;
                    if (bj < 0 || bj >= NBY) continue;

                    const int i0 = bi * bx + 1;
                    const int i1 = std::min((int)rows, i0 + bx);
                    const int j0 = bj * by + 1;
                    const int j1 = std::min((int)cols, j0 + by);

                    // ===== Task zavisnosti po blokovima =====
                    if (bi == 0 && bj == 0) {
                        #pragma omp task firstprivate(i0,i1,j0,j1) \
                                         shared(S,X,Y,rows,cols,match_score,mismatch_score,gap_penalty) \
                                         depend(inout: token(0,0)) \
                                         in_reduction(+:visited)
                        {
                            for (int i = i0; i < i1; ++i)
                                for (int j = j0; j < j1; ++j) {
                                    float m   = S[i-1][j-1] + (X[i-1]==Y[j-1]?match_score:mismatch_score);
                                    float del = S[i-1][j]   + gap_penalty;
                                    float ins = S[i][j-1]   + gap_penalty;
                                    S[i][j] = std::max({m, del, ins});
                                    ++visited;
                                }
                        }
                    } 
                    else if (bi == 0) {
                        #pragma omp task firstprivate(i0,i1,j0,j1,bj) \
                                         shared(S,X,Y,rows,cols,match_score,mismatch_score,gap_penalty) \
                                         depend(in: token(0,bj-1)) depend(inout: token(0,bj)) \
                                         in_reduction(+:visited)
                        {
                            for (int i = i0; i < i1; ++i)
                                for (int j = j0; j < j1; ++j) {
                                    float m   = S[i-1][j-1] + (X[i-1]==Y[j-1]?match_score:mismatch_score);
                                    float del = S[i-1][j]   + gap_penalty;
                                    float ins = S[i][j-1]   + gap_penalty;
                                    S[i][j] = std::max({m, del, ins});
                                    ++visited;
                                }
                        }
                    } 
                    else if (bj == 0) {
                        #pragma omp task firstprivate(i0,i1,j0,j1,bi) \
                                         shared(S,X,Y,rows,cols,match_score,mismatch_score,gap_penalty) \
                                         depend(in: token(bi-1,0)) depend(inout: token(bi,0)) \
                                         in_reduction(+:visited)
                        {
                            for (int i = i0; i < i1; ++i)
                                for (int j = j0; j < j1; ++j) {
                                    float m   = S[i-1][j-1] + (X[i-1]==Y[j-1]?match_score:mismatch_score);
                                    float del = S[i-1][j]   + gap_penalty;
                                    float ins = S[i][j-1]   + gap_penalty;
                                    S[i][j] = std::max({m, del, ins});
                                    ++visited;
                                }
                        }
                    } 
                    else {
                        #pragma omp task firstprivate(i0,i1,j0,j1,bi,bj) \
                                         shared(S,X,Y,rows,cols,match_score,mismatch_score,gap_penalty) \
                                         depend(in: token(bi-1,bj), token(bi,bj-1)) \
                                         depend(inout: token(bi,bj)) \
                                         in_reduction(+:visited)
                        {
                            for (int i = i0; i < i1; ++i)
                                for (int j = j0; j < j1; ++j) {
                                    float m   = S[i-1][j-1] + (X[i-1]==Y[j-1]?match_score:mismatch_score);
                                    float del = S[i-1][j]   + gap_penalty;
                                    float ins = S[i][j-1]   + gap_penalty;
                                    S[i][j] = std::max({m, del, ins});
                                    ++visited;
                                }
                        }
                    }
                }
            }
        } // taskgroup
    } // parallel/single

    return visited;
}
