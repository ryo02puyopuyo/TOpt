/*
	TOpt: An Efficient Quantum Compiler that Reduces the T Count
	Copyright (C) 2018  Luke Heyfron

	This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
// ステップ1: M4RI を他の全て（C++標準ライブラリやプロジェクトヘッダ）よりも先にインクルードする
extern "C" {
    #include <m4ri/m4ri.h>
}
#include <chrono>

// ステップ2: C++標準ライブラリ
#include <iostream>
#include <cmath>
#include <set>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <map>
#include <random>    // random_device, mt19937 用
#include <algorithm> // shuffle 用
#include <iomanip>
// (他にもあればここに追加)

// ステップ3: プロジェクト固有のヘッダ
#include "GateSynthesisMatrix.h"
#include "LCL/LCL_Mat_GF2.h"
#include "LCL/Core/LCL_ConsoleOut.h"
#include "LCL/LCL_Int.h"

// ステップ4: using namespace は全ての #include の「後」に置く
using namespace std;
using namespace LCL_ConsoleOut;

bool** M4RI_direct_nullspace(mzd_t* A_in, int& out_d);

namespace {
struct CleanupTrace {
    std::vector<int> zeroed_cols;
    std::vector<std::pair<int, int>> swaps;
};

inline int aa_table_index(int n, int i, int j) {
    int lo = (i < j) ? i : j;
    int hi = (i < j) ? j : i;
    return lo * n - lo * (lo + 1) / 2 + hi - lo - 1;
}

void build_aa_table(mzd_t* AA_tab, mzd_t* A_m4ri, int n, int active_cols) {
    if (!AA_tab) return;
    mzd_t* AA_win = mzd_init_window(AA_tab, 0, 0, AA_tab->nrows, active_cols);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            int idx = aa_table_index(n, i, j);
            mzd_row_clear_offset(AA_win, idx, 0);
            for (int w = 0; w < AA_win->width; ++w) {
                mzd_row(AA_win, idx)[w] = mzd_row(A_m4ri, i)[w] & mzd_row(A_m4ri, j)[w];
            }
        }
    }
    mzd_free_window(AA_win);
}

void cleanup_with_trace(bool** A, int n, int m, int& mp, CleanupTrace& trace) {
    std::vector<char> was_zeroed(m, 0);

    for (int j1 = 0; j1 < (m - 1); ++j1) {
        for (int j2 = (j1 + 1); j2 < m; ++j2) {
            bool same = 1;
            for (int i = 0; same && (i < n); ++i) {
                same *= (A[i][j1] == A[i][j2]);
            }
            if (same) {
                if (!was_zeroed[j1]) {
                    was_zeroed[j1] = 1;
                    trace.zeroed_cols.push_back(j1);
                }
                if (!was_zeroed[j2]) {
                    was_zeroed[j2] = 1;
                    trace.zeroed_cols.push_back(j2);
                }
                for (int i = 0; i < n; ++i) {
                    A[i][j1] = 0;
                    A[i][j2] = 0;
                }
            }
        }
    }

    int j_end = (m - 1);
    int non_zero_count = 0;
    for (int j = 0; j < m; ++j) {
        int sum = 0;
        for (int i = 0; i < n; ++i) {
            sum += A[i][j];
        }
        if (!sum) {
            bool found = false;
            while ((!found) && (j_end > j)) {
                for (int i = 0; (!found) && (i < n); ++i) {
                    found = A[i][j_end];
                }
                if (!found) j_end--;
            }
            if (found) {
                LCL_Mat_GF2::swapcol(A, n, m, j, j_end);
                trace.swaps.push_back(std::make_pair(j, j_end));
                non_zero_count++;
            }
        } else {
            non_zero_count++;
        }
    }
    mp = non_zero_count;
}

void apply_cleanup_trace_to_aa_table(mzd_t* AA_tab, int aa_rows, const CleanupTrace& trace) {
    if (!AA_tab) return;

    for (size_t idx = 0; idx < trace.zeroed_cols.size(); ++idx) {
        int col = trace.zeroed_cols[idx];
        for (int row = 0; row < aa_rows; ++row) {
            mzd_write_bit(AA_tab, row, col, 0);
        }
    }

    for (size_t idx = 0; idx < trace.swaps.size(); ++idx) {
        mzd_col_swap(AA_tab, trace.swaps[idx].first, trace.swaps[idx].second);
    }
}

void update_aa_table_after_hit(mzd_t* AA_tab, mzd_t* A_m4ri_full, int n, int active_cols,
                               const CleanupTrace& trace, const std::vector<int>& changed_rows) {
    if (!AA_tab) return;

    int aa_rows = n * (n - 1) / 2;
    apply_cleanup_trace_to_aa_table(AA_tab, aa_rows, trace);

    if (changed_rows.empty()) return;

    std::vector<char> touched(n, 0);
    for (size_t idx = 0; idx < changed_rows.size(); ++idx) {
        touched[changed_rows[idx]] = 1;
    }

    mzd_t* A_win = mzd_init_window(A_m4ri_full, 0, 0, n, active_cols);
    mzd_t* AA_win = mzd_init_window(AA_tab, 0, 0, aa_rows, active_cols);

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (!(touched[i] || touched[j])) continue;
            int idx = aa_table_index(n, i, j);
            mzd_row_clear_offset(AA_win, idx, 0);
            for (int w = 0; w < AA_win->width; ++w) {
                mzd_row(AA_win, idx)[w] = mzd_row(A_win, i)[w] & mzd_row(A_win, j)[w];
            }
        }
    }

    mzd_free_window(AA_win);
    mzd_free_window(A_win);
}

struct DynamicBasisRow {
    std::vector<unsigned char> bits;
    int pivot;
    int source_row_id;
};

struct DynamicChiState {
    std::vector<std::vector<unsigned char>> chi_rows;
    std::vector<unsigned char> row_active;
    std::vector<DynamicBasisRow> basis;
    std::vector<int> row_to_basis_slot;
    std::vector<int> pivot_to_basis_slot;
    std::vector<unsigned char> x_prev;
};

struct DynamicCandidate {
    int c1;
    int c2;
    int x_weight;
    std::vector<unsigned char> x_bits;
};

struct PackedDynamicBasisRow {
    std::vector<uint64_t> bits;
    int pivot;
    int source_row_id;
};

struct PackedChiState {
    std::vector<std::vector<uint64_t>> chi_rows;
    std::vector<unsigned char> row_active;
    std::vector<PackedDynamicBasisRow> basis;
    std::vector<int> row_to_basis_slot;
    std::vector<int> pivot_to_basis_slot;
    std::vector<unsigned char> x_prev;
    int words;
};

inline int dynamic_row_id(int alpha, int beta, int gamma, int n) {
    return alpha * n * n + beta * n + gamma;
}

inline unsigned char dynamic_chi_row_active(const std::vector<unsigned char>& x_bits,
                                            int alpha, int beta, int gamma) {
    return (unsigned char)(x_bits[alpha] | x_bits[beta] | x_bits[gamma]);
}

inline int packed_word_count(int m) {
    return (m + 63) / 64;
}

inline int packed_get_bit(const std::vector<uint64_t>& row, int col) {
    return (int)((row[col >> 6] >> (col & 63)) & 1ULL);
}

inline void packed_set_bit(std::vector<uint64_t>& row, int col) {
    row[col >> 6] |= (1ULL << (col & 63));
}

void packed_xor_row(std::vector<uint64_t>& dst, const std::vector<uint64_t>& src) {
    for (int w = 0; w < (int)dst.size(); ++w) dst[w] ^= src[w];
}

bool packed_is_zero_row(const std::vector<uint64_t>& row) {
    for (int w = 0; w < (int)row.size(); ++w) {
        if (row[w]) return false;
    }
    return true;
}

int packed_find_leftmost_one(const std::vector<uint64_t>& row, int m) {
    for (int w = 0; w < (int)row.size(); ++w) {
        uint64_t x = row[w];
        if (x) {
            int bit = __builtin_ctzll(x);
            int col = w * 64 + bit;
            return (col < m) ? col : -1;
        }
    }
    return -1;
}

int dynamic_find_leftmost_one(const std::vector<unsigned char>& row) {
    for (int c = 0; c < (int)row.size(); ++c) {
        if (row[c]) return c;
    }
    return -1;
}

bool dynamic_is_zero_row(const std::vector<unsigned char>& row) {
    for (int c = 0; c < (int)row.size(); ++c) {
        if (row[c]) return false;
    }
    return true;
}

void dynamic_xor_row(std::vector<unsigned char>& dst, const std::vector<unsigned char>& src) {
    for (int c = 0; c < (int)dst.size(); ++c) dst[c] ^= src[c];
}

void build_chi_row_bool(bool** A, const std::vector<unsigned char>& x_bits, int m,
                        int alpha, int beta, int gamma, std::vector<unsigned char>& out_row) {
    unsigned char x_a = x_bits[alpha];
    unsigned char x_b = x_bits[beta];
    unsigned char x_c = x_bits[gamma];
    unsigned char term_const = x_a & x_b & x_c;

    for (int j = 0; j < m; ++j) {
        int val = term_const
            ^ ((x_a & x_b) & A[gamma][j])
            ^ ((x_b & x_c) & A[alpha][j])
            ^ ((x_c & x_a) & A[beta][j])
            ^ (x_a & A[beta][j] & A[gamma][j])
            ^ (x_b & A[gamma][j] & A[alpha][j])
            ^ (x_c & A[alpha][j] & A[beta][j]);
        out_row[j] = (unsigned char)(val & 1);
    }
}

void build_chi_row_packed(bool** A, const std::vector<unsigned char>& x_bits, int m,
                          int alpha, int beta, int gamma, std::vector<uint64_t>& out_row) {
    std::fill(out_row.begin(), out_row.end(), 0ULL);
    unsigned char x_a = x_bits[alpha];
    unsigned char x_b = x_bits[beta];
    unsigned char x_c = x_bits[gamma];
    unsigned char term_const = x_a & x_b & x_c;

    for (int j = 0; j < m; ++j) {
        int val = term_const
            ^ ((x_a & x_b) & A[gamma][j])
            ^ ((x_b & x_c) & A[alpha][j])
            ^ ((x_c & x_a) & A[beta][j])
            ^ (x_a & A[beta][j] & A[gamma][j])
            ^ (x_b & A[gamma][j] & A[alpha][j])
            ^ (x_c & A[alpha][j] & A[beta][j]);
        if (val & 1) packed_set_bit(out_row, j);
    }
}

void dynamic_add_row_to_basis(const std::vector<unsigned char>& row_in, int source_row_id, int m,
                              std::vector<DynamicBasisRow>& basis,
                              std::vector<int>& row_to_basis_slot,
                              std::vector<int>& pivot_to_basis_slot) {
    std::vector<unsigned char> row = row_in;
    for (int bi = 0; bi < (int)basis.size(); ++bi) {
        int pivot = basis[bi].pivot;
        if (pivot >= 0 && row[pivot]) dynamic_xor_row(row, basis[bi].bits);
    }

    int pivot = dynamic_find_leftmost_one(row);
    if (pivot < 0) return;

    for (int bi = 0; bi < (int)basis.size(); ++bi) {
        if (basis[bi].bits[pivot]) dynamic_xor_row(basis[bi].bits, row);
    }

    DynamicBasisRow new_row;
    new_row.bits = row;
    new_row.pivot = pivot;
    new_row.source_row_id = source_row_id;
    basis.push_back(new_row);
    std::sort(basis.begin(), basis.end(), [](const DynamicBasisRow& a, const DynamicBasisRow& b) {
        return a.pivot < b.pivot;
    });

    std::fill(row_to_basis_slot.begin(), row_to_basis_slot.end(), -1);
    std::fill(pivot_to_basis_slot.begin(), pivot_to_basis_slot.end(), -1);
    for (int bi = 0; bi < (int)basis.size(); ++bi) {
        if (basis[bi].source_row_id >= 0 && basis[bi].source_row_id < (int)row_to_basis_slot.size()) {
            row_to_basis_slot[basis[bi].source_row_id] = bi;
        }
        if (basis[bi].pivot >= 0 && basis[bi].pivot < (int)pivot_to_basis_slot.size()) {
            pivot_to_basis_slot[basis[bi].pivot] = bi;
        }
    }
}

void dynamic_rebuild_basis_from_state(DynamicChiState& state, int m) {
    state.basis.clear();
    std::fill(state.row_to_basis_slot.begin(), state.row_to_basis_slot.end(), -1);
    std::fill(state.pivot_to_basis_slot.begin(), state.pivot_to_basis_slot.end(), -1);
    for (int r = 0; r < (int)state.chi_rows.size(); ++r) {
        if (!state.row_active[r]) continue;
        if (dynamic_is_zero_row(state.chi_rows[r])) continue;
        dynamic_add_row_to_basis(state.chi_rows[r], r, m, state.basis, state.row_to_basis_slot, state.pivot_to_basis_slot);
    }
}

void dynamic_reindex_basis(std::vector<DynamicBasisRow>& basis,
                           std::vector<int>& row_to_basis_slot,
                           std::vector<int>& pivot_to_basis_slot) {
    std::sort(basis.begin(), basis.end(), [](const DynamicBasisRow& a, const DynamicBasisRow& b) {
        return a.pivot < b.pivot;
    });

    std::fill(row_to_basis_slot.begin(), row_to_basis_slot.end(), -1);
    std::fill(pivot_to_basis_slot.begin(), pivot_to_basis_slot.end(), -1);
    for (int bi = 0; bi < (int)basis.size(); ++bi) {
        if (basis[bi].source_row_id >= 0 && basis[bi].source_row_id < (int)row_to_basis_slot.size()) {
            row_to_basis_slot[basis[bi].source_row_id] = bi;
        }
        if (basis[bi].pivot >= 0 && basis[bi].pivot < (int)pivot_to_basis_slot.size()) {
            pivot_to_basis_slot[basis[bi].pivot] = bi;
        }
    }
}

bool dynamic_membership_test(const std::vector<unsigned char>& target,
                             const std::vector<DynamicBasisRow>& basis) {
    std::vector<unsigned char> row = target;
    for (int bi = 0; bi < (int)basis.size(); ++bi) {
        int pivot = basis[bi].pivot;
        if (pivot >= 0 && row[pivot]) dynamic_xor_row(row, basis[bi].bits);
    }
    return dynamic_is_zero_row(row);
}

void packed_add_row_to_basis(const std::vector<uint64_t>& row_in, int source_row_id, int m,
                             std::vector<PackedDynamicBasisRow>& basis,
                             std::vector<int>& row_to_basis_slot,
                             std::vector<int>& pivot_to_basis_slot) {
    std::vector<uint64_t> row = row_in;
    for (int bi = 0; bi < (int)basis.size(); ++bi) {
        int pivot = basis[bi].pivot;
        if (pivot >= 0 && packed_get_bit(row, pivot)) packed_xor_row(row, basis[bi].bits);
    }

    int pivot = packed_find_leftmost_one(row, m);
    if (pivot < 0) return;

    for (int bi = 0; bi < (int)basis.size(); ++bi) {
        if (packed_get_bit(basis[bi].bits, pivot)) packed_xor_row(basis[bi].bits, row);
    }

    PackedDynamicBasisRow new_row;
    new_row.bits = row;
    new_row.pivot = pivot;
    new_row.source_row_id = source_row_id;
    basis.push_back(new_row);
    std::sort(basis.begin(), basis.end(), [](const PackedDynamicBasisRow& a, const PackedDynamicBasisRow& b) {
        return a.pivot < b.pivot;
    });

    std::fill(row_to_basis_slot.begin(), row_to_basis_slot.end(), -1);
    std::fill(pivot_to_basis_slot.begin(), pivot_to_basis_slot.end(), -1);
    for (int bi = 0; bi < (int)basis.size(); ++bi) {
        if (basis[bi].source_row_id >= 0 && basis[bi].source_row_id < (int)row_to_basis_slot.size()) {
            row_to_basis_slot[basis[bi].source_row_id] = bi;
        }
        if (basis[bi].pivot >= 0 && basis[bi].pivot < (int)pivot_to_basis_slot.size()) {
            pivot_to_basis_slot[basis[bi].pivot] = bi;
        }
    }
}

void packed_rebuild_basis_from_state(PackedChiState& state, int m) {
    state.basis.clear();
    std::fill(state.row_to_basis_slot.begin(), state.row_to_basis_slot.end(), -1);
    std::fill(state.pivot_to_basis_slot.begin(), state.pivot_to_basis_slot.end(), -1);
    for (int r = 0; r < (int)state.chi_rows.size(); ++r) {
        if (!state.row_active[r]) continue;
        if (packed_is_zero_row(state.chi_rows[r])) continue;
        packed_add_row_to_basis(state.chi_rows[r], r, m, state.basis, state.row_to_basis_slot, state.pivot_to_basis_slot);
    }
}

bool packed_membership_test_pair(int c1, int c2, int m, const std::vector<PackedDynamicBasisRow>& basis) {
    std::vector<uint64_t> row(packed_word_count(m), 0ULL);
    packed_set_bit(row, c1);
    packed_set_bit(row, c2);
    for (int bi = 0; bi < (int)basis.size(); ++bi) {
        int pivot = basis[bi].pivot;
        if (pivot >= 0 && packed_get_bit(row, pivot)) packed_xor_row(row, basis[bi].bits);
    }
    return packed_is_zero_row(row);
}

void packed_rebuild_basis_from_state_m4ri(PackedChiState& state, int m) {
    std::vector<int> active_rows;
    active_rows.reserve(state.chi_rows.size());
    for (int r = 0; r < (int)state.chi_rows.size(); ++r) {
        if (!state.row_active[r]) continue;
        if (packed_is_zero_row(state.chi_rows[r])) continue;
        active_rows.push_back(r);
    }

    state.basis.clear();
    std::fill(state.row_to_basis_slot.begin(), state.row_to_basis_slot.end(), -1);
    std::fill(state.pivot_to_basis_slot.begin(), state.pivot_to_basis_slot.end(), -1);
    if (active_rows.empty()) return;

    mzd_t* M = mzd_init((rci_t)active_rows.size(), m);
    for (int rr = 0; rr < (int)active_rows.size(); ++rr) {
        int src_r = active_rows[rr];
        for (int c = 0; c < m; ++c) {
            if (packed_get_bit(state.chi_rows[src_r], c)) mzd_write_bit(M, rr, c, 1);
        }
    }

    rci_t rank = mzd_echelonize(M, 1);
    for (rci_t rr = 0; rr < rank; ++rr) {
        std::vector<uint64_t> row(state.words, 0ULL);
        for (int c = 0; c < m; ++c) {
            if (mzd_read_bit(M, rr, c)) packed_set_bit(row, c);
        }
        int pivot = packed_find_leftmost_one(row, m);
        if (pivot < 0) continue;

        PackedDynamicBasisRow br;
        br.bits = row;
        br.pivot = pivot;
        br.source_row_id = active_rows[rr];
        state.basis.push_back(br);
    }

    std::sort(state.basis.begin(), state.basis.end(), [](const PackedDynamicBasisRow& a, const PackedDynamicBasisRow& b) {
        return a.pivot < b.pivot;
    });
    for (int bi = 0; bi < (int)state.basis.size(); ++bi) {
        if (state.basis[bi].source_row_id >= 0 && state.basis[bi].source_row_id < (int)state.row_to_basis_slot.size()) {
            state.row_to_basis_slot[state.basis[bi].source_row_id] = bi;
        }
        if (state.basis[bi].pivot >= 0 && state.basis[bi].pivot < (int)state.pivot_to_basis_slot.size()) {
            state.pivot_to_basis_slot[state.basis[bi].pivot] = bi;
        }
    }

    mzd_free(M);
}

std::vector<std::vector<unsigned char>> packed_nullspace_basis(
    const std::vector<std::vector<uint64_t>>& chi_rows,
    const std::vector<unsigned char>& row_active,
    int m) {
    std::vector<std::vector<uint64_t>> mat;
    mat.reserve(chi_rows.size());
    for (int r = 0; r < (int)chi_rows.size(); ++r) {
        if (row_active[r]) mat.push_back(chi_rows[r]);
    }

    int rows = (int)mat.size();
    int rank = 0;
    std::vector<int> pivot_cols;
    pivot_cols.reserve(std::min(rows, m));

    for (int col = 0; col < m && rank < rows; ++col) {
        int pivot_row = -1;
        for (int r = rank; r < rows; ++r) {
            if (packed_get_bit(mat[r], col)) {
                pivot_row = r;
                break;
            }
        }
        if (pivot_row < 0) continue;
        if (pivot_row != rank) std::swap(mat[pivot_row], mat[rank]);

        for (int r = 0; r < rows; ++r) {
            if (r != rank && packed_get_bit(mat[r], col)) packed_xor_row(mat[r], mat[rank]);
        }
        pivot_cols.push_back(col);
        rank++;
    }

    std::vector<unsigned char> is_pivot_col(m, 0);
    for (int i = 0; i < (int)pivot_cols.size(); ++i) is_pivot_col[pivot_cols[i]] = 1;

    std::vector<std::vector<unsigned char>> basis;
    for (int free_col = 0; free_col < m; ++free_col) {
        if (is_pivot_col[free_col]) continue;
        std::vector<unsigned char> vec(m, 0);
        vec[free_col] = 1;
        for (int i = 0; i < rank; ++i) {
            int pivot = pivot_cols[i];
            vec[pivot] = (unsigned char)packed_get_bit(mat[i], free_col);
        }
        basis.push_back(vec);
    }
    return basis;
}

std::vector<std::vector<unsigned char>> dynamic_nullspace_basis(
    const std::vector<std::vector<unsigned char>>& chi_rows,
    const std::vector<unsigned char>& row_active,
    int m) {
    std::vector<std::vector<unsigned char>> mat;
    mat.reserve(chi_rows.size());
    for (int r = 0; r < (int)chi_rows.size(); ++r) {
        if (row_active[r]) mat.push_back(chi_rows[r]);
    }

    int rows = (int)mat.size();
    int rank = 0;
    std::vector<int> pivot_cols;
    pivot_cols.reserve(std::min(rows, m));

    for (int col = 0; col < m && rank < rows; ++col) {
        int pivot_row = -1;
        for (int r = rank; r < rows; ++r) {
            if (mat[r][col]) {
                pivot_row = r;
                break;
            }
        }
        if (pivot_row < 0) continue;
        if (pivot_row != rank) std::swap(mat[pivot_row], mat[rank]);

        for (int r = 0; r < rows; ++r) {
            if (r != rank && mat[r][col]) dynamic_xor_row(mat[r], mat[rank]);
        }
        pivot_cols.push_back(col);
        rank++;
    }

    std::vector<unsigned char> is_pivot_col(m, 0);
    for (int i = 0; i < (int)pivot_cols.size(); ++i) is_pivot_col[pivot_cols[i]] = 1;

    std::vector<std::vector<unsigned char>> basis;
    for (int free_col = 0; free_col < m; ++free_col) {
        if (is_pivot_col[free_col]) continue;
        std::vector<unsigned char> vec(m, 0);
        vec[free_col] = 1;
        for (int i = 0; i < rank; ++i) {
            int pivot = pivot_cols[i];
            vec[pivot] = mat[i][free_col];
        }
        basis.push_back(vec);
    }
    return basis;
}

void dynamic_collect_affected_rows(const std::vector<unsigned char>& x_prev,
                                   const std::vector<unsigned char>& x_cur,
                                   int n,
                                   std::vector<int>& affected_rows) {
    std::vector<unsigned char> marked(n * n * n, 0);
    for (int k = 0; k < n; ++k) {
        if ((x_prev[k] ^ x_cur[k]) == 0) continue;
        for (int beta = 0; beta < n; ++beta) {
            for (int gamma = 0; gamma < n; ++gamma) {
                marked[dynamic_row_id(k, beta, gamma, n)] = 1;
                marked[dynamic_row_id(beta, k, gamma, n)] = 1;
                marked[dynamic_row_id(beta, gamma, k, n)] = 1;
            }
        }
    }
    affected_rows.clear();
    for (int r = 0; r < (int)marked.size(); ++r) {
        if (marked[r]) affected_rows.push_back(r);
    }
}

void dynamic_build_full_chi_state(DynamicChiState& state, bool** A, const std::vector<unsigned char>& x_bits, int n, int m) {
    int row_count = n * n * n;
    state.chi_rows.assign(row_count, std::vector<unsigned char>(m, 0));
    state.row_active.assign(row_count, 1);
    state.row_to_basis_slot.assign(row_count, -1);
    state.pivot_to_basis_slot.assign(m, -1);
    state.x_prev = x_bits;

    for (int alpha = 0; alpha < n; ++alpha) {
        for (int beta = 0; beta < n; ++beta) {
            for (int gamma = 0; gamma < n; ++gamma) {
                int rid = dynamic_row_id(alpha, beta, gamma, n);
                build_chi_row_bool(A, x_bits, m, alpha, beta, gamma, state.chi_rows[rid]);
            }
        }
    }
    dynamic_rebuild_basis_from_state(state, m);
}

void dynamic_build_compressed_chi_state(DynamicChiState& state, bool** A, const std::vector<unsigned char>& x_bits, int n, int m) {
    int row_count = n * n * n;
    state.chi_rows.assign(row_count, std::vector<unsigned char>(m, 0));
    state.row_active.assign(row_count, 0);
    state.row_to_basis_slot.assign(row_count, -1);
    state.pivot_to_basis_slot.assign(m, -1);
    state.x_prev = x_bits;

    for (int alpha = 0; alpha < n; ++alpha) {
        for (int beta = 0; beta < n; ++beta) {
            for (int gamma = 0; gamma < n; ++gamma) {
                int rid = dynamic_row_id(alpha, beta, gamma, n);
                unsigned char active = dynamic_chi_row_active(x_bits, alpha, beta, gamma);
                state.row_active[rid] = active;
                if (!active) continue;
                build_chi_row_bool(A, x_bits, m, alpha, beta, gamma, state.chi_rows[rid]);
            }
        }
    }
    dynamic_rebuild_basis_from_state(state, m);
}

void packed_build_compressed_chi_state(PackedChiState& state, bool** A, const std::vector<unsigned char>& x_bits, int n, int m) {
    int row_count = n * n * n;
    state.words = packed_word_count(m);
    state.chi_rows.assign(row_count, std::vector<uint64_t>(state.words, 0ULL));
    state.row_active.assign(row_count, 0);
    state.row_to_basis_slot.assign(row_count, -1);
    state.pivot_to_basis_slot.assign(m, -1);
    state.x_prev = x_bits;

    for (int alpha = 0; alpha < n; ++alpha) {
        for (int beta = 0; beta < n; ++beta) {
            for (int gamma = 0; gamma < n; ++gamma) {
                int rid = dynamic_row_id(alpha, beta, gamma, n);
                unsigned char active = dynamic_chi_row_active(x_bits, alpha, beta, gamma);
                state.row_active[rid] = active;
                if (!active) continue;
                build_chi_row_packed(A, x_bits, m, alpha, beta, gamma, state.chi_rows[rid]);
            }
        }
    }
    packed_rebuild_basis_from_state(state, m);
}

void packed_build_compressed_chi_state_m4ri(PackedChiState& state, bool** A, const std::vector<unsigned char>& x_bits, int n, int m) {
    int row_count = n * n * n;
    state.words = packed_word_count(m);
    state.chi_rows.assign(row_count, std::vector<uint64_t>(state.words, 0ULL));
    state.row_active.assign(row_count, 0);
    state.row_to_basis_slot.assign(row_count, -1);
    state.pivot_to_basis_slot.assign(m, -1);
    state.x_prev = x_bits;

    for (int alpha = 0; alpha < n; ++alpha) {
        for (int beta = 0; beta < n; ++beta) {
            for (int gamma = 0; gamma < n; ++gamma) {
                int rid = dynamic_row_id(alpha, beta, gamma, n);
                unsigned char active = dynamic_chi_row_active(x_bits, alpha, beta, gamma);
                state.row_active[rid] = active;
                if (!active) continue;
                build_chi_row_packed(A, x_bits, m, alpha, beta, gamma, state.chi_rows[rid]);
            }
        }
    }
    packed_rebuild_basis_from_state_m4ri(state, m);
}

void dynamic_update_chi_state(DynamicChiState& state, bool** A, const std::vector<unsigned char>& x_bits, int n, int m,
                              std::vector<int>& affected_rows) {
    dynamic_collect_affected_rows(state.x_prev, x_bits, n, affected_rows);
    for (int idx = 0; idx < (int)affected_rows.size(); ++idx) {
        int rid = affected_rows[idx];
        int alpha = rid / (n * n);
        int beta = (rid / n) % n;
        int gamma = rid % n;
        build_chi_row_bool(A, x_bits, m, alpha, beta, gamma, state.chi_rows[rid]);
    }
    state.x_prev = x_bits;
    dynamic_rebuild_basis_from_state(state, m);
}

bool dynamic_update_chi_state_local_repair(DynamicChiState& state, bool** A, const std::vector<unsigned char>& x_bits,
                                           int n, int m, std::vector<int>& affected_rows, bool& rebuilt_basis) {
    dynamic_collect_affected_rows(state.x_prev, x_bits, n, affected_rows);
    bool basis_touched = false;
    for (int idx = 0; idx < (int)affected_rows.size(); ++idx) {
        int rid = affected_rows[idx];
        if (state.row_to_basis_slot[rid] != -1) {
            basis_touched = true;
        }
        int alpha = rid / (n * n);
        int beta = (rid / n) % n;
        int gamma = rid % n;
        build_chi_row_bool(A, x_bits, m, alpha, beta, gamma, state.chi_rows[rid]);
    }

    state.x_prev = x_bits;
    rebuilt_basis = basis_touched;
    if (basis_touched) {
        dynamic_rebuild_basis_from_state(state, m);
        return true;
    }

    for (int idx = 0; idx < (int)affected_rows.size(); ++idx) {
        int rid = affected_rows[idx];
        if (dynamic_is_zero_row(state.chi_rows[rid])) continue;
        dynamic_add_row_to_basis(state.chi_rows[rid], rid, m, state.basis, state.row_to_basis_slot, state.pivot_to_basis_slot);
    }
    return true;
}

bool dynamic_update_chi_state_local_repair_compressed(DynamicChiState& state, bool** A, const std::vector<unsigned char>& x_bits,
                                                      int n, int m, std::vector<int>& affected_rows, bool& rebuilt_basis) {
    dynamic_collect_affected_rows(state.x_prev, x_bits, n, affected_rows);
    bool basis_touched = false;
    for (int idx = 0; idx < (int)affected_rows.size(); ++idx) {
        int rid = affected_rows[idx];
        if (state.row_to_basis_slot[rid] != -1) basis_touched = true;

        int alpha = rid / (n * n);
        int beta = (rid / n) % n;
        int gamma = rid % n;
        unsigned char active = dynamic_chi_row_active(x_bits, alpha, beta, gamma);
        state.row_active[rid] = active;
        if (active) build_chi_row_bool(A, x_bits, m, alpha, beta, gamma, state.chi_rows[rid]);
    }

    state.x_prev = x_bits;
    rebuilt_basis = basis_touched;
    if (basis_touched) {
        dynamic_rebuild_basis_from_state(state, m);
        return true;
    }

    for (int idx = 0; idx < (int)affected_rows.size(); ++idx) {
        int rid = affected_rows[idx];
        if (!state.row_active[rid]) continue;
        if (dynamic_is_zero_row(state.chi_rows[rid])) continue;
        dynamic_add_row_to_basis(state.chi_rows[rid], rid, m, state.basis, state.row_to_basis_slot, state.pivot_to_basis_slot);
    }
    return true;
}

bool packed_update_chi_state_local_repair_compressed(PackedChiState& state, bool** A, const std::vector<unsigned char>& x_bits,
                                                     int n, int m, std::vector<int>& affected_rows, bool& rebuilt_basis) {
    dynamic_collect_affected_rows(state.x_prev, x_bits, n, affected_rows);
    bool basis_touched = false;
    for (int idx = 0; idx < (int)affected_rows.size(); ++idx) {
        int rid = affected_rows[idx];
        if (state.row_to_basis_slot[rid] != -1) basis_touched = true;

        int alpha = rid / (n * n);
        int beta = (rid / n) % n;
        int gamma = rid % n;
        unsigned char active = dynamic_chi_row_active(x_bits, alpha, beta, gamma);
        state.row_active[rid] = active;
        if (active) build_chi_row_packed(A, x_bits, m, alpha, beta, gamma, state.chi_rows[rid]);
    }

    state.x_prev = x_bits;
    rebuilt_basis = basis_touched;
    if (basis_touched) {
        packed_rebuild_basis_from_state(state, m);
        return true;
    }

    for (int idx = 0; idx < (int)affected_rows.size(); ++idx) {
        int rid = affected_rows[idx];
        if (!state.row_active[rid]) continue;
        if (packed_is_zero_row(state.chi_rows[rid])) continue;
        packed_add_row_to_basis(state.chi_rows[rid], rid, m, state.basis, state.row_to_basis_slot, state.pivot_to_basis_slot);
    }
    return true;
}

bool packed_update_chi_state_local_repair_compressed_m4ri(PackedChiState& state, bool** A, const std::vector<unsigned char>& x_bits,
                                                          int n, int m, std::vector<int>& affected_rows, bool& rebuilt_basis) {
    dynamic_collect_affected_rows(state.x_prev, x_bits, n, affected_rows);
    bool basis_touched = false;
    for (int idx = 0; idx < (int)affected_rows.size(); ++idx) {
        int rid = affected_rows[idx];
        if (state.row_to_basis_slot[rid] != -1) basis_touched = true;

        int alpha = rid / (n * n);
        int beta = (rid / n) % n;
        int gamma = rid % n;
        unsigned char active = dynamic_chi_row_active(x_bits, alpha, beta, gamma);
        state.row_active[rid] = active;
        if (active) build_chi_row_packed(A, x_bits, m, alpha, beta, gamma, state.chi_rows[rid]);
    }

    state.x_prev = x_bits;
    rebuilt_basis = basis_touched;
    if (basis_touched) {
        packed_rebuild_basis_from_state_m4ri(state, m);
        return true;
    }

    for (int idx = 0; idx < (int)affected_rows.size(); ++idx) {
        int rid = affected_rows[idx];
        if (!state.row_active[rid]) continue;
        if (packed_is_zero_row(state.chi_rows[rid])) continue;
        packed_add_row_to_basis(state.chi_rows[rid], rid, m, state.basis, state.row_to_basis_slot, state.pivot_to_basis_slot);
    }
    return true;
}

bool** packed_active_rows_to_m4ri_nullspace(const PackedChiState& state, int m, int& d_ns) {
    int rows = 0;
    for (int r = 0; r < (int)state.chi_rows.size(); ++r) if (state.row_active[r]) rows++;
    mzd_t* M = mzd_init(rows, m);
    int rr = 0;
    for (int r = 0; r < (int)state.chi_rows.size(); ++r) {
        if (!state.row_active[r]) continue;
        for (int c = 0; c < m; ++c) {
            if (packed_get_bit(state.chi_rows[r], c)) mzd_write_bit(M, rr, c, 1);
        }
        rr++;
    }
    bool** NS = M4RI_direct_nullspace(M, d_ns);
    mzd_free(M);
    return NS;
}

bool dynamic_try_single_basis_repair(DynamicChiState& state,
                                     const std::vector<int>& affected_rows,
                                     int m,
                                     int& repaired_basis_rows) {
    std::vector<int> affected_basis_rows;
    affected_basis_rows.reserve(affected_rows.size());
    for (int idx = 0; idx < (int)affected_rows.size(); ++idx) {
        int rid = affected_rows[idx];
        if (state.row_to_basis_slot[rid] != -1) affected_basis_rows.push_back(rid);
    }

    repaired_basis_rows = (int)affected_basis_rows.size();
    if ((int)affected_basis_rows.size() != 1) return false;

    const int old_basis_size = (int)state.basis.size();
    const int bad_row_id = affected_basis_rows[0];
    const int bad_slot = state.row_to_basis_slot[bad_row_id];
    if (bad_slot < 0 || bad_slot >= old_basis_size) return false;

    std::vector<DynamicBasisRow> tmp_basis;
    tmp_basis.reserve(old_basis_size);
    for (int bi = 0; bi < old_basis_size; ++bi) {
        if (bi == bad_slot) continue;
        tmp_basis.push_back(state.basis[bi]);
    }

    std::vector<int> tmp_row_to_basis_slot(state.row_to_basis_slot.size(), -1);
    std::vector<int> tmp_pivot_to_basis_slot(state.pivot_to_basis_slot.size(), -1);
    dynamic_reindex_basis(tmp_basis, tmp_row_to_basis_slot, tmp_pivot_to_basis_slot);

    for (int idx = 0; idx < (int)affected_rows.size(); ++idx) {
        int rid = affected_rows[idx];
        if (dynamic_is_zero_row(state.chi_rows[rid])) continue;
        dynamic_add_row_to_basis(state.chi_rows[rid], rid, m, tmp_basis, tmp_row_to_basis_slot, tmp_pivot_to_basis_slot);
    }

    if ((int)tmp_basis.size() != old_basis_size) return false;

    state.basis.swap(tmp_basis);
    state.row_to_basis_slot.swap(tmp_row_to_basis_slot);
    state.pivot_to_basis_slot.swap(tmp_pivot_to_basis_slot);
    return true;
}

bool dynamic_update_chi_state_local_repair_k1(DynamicChiState& state, bool** A, const std::vector<unsigned char>& x_bits,
                                              int n, int m, std::vector<int>& affected_rows,
                                              bool& rebuilt_basis, bool& repaired_basis, int& repaired_basis_rows) {
    dynamic_collect_affected_rows(state.x_prev, x_bits, n, affected_rows);
    bool basis_touched = false;
    for (int idx = 0; idx < (int)affected_rows.size(); ++idx) {
        int rid = affected_rows[idx];
        if (state.row_to_basis_slot[rid] != -1) basis_touched = true;
        int alpha = rid / (n * n);
        int beta = (rid / n) % n;
        int gamma = rid % n;
        build_chi_row_bool(A, x_bits, m, alpha, beta, gamma, state.chi_rows[rid]);
    }

    state.x_prev = x_bits;
    rebuilt_basis = false;
    repaired_basis = false;
    repaired_basis_rows = 0;

    if (dynamic_try_single_basis_repair(state, affected_rows, m, repaired_basis_rows)) {
        repaired_basis = true;
        return true;
    }

    rebuilt_basis = basis_touched;
    if (basis_touched) {
        dynamic_rebuild_basis_from_state(state, m);
        return true;
    }

    for (int idx = 0; idx < (int)affected_rows.size(); ++idx) {
        int rid = affected_rows[idx];
        if (dynamic_is_zero_row(state.chi_rows[rid])) continue;
        dynamic_add_row_to_basis(state.chi_rows[rid], rid, m, state.basis, state.row_to_basis_slot, state.pivot_to_basis_slot);
    }
    return true;
}
} // namespace
// --- M4RI 統合ラッパー関数 ---
// (これらの関数は LempelX / LempelX2 から呼び出されます)

static mzd_t* convert_to_mzd(bool const** A, int n, int m) {
    mzd_t* M = mzd_init(n, m);
    for (rci_t r = 0; r < n; ++r) {
        if (A[r]) {
            for (rci_t c = 0; c < m; ++c) {
                if (A[r][c]) mzd_write_bit(M, r, c, 1);
            }
        }
    }
    return M;
}

static bool** convert_from_mzd(mzd_t* X) {
    if (X == NULL) return NULL;
    int out_rows = X->nrows;
    int out_cols = X->ncols;
    bool** out = LCL_Mat_GF2::construct(out_rows, out_cols);
    if (!out) return NULL;
    for (rci_t r = 0; r < out_rows; ++r) {
        for (rci_t c = 0; c < out_cols; ++c) {
            out[r][c] = (bool)mzd_read_bit(X, r, c);
        }
    }
    return out;
}

/**
 * @brief M4RI の nullspace (mzd_kernel_left_pluq) を呼び出すラッパー関数。
 * LCL_Mat_GF2::nullspace と同じインターフェースを提供します。
 */
bool** GateSynthesisMatrix::M4RI_wrapper_for_nullspace(bool const** A, int n, int m, int& out_d) {
    // 1. bool** から mzd_t* へ変換 (これがコピーの代わりにもなります)
    mzd_t* A_m4ri = convert_to_mzd(A, n, m);
    
    // 2. M4RI のカーネル（零空間）計算関数を呼び出す
    // 
    // これがM4RIによるガウスの消去法（PLUQ分解経由）の実行部分です。
    mzd_t* X_m4ri = mzd_kernel_left_pluq(A_m4ri, 0); // 0 = デフォルトのカットオフ

    bool** NS = NULL; // LCL 形式の出力行列

    // 3. 結果を処理する
    if (X_m4ri == NULL) {
        // カーネルが自明（ランク = m）の場合
        out_d = 0;
        NS = NULL;
    } else {
        // カーネルが見つかった場合
        out_d = X_m4ri->ncols; // カーネルの次元（自由度）が d となります
        
        // 4. mzd_t* から bool** へ結果を変換する
        NS = convert_from_mzd(X_m4ri);
    }

    // 5. M4RIで確保したメモリを解放する
    mzd_free(A_m4ri); // 入力行列のコピー
    if (X_m4ri) {
        mzd_free(X_m4ri); // 結果の行列
    }
    
    return NS;
}

/**
 * @brief 新規追加: 既に mzd_t* 形式の行列を受け取り、Nullspaceを計算して LCL形式で返す
 * Chi行列計算後の再変換コストを防ぐために使用します。
 A_in はここで破壊される可能性があるため
 */
bool** M4RI_direct_nullspace(mzd_t* A_in, int& out_d) {

    mzd_t* X_m4ri = mzd_kernel_left_pluq(A_in, 0);

    bool** NS = NULL;
    if (X_m4ri == NULL) {
        out_d = 0;
    } else {
        out_d = X_m4ri->ncols;
        NS = convert_from_mzd(X_m4ri);
    }

    if (X_m4ri) mzd_free(X_m4ri);
    return NS;
}

/**
 * @brief 新規追加: Chi行列を直接 M4RI (mzd_t) に書き込む関数
 * 入力 A, x はアクセスの利便性から LCL(bool**) のままにし、出力 Aext のみ mzd_t にします。
 */
/*void GateSynthesisMatrix::Chi_M4RI(bool** A, bool** x, int n, int m, mzd_t* Aext) {
    int i = 0;
    int temp = 0;
    
    for(int alpha = 0; alpha < n; alpha++) {
        bool x_alpha = x[alpha][0];
        for(int beta = 0; beta < n; beta++) {
            bool x_beta = x[beta][0];
            for(int gamma = 0; gamma < n; gamma++) {
                bool x_gamma = x[gamma][0];
                
                bool term_const = x_alpha && x_beta && x_gamma;
                
                for(int j = 0; j < m; j++) {
                    
                    bool A_a = A[alpha][j];
                    bool A_b = A[beta][j];
                    bool A_c = A[gamma][j];

                    int val = term_const
                            ^ (x_alpha & x_beta & A_c)
                            ^ (x_beta & x_gamma & A_a)
                            ^ (x_gamma & x_alpha & A_b)
                            ^ (x_alpha & A_b & A_c)
                            ^ (x_beta & A_c & A_a)
                            ^ (x_gamma & A_a & A_b);
                    
                    if (val & 1) {
                        mzd_write_bit(Aext, i, j, 1);
                    }
                }
                i++;
            }
        }
    }
}
    */

// 高速化版: 行単位でビット演算を行う
void GateSynthesisMatrix::Chi_M4RI(mzd_t* A, bool** x, int n, int m, mzd_t* Aext) {
    // Aext は呼び出し元でゼロクリアされている前提ですが、念の為ここでもクリア可
    // ここでは上書き合成していくので、ゼロクリアは呼び出し元(LempelX2)の責任とします。

    int row_idx = 0;
    
    // M4RIの行データへのポインタを使ってワード単位で演算するための準備
    // A->width は "ワード(64bit)の数" です
    int width = A->width; 

    for(int alpha = 0; alpha < n; alpha++) {
        bool x_a = x[alpha][0];
        
        for(int beta = 0; beta < n; beta++) {
            bool x_b = x[beta][0];
            
            for(int gamma = 0; gamma < n; gamma++) {
                bool x_c = x[gamma][0];
                
                // 定数項: x_a * x_b * x_c
                // これが1なら、その行のすべての要素に1を足す（=ビット反転）
                bool term_const = x_a && x_b && x_c;
                
                // 行単位の計算
                // Aext[row_idx] = term_const 
                //               ^ (x_a*x_b)*A[gamma] ^ (x_b*x_c)*A[alpha] ^ (x_c*x_a)*A[beta]
                //               ^ x_a*(A[beta]&A[gamma]) ^ ...

                // 各ワード(64bit塊)ごとにループ
                for(int w = 0; w < width; w++) {
                    word res = 0;

                    // 1. 定数項 (全ビット1 または 0)
                    if(term_const) res = ~res; // 全ビット反転(111...111)

                    // 2. 線形項 (x_a*x_b * A[gamma] など)
                    if(x_a && x_b) res ^= mzd_row(A, gamma)[w];
                    if(x_b && x_c) res ^= mzd_row(A, alpha)[w];
                    if(x_c && x_a) res ^= mzd_row(A, beta)[w];

                    // 3. 2次項 (x_a * A[beta] * A[gamma] など)
                    // ワード同士の AND をとってから XOR
                    if(x_a) res ^= (mzd_row(A, beta)[w] & mzd_row(A, gamma)[w]);
                    if(x_b) res ^= (mzd_row(A, gamma)[w] & mzd_row(A, alpha)[w]);
                    if(x_c) res ^= (mzd_row(A, alpha)[w] & mzd_row(A, beta)[w]);

                    // 結果を書き込み
                    mzd_row(Aext, row_idx)[w] = res;
                }
                
                row_idx++;
            }
        }
    }
}

bool** GateSynthesisMatrix::from_signature(bool*** S, int n, int& mp) {
    int N = (int)pow(2,n);
    bool* a = new bool[N];
    for(int i = 0; i < N; i++) a[i] = 0;

    int I;
    for(int i = 0; i < n; i++) {
        if(S[i][i][i]) {
            I = (int)pow(2,i);
            a[I] = !a[I];
        }
        if(i<(n-1)) {
            for(int j = (i+1); j < n; j++) {
                if(S[i][j][j]) {
                    I = (int)pow(2,i);
                    a[I] = !a[I];
                    I = (int)pow(2,j);
                    a[I] = !a[I];
                    I = (int)pow(2,i)+(int)pow(2,j);
                    a[I] = !a[I];
                }
                if(j<(n-1)) {
                    for(int k = (j+1); k < n; k++) {
                        if(S[i][j][k]) {
                            I = (int)pow(2,i);
                            a[I] = !a[I];
                            I = (int)pow(2,j);
                            a[I] = !a[I];
                            I = (int)pow(2,k);
                            a[I] = !a[I];
                            I = (int)pow(2,i)+(int)pow(2,j);
                            a[I] = !a[I];
                            I = (int)pow(2,i)+(int)pow(2,k);
                            a[I] = !a[I];
                            I = (int)pow(2,j)+(int)pow(2,k);
                            a[I] = !a[I];
                            I = (int)pow(2,i)+(int)pow(2,j)+(int)pow(2,k);
                            a[I] = !a[I];
                        }
                    }
                }
            }
        }
    }
    int m = 0;
    for(int i = 0; i < N; i++) {
        if(a[i]) {
            m++;
        }
    }
    bool** out = LCL_Mat_GF2::construct(n,m);
    int j = 0;
    for(int I = 0; I < N; I++) {
        if(a[I]) {
            int Ip = I;
            int i = 0;
            while(Ip>0) {
                out[i][j] = (Ip%2);
                Ip /= 2;
                i++;
            }
            j++;
        }
    }

    delete [] a;

    mp = m;
    return out;
}

void GateSynthesisMatrix::cleanup(bool** A, int n, int m, int& mp) {
    for(int j1 = 0; j1 < (m-1); j1++) {
        for(int j2 = (j1+1); j2 < m; j2++) {
            bool same = 1;
            for(int i = 0; same&&(i < n); i++) {
                same *= (A[i][j1]==A[i][j2]);
            }
            if(same) {
                for(int i = 0; i < n; i++) {
                    A[i][j1]=0;
                    A[i][j2]=0;
                }
            }
        }
    }

    int j_end = (m-1);
    int non_zero_count = 0;
    for(int j = 0; j < m; j++) {
        int sum = 0;
        for(int i = 0; i < n; i++) {
            sum += A[i][j];
        }
        if(!sum) {
            bool found = false;
            while((!found)&&(j_end>j)) {
                for(int i = 0; (!found)&&(i < n); i++) {
                    found = A[i][j_end];
                }
                if(!found) j_end--;
            }
            if(found) {
                LCL_Mat_GF2::swapcol(A,n,m,j,j_end);
                non_zero_count++;
            }
        } else {
            non_zero_count++;
        }
    }
    mp = non_zero_count;
}

void GateSynthesisMatrix::LempelX(bool** A, int n, int m, int& omp) {
    int this_m = m;
    bool** x = LCL_Mat_GF2::construct(n,1);
    bool** nv = LCL_Mat_GF2::construct(1,m+1);
    bool** xnv = LCL_Mat_GF2::construct(n,m+1);
    bool** A_xnv = LCL_Mat_GF2::construct(n,m+1);
    int n_ext = n+((n*(n-1)*(n-2))/6);
    bool** A_ext = LCL_Mat_GF2::construct(n_ext,m+1);

    //時間計測用
    std::cout << "in LemopelX" <<endl;
    std::chrono::microseconds g_total_nullspace_duration(0);

    /*bool __x[n][1];
    bool __nv[1][m+1];
    bool __xnv[n][m+1];
    bool __A_xnv[n][m+1];
    bool __A_ext[n_ext][m+1];

    bool* _x = (bool*)&__x;
    bool* _nv = (bool*)&__nv;
    bool* _xnv = (bool*)&__xnv;
    bool* _A_xnv = (bool*)&__A_xnv;
    bool* _A_ext = (bool*)&__A_ext;

    bool** x = (bool**)&_x;
    bool** nv = (bool**)&_nv;
    bool** xnv = (bool**)&_xnv;
    bool** A_xnv = (bool**)&_A_xnv;
    bool** A_ext = (bool**)&_A_ext;*/

    bool found = 1;
    int round = 0;
    while(found&&(round<m)) {
        cout << "A for round " << round << ":" << endl;
        LCL_Mat_GF2::print((const bool**)A,n,this_m);
        found = 0;
        LOut(); cout << "Round = " << round << endl;
        LCL_Mat_GF2::copy((const bool**)A,n,this_m,A_ext);
		int col_perm[this_m]; for(int i = 0; i < this_m; i++) col_perm[i]=i;
		LCL_Int::randperm(col_perm,this_m);
        for(int j1_ind = 0; (!found)&&(j1_ind < (this_m-1)); j1_ind++) {
            for(int j2_ind = (j1_ind+1); (!found)&&(j2_ind < this_m); j2_ind++) {
				int j1 = col_perm[j1_ind];
				int j2 = col_perm[j2_ind];
                for(int i = 0; i < n; i++) {
                    x[i][0] = (A[i][j1] + A[i][j2])%2;
                }
                int I = 0;
                for(int a = 0; a < (n-2); a++) {
                    for(int b = (a+1); b < (n-1); b++) {
                        for(int c = (b+1); c < n; c++) {
                            for(int j = 0; j < this_m; j++) {
                                A_ext[n+I][j] = (x[a][0]*A[b][j]*A[c][j] + x[b][0]*A[c][j]*A[a][j] + x[c][0]*A[a][j]*A[b][j])%2;
                            }
                            I++;
                        }
                    }
                }
                int d=-1;
                int m_NS = this_m;


                //零空間計算
                auto start_ns = std::chrono::high_resolution_clock::now();
                bool** NS = LCL_Mat_GF2::nullspace((const bool**)A_ext,n_ext,m_NS,d);
                auto end_ns = std::chrono::high_resolution_clock::now();
                auto duration_this_call = std::chrono::duration_cast<std::chrono::microseconds>(end_ns - start_ns);
                g_total_nullspace_duration += duration_this_call;

                //bool** NS = ::M4RI_wrapper_for_nullspace((const bool**)A_ext,n_ext,m_NS,d);
                found = 0;
                int nsv = -1;
                for(int h = 0; (!found)&&(h < d); h++) {
                    found = (NS[j1][h]+NS[j2][h])%2;
                    if(found) {
                        nsv = h;
                    }
                }
                if(found) {
                    int weight_nv = 0;
                    for(int h = 0; h < this_m; h++) {
                        nv[0][h] = NS[h][nsv];
                        weight_nv += nv[0][h];
                    }
                    if(weight_nv%2) {
                        nv[0][this_m]=1;
                        for(int h = 0; h < n; h++) A[h][this_m]=0;
                        this_m++;
                    }
                    LCL_Mat_GF2::times((const bool**)x,(const bool**)nv,n,1,this_m,xnv);
                    LCL_Mat_GF2::add((const bool**)A,(const bool**)xnv,n,this_m,A_xnv);
                    LCL_Mat_GF2::copy((const bool**)A_xnv,n,this_m,A);
                    int mp=-1;
                    GateSynthesisMatrix::cleanup(A,n,this_m,mp);
                    this_m = mp;
                }
                if(NS) LCL_Mat_GF2::destruct(NS,m_NS,d);
            }
        }
        round++;
        cout << "Round end" << endl;
    }
    cout << "OUT OF LOOP" << endl;

    /*LCL_Mat_GF2::print((const bool**)A,n,this_m,"A\n");
    LCL_Mat_GF2::print((const bool**)x,n,1,"x\n");
    LCL_Mat_GF2::print((const bool**)A_ext,n_ext,m+1,"A_ext\n");
    LCL_Mat_GF2::print((const bool**)nv,1,m+1,"nv\n");
    LCL_Mat_GF2::print((const bool**)xnv,n,m+1,"xnv\n");
    LCL_Mat_GF2::print((const bool**)A_xnv,n,m+1,"A_xnv\n");*/

    cout << "Destroying x" << endl;
    LCL_Mat_GF2::destruct(x,n,1);
    cout << "Destroying A_ext" << endl;
    LCL_Mat_GF2::destruct(A_ext,n_ext,m+1);
    cout << "Destroying nv" << endl;
    LCL_Mat_GF2::destruct(nv,1,m+1);
    cout << "Destroying xnv" << endl;
    LCL_Mat_GF2::destruct(xnv,n,m+1);
    cout << "Destroying A_xnv" << endl;
    LCL_Mat_GF2::destruct(A_xnv,n,m+1);
    omp = this_m;
    cout << "END OF LEMPELX" << endl;
}

void GateSynthesisMatrix::Chi(bool** A, bool** x, int n, int m, bool** Aext) {
    int i = 0;
    int temp = 0;
    for(int alpha = 0; alpha < n; alpha++) {
        for(int beta = 0; beta < n; beta++) {
            for(int gamma = 0; gamma < n; gamma++) {
                for(int j = 0; j < m; j++) {
                    temp = x[alpha][0]*x[beta][0]*x[gamma][0]
                        + x[alpha][0]*x[beta][0]*A[gamma][j]
                        + x[beta][0]*x[gamma][0]*A[alpha][j]
                        + x[gamma][0]*x[alpha][0]*A[beta][j]
                        + x[alpha][0]*A[beta][j]*A[gamma][j]
                        + x[beta][0]*A[gamma][j]*A[alpha][j]
                        + x[gamma][0]*A[alpha][j]*A[beta][j];
                    Aext[i][j] = (temp%2);
                }
                i++;
            }
        }
    }
}


void GateSynthesisMatrix::ChiPrime(bool** A, bool** x, int n, int m, bool** Aext) {
    int i = 0;
    int temp = 0;
    for(int alpha = 0; alpha < n; alpha++) {
        for(int beta = 0; beta < n; beta++) {
            for(int gamma = 0; gamma < n; gamma++) {
                for(int j = 0; j < m; j++) {
                    temp = x[alpha][0]*x[beta][0]*A[gamma][j]
                        + x[beta][0]*x[gamma][0]*A[alpha][j]
                        + x[gamma][0]*x[alpha][0]*A[beta][j]
                        + x[alpha][0]*A[beta][j]*A[gamma][j]
                        + x[beta][0]*A[gamma][j]*A[alpha][j]
                        + x[gamma][0]*A[alpha][j]*A[beta][j];
                    Aext[i][j] = (temp%2);
                }
                i++;
            }
        }
    }
}

// ハミング距離とペア情報を保持する構造体
struct ColumnPair {
    int j1;
    int j2;
    int distance;
    // 距離が小さい順にソートするための比較演算子
    bool operator<(const ColumnPair& other) const {
        return distance < other.distance;
    }
};

void GateSynthesisMatrix::LempelX3(bool** A, int n, int m, int& omp) {
    auto start_total = std::chrono::high_resolution_clock::now();
    
    std::cout << endl << "in LempelX3 (Hamming Distance Version)" << endl;
    int this_m = m;
    
    // 作業用領域の確保
    bool** x = LCL_Mat_GF2::construct(n, 1);
    int n_chi_A = n * n * n;
    bool** chi_A = LCL_Mat_GF2::construct(n_chi_A, m + 1);
    bool** Anew = LCL_Mat_GF2::construct(n, m + 1);
    bool** Abest = LCL_Mat_GF2::construct(n, m + 1);
    
    LCL_Mat_GF2::copy((const bool**)A, n, m, Abest);
    int m_best = m;

    std::chrono::microseconds total_ns_duration(0);
    std::chrono::microseconds total_chi_duration(0);

    bool found = true;
    int round = 0;

    while (found && (round < m)) {
        found = false;
        LOut(); cout << "Round = " << round << " | Current Columns: " << this_m << endl;

        // --- ステップ1: 全ペアのハミング距離を計算 ---
        std::vector<ColumnPair> pair_list;
        pair_list.reserve(this_m * (this_m - 1) / 2);

        for (int j1 = 0; j1 < this_m - 1; j1++) {
            for (int j2 = j1 + 1; j2 < this_m; j2++) {
                int dist = 0;
                for (int i = 0; i < n; i++) {
                    if (A[i][j1] != A[i][j2]) dist++;
                }
                pair_list.push_back({j1, j2, dist});
            }
        }

        // --- ステップ2: ハミング距離が小さい順にソート ---
        std::sort(pair_list.begin(), pair_list.end());

        // --- ステップ3: ソートされた順にペアを試行 ---
        for (const auto& pair : pair_list) {
            if (found) break; // 削減が見つかったら次のラウンドへ

            int c1 = pair.j1;
            int c2 = pair.j2;

            // x = A[c1] + A[c2] (mod 2)
            for (int i = 0; i < n; i++) {
                x[i][0] = (A[i][c1] + A[i][c2]) % 2;
            }

            // Chi行列の計算
            auto s_chi = std::chrono::high_resolution_clock::now();
            GateSynthesisMatrix::Chi(A, x, n, this_m, chi_A);
            total_chi_duration += std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - s_chi);

            // Nullspaceの計算
            int d = 0;
            auto s_ns = std::chrono::high_resolution_clock::now();
            bool** NS = LCL_Mat_GF2::nullspace((const bool**)chi_A, n_chi_A, this_m, d);
            total_ns_duration += std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - s_ns);

            // 有効な解があるか確認
            int nsv = -1;
            for (int h = 0; (!found) && (h < d); h++) {
                if ((NS[c1][h] + NS[c2][h]) % 2 == 1) {
                    found = true;
                    nsv = h;
                }
            }

            if (found) {
                // 行列 A の更新
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < this_m; j++) {
                        Anew[i][j] = (A[i][j] + x[i][0] * NS[j][nsv]) % 2;
                    }
                }
                
                int mp;
                GateSynthesisMatrix::cleanup(Anew, n, this_m, mp);
                
                if (mp < m_best) {
                    LCL_Mat_GF2::copy((const bool**)Anew, n, mp, Abest);
                    m_best = mp;
                    cout << "  [HIT] Reduced to m = " << mp << " using Dist = " << pair.distance << endl;
                }
            }

            if (d > 0 && NS != NULL) {
                LCL_Mat_GF2::destruct(NS, this_m, d);
            }
        }

        // ベストな状態を A に反映して次のラウンドへ
        LCL_Mat_GF2::copy((const bool**)Abest, n, m_best, A);
        this_m = m_best;
        round++;
    }

    // 後処理
    LCL_Mat_GF2::destruct(x, n, 1);
    LCL_Mat_GF2::destruct(chi_A, n_chi_A, m + 1);
    LCL_Mat_GF2::destruct(Anew, n, m + 1);
    LCL_Mat_GF2::destruct(Abest, n, m + 1);
    
    omp = this_m;
    
    auto end_total = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
    
    cout << "Total LempelX2 (Hamming) time: " << duration.count() << " ms" << endl;
    cout << "Final Nullspace total time: " << std::chrono::duration_cast<std::chrono::milliseconds>(total_ns_duration).count() << " ms" << endl;
}
void GateSynthesisMatrix::LempelX2(bool** A, int n, int m, int& omp) {

    auto start_total = std::chrono::high_resolution_clock::now();
    auto end_total = std::chrono::high_resolution_clock::now();
    auto duration_total = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
    //M4RI版
    if(false){
        bool** A_copy = LCL_Mat_GF2::construct(n, m + 1);
        LCL_Mat_GF2::copy((const bool**)A, n, m, A_copy);
        int omp_m4ri = m;
        //時間計測
        start_total = std::chrono::high_resolution_clock::now();
        LempelX2_M4RI_DetailedStats(A_copy, n, m, omp_m4ri);
        end_total = std::chrono::high_resolution_clock::now();
        duration_total = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
        std::cout << "Total LempelX2_M4RI time: " << duration_total.count() << " ms" << endl;
        std::cout << "M4RI result size: " << omp_m4ri << " columns" << endl;
        LCL_Mat_GF2::destruct(A_copy, n, m + 1);
        return;
    }
    //Hamming距離版
    if(false){
        std::cout << "\n--- Testing Hamming Distance Version ---" << endl;
        // 元の行列Aから新しいコピーを作成 (公平な比較のため)
        bool** A_copy = LCL_Mat_GF2::construct(n, m + 1);
        LCL_Mat_GF2::copy((const bool**)A, n, m, A_copy);
        
        int omp_hamming = m;

        // 時間計測開始
        start_total = std::chrono::high_resolution_clock::now();
        GateSynthesisMatrix::LempelX2_M4RI_Hamming(A_copy, n, m, omp_hamming);
        end_total = std::chrono::high_resolution_clock::now();

        //消すかも
        LCL_Mat_GF2::copy((const bool**)A_copy, n, omp_hamming, A);
        omp = omp_hamming;
        //
        
        duration_total = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
        std::cout << "Total LempelX2_M4RI_Hamming time: " << duration_total.count() << " ms" << endl;
        std::cout << "Hamming result size: " << omp_hamming << " columns" << endl;

        // メモリ解放
        LCL_Mat_GF2::destruct(A_copy, n, m + 1);
        return;
    }
    //ランダム＋ビームサーチ版
    /*
    if(false){
        std::cout << "\n--- Testing randaom Beam Search Version (K=5, N=50) ---" << endl;
        bool** A_copy = LCL_Mat_GF2::construct(n, m + 1);
        LCL_Mat_GF2::copy((const bool**)A, n, m, A_copy);
        
        int omp_beam = m;
        // ビーム幅などのパラメータ（必要に応じて引数化してください）
        int beam_width = 5; 

        start_total = std::chrono::high_resolution_clock::now();
        
        // ビームサーチ版の関数を呼び出し
        GateSynthesisMatrix::LempelX2_M4RI_RandomBeamSearch(A_copy, n, m, omp_beam);
        
        end_total = std::chrono::high_resolution_clock::now();
        duration_total = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
        
        std::cout << "Total LempelX2_M4RI_randam_BeamSearch time: " << duration_total.count() << " ms" << endl;
        std::cout << "random_Beam Search result size: " << omp_beam << " columns" << endl;

        LCL_Mat_GF2::destruct(A_copy, n, m + 1);
    }
    */
    //ハミング距離＋ビームサーチ版
    if(false){
        std::cout << "\n--- Testing raoundrobinhammingBeamSearch Version (K=5, N=50) ---" << endl;
        bool** A_copy = LCL_Mat_GF2::construct(n, m + 1);
        LCL_Mat_GF2::copy((const bool**)A, n, m, A_copy);
        
        int omp_beam = m;
        // ビーム幅などのパラメータ（必要に応じて引数化してください）
        int beam_width = 5; 

        start_total = std::chrono::high_resolution_clock::now();
        
        // ビームサーチ版の関数を呼び出し
        GateSynthesisMatrix::LempelX2_M4RI_BeamSearch(A_copy, n, m, omp_beam);
        
        end_total = std::chrono::high_resolution_clock::now();
        duration_total = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
        
        std::cout << "Total LempelX2_M4RI_hammingBeamSearch time: " << duration_total.count() << " ms" << endl;
        std::cout << "hamming Beam Search result size: " << omp_beam << " columns" << endl;

        LCL_Mat_GF2::destruct(A_copy, n, m + 1);
    }

    if(false){
        std::cout << "\n--- Testing sequentialBeamSearch Version (K=5, N=50) ---" << endl;
        bool** A_copy = LCL_Mat_GF2::construct(n, m + 1);
        LCL_Mat_GF2::copy((const bool**)A, n, m, A_copy);
        
        int omp_beam = m;
        // ビーム幅などのパラメータ（必要に応じて引数化してください）
        int beam_width = 5; 

        start_total = std::chrono::high_resolution_clock::now();
        
        // ビームサーチ版の関数を呼び出し
        GateSynthesisMatrix::LempelX2_M4RI_SequentialBeamSearch(A_copy, n, m, omp_beam);
        
        end_total = std::chrono::high_resolution_clock::now();
        duration_total = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
        
        std::cout << "Total LempelX2_M4RI_hammingBeamSearch time: " << duration_total.count() << " ms" << endl;
        std::cout << "hamming Beam Search result size: " << omp_beam << " columns" << endl;

        LCL_Mat_GF2::destruct(A_copy, n, m + 1);
        return;
    }
    
    
    //オリジナル版
    start_total = std::chrono::high_resolution_clock::now();
    std::cout << endl<<"in LempelX2" <<endl;
    int this_m = m;
    bool** x = LCL_Mat_GF2::construct(n,1);
    bool** y = LCL_Mat_GF2::construct(m+1,1);
    bool** xyT = LCL_Mat_GF2::construct(n,m+1);
    bool** A_xyT = LCL_Mat_GF2::construct(n,m+1);
    int n_chi_A = n*n*n;
    bool** chi_A = LCL_Mat_GF2::construct(n_chi_A,m+1);

    std::chrono::microseconds g_total_nullspace_duration(0);
    std::chrono::microseconds total_chi_duration(0);

    // Column randomizers
    int* r_j1 = new int[m];
    int* r_j2 = new int[m];

    // Provisional Anew
    bool** Anew = LCL_Mat_GF2::construct(n,m+1);
    LCL_Mat_GF2::copy((const bool**)A,n,m,Anew);
    //int m_new = m;

    // Current best A
    bool** Abest = LCL_Mat_GF2::construct(n,m+1);
    LCL_Mat_GF2::copy((const bool**)A,n,m,Abest);
    int m_best = m;

    bool found = 1;
    int round = 0;
    while(found&&(round<m)) {
        found = 0;
        LOut(); cout << "Round = " << round << endl;
        //後で消す
        std::cout << "Round " << round << " | Current Columns (m): " << this_m << std::endl;
        //LCL_Mat_GF2::print(A,n,this_m,"A: ");
        LCL_Int::randperm(r_j1,this_m-1);
        for(int j1 = 0; (!found)&&(j1 < (this_m-1)); j1++) {
            int this_col_1 = r_j1[j1];
            LCL_Int::randperm(r_j2,this_m-1-j1,j1+1);
            for(int j2 = 0; (!found)&&(j2 < (this_m-1-j1)); j2++) {
                int this_col_2 = r_j2[j2];

                for(int i = 0; i < n; i++) {
                    x[i][0] = (A[i][this_col_1] + A[i][this_col_2])%2;
                }

                // Check chi.
                auto start_chi = std::chrono::high_resolution_clock::now();
                GateSynthesisMatrix::Chi(A,x,n,this_m,chi_A);
                auto end_chi = std::chrono::high_resolution_clock::now();
                auto duration_this_call = std::chrono::duration_cast<std::chrono::microseconds>(end_chi - start_chi);
                total_chi_duration += duration_this_call;

                //LCL_Mat_GF2::print(chi_A,n*n*n,this_m);
                int d = 0;
                //時間計測
                auto start_ns = std::chrono::high_resolution_clock::now();
                bool** NS = LCL_Mat_GF2::nullspace((const bool**)chi_A,n_chi_A,this_m,d);
                //bool** NS = ::M4RI_wrapper_for_nullspace((const bool**)chi_A,n_chi_A,this_m,d);
                auto end_ns = std::chrono::high_resolution_clock::now();
                duration_this_call = std::chrono::duration_cast<std::chrono::microseconds>(end_ns - start_ns);
                g_total_nullspace_duration += duration_this_call;

                //LCL_Mat_GF2::print(NS,this_m,d,"NS: ");
                found = 0;
                int nsv = -1;
                for(int h = 0; (!found)&&(h < d); h++) {
                    found = (NS[this_col_1][h]+NS[this_col_2][h])%2;
                    if(found) {
                        nsv = h;
                        //cout << "j1 = " << this_col_1 << ", j2 = " << this_col_2 << endl;
                    }
                }
                if(found) {
                    //for(int j = 0; j < this_m; j++) y[j][0] = NS[j][nsv];
                    for(int i = 0; i < n; i++) {
                        for(int j = 0; j < this_m; j++) {
                            int temp = (A[i][j] + x[i][0]*NS[j][nsv])%2;
                            Anew[i][j] = temp;
                        }
                    }
                    int mp;
                    GateSynthesisMatrix::cleanup(Anew,n,this_m,mp);
                    //this_m = mp;
                    if(mp<m_best) {
                        LCL_Mat_GF2::copy((const bool**)Anew,n,mp,Abest);
                        m_best = mp;
                    }
                } else {
                    if(false/*Conditions for chi prime*/) {
                        // Do chi prime
                    } else {
                        // Do single column
                    }
                }
                if(d) LCL_Mat_GF2::destruct(NS,this_m,d);
            }
        }
        LCL_Mat_GF2::copy((const bool**)Abest,n,m_best,A);
        this_m = m_best;
        round++;
    }

    std::cout << "nullspace time: "
          << std::chrono::duration_cast<std::chrono::milliseconds>(g_total_nullspace_duration).count()
          << " ms" 
          << ", chi time: "
          << std::chrono::duration_cast<std::chrono::milliseconds>(total_chi_duration).count()
          << " ms"
          << std::endl;

    delete [] r_j1; r_j1 = NULL;
    delete [] r_j2; r_j2 = NULL;
    LCL_Mat_GF2::destruct(x,n,1);
    LCL_Mat_GF2::destruct(chi_A,n_chi_A,m+1);
    LCL_Mat_GF2::destruct(y,m+1,1);
    LCL_Mat_GF2::destruct(xyT,n,m+1);
    LCL_Mat_GF2::destruct(A_xyT,n,m+1);
    LCL_Mat_GF2::destruct(Anew,n,m+1);
    LCL_Mat_GF2::destruct(Abest,n,m+1);
    omp = this_m;

    //
    end_total = std::chrono::high_resolution_clock::now();
    duration_total = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
    std::cout << "Total LempelX2 time: " << duration_total.count() << " ms" << endl;
}

void GateSynthesisMatrix::LempelX2_M4RI(bool** A, int n, int m, int& omp) {
    auto start_total = std::chrono::high_resolution_clock::now();
    std::cout << "in LempelX2_M4RIaaa" << endl;
    int this_m = m;
    int initial_total_m = m; // ★追加: 初期の列数を保持

    // 小さい行列は LCL のままで
    bool** x = LCL_Mat_GF2::construct(n, 1);
    bool** y = LCL_Mat_GF2::construct(m + 1, 1);
    bool** xyT = LCL_Mat_GF2::construct(n, m + 1);
    bool** A_xyT = LCL_Mat_GF2::construct(n, m + 1);
    
    int n_chi_A = n * n * n;
    
    // chi_A は mzd_t* (M4RI型) で作成
    mzd_t* chi_A = mzd_init(n_chi_A, m + 1);
    mzd_t* A_m4ri = convert_to_mzd((bool const**)A, n, m + 1);

    std::chrono::microseconds g_total_nullspace_duration(0);
    std::chrono::microseconds total_chi_duration(0);

    int* r_j1 = new int[m];
    int* r_j2 = new int[m];

    bool** Anew = LCL_Mat_GF2::construct(n, m + 1);
    LCL_Mat_GF2::copy((const bool**)A, n, m, Anew);

    bool** Abest = LCL_Mat_GF2::construct(n, m + 1);
    LCL_Mat_GF2::copy((const bool**)A, n, m, Abest);
    int m_best = m;

    bool found = 1;
    int round = 0;

    //零空間を求めた回数
    long long total_attempts = 0;
    //成功した回数
    long long successful_attempts = 0;

    while(found && (round < m)) {
        found = 0;
        // ★ラウンドごとの統計用カウンタ
        long long round_attempts = 0;
        long long round_successful_attempts = 0;

        int m_at_round_start = this_m; // ★ ラウンド開始時の列数を記録

        LOut(); cout << "Round = " << round << " (Start m: " << this_m << ")" << endl;
        
        LCL_Int::randperm(r_j1, this_m - 1);
        for(int j1 = 0; (!found) && (j1 < (this_m - 1)); j1++) {
            int this_col_1 = r_j1[j1];
            LCL_Int::randperm(r_j2, this_m - 1 - j1, j1 + 1);
            for(int j2 = 0; (!found) && (j2 < (this_m - 1 - j1)); j2++) {
                int this_col_2 = r_j2[j2];

                for(int i = 0; i < n; i++) {
                    x[i][0] = (A[i][this_col_1] + A[i][this_col_2]) % 2;
                }

                // mzd_t は再利用時にゼロクリアが必要 (簡易実装)
                for(int r=0; r<chi_A->nrows; r++) {
                      for(int c=0; c<chi_A->width; c++) mzd_row(chi_A, r)[c] = 0;
                }

                // Chi行列の作成 (M4RI版)
                auto start_chi = std::chrono::high_resolution_clock::now();
                GateSynthesisMatrix::Chi_M4RI(A_m4ri, x, n, this_m, chi_A);
                auto end_chi = std::chrono::high_resolution_clock::now();
                auto duration_this_call = std::chrono::duration_cast<std::chrono::microseconds>(end_chi - start_chi);
                total_chi_duration += duration_this_call;

                int d = 0;

                auto start_ns = std::chrono::high_resolution_clock::now();
                bool** NS = M4RI_direct_nullspace(chi_A, d);

                auto end_ns = std::chrono::high_resolution_clock::now();
                duration_this_call = std::chrono::duration_cast<std::chrono::microseconds>(end_ns - start_ns);
                g_total_nullspace_duration += duration_this_call;

                // found = 0; // ループ条件ですでに0になっているので不要だが念のため
                int nsv = -1;

                round_attempts++;
                //null空間を求めた回数をカウント
                total_attempts++;

                for(int h = 0; (!found) && (h < d); h++) {
                    found = (NS[this_col_1][h] + NS[this_col_2][h]) % 2;
                    if(found) {
                        nsv = h;
                        //null空間が成功した回数をカウント
                        successful_attempts++;
                        round_successful_attempts++;
                    }
                }
                
                if(found) {
                    for(int i = 0; i < n; i++) {
                        for(int j = 0; j < this_m; j++) {
                            int temp = (A[i][j] + x[i][0] * NS[j][nsv]) % 2;
                            Anew[i][j] = temp;
                        }
                    }

                    int mp;
                    // int before_cleanup_m = this_m; // 既に m_at_round_start があるのでここは計算用
                    GateSynthesisMatrix::cleanup(Anew, n, this_m, mp);
                    
                    if(mp < m_best) {
                        // ★ ヒット時の瞬間的な削除数を表示
                        // std::cout << "  [HIT!] Deleted: " << (this_m - mp) << " cols" << std::endl;
                        LCL_Mat_GF2::copy((const bool**)Anew, n, mp, Abest);
                        m_best = mp;
                    }else{
                        std::cout << "No improvement: " << mp << " >= " << m_best << std::endl;
                    }
                }
                
                if(NS) LCL_Mat_GF2::destruct(NS, this_m, d);
            }
        }

        // ★ラウンドごとの集計表示
        if (round_attempts > 0) {
            double hit_rate = (double)round_successful_attempts / round_attempts * 100.0;
            int deleted_this_round = m_at_round_start - m_best; // このラウンドでの削減数
            
            std::cout << ">> Round " << round << " Summary:" << std::endl;
            std::cout << "   Hit Rate: " << hit_rate << "% (" << round_successful_attempts << "/" << round_attempts << ")" << std::endl;
            std::cout << "   Reduced : " << deleted_this_round << " cols (" << m_at_round_start << " -> " << m_best << ")" << std::endl;
        }

        LCL_Mat_GF2::copy((const bool**)Abest, n, m_best, A);

        // ★追加: A_m4ri もベストな解で同期 (全体を書き直す)
        for(int r=0; r<n; r++) {
            for(int c=0; c<m_best; c++) {
                mzd_write_bit(A_m4ri, r, c, A[r][c]);
            }
            // 残りの列は0にしておく（次のラウンドのため）
            for(int c=m_best; c<m+1; c++) {
                mzd_write_bit(A_m4ri, r, c, 0);
            }
        }
        this_m = m_best;
        round++;
    }

    // ★ 最終結果の表示
    auto end_total = std::chrono::high_resolution_clock::now();
    auto total_dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
    std::cout << "\n=== Optimization Summary ===" << std::endl;
    std::cout << "Algorithm       : LempelX2_M4RI" << std::endl;
    std::cout << "Initial T-count : " << initial_total_m << std::endl;
    std::cout << "Final T-count   : " << this_m << std::endl;
    std::cout << "Total Reduced   : " << (initial_total_m - this_m) << " gates" << std::endl;
    std::cout << "Execution Time  : " << total_dur.count() << " ms" << std::endl;
    std::cout << "Chi calculation : " << total_chi_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Nullspace calc  : " << g_total_nullspace_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Total attempts  : " << total_attempts << " (Success: " << successful_attempts << ")" << std::endl;
    std::cout << "============================" << std::endl;

    delete [] r_j1; 
    delete [] r_j2;
    LCL_Mat_GF2::destruct(x, n, 1);
    
    mzd_free(chi_A); // mzd_freeで解放
    mzd_free(A_m4ri); // ★解放忘れずに
    
    LCL_Mat_GF2::destruct(y, m + 1, 1);
    LCL_Mat_GF2::destruct(xyT, n, m + 1);
    LCL_Mat_GF2::destruct(A_xyT, n, m + 1);
    LCL_Mat_GF2::destruct(Anew, n, m + 1);
    LCL_Mat_GF2::destruct(Abest, n, m + 1);
    omp = this_m;
    cout << "END OF LEMPELX2" << endl;
}

void GateSynthesisMatrix::LempelX2_M4RI_DetailedStats(bool** A, int n, int m, int& omp) {
    std::cout << "--- Starting LempelX2_M4RI_DetailedStats ---" << std::endl;
    int this_m = m;
    int initial_total_m = m; 

    // テンポラリ行列の構築
    bool** x = LCL_Mat_GF2::construct(n, 1);
    bool** Anew = LCL_Mat_GF2::construct(n, m + 1);
    bool** Abest = LCL_Mat_GF2::construct(n, m + 1);
    LCL_Mat_GF2::copy((const bool**)A, n, m, Abest);
    int m_best = m;

    // M4RI構造の初期化
    int n_chi_A = n * n * n;
    mzd_t* chi_A = mzd_init(n_chi_A, m + 1);
    mzd_t* A_m4ri = convert_to_mzd((bool const**)A, n, m + 1);

    // 統計・時間計測用
    std::chrono::microseconds total_ns_duration(0);
    std::chrono::microseconds total_chi_duration(0);
    long long total_attempts = 0;
    long long total_successful_finds = 0;

    int* r_j1 = new int[m];
    int* r_j2 = new int[m];

    bool improved_in_round = true;
    int round = 0;

    while (improved_in_round && (round < m)) {
        improved_in_round = false;
        int m_at_round_start = this_m;
        long long round_attempts = 0;

        std::cout << "Round " << round << " (Current m: " << this_m << ")" << std::endl;
        
        LCL_Int::randperm(r_j1, this_m - 1);
        for (int i_idx = 0; (!improved_in_round) && (i_idx < (this_m - 1)); i_idx++) {
            int c1 = r_j1[i_idx];
            LCL_Int::randperm(r_j2, this_m - 1 - i_idx, i_idx + 1);
            
            for (int j_idx = 0; (!improved_in_round) && (j_idx < (this_m - 1 - i_idx)); j_idx++) {
                int c2 = r_j2[j_idx];

                // 1. xベクトル作成 & ハミング距離計算
                int dist = 0;
                for (int k = 0; k < n; k++) {
                    x[k][0] = (A[k][c1] + A[k][c2]) % 2;
                    if (x[k][0]) dist++;
                }

                // 2. Chi行列構築 (M4RI)
                for (int r = 0; r < chi_A->nrows; r++) {
                    mzd_row_clear_offset(chi_A, r, 0); // 高速な行クリア
                }
                
                auto t1 = std::chrono::high_resolution_clock::now();
                GateSynthesisMatrix::Chi_M4RI(A_m4ri, x, n, this_m, chi_A);
                auto t2 = std::chrono::high_resolution_clock::now();
                total_chi_duration += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);

                // 3. 零空間計算
                int d = 0;
                auto t3 = std::chrono::high_resolution_clock::now();
                bool** NS = M4RI_direct_nullspace(chi_A, d);
                auto t4 = std::chrono::high_resolution_clock::now();
                total_ns_duration += std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3);

                round_attempts++;
                total_attempts++;

                // 4. 有効な削減ベクトルの探索
                int nsv = -1;
                bool found_vec = false;
                for (int h = 0; h < d; h++) {
                    if ((NS[c1][h] + NS[c2][h]) % 2 == 1) {
                        nsv = h;
                        found_vec = true;
                        break;
                    }
                }

                if (found_vec) {
                    // 行列更新
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < this_m; j++) {
                            Anew[i][j] = (A[i][j] + x[i][0] * NS[j][nsv]) % 2;
                        }
                    }

                    int mp;
                    GateSynthesisMatrix::cleanup(Anew, n, this_m, mp);

                    if (mp < m_best) {
                        // ★ 削減成功時の表示
                        std::cout << "  [HIT!] Dist=" << dist 
                                  << " | Reduced: " << (m_best - mp) << " cols (" 
                                  << m_best << " -> " << mp << ")" << std::endl;

                        LCL_Mat_GF2::copy((const bool**)Anew, n, mp, Abest);
                        m_best = mp;
                        improved_in_round = true;
                        total_successful_finds++;
                    }
                }
                if (NS) LCL_Mat_GF2::destruct(NS, this_m, d);
            }
        }

        if (improved_in_round) {
            // A と A_m4ri を最新の状態に同期
            LCL_Mat_GF2::copy((const bool**)Abest, n, m_best, A);
            for (int r = 0; r < n; r++) {
                for (int c = 0; c < m_best; c++) mzd_write_bit(A_m4ri, r, c, A[r][c]);
                for (int c = m_best; c < m + 1; c++) mzd_write_bit(A_m4ri, r, c, 0);
            }
            this_m = m_best;
        }
        round++;
    }

    // --- 最終レポート ---
    std::cout << "\n====================================" << std::endl;
    std::cout << "Final Results for " << __func__ << std::endl;
    std::cout << "  Initial T-count : " << initial_total_m << std::endl;
    std::cout << "  Final T-count   : " << this_m << std::endl;
    std::cout << "  Total Reduction : " << (initial_total_m - this_m) << " gates" << std::endl;
    std::cout << "  Total Attempts  : " << total_attempts << std::endl;
    std::cout << "  Chi Time        : " << std::chrono::duration_cast<std::chrono::milliseconds>(total_chi_duration).count() << " ms" << std::endl;
    std::cout << "  Nullspace Time  : " << std::chrono::duration_cast<std::chrono::milliseconds>(total_ns_duration).count() << " ms" << std::endl;
    std::cout << "====================================\n" << std::endl;

    // 後片付け
    delete[] r_j1; delete[] r_j2;
    LCL_Mat_GF2::destruct(x, n, 1);
    LCL_Mat_GF2::destruct(Anew, n, m + 1);
    LCL_Mat_GF2::destruct(Abest, n, m + 1);
    mzd_free(chi_A);
    mzd_free(A_m4ri);
    omp = this_m;
}

// ペアと距離を管理する構造体
struct ColPair {
    int c1;
    int c2;
    int dist;
    // 距離が小さい順にソートするための演算子オーバーロード
    bool operator<(const ColPair& other) const {
        return dist < other.dist;
    }
};

/*
void GateSynthesisMatrix::LempelX2_M4RI_Hamming(bool** A, int n, int m, int& omp) {
    std::cout << "in LempelX2_M4RI_Hamming" << endl;
    int this_m = m;
    int initial_total_m = m;

    // 小さい行列は LCL のままで
    bool** x = LCL_Mat_GF2::construct(n, 1);
    // y, xyT, A_xyT は使用していないようなので省略しても良いですが、念のため残すならそのままで
    
    int n_chi_A = n * n * n;
    mzd_t* chi_A = mzd_init(n_chi_A, m + 1);
    mzd_t* A_m4ri = convert_to_mzd((bool const**)A, n, m + 1);

    std::chrono::microseconds g_total_nullspace_duration(0);
    std::chrono::microseconds total_chi_duration(0);
    std::chrono::microseconds total_hamming_duration(0);

    // ランダム用配列は不要になったので削除またはコメントアウト
    // int* r_j1 = new int[m];
    // int* r_j2 = new int[m];

    bool** Anew = LCL_Mat_GF2::construct(n, m + 1);
    LCL_Mat_GF2::copy((const bool**)A, n, m, Anew);

    bool** Abest = LCL_Mat_GF2::construct(n, m + 1);
    LCL_Mat_GF2::copy((const bool**)A, n, m, Abest);
    int m_best = m;

    bool found = 1;
    int round = 0;

    long long total_attempts = 0;
    long long successful_attempts = 0;

    // ★ ペア候補を格納するベクター（メモリ確保の効率化のため外で宣言）
    std::vector<ColPair> candidates;

    while(found && (round < m)) {
        found = 0;
        long long round_attempts = 0;
        long long round_successful_attempts = 0;
        int m_at_round_start = this_m;

        LOut(); cout << "Round = " << round << " (Start m: " << this_m << ")" << endl;

        // ---------------------------------------------------------
        // ★ 変更点: ランダムではなく、ハミング距離を計算してソートする
        // ---------------------------------------------------------
        auto start_hamming = std::chrono::high_resolution_clock::now();
        candidates.clear();
        candidates.reserve(this_m * (this_m - 1) / 2); // メモリ予約

        // 統計用マップ（距離ごとの出現回数）
        std::map<int, int> dist_stats;

        // 全ペアのハミング距離を計算
        for (int j1 = 0; j1 < this_m; ++j1) {
            for (int j2 = j1 + 1; j2 < this_m; ++j2) {
                int dist = 0;
                // n行分の排他的論理和(不一致)をカウント
                for (int k = 0; k < n; ++k) {
                    if (A[k][j1] != A[k][j2]) {
                        dist++;
                    }
                }
                candidates.push_back({j1, j2, dist});
                dist_stats[dist]++; // ★ ここを追加してください
            }
        }

        // 距離が小さい順にソート (昇順)
        std::sort(candidates.begin(), candidates.end());

        auto end_hamming = std::chrono::high_resolution_clock::now();
        
        total_hamming_duration += std::chrono::duration_cast<std::chrono::microseconds>(end_hamming - start_hamming);

        // ★ ハミング距離の一覧（統計情報）を表示
        std::cout << "--- Hamming Distance Stats ---" << std::endl;
        for (std::map<int, int>::const_iterator it = dist_stats.begin(); it != dist_stats.end(); ++it) {
            // it->first が dist（ハミング距離）、it->second が count（出現回数）に対応します
            std::cout << " Dist " << it->first << ": " << it->second << " pairs" << std::endl;
        }
        std::cout << "-------------------------------" << std::endl;

        // ソートされたペア順に試行
        for (const auto& pair : candidates) {
            if (found) break; // 見つかったら次のラウンドへ

            int this_col_1 = pair.c1;
            int this_col_2 = pair.c2;

            // -----------------------------------------------------
            // 以下、既存の処理
            // -----------------------------------------------------

            for(int i = 0; i < n; i++) {
                x[i][0] = (A[i][this_col_1] + A[i][this_col_2]) % 2;
            }

            for(int r=0; r<chi_A->nrows; r++) {
                    for(int c=0; c<chi_A->width; c++) mzd_row(chi_A, r)[c] = 0;
            }

            auto start_chi = std::chrono::high_resolution_clock::now();
            GateSynthesisMatrix::Chi_M4RI(A_m4ri, x, n, this_m, chi_A);
            auto end_chi = std::chrono::high_resolution_clock::now();
            total_chi_duration += std::chrono::duration_cast<std::chrono::microseconds>(end_chi - start_chi);

            int d = 0;
            auto start_ns = std::chrono::high_resolution_clock::now();
            bool** NS = M4RI_direct_nullspace(chi_A, d);
            auto end_ns = std::chrono::high_resolution_clock::now();
            g_total_nullspace_duration += std::chrono::duration_cast<std::chrono::microseconds>(end_ns - start_ns);

            round_attempts++;
            total_attempts++;

            int nsv = -1;
            for(int h = 0; (!found) && (h < d); h++) {
                found = (NS[this_col_1][h] + NS[this_col_2][h]) % 2;
                if(found) {
                    nsv = h;
                    successful_attempts++;
                    round_successful_attempts++;
                }
            }
            
            if(found) {
                for(int i = 0; i < n; i++) {
                    for(int j = 0; j < this_m; j++) {
                        int temp = (A[i][j] + x[i][0] * NS[j][nsv]) % 2;
                        Anew[i][j] = temp;
                    }
                }
                int mp;
                int before_cleanup_m = this_m;
                GateSynthesisMatrix::cleanup(Anew, n, this_m, mp);
                
                if(mp < m_best) {
                    std::cout << "  [HIT!] Dist=" << pair.dist << " | Deleted: " << (before_cleanup_m - mp) << " cols" << std::endl;
                    LCL_Mat_GF2::copy((const bool**)Anew, n, mp, Abest);
                    m_best = mp;
                } else {
                    std::cout << "No improvement: " << mp << " >= " << m_best << std::endl;
                }
            }
            
            if(NS) LCL_Mat_GF2::destruct(NS, this_m, d);
        }

        // ラウンド結果表示
        if (round_attempts > 0) {
            double hit_rate = (double)round_successful_attempts / round_attempts * 100.0;
            std::cout << ">> Round " << round << " Result: Hit Rate = " << hit_rate << "% "
                      << "(" << round_successful_attempts << "/" << round_attempts << ")" << std::endl;
        }

        LCL_Mat_GF2::copy((const bool**)Abest, n, m_best, A);
        
        // A_m4ri の同期
        for(int r=0; r<n; r++) {
            for(int c=0; c<m_best; c++) {
                mzd_write_bit(A_m4ri, r, c, A[r][c]);
            }
            for(int c=m_best; c<m+1; c++) {
                mzd_write_bit(A_m4ri, r, c, 0);
            }
        }
        this_m = m_best;
        round++;
    }

// ★ 最終的な時間の表示
    std::cout << "\n=== Final Performance Statistics ===" << std::endl;
    std::cout << "Hamming time:   " << std::chrono::duration_cast<std::chrono::milliseconds>(total_hamming_duration).count() << " ms" << std::endl;
    std::cout << "Chi time:       " << std::chrono::duration_cast<std::chrono::milliseconds>(total_chi_duration).count() << " ms" << std::endl;
    std::cout << "Nullspace time: " << std::chrono::duration_cast<std::chrono::milliseconds>(g_total_nullspace_duration).count() << " ms" << std::endl;
    std::cout << "Total attempts: " << total_attempts << " (Success: " << successful_attempts << ")" << std::endl;

    // delete [] r_j1; 
    // delete [] r_j2;
    LCL_Mat_GF2::destruct(x, n, 1);
    mzd_free(chi_A);
    mzd_free(A_m4ri);
    LCL_Mat_GF2::destruct(Anew, n, m + 1);
    LCL_Mat_GF2::destruct(Abest, n, m + 1);
    omp = this_m;
    cout << "END OF LEMPELX2" << endl;
}
    */

void GateSynthesisMatrix::LempelX2_M4RI_Hamming(bool** A, int n, int m, int& omp) {
    auto start_total = std::chrono::high_resolution_clock::now();
    std::cout << "in LempelX2_M4RI_Hamming　安定版　m4ri" << endl;
    int this_m = m;
    int initial_total_m = m;

    bool** x = LCL_Mat_GF2::construct(n, 1);
    
    int n_chi_A = n * n * n;
    // 最大サイズ (m + 1) でメモリを確保
    mzd_t* chi_A_full = mzd_init(n_chi_A, m + 1);
    mzd_t* A_m4ri_full = convert_to_mzd((bool const**)A, n, m + 1);

    std::chrono::microseconds g_total_nullspace_duration(0);
    std::chrono::microseconds total_chi_duration(0);
    std::chrono::microseconds total_hamming_duration(0);

    bool** Anew = LCL_Mat_GF2::construct(n, m + 1);
    LCL_Mat_GF2::copy((const bool**)A, n, m, Anew);

    bool** Abest = LCL_Mat_GF2::construct(n, m + 1);
    LCL_Mat_GF2::copy((const bool**)A, n, m, Abest);
    int m_best = m;

    bool found = 1;
    int round = 0;
    long long total_attempts = 0;
    long long successful_attempts = 0;

    std::vector<ColPair> candidates;

    while(found && (round < m)) {
        found = 0;
        long long round_attempts = 0;
        long long round_successful_attempts = 0;

        LOut() << "Round = " << round << " (Start m: " << this_m << ")" << endl;

        // --- ハミング距離の計算とソート ---
        auto start_hamming = std::chrono::high_resolution_clock::now();
        candidates.clear();
        candidates.reserve(this_m * (this_m - 1) / 2);
        std::map<int, int> dist_stats;

        for (int j1 = 0; j1 < this_m; ++j1) {
            for (int j2 = j1 + 1; j2 < this_m; ++j2) {
                int dist = 0;
                for (int k = 0; k < n; ++k) {
                    if (A[k][j1] != A[k][j2]) dist++;
                }
                candidates.push_back({j1, j2, dist});
                dist_stats[dist]++;
            }
        }
        std::sort(candidates.begin(), candidates.end());
        auto end_hamming = std::chrono::high_resolution_clock::now();
        total_hamming_duration += std::chrono::duration_cast<std::chrono::microseconds>(end_hamming - start_hamming);

        std::cout << "--- Hamming Distance Stats ---" << std::endl;
    for (std::map<int, int>::const_iterator it = dist_stats.begin(); it != dist_stats.end(); ++it) {
        std::cout << " Dist " << it->first << ": " << it->second << " pairs" << std::endl;
    }
        std::cout << "-------------------------------" << std::endl;

        // ★ mzd_init_window による次元管理
        // 0列目から this_m 列目まで（exclusive）のウィンドウを作成
        mzd_t* A_m4ri = mzd_init_window(A_m4ri_full, 0, 0, n, this_m);
        mzd_t* chi_A = mzd_init_window(chi_A_full, 0, 0, n_chi_A, this_m);

        for (const auto& pair : candidates) {
            if (found) break;

            int this_col_1 = pair.c1;
            int this_col_2 = pair.c2;

            for(int i = 0; i < n; i++) {
                x[i][0] = (A[i][this_col_1] + A[i][this_col_2]) % 2;
            }

            // ウィンドウに対してゼロクリア（内部のデータ chi_A_full もクリアされる）
            mzd_set_ui(chi_A, 0);

            auto start_chi = std::chrono::high_resolution_clock::now();
            GateSynthesisMatrix::Chi_M4RI(A_m4ri, x, n, this_m, chi_A);
            auto end_chi = std::chrono::high_resolution_clock::now();
            total_chi_duration += std::chrono::duration_cast<std::chrono::microseconds>(end_chi - start_chi);

            int d = 0;
            auto start_ns = std::chrono::high_resolution_clock::now();
            // ウィンドウを渡すことで、内部の M4RI 関数は this_m 次元として処理する
            bool** NS = M4RI_direct_nullspace(chi_A, d);
            auto end_ns = std::chrono::high_resolution_clock::now();
            g_total_nullspace_duration += std::chrono::duration_cast<std::chrono::microseconds>(end_ns - start_ns);

            round_attempts++;
            total_attempts++;

            int nsv = -1;
            for(int h = 0; (!found) && (h < d); h++) {
                found = (NS[this_col_1][h] + NS[this_col_2][h]) % 2;
                if(found) {
                    nsv = h;
                    successful_attempts++;
                    round_successful_attempts++;
                }
            }
            
            if(found) {
                for(int i = 0; i < n; i++) {
                    for(int j = 0; j < this_m; j++) {
                        Anew[i][j] = (A[i][j] + x[i][0] * NS[j][nsv]) % 2;
                    }
                }
                int mp;
                int before_cleanup_m = this_m;
                GateSynthesisMatrix::cleanup(Anew, n, this_m, mp);
                
                if(mp < m_best) {
                    std::cout << "  [HIT!] Dist=" << pair.dist << " | Deleted: " << (before_cleanup_m - mp) << " cols" << std::endl;
                    LCL_Mat_GF2::copy((const bool**)Anew, n, mp, Abest);
                    m_best = mp;
                } else {
                    std::cout << "No improvement: " << mp << " >= " << m_best << std::endl;
                }
            }
            if(NS) LCL_Mat_GF2::destruct(NS, this_m, d);
        }

        // ラウンド終了後にウィンドウを解放（実データは破棄されない）
        mzd_free_window(A_m4ri);
        mzd_free_window(chi_A);

        if (round_attempts > 0) {
            double hit_rate = (double)round_successful_attempts / round_attempts * 100.0;
            std::cout << ">> Round " << round << " Result: Hit Rate = " << hit_rate << "% "
                      << "(" << round_successful_attempts << "/" << round_attempts << ")" << std::endl;
        }

        // A (bool**) と A_m4ri_full (mzd_t) の同期
        LCL_Mat_GF2::copy((const bool**)Abest, n, m_best, A);
        for(int r=0; r<n; r++) {
            for(int c=0; c<m_best; c++) {
                mzd_write_bit(A_m4ri_full, r, c, A[r][c]);
            }
            // 余剰次元のクリア（任意ですが整合性のため）
            for(int c=m_best; c<m+1; c++) {
                mzd_write_bit(A_m4ri_full, r, c, 0);
            }
        }
        this_m = m_best;
        round++;
    }

    std::cout << "\n=== Final Performance Statistics ===" << std::endl;
    std::cout << "Hamming time:   " << std::chrono::duration_cast<std::chrono::milliseconds>(total_hamming_duration).count() << " ms" << std::endl;
    std::cout << "Chi time:       " << std::chrono::duration_cast<std::chrono::milliseconds>(total_chi_duration).count() << " ms" << std::endl;
    std::cout << "Nullspace time: " << std::chrono::duration_cast<std::chrono::milliseconds>(g_total_nullspace_duration).count() << " ms" << std::endl;
    std::cout << "Total attempts: " << total_attempts << " (Success: " << successful_attempts << ")" << std::endl;

    LCL_Mat_GF2::destruct(x, n, 1);
    mzd_free(chi_A_full);
    mzd_free(A_m4ri_full);
    LCL_Mat_GF2::destruct(Anew, n, m + 1);
    LCL_Mat_GF2::destruct(Abest, n, m + 1);
    omp = this_m;
    auto end_total = std::chrono::high_resolution_clock::now();
    auto total_dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
    std::cout << "\n=== Optimization Summary ===" << std::endl;
    std::cout << "Algorithm       : LempelX2_M4RI_Hamming" << std::endl;
    std::cout << "Initial T-count : " << initial_total_m << std::endl;
    std::cout << "Final T-count   : " << omp << std::endl;
    std::cout << "Total Reduced   : " << (initial_total_m - omp) << " gates" << std::endl;
    std::cout << "Execution Time  : " << total_dur.count() << " ms" << std::endl;
    std::cout << "Chi calculation : " << total_chi_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Nullspace calc  : " << g_total_nullspace_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "============================" << std::endl;
    cout << "END OF LEMPELX2 安定版" << endl;
}



// --- 履歴・状態管理用の構造体 ---
struct NodeHistory {
    int id, pid, round, c1, c2, dist, m, deltaM;
};

struct TODDState {
    bool** A;
    int m;
    int nodeID, parentID, round, c1, c2, dist, deltaM;

    TODDState(bool** matrix, int cols, int id, int pid, int r, int _c1, int _c2, int d, int dm) 
        : A(matrix), m(cols), nodeID(id), parentID(pid), round(r), c1(_c1), c2(_c2), dist(d), deltaM(dm) {}

    bool operator<(const TODDState& other) const {
        if (m != other.m) return m < other.m;
        return dist < other.dist; // タイブレーク：ハミング距離が短い方を優先
    }
};

// 各親ノードの探索進捗を管理するトラッカー
struct BeamParentTracker {
    int beamIdx;
    TODDState* state;
    std::vector<ColPair> pairs; 
    size_t nextPairIdx;
    mzd_t* A_m4ri;

    BeamParentTracker(int b, TODDState* s, std::vector<ColPair> p, mzd_t* m)
        : beamIdx(b), state(s), pairs(p), nextPairIdx(0), A_m4ri(m) {}
};

// --- 補助関数 ---

// 行列のディープコピー用
static bool** copy_matrix_local(bool** src, int n, int m, int max_m) {
    bool** dst = LCL_Mat_GF2::construct(n, max_m);
    LCL_Mat_GF2::copy((const bool**)src, n, m, dst);
    return dst;
}

// 最終的な削減パス（木のたどり方）を出力
static void print_optimization_path(int bestNodeID, const std::vector<NodeHistory>& history, int initial_m) {
    if (bestNodeID <= 0) return;
    std::vector<NodeHistory> path;
    int currentID = bestNodeID;

    while (currentID != 0) {
        bool found = false;
        for (const auto& h : history) {
            if (h.id == currentID) {
                path.push_back(h);
                currentID = h.pid;
                found = true;
                break;
            }
        }
        if (!found) break;
    }
    std::reverse(path.begin(), path.end());

    std::cout << "\n============================================================" << std::endl;
    std::cout << "Final Adopted Path (Diversity-Aware Beam Search)" << std::endl;
    std::cout << "============================================================" << std::endl;
    printf("%-8s %-15s %-10s %-10s %s\n", "Round", "Pair(c1,c2)", "Dist", "Delta", "Current m");
    std::cout << "------------------------------------------------------------" << std::endl;
    printf("%-8s %-15s %-10s %-10s %d\n", "Start", "-", "-", "-", initial_m);

    for (const auto& h : path) {
        if (h.id == 1 && h.round == 0 && h.deltaM == 0) continue; // Root dummy skip if any
        std::string pStr = "(" + std::to_string(h.c1) + "," + std::to_string(h.c2) + ")";
        printf("%-8d %-15s %-10d %-10d %d\n", h.round, pStr.c_str(), h.dist, h.deltaM, h.m);
    }
    std::cout << "============================================================\n" << std::endl;
}

// --- メイン関数：ラウンドロビン方式ビームサーチ ---
void GateSynthesisMatrix::LempelX2_M4RI_BeamSearch(bool** A_init, int n, int m_init, int& omp) {
    auto start_total = std::chrono::high_resolution_clock::now();
    int initial_m = m_init;
    int K = 5; int N = 50; int max_m = m_init + 1;
    int nodeCounter = 0;
    std::vector<NodeHistory> history;
    std::vector<TODDState> current_beam;

    // 初期ノード ID=1
    current_beam.push_back(TODDState(copy_matrix_local(A_init, n, m_init, max_m), m_init, ++nodeCounter, 0, -1, -1, -1, 0, 0));
    history.push_back({nodeCounter, 0, -1, -1, -1, 0, m_init, 0});

    bool improved = true; int round = 0;
    bool** Anew = LCL_Mat_GF2::construct(n, max_m);
    bool** x_vec = LCL_Mat_GF2::construct(n, 1);

    while (improved && round < m_init) {
        improved = false;
        std::vector<TODDState> next_candidates;
        std::cout << "--- Round " << round << " | Beam Size: " << current_beam.size() << " (Round-Robin Mode) ---" << std::endl;

        // 1. 各親ノードのトラッカー準備
        std::vector<BeamParentTracker> trackers;
        for (size_t i = 0; i < current_beam.size(); ++i) {
            TODDState& s = current_beam[i];
            std::vector<ColPair> p_list;
            for (int j1 = 0; j1 < s.m; ++j1) {
                for (int j2 = j1 + 1; j2 < s.m; ++j2) {
                    int d_val = 0;
                    for (int k = 0; k < n; ++k) if (s.A[k][j1] != s.A[k][j2]) d_val++;
                    p_list.push_back({j1, j2, d_val});
                }
            }
            std::sort(p_list.begin(), p_list.end()); // ハミング距離順
            trackers.emplace_back((int)i, &s, p_list, convert_to_mzd((bool const**)s.A, n, max_m));
        }

        mzd_t* chi_A = mzd_init(n * n * n, max_m);
        bool anyWorkLeft = true;

        // インターリーブ探索：各親から交互に1ペアずつ試行
        while (next_candidates.size() < (size_t)N && anyWorkLeft) {
            anyWorkLeft = false;
            for (auto& trk : trackers) {
                if (trk.nextPairIdx < trk.pairs.size()) {
                    anyWorkLeft = true;
                    ColPair& cp = trk.pairs[trk.nextPairIdx++];

                    for(int i = 0; i < n; i++) x_vec[i][0] = (trk.state->A[i][cp.c1] + trk.state->A[i][cp.c2]) % 2;
                    for(int r=0; r<chi_A->nrows; r++) mzd_row_clear_offset(chi_A, r, 0);
                    
                    GateSynthesisMatrix::Chi_M4RI(trk.A_m4ri, x_vec, n, trk.state->m, chi_A);
                    int d_ns = 0;
                    bool** NS = M4RI_direct_nullspace(chi_A, d_ns);
                    
                    bool found_vec = false;
                    int nsv = -1;
                    for(int h = 0; h < d_ns; h++) {
                        if ((NS[cp.c1][h] + NS[cp.c2][h]) % 2 == 1) { nsv = h; found_vec = true; break; }
                    }

                    if (found_vec) {
                        for(int i = 0; i < n; i++) {
                            for(int j = 0; j < trk.state->m; j++) Anew[i][j] = (trk.state->A[i][j] + x_vec[i][0] * NS[j][nsv]) % 2;
                        }
                        int mp;
                        GateSynthesisMatrix::cleanup(Anew, n, trk.state->m, mp);
                        if (mp < trk.state->m) {
                            int currentID = ++nodeCounter;
                            next_candidates.push_back(TODDState(copy_matrix_local(Anew, n, mp, max_m), mp, currentID, trk.state->nodeID, round, cp.c1, cp.c2, cp.dist, trk.state->m - mp));
                            history.push_back({currentID, trk.state->nodeID, round, cp.c1, cp.c2, cp.dist, mp, trk.state->m - mp});
                            improved = true;
                            if (next_candidates.size() >= (size_t)N) break;
                        }
                    }
                    if (NS) LCL_Mat_GF2::destruct(NS, trk.state->m, d_ns);
                    if (next_candidates.size() >= (size_t)N) break;
                }
            }
        }

        // 3. 次世代選別（Pruning）
        if (!next_candidates.empty()) {
            std::sort(next_candidates.begin(), next_candidates.end());
            
            for (auto& trk : trackers) {
                LCL_Mat_GF2::destruct(trk.state->A, n, max_m);
                mzd_free(trk.A_m4ri);
            }
            current_beam.clear();

            for (size_t i = 0; i < next_candidates.size(); ++i) {
                if (i < (size_t)K) current_beam.push_back(next_candidates[i]);
                else LCL_Mat_GF2::destruct(next_candidates[i].A, n, max_m);
            }
            std::cout << "  Round " << round << " Finished. Best m: " << current_beam[0].m << " (Found " << next_candidates.size() << " paths)" << std::endl;
        } else {
            for (auto& trk : trackers) mzd_free(trk.A_m4ri);
        }
        mzd_free(chi_A);
        round++;
    }

    if (!current_beam.empty()) {
        print_optimization_path(current_beam[0].nodeID, history, m_init);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < current_beam[0].m; ++j) A_init[i][j] = current_beam[0].A[i][j];
        }
        omp = current_beam[0].m;
        for (size_t i = 0; i < current_beam.size(); ++i) LCL_Mat_GF2::destruct(current_beam[i].A, n, max_m);
    }
    LCL_Mat_GF2::destruct(Anew, n, max_m);
    LCL_Mat_GF2::destruct(x_vec, n, 1);
    auto end_total = std::chrono::high_resolution_clock::now();
    auto total_dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
    std::cout << "\n=== Optimization Summary ===" << std::endl;
    std::cout << "Algorithm       : LempelX2_M4RI_BeamSearch" << std::endl;
    std::cout << "Initial T-count : " << initial_m << std::endl;
    std::cout << "Final T-count   : " << omp << std::endl;
    std::cout << "Total Reduced   : " << (initial_m - omp) << " gates" << std::endl;
    std::cout << "Execution Time  : " << total_dur.count() << " ms" << std::endl;
    std::cout << "============================" << std::endl;
}

// 2. Random Beam Search版（多様性重視の構成に準拠）
void GateSynthesisMatrix::LempelX2_M4RI_RandomBeamSearch(bool** A_init, int n, int m_init, int& omp) {
    auto start_total = std::chrono::high_resolution_clock::now();
    int initial_m = m_init;
    int K = 5; int N = 50; int max_m = m_init + 1;
    int nodeCounter = 0;
    std::vector<NodeHistory> history;
    std::vector<TODDState> current_beam;

    current_beam.push_back(TODDState(copy_matrix_local(A_init, n, m_init, max_m), m_init, ++nodeCounter, 0, -1, -1, -1, 0, 0));
    history.push_back({nodeCounter, 0, -1, -1, -1, 0, m_init, 0});

    std::random_device rd; std::mt19937 g(rd());
    bool improved = true; int round = 0;
    bool** Anew = LCL_Mat_GF2::construct(n, max_m);
    bool** x_vec = LCL_Mat_GF2::construct(n, 1);

    while (improved && round < m_init) {
        improved = false;
        std::vector<TODDState> next_candidates;
        std::cout << "--- Round " << round << " | Beam Size: " << current_beam.size() << " (Random Mode) ---" << std::endl;

        bool limit_reached = false;
        for (size_t b = 0; b < current_beam.size() && !limit_reached; ++b) {
            TODDState& state = current_beam[b];
            std::vector<ColPair> candidates;
            for (int j1 = 0; j1 < state.m; ++j1) {
                for (int j2 = j1 + 1; j2 < state.m; ++j2) { candidates.push_back({j1, j2, 0}); }
            }
            std::shuffle(candidates.begin(), candidates.end(), g);

            mzd_t* A_m4ri = convert_to_mzd((bool const**)state.A, n, max_m);
            mzd_t* chi_A = mzd_init(n * n * n, max_m);

            for (size_t p = 0; p < candidates.size(); ++p) {
                int c1 = candidates[p].c1; int c2 = candidates[p].c2;
                int dist = 0;
                for(int i = 0; i < n; i++) {
                    x_vec[i][0] = (state.A[i][c1] + state.A[i][c2]) % 2;
                    if (x_vec[i][0]) dist++;
                }

                for(int r=0; r<chi_A->nrows; r++) mzd_row_clear_offset(chi_A, r, 0);
                GateSynthesisMatrix::Chi_M4RI(A_m4ri, x_vec, n, state.m, chi_A);
                int d_ns = 0; bool** NS = M4RI_direct_nullspace(chi_A, d_ns);
                int nsv = -1; bool found = false;
                for(int h = 0; (!found) && (h < d_ns); h++) {
                    if ((NS[c1][h] + NS[c2][h]) % 2 == 1) { nsv = h; found = true; }
                }

                if (found) {
                    for(int i = 0; i < n; i++) { for(int j = 0; j < state.m; j++) Anew[i][j] = (state.A[i][j] + x_vec[i][0] * NS[j][nsv]) % 2; }
                    int mp; GateSynthesisMatrix::cleanup(Anew, n, state.m, mp);
                    if (mp < state.m) {
                        int currentID = ++nodeCounter;
                        next_candidates.push_back(TODDState(copy_matrix_local(Anew, n, mp, max_m), mp, currentID, state.nodeID, round, c1, c2, dist, state.m - mp));
                        history.push_back({currentID, state.nodeID, round, c1, c2, dist, mp, state.m - mp});
                        improved = true;
                        if (next_candidates.size() >= (size_t)N) limit_reached = true;
                    }
                }
                if (NS) LCL_Mat_GF2::destruct(NS, state.m, d_ns);
                if (limit_reached) break; 
            }
            mzd_free(chi_A); mzd_free(A_m4ri);
        }

        if (!next_candidates.empty()) {
            std::sort(next_candidates.begin(), next_candidates.end());
            for (size_t i = 0; i < current_beam.size(); ++i) LCL_Mat_GF2::destruct(current_beam[i].A, n, max_m);
            current_beam.clear();
            for (size_t i = 0; i < next_candidates.size(); ++i) {
                if (i < (size_t)K) current_beam.push_back(next_candidates[i]);
                else LCL_Mat_GF2::destruct(next_candidates[i].A, n, max_m);
            }
            std::cout << "  Round " << round << " Best m: " << current_beam[0].m << std::endl;
        }
        round++;
    }

    if (!current_beam.empty()) {
        print_optimization_path(current_beam[0].nodeID, history, m_init);
        for (int i = 0; i < n; ++i) { for (int j = 0; j < current_beam[0].m; ++j) A_init[i][j] = current_beam[0].A[i][j]; }
        omp = current_beam[0].m;
        for (size_t i = 0; i < current_beam.size(); ++i) LCL_Mat_GF2::destruct(current_beam[i].A, n, max_m);
    }
    LCL_Mat_GF2::destruct(Anew, n, max_m);
    LCL_Mat_GF2::destruct(x_vec, n, 1);
    auto end_total = std::chrono::high_resolution_clock::now();
    auto total_dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
    std::cout << "\n=== Optimization Summary ===" << std::endl;
    std::cout << "Algorithm       : LempelX2_M4RI_RandomBeamSearch" << std::endl;
    std::cout << "Initial T-count : " << initial_m << std::endl;
    std::cout << "Final T-count   : " << omp << std::endl;
    std::cout << "Total Reduced   : " << (initial_m - omp) << " gates" << std::endl;
    std::cout << "Execution Time  : " << total_dur.count() << " ms" << std::endl;
    std::cout << "============================" << std::endl;
}


void GateSynthesisMatrix::LempelX2_M4RI_SequentialBeamSearch(bool** A_init, int n, int m_init, int& omp) {
    auto start_total = std::chrono::high_resolution_clock::now();
    int initial_m = m_init;
    int K = 5;  // ビーム幅
    int N = 50; // 候補数上限（これに達したら次の親ノードには行かずに打ち切り）
    int max_m = m_init + 1;
    int nodeCounter = 0;
    std::vector<NodeHistory> history;
    std::vector<TODDState> current_beam;

    // 初期状態を登録
    current_beam.push_back(TODDState(copy_matrix_local(A_init, n, m_init, max_m), m_init, ++nodeCounter, 0, -1, -1, -1, 0, 0));
    history.push_back({nodeCounter, 0, -1, -1, -1, 0, m_init, 0});

    bool improved = true;
    int round = 0;
    bool** Anew = LCL_Mat_GF2::construct(n, max_m);
    bool** x_vec = LCL_Mat_GF2::construct(n, 1);

    while (improved && round < m_init) {
        improved = false;
        std::vector<TODDState> next_candidates;
        std::cout << "--- Round " << round << " | Beam Size: " << current_beam.size() << " (Sequential Mode) ---" << std::endl;

        bool limit_reached = false;

        // 【深さ優先的な挙動】親ノードを順番に走査し、合計がNに達したら即終了
        for (size_t b = 0; b < current_beam.size() && !limit_reached; ++b) {
            TODDState& state = current_beam[b];
            
            // カラムペアと距離の算出
            std::vector<ColPair> candidates;
            for (int j1 = 0; j1 < state.m; ++j1) {
                for (int j2 = j1 + 1; j2 < state.m; ++j2) {
                    int d_val = 0;
                    for (int k = 0; k < n; ++k) if (state.A[k][j1] != state.A[k][j2]) d_val++;
                    candidates.push_back({j1, j2, d_val});
                }
            }
            std::sort(candidates.begin(), candidates.end()); // ハミング距離順

            mzd_t* A_m4ri = convert_to_mzd((bool const**)state.A, n, max_m);
            mzd_t* chi_A = mzd_init(n * n * n, max_m);

            for (size_t p = 0; p < candidates.size(); ++p) {
                int c1 = candidates[p].c1; int c2 = candidates[p].c2; int cur_dist = candidates[p].dist;

                for(int i = 0; i < n; i++) x_vec[i][0] = (state.A[i][c1] + state.A[i][c2]) % 2;
                for(int r=0; r<chi_A->nrows; r++) mzd_row_clear_offset(chi_A, r, 0);

                GateSynthesisMatrix::Chi_M4RI(A_m4ri, x_vec, n, state.m, chi_A);
                int d_ns = 0;
                bool** NS = M4RI_direct_nullspace(chi_A, d_ns);
                
                int nsv = -1; bool found_vec = false;
                for(int h = 0; h < d_ns; h++) {
                    if ((NS[c1][h] + NS[c2][h]) % 2 == 1) { nsv = h; found_vec = true; break; }
                }

                if (found_vec) {
                    for(int i = 0; i < n; i++) {
                        for(int j = 0; j < state.m; j++) Anew[i][j] = (state.A[i][j] + x_vec[i][0] * NS[j][nsv]) % 2;
                    }
                    int mp;
                    GateSynthesisMatrix::cleanup(Anew, n, state.m, mp);
                    if (mp < state.m) {
                        int currentID = ++nodeCounter;
                        // 削減案が見つかるたびに表示
                        std::cout << "  [HIT!] Beam[" << b << "] Dist=" << cur_dist 
                                  << " | Deleted: " << (state.m - mp) << " cols" << std::endl;

                        next_candidates.push_back(TODDState(copy_matrix_local(Anew, n, mp, max_m), mp, currentID, state.nodeID, round, c1, c2, cur_dist, state.m - mp));
                        history.push_back({currentID, state.nodeID, round, c1, c2, cur_dist, mp, state.m - mp});
                        improved = true;

                        // 合計が50個に達したら、現在の親の探索も、次の親の探索も打ち切り
                        if (next_candidates.size() >= (size_t)N) {
                            limit_reached = true;
                            break; 
                        }
                    }
                }
                if (NS) LCL_Mat_GF2::destruct(NS, state.m, d_ns);
            }
            mzd_free(chi_A);
            mzd_free(A_m4ri);
        }

        // 次世代の選別
        if (!next_candidates.empty()) {
            std::sort(next_candidates.begin(), next_candidates.end());
            
            // 古いビームを解放
            for (size_t i = 0; i < current_beam.size(); ++i) LCL_Mat_GF2::destruct(current_beam[i].A, n, max_m);
            current_beam.clear();

            // 新しいエリート5つを採用
            for (size_t i = 0; i < next_candidates.size(); ++i) {
                if (i < (size_t)K) current_beam.push_back(next_candidates[i]);
                else LCL_Mat_GF2::destruct(next_candidates[i].A, n, max_m);
            }
            std::cout << "  Round " << round << " Finished. Best m: " << current_beam[0].m << " (Found " << next_candidates.size() << " candidates)" << std::endl;
        }
        round++;
    }

    // 最終報告
    if (!current_beam.empty()) {
        print_optimization_path(current_beam[0].nodeID, history, m_init);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < current_beam[0].m; ++j) A_init[i][j] = current_beam[0].A[i][j];
        }
        omp = current_beam[0].m;
        for (size_t i = 0; i < current_beam.size(); ++i) LCL_Mat_GF2::destruct(current_beam[i].A, n, max_m);
    }

    LCL_Mat_GF2::destruct(Anew, n, max_m);
    LCL_Mat_GF2::destruct(x_vec, n, 1);
    auto end_total = std::chrono::high_resolution_clock::now();
    auto total_dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
    std::cout << "\n=== Optimization Summary ===" << std::endl;
    std::cout << "Algorithm       : LempelX2_M4RI_SequentialBeamSearch" << std::endl;
    std::cout << "Initial T-count : " << initial_m << std::endl;
    std::cout << "Final T-count   : " << omp << std::endl;
    std::cout << "Total Reduced   : " << (initial_m - omp) << " gates" << std::endl;
    std::cout << "Execution Time  : " << total_dur.count() << " ms" << std::endl;
    std::cout << "============================" << std::endl;
}

void GateSynthesisMatrix::SparsifyAndTrack_Bool(bool** A, int n, int m, std::vector<std::pair<int, int>>& cnot_history) {
    int pivot_row = 0;
    for (int j = 0; j < m && pivot_row < n; ++j) {
        int sel = -1;
        for (int i = pivot_row; i < n; ++i) {
            if (A[i][j]) { sel = i; break; }
        }
        if (sel != -1) {
            if (sel != pivot_row) {
                for (int c = 0; c < m; ++c) {
                    bool temp = A[sel][c];
                    A[sel][c] = A[pivot_row][c];
                    A[pivot_row][c] = temp;
                }
                cnot_history.push_back({sel, pivot_row});
                cnot_history.push_back({pivot_row, sel});
                cnot_history.push_back({sel, pivot_row});
            }
            for (int i = 0; i < n; ++i) {
                if (i != pivot_row && A[i][j]) {
                    for (int c = 0; c < m; ++c) {
                        A[i][c] = (A[i][c] ^ A[pivot_row][c]);
                    }
                    cnot_history.push_back({pivot_row, i});
                }
            }
            pivot_row++;
        }
    }
}

void GateSynthesisMatrix::LempelX2_M4RI_Hamming_Preprocess(bool** A, int n, int m, int& omp) {
    auto start_total = std::chrono::high_resolution_clock::now();
    std::cout << "\n[Improved TODD] Starting with Row Preprocessing + M4RI Hamming Search" << std::endl;
    int this_m = m; int initial_m = m; int n_chi_A = n * n * n;

    auto print_density_info = [&](const std::string& prefix, bool** M, int r_max, int c_max) {
        long long ones = 0;
        long long total = (long long)r_max * c_max;
        for (int r = 0; r < r_max; ++r) for (int c = 0; c < c_max; ++c) if (M[r][c]) ones++;
        double percent = (total > 0) ? (100.0 * ones / total) : 0.0;
        std::cout << prefix << " Density: " << ones << " ones (" << ((double)((long long)(percent * 100)) / 100.0) << "%)." << std::endl;
    };

    std::cout << "[Initial] Matrix size: " << n << "x" << m << std::endl;
    print_density_info("[Initial]", A, n, m);
    std::vector<std::pair<int, int>> pre_cnots;
    SparsifyAndTrack_Bool(A, n, m, pre_cnots);
    print_density_info("[Post-Sparsify]", A, n, m);

    mzd_t* A_m4ri_full = convert_to_mzd((bool const**)A, n, m + 1);

    bool** x_vec = LCL_Mat_GF2::construct(n, 1);
    bool** Anew = LCL_Mat_GF2::construct(n, m + 1);
    bool** Abest = LCL_Mat_GF2::construct(n, m + 1);
    LCL_Mat_GF2::copy((const bool**)A, n, m, Abest);
    int m_best = m;
    mzd_t* chi_A_full = mzd_init(n_chi_A, m + 1);
    std::chrono::microseconds total_ns_duration{0};
    std::chrono::microseconds total_chi_duration{0};
    bool found = true; int round = 0;

    while (found && (round < m)) {
        found = false;
        std::cout << "--- Round " << round << " | Current Columns: " << this_m << " ---" << std::endl;
        auto round_start = std::chrono::high_resolution_clock::now();

        std::vector<ColPair> candidates;
        for (int j1 = 0; j1 < this_m; ++j1) {
            for (int j2 = j1 + 1; j2 < this_m; ++j2) {
                int dist = 0;
                for (int k = 0; k < n; ++k) if (A[k][j1] != A[k][j2]) dist++;
                candidates.push_back({j1, j2, dist});
            }
        }
        std::sort(candidates.begin(), candidates.end());

        mzd_t* A_m4ri_win = mzd_init_window(A_m4ri_full, 0, 0, n, this_m);
        mzd_t* chi_A_win = mzd_init_window(chi_A_full, 0, 0, n_chi_A, this_m);

        for (const auto& pair : candidates) {
            if (found) break;
            for (int i = 0; i < n; i++) x_vec[i][0] = (A[i][pair.c1] + A[i][pair.c2]) % 2;
            mzd_set_ui(chi_A_win, 0);
            
            auto s_chi = std::chrono::high_resolution_clock::now();
            GateSynthesisMatrix::Chi_M4RI(A_m4ri_win, x_vec, n, this_m, chi_A_win);
            total_chi_duration += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - s_chi);

            int d_ns = 0;
            auto s_ns = std::chrono::high_resolution_clock::now();
            bool** NS = M4RI_direct_nullspace(chi_A_win, d_ns);
            total_ns_duration += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - s_ns);

            for (int h = 0; h < d_ns; h++) {
                if ((NS[pair.c1][h] + NS[pair.c2][h]) % 2 == 1) {
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < this_m; j++) Anew[i][j] = (A[i][j] + x_vec[i][0] * NS[j][h]) % 2;
                    }
                    int mp;
                    GateSynthesisMatrix::cleanup(Anew, n, this_m, mp);
                    if (mp < this_m) {
                        std::cout << "  [HIT!] Pair(" << pair.c1 << "," << pair.c2 << ") Dist=" << pair.dist 
                                  << " | " << this_m << " -> " << mp << " columns" << std::endl;
                        LCL_Mat_GF2::copy((const bool**)Anew, n, mp, Abest);
                        m_best = mp;
                        found = true;
                        break;
                    }
                }
            }
            if (NS) LCL_Mat_GF2::destruct(NS, this_m, d_ns);
        }

        mzd_free_window(A_m4ri_win);
        mzd_free_window(chi_A_win);

        if (found) {
            LCL_Mat_GF2::copy((const bool**)Abest, n, m_best, A);
            for (int r = 0; r < n; r++) {
                for (int c = 0; c < m_best; c++) mzd_write_bit(A_m4ri_full, r, c, A[r][c]);
                for (int c = m_best; c < m + 1; c++) mzd_write_bit(A_m4ri_full, r, c, 0);
            }
            this_m = m_best;
        }
        auto round_end = std::chrono::high_resolution_clock::now();
        std::cout << "Round " << round << " Finished in " << std::chrono::duration_cast<std::chrono::milliseconds>(round_end - round_start).count() << " ms." << std::endl;
        round++;
    }

    omp = this_m;
    auto end_total = std::chrono::high_resolution_clock::now();
    auto total_dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);

    std::cout << "\n=== Optimization Summary ===" << std::endl;
    std::cout << "Initial T-count : " << initial_m << std::endl;
    std::cout << "Final T-count   : " << omp << std::endl;
    std::cout << "Total Reduced   : " << (initial_m - omp) << " gates" << std::endl;
    std::cout << "Execution Time  : " << total_dur.count() << " ms" << std::endl;
    std::cout << "Chi calculation : " << total_chi_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Nullspace calc  : " << total_ns_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Preprocessing CNOTs: " << pre_cnots.size() << std::endl;
    std::cout << "============================" << std::endl;

    LCL_Mat_GF2::destruct(x_vec, n, 1);
    LCL_Mat_GF2::destruct(Anew, n, m + 1);
    LCL_Mat_GF2::destruct(Abest, n, m + 1);
    mzd_free(A_m4ri_full);
    mzd_free(chi_A_full);
}
void GateSynthesisMatrix::GreedyWeightReduction_Bool(bool** A, int n, int m, std::vector<std::pair<int, int>>& cnot_history) {
    if (n < 2) return;

    // Pack A into uint64_t for fast XOR and popcount
    int m_packed = (m + 63) / 64;
    std::vector<std::vector<uint64_t>> A_p(n, std::vector<uint64_t>(m_packed, 0));
    for (int i = 0; i < n; ++i) {
        for (int c = 0; c < m; ++c) {
            if (A[i][c]) A_p[i][c / 64] |= (1ULL << (c % 64));
        }
    }

    auto get_weight = [&](int i) {
        int w = 0;
        for (int p = 0; p < m_packed; ++p) w += __builtin_popcountll(A_p[i][p]);
        return w;
    };

    auto get_total_weight = [&]() {
        long long total = 0;
        for (int i = 0; i < n; ++i) total += get_weight(i);
        return total;
    };

    long long initial_total_weight = get_total_weight();

    bool improved = true;
    while (improved) {
        improved = false;
        int max_reduction = 0;
        int best_i = -1;
        std::vector<int> best_sources;

        std::vector<int> current_weights(n);
        for (int i = 0; i < n; ++i) current_weights[i] = get_weight(i);

        // 1-step (Row i ^ Row j)
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                int nw = 0;
                for (int p = 0; p < m_packed; ++p) nw += __builtin_popcountll(A_p[i][p] ^ A_p[j][p]);
                int red = current_weights[i] - nw;
                if (red > max_reduction) {
                    max_reduction = red; best_i = i; best_sources = {j};
                }
            }
        }

        // 2-step (Row i ^ Row j ^ Row k)
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                for (int k = j + 1; k < n; ++k) {
                    if (i == k) continue;
                    int nw = 0;
                    for (int p = 0; p < m_packed; ++p) nw += __builtin_popcountll(A_p[i][p] ^ A_p[j][p] ^ A_p[k][p]);
                    int red = current_weights[i] - nw;
                    if (red > max_reduction) {
                        max_reduction = red; best_i = i; best_sources = {j, k};
                    }
                }
            }
        }

        // 3-step (Row i ^ Row j ^ Row k ^ Row l)
        if (n >= 4 && n <= 100) {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    if (i == j) continue;
                    for (int k = j + 1; k < n; ++k) {
                        if (i == k) continue;
                        for (int l = k + 1; l < n; ++l) {
                            if (i == l) continue;
                            int nw = 0;
                            for (int p = 0; p < m_packed; ++p) nw += __builtin_popcountll(A_p[i][p] ^ A_p[j][p] ^ A_p[k][p] ^ A_p[l][p]);
                            int red = current_weights[i] - nw;
                            if (red > max_reduction) {
                                max_reduction = red; best_i = i; best_sources = {j, k, l};
                            }
                        }
                    }
                }
            }
        }

        if (max_reduction > 0) {
            std::cout << "  [Multi-step Greedy] XORing row " << best_i << " with " << best_sources.size() << " rows. Reduction: " << max_reduction << std::endl;
            for (int src : best_sources) {
                for (int p = 0; p < m_packed; ++p) A_p[best_i][p] ^= A_p[src][p];
                cnot_history.push_back({src, best_i});
            }
            improved = true;
        }
    }

    // Unpack A_p back to A
    for (int i = 0; i < n; ++i) {
        for (int c = 0; c < m; ++c) {
            A[i][c] = (A_p[i][c / 64] >> (c % 64)) & 1ULL;
        }
    }

    long long final_total_weight = get_total_weight();
    std::cout << "[Greedy TODD] Preprocessing complete. Total weight: " << initial_total_weight << " -> " << final_total_weight 
              << " (" << (initial_total_weight - final_total_weight) << " ones removed)" << std::endl;
}

void GateSynthesisMatrix::LempelX2_M4RI_GreedyPreprocess(bool** A, int n, int m, int& omp) {
    auto start_total = std::chrono::high_resolution_clock::now();
    std::cout << "\n[Greedy TODD] Starting with Greedy Weight Reduction + M4RI Hamming Search" << std::endl;
    int this_m = m; int initial_m = m; int n_chi_A = n * n * n;

    auto print_density_info = [&](const std::string& prefix, bool** M, int r_max, int c_max) {
        long long ones = 0;
        long long total = (long long)r_max * c_max;
        for (int r = 0; r < r_max; ++r) for (int c = 0; c < c_max; ++c) if (M[r][c]) ones++;
        double percent = (total > 0) ? (100.0 * ones / total) : 0.0;
        std::cout << prefix << " Density: " << ones << " ones (" << ((double)((long long)(percent * 100)) / 100.0) << "%)." << std::endl;
    };

    std::cout << "[Initial] Matrix size: " << n << "x" << m << std::endl;
    print_density_info("[Initial]", A, n, m);
    std::vector<std::pair<int, int>> pre_cnots;
    GreedyWeightReduction_Bool(A, n, m, pre_cnots);
    print_density_info("[Post-Greedy]", A, n, m);

    mzd_t* A_m4ri_full = convert_to_mzd((bool const**)A, n, m + 1);

    bool** x_vec = LCL_Mat_GF2::construct(n, 1);
    bool** Anew = LCL_Mat_GF2::construct(n, m + 1);
    bool** Abest = LCL_Mat_GF2::construct(n, m + 1);
    LCL_Mat_GF2::copy((const bool**)A, n, m, Abest);
    int m_best = m;
    mzd_t* chi_A_full = mzd_init(n_chi_A, m + 1);
    std::chrono::microseconds total_ns_duration{0};
    std::chrono::microseconds total_chi_duration{0};
    bool found = true; int round = 0;

    while (found && (round < m)) {
        found = false;
        std::cout << "--- Round " << round << " | Current Columns: " << this_m << " ---" << std::endl;
        auto round_start = std::chrono::high_resolution_clock::now();

        std::vector<ColPair> candidates;
        for (int j1 = 0; j1 < this_m; ++j1) {
            for (int j2 = j1 + 1; j2 < this_m; ++j2) {
                int dist = 0;
                for (int k = 0; k < n; ++k) if (A[k][j1] != A[k][j2]) dist++;
                candidates.push_back({j1, j2, dist});
            }
        }
        std::sort(candidates.begin(), candidates.end());

        mzd_t* A_m4ri_win = mzd_init_window(A_m4ri_full, 0, 0, n, this_m);
        mzd_t* chi_A_win = mzd_init_window(chi_A_full, 0, 0, n_chi_A, this_m);

        for (const auto& pair : candidates) {
            if (found) break;
            for (int i = 0; i < n; i++) x_vec[i][0] = (A[i][pair.c1] + A[i][pair.c2]) % 2;
            mzd_set_ui(chi_A_win, 0);
            
            auto s_chi = std::chrono::high_resolution_clock::now();
            GateSynthesisMatrix::Chi_M4RI(A_m4ri_win, x_vec, n, this_m, chi_A_win);
            total_chi_duration += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - s_chi);

            int d_ns = 0;
            auto s_ns = std::chrono::high_resolution_clock::now();
            bool** NS = M4RI_direct_nullspace(chi_A_win, d_ns);
            total_ns_duration += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - s_ns);

            for (int h = 0; h < d_ns; h++) {
                if ((NS[pair.c1][h] + NS[pair.c2][h]) % 2 == 1) {
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < this_m; j++) Anew[i][j] = (A[i][j] + x_vec[i][0] * NS[j][h]) % 2;
                    }
                    int mp;
                    GateSynthesisMatrix::cleanup(Anew, n, this_m, mp);
                    if (mp < this_m) {
                        std::cout << "  [HIT!] Pair(" << pair.c1 << "," << pair.c2 << ") Dist=" << pair.dist 
                                  << " | " << this_m << " -> " << mp << " columns" << std::endl;
                        LCL_Mat_GF2::copy((const bool**)Anew, n, mp, Abest);
                        m_best = mp;
                        found = true;
                        break;
                    }
                }
            }
            if (NS) LCL_Mat_GF2::destruct(NS, this_m, d_ns);
        }

        mzd_free_window(A_m4ri_win);
        mzd_free_window(chi_A_win);

        if (found) {
            LCL_Mat_GF2::copy((const bool**)Abest, n, m_best, A);
            for (int r = 0; r < n; r++) {
                for (int c = 0; c < m_best; c++) mzd_write_bit(A_m4ri_full, r, c, A[r][c]);
                for (int c = m_best; c < m + 1; c++) mzd_write_bit(A_m4ri_full, r, c, 0);
            }
            this_m = m_best;
        }
        auto round_end = std::chrono::high_resolution_clock::now();
        std::cout << "Round " << round << " Finished in " << std::chrono::duration_cast<std::chrono::milliseconds>(round_end - round_start).count() << " ms." << std::endl;
        round++;
    }

    omp = this_m;
    auto end_total = std::chrono::high_resolution_clock::now();
    auto total_dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);

    std::cout << "\n=== Optimization Summary ===" << std::endl;
    std::cout << "Algorithm       : LempelX2_M4RI_GreedyPreprocess" << std::endl;
    std::cout << "Initial T-count : " << initial_m << std::endl;
    std::cout << "Final T-count   : " << omp << std::endl;
    std::cout << "Total Reduced   : " << (initial_m - omp) << " gates" << std::endl;
    std::cout << "Execution Time  : " << total_dur.count() << " ms" << std::endl;
    std::cout << "Chi calculation : " << total_chi_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Nullspace calc  : " << total_ns_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Preprocessing CNOTs: " << pre_cnots.size() << std::endl;
    std::cout << "============================" << std::endl;

    LCL_Mat_GF2::destruct(x_vec, n, 1);
    LCL_Mat_GF2::destruct(Anew, n, m + 1);
    LCL_Mat_GF2::destruct(Abest, n, m + 1);
    mzd_free(A_m4ri_full);
}

struct FastColPair {
    int c1;
    int c2;
    int dist;
    std::vector<uint64_t> x_val;
    bool operator<(const FastColPair& o) const {
        if (dist != o.dist) return dist < o.dist;
        for (size_t i = 0; i < x_val.size(); ++i) {
            if (x_val[i] != o.x_val[i]) return x_val[i] < o.x_val[i];
        }
        if (c1 != o.c1) return c1 < o.c1;
        return c2 < o.c2;
    }
};

void GateSynthesisMatrix::LempelX2_M4RI_GreedyPreprocess_Fast(bool** A, int n, int m, int& omp) {
    auto start_total = std::chrono::high_resolution_clock::now();
    std::cout << "\n[Greedy TODD] Starting with Fast Hamming+Lazy Chi+Packed M4RI" << std::endl;
    int this_m = m; int initial_m = m;

    auto print_density_info = [&](const std::string& prefix, bool** M, int r_max, int c_max) {
        long long ones = 0;
        long long total = (long long)r_max * c_max;
        for (int r = 0; r < r_max; ++r) for (int c = 0; c < c_max; ++c) if (M[r][c]) ones++;
        double percent = (total > 0) ? (100.0 * ones / total) : 0.0;
        std::cout << prefix << " Density: " << ones << " ones (" << ((double)((long long)(percent * 100)) / 100.0) << "%)." << std::endl;
    };

    std::cout << "[Initial] Matrix size: " << n << "x" << m << std::endl;
    print_density_info("[Initial]", A, n, m);

    mzd_t* A_m4ri_full = convert_to_mzd((bool const**)A, n, m + 1);

    bool** x_vec = LCL_Mat_GF2::construct(n, 1);
    bool** Anew = LCL_Mat_GF2::construct(n, m + 1);
    bool** Abest = LCL_Mat_GF2::construct(n, m + 1);
    LCL_Mat_GF2::copy((const bool**)A, n, m, Abest);
    int m_best = m;
    
    std::chrono::microseconds total_ns_duration{0};
    std::chrono::microseconds total_chi_duration{0};
    bool found = true; int round = 0;

    int num_words = (n + 63) / 64;

    while (found && (round < m)) {
        found = false;
        std::cout << "--- Round " << round << " | Current Columns: " << this_m << " ---" << std::endl;
        auto round_start = std::chrono::high_resolution_clock::now();

        std::vector<std::vector<uint64_t>> cols(this_m, std::vector<uint64_t>(num_words, 0));
        for (int j = 0; j < this_m; ++j) {
            for (int i = 0; i < n; ++i) {
                if (A[i][j]) cols[j][i / 64] |= (1ULL << (i % 64));
            }
        }

        std::vector<FastColPair> candidates;
        for (int j1 = 0; j1 < this_m; ++j1) {
            for (int j2 = j1 + 1; j2 < this_m; ++j2) {
                int dist = 0;
                std::vector<uint64_t> x_val(num_words, 0);
                for (int w = 0; w < num_words; ++w) {
                    x_val[w] = cols[j1][w] ^ cols[j2][w];
                    dist += __builtin_popcountll(x_val[w]);
                }
                candidates.push_back({j1, j2, dist, x_val});
            }
        }
        std::sort(candidates.begin(), candidates.end());

        mzd_t* A_m4ri_win = mzd_init_window(A_m4ri_full, 0, 0, n, this_m);

        bool has_cached_ns = false;
        std::vector<uint64_t> last_x_val;
        bool** NS_cached = nullptr;
        int d_ns_cached = 0;

        for (const auto& pair : candidates) {
            if (found) break;
            
            bool same_x = (has_cached_ns && pair.x_val == last_x_val);
            
            if (!same_x) {
                if (NS_cached) {
                    LCL_Mat_GF2::destruct(NS_cached, this_m, d_ns_cached);
                    NS_cached = nullptr;
                }
                
                for (int i = 0; i < n; i++) x_vec[i][0] = (A[i][pair.c1] + A[i][pair.c2]) % 2;
                
                int non_zero_rows = 0;
                for(int alpha = 0; alpha < n; alpha++) {
                    bool x_a = x_vec[alpha][0];
                    for(int beta = 0; beta < n; beta++) {
                        bool x_b = x_vec[beta][0];
                        for(int gamma = 0; gamma < n; gamma++) {
                            bool x_c = x_vec[gamma][0];
                            if (x_a || x_b || x_c) non_zero_rows++;
                        }
                    }
                }
                if (non_zero_rows == 0) non_zero_rows = 1;
                
                mzd_t* chi_A_packed = mzd_init(non_zero_rows, this_m);
                
                auto s_chi = std::chrono::high_resolution_clock::now();
                int row_idx = 0;
                for(int alpha = 0; alpha < n; alpha++) {
                    bool x_a = x_vec[alpha][0];
                    for(int beta = 0; beta < n; beta++) {
                        bool x_b = x_vec[beta][0];
                        for(int gamma = 0; gamma < n; gamma++) {
                            bool x_c = x_vec[gamma][0];
                            if (!(x_a || x_b || x_c)) continue;
                            
                            bool term_const = x_a && x_b && x_c;
                            for(int w = 0; w < chi_A_packed->width; w++) {
                                word res = 0;
                                if(term_const) res = ~res;
                                if(x_a && x_b) res ^= mzd_row(A_m4ri_win, gamma)[w];
                                if(x_b && x_c) res ^= mzd_row(A_m4ri_win, alpha)[w];
                                if(x_c && x_a) res ^= mzd_row(A_m4ri_win, beta)[w];
                                if(x_a) res ^= (mzd_row(A_m4ri_win, beta)[w] & mzd_row(A_m4ri_win, gamma)[w]);
                                if(x_b) res ^= (mzd_row(A_m4ri_win, gamma)[w] & mzd_row(A_m4ri_win, alpha)[w]);
                                if(x_c) res ^= (mzd_row(A_m4ri_win, alpha)[w] & mzd_row(A_m4ri_win, beta)[w]);
                                mzd_row(chi_A_packed, row_idx)[w] = res;
                            }
                            row_idx++;
                        }
                    }
                }
                total_chi_duration += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - s_chi);

                d_ns_cached = 0;
                auto s_ns = std::chrono::high_resolution_clock::now();
                NS_cached = M4RI_direct_nullspace(chi_A_packed, d_ns_cached);
                total_ns_duration += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - s_ns);
                
                mzd_free(chi_A_packed);
                
                last_x_val = pair.x_val;
                has_cached_ns = true;
            }

            for (int h = 0; h < d_ns_cached; h++) {
                if ((NS_cached[pair.c1][h] + NS_cached[pair.c2][h]) % 2 == 1) {
                    for (int i = 0; i < n; i++) x_vec[i][0] = (A[i][pair.c1] + A[i][pair.c2]) % 2;

                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < this_m; j++) Anew[i][j] = (A[i][j] + x_vec[i][0] * NS_cached[j][h]) % 2;
                    }
                    int mp;
                    GateSynthesisMatrix::cleanup(Anew, n, this_m, mp);
                    if (mp < this_m) {
                        std::cout << "  [HIT!] Pair(" << pair.c1 << "," << pair.c2 << ") Dist=" << pair.dist 
                                  << " | " << this_m << " -> " << mp << " columns" << std::endl;
                        LCL_Mat_GF2::copy((const bool**)Anew, n, mp, Abest);
                        m_best = mp;
                        found = true;
                        break;
                    }
                }
            }
        }
        
        if (NS_cached) {
            LCL_Mat_GF2::destruct(NS_cached, this_m, d_ns_cached);
            NS_cached = nullptr;
        }

        mzd_free_window(A_m4ri_win);

        if (found) {
            LCL_Mat_GF2::copy((const bool**)Abest, n, m_best, A);
            for (int r = 0; r < n; r++) {
                for (int c = 0; c < m_best; c++) mzd_write_bit(A_m4ri_full, r, c, A[r][c]);
                for (int c = m_best; c < m + 1; c++) mzd_write_bit(A_m4ri_full, r, c, 0);
            }
            this_m = m_best;
        }
        auto round_end = std::chrono::high_resolution_clock::now();
        std::cout << "Round " << round << " Finished in " << std::chrono::duration_cast<std::chrono::milliseconds>(round_end - round_start).count() << " ms." << std::endl;
        round++;
    }

    omp = this_m;
    auto end_total = std::chrono::high_resolution_clock::now();
    auto total_dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);

    std::cout << "\n=== Optimization Summary ===" << std::endl;
    std::cout << "Algorithm       : LempelX2_M4RI_GreedyPreprocess_Fast" << std::endl;
    std::cout << "Initial T-count : " << initial_m << std::endl;
    std::cout << "Final T-count   : " << omp << std::endl;
    std::cout << "Total Reduced   : " << (initial_m - omp) << " gates" << std::endl;
    std::cout << "Execution Time  : " << total_dur.count() << " ms" << std::endl;
    std::cout << "Chi calculation : " << total_chi_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Nullspace calc  : " << total_ns_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "============================" << std::endl;

    LCL_Mat_GF2::destruct(x_vec, n, 1);
    LCL_Mat_GF2::destruct(Anew, n, m + 1);
    LCL_Mat_GF2::destruct(Abest, n, m + 1);
}

void GateSynthesisMatrix::LempelX2_M4RI_Experimental(bool** A, int n, int m, int& omp, bool use_packed_chi, bool use_aa_table, bool use_memoization, bool use_random_sketch, int sketch_margin) {
    auto start_total = std::chrono::high_resolution_clock::now();
    std::cout << "\n[Experimental TODD] Flags - Packed: " << use_packed_chi << ", AA Table: " << use_aa_table << ", Memo: " << use_memoization << ", Sketch: " << use_random_sketch << " (margin=" << sketch_margin << ")" << std::endl;
    int this_m = m; int initial_m = m;
    
    mzd_t* A_m4ri_full = convert_to_mzd((bool const**)A, n, m + 1);

    bool** x_vec = LCL_Mat_GF2::construct(n, 1);
    bool** Anew = LCL_Mat_GF2::construct(n, m + 1);
    bool** Abest = LCL_Mat_GF2::construct(n, m + 1);
    LCL_Mat_GF2::copy((const bool**)A, n, m, Abest);
    int m_best = m;
    
    std::chrono::microseconds total_ns_duration{0};
    std::chrono::microseconds total_chi_duration{0};
    std::chrono::microseconds total_aa_tb_duration{0};
    bool found = true; int round = 0;
    long long total_sketch_rejections = 0;
    long long total_sketch_pass = 0;
    long long total_sketch_skipped = 0;  // chi too small for sketch
    long long total_pairs_tested = 0;
    long long total_full_ns_empty = 0;   // track how often full nullspace is empty (for observability)
    std::chrono::microseconds total_sketch_duration{0};
    std::mt19937 rng(42); // Fixed seed for reproducibility

    int num_words = (n + 63) / 64;

    mzd_t* AA_tab = nullptr;
    if (use_aa_table) {
        AA_tab = mzd_init(n * (n - 1) / 2, m + 1);
        auto s_aa = std::chrono::high_resolution_clock::now();
        build_aa_table(AA_tab, A_m4ri_full, n, this_m);
        total_aa_tb_duration += std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - s_aa);
    }

    while (found && (round < m)) {
        found = false;
        std::cout << "--- Round " << round << " | Current Columns: " << this_m << " ---" << std::endl;
        auto round_start = std::chrono::high_resolution_clock::now();
        CleanupTrace best_cleanup_trace;
        std::vector<int> best_changed_rows;

        std::vector<std::vector<uint64_t>> cols(this_m, std::vector<uint64_t>(num_words, 0));
        for (int j = 0; j < this_m; ++j) {
            for (int i = 0; i < n; ++i) {
                if (A[i][j]) cols[j][i / 64] |= (1ULL << (i % 64));
            }
        }

        std::vector<FastColPair> candidates;
        for (int j1 = 0; j1 < this_m; ++j1) {
            for (int j2 = j1 + 1; j2 < this_m; ++j2) {
                int dist = 0;
                std::vector<uint64_t> x_val(num_words, 0);
                for (int w = 0; w < num_words; ++w) {
                    x_val[w] = cols[j1][w] ^ cols[j2][w];
                    dist += __builtin_popcountll(x_val[w]);
                }
                candidates.push_back({j1, j2, dist, x_val});
            }
        }
        
        if (use_memoization) {
            std::sort(candidates.begin(), candidates.end());
        } else {
            std::sort(candidates.begin(), candidates.end(), [](const FastColPair& a, const FastColPair& b) {
                if (a.dist != b.dist) return a.dist < b.dist;
                if (a.c1 != b.c1) return a.c1 < b.c1;
                return a.c2 < b.c2;
            });
        }

        mzd_t* A_m4ri_win = mzd_init_window(A_m4ri_full, 0, 0, n, this_m);

        auto get_AA = [&](int i, int j, int w, mzd_t* AA_win) -> word {
            if (i == j) return mzd_row(A_m4ri_win, i)[w];
            int idx = aa_table_index(n, i, j);
            return mzd_row(AA_win, idx)[w];
        };

        bool has_cached_ns = false;
        std::vector<uint64_t> last_x_val;
        bool** NS_cached = nullptr;
        int d_ns_cached = 0;

        for (const auto& pair : candidates) {
            if (found) break;
            
            bool same_x = (use_memoization && has_cached_ns && pair.x_val == last_x_val);
            
                if (!same_x) {
                if (NS_cached) {
                    LCL_Mat_GF2::destruct(NS_cached, this_m, d_ns_cached);
                    NS_cached = nullptr;
                }
                
                for (int i = 0; i < n; i++) x_vec[i][0] = (A[i][pair.c1] + A[i][pair.c2]) % 2;
                
                int non_zero_rows = 0;
                if (use_packed_chi) {
                    for(int alpha = 0; alpha < n; alpha++) {
                        bool x_a = x_vec[alpha][0];
                        for(int beta = 0; beta < n; beta++) {
                            bool x_b = x_vec[beta][0];
                            for(int gamma = 0; gamma < n; gamma++) {
                                bool x_c = x_vec[gamma][0];
                                if (x_a || x_b || x_c) non_zero_rows++;
                            }
                        }
                    }
                    if (non_zero_rows == 0) non_zero_rows = 1;
                } else {
                    non_zero_rows = n * n * n;
                }
                
                mzd_t* chi_A_experimental = mzd_init(non_zero_rows, this_m);
                
                auto s_chi = std::chrono::high_resolution_clock::now();
                int row_idx = 0;
                
                mzd_t* AA_win = nullptr;
                if (use_aa_table) {
                    AA_win = mzd_init_window(AA_tab, 0, 0, AA_tab->nrows, this_m);
                }

                for(int alpha = 0; alpha < n; alpha++) {
                    bool x_a = x_vec[alpha][0];
                    for(int beta = 0; beta < n; beta++) {
                        bool x_b = x_vec[beta][0];
                        for(int gamma = 0; gamma < n; gamma++) {
                            bool x_c = x_vec[gamma][0];
                            if (use_packed_chi && !(x_a || x_b || x_c)) continue;
                            
                            bool term_const = x_a && x_b && x_c;
                            for(int w = 0; w < chi_A_experimental->width; w++) {
                                word res = 0;
                                if(term_const) res = ~res;
                                
                                if (use_aa_table) {
                                    if(x_a) res ^= get_AA(beta, gamma, w, AA_win);
                                    if(x_b) res ^= get_AA(gamma, alpha, w, AA_win);
                                    if(x_c) res ^= get_AA(alpha, beta, w, AA_win);
                                    
                                    if(x_a && x_b) res ^= mzd_row(A_m4ri_win, gamma)[w];
                                    if(x_b && x_c) res ^= mzd_row(A_m4ri_win, alpha)[w];
                                    if(x_c && x_a) res ^= mzd_row(A_m4ri_win, beta)[w];
                                } else {
                                    if(x_a && x_b) res ^= mzd_row(A_m4ri_win, gamma)[w];
                                    if(x_b && x_c) res ^= mzd_row(A_m4ri_win, alpha)[w];
                                    if(x_c && x_a) res ^= mzd_row(A_m4ri_win, beta)[w];
                                    if(x_a) res ^= (mzd_row(A_m4ri_win, beta)[w] & mzd_row(A_m4ri_win, gamma)[w]);
                                    if(x_b) res ^= (mzd_row(A_m4ri_win, gamma)[w] & mzd_row(A_m4ri_win, alpha)[w]);
                                    if(x_c) res ^= (mzd_row(A_m4ri_win, alpha)[w] & mzd_row(A_m4ri_win, beta)[w]);
                                }
                                mzd_row(chi_A_experimental, row_idx)[w] = res;
                            }
                            row_idx++;
                        }
                    }
                }
                
                if (use_aa_table) {
                    mzd_free_window(AA_win);
                }
                
                total_chi_duration += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - s_chi);

                // === Random Sketch Early Rejection (pair-condition aware) ===
                // Since chi nullspace is always non-empty, empty-check is useless.
                // Instead, check if sketch nullspace contains v with v[c1]⊕v[c2]=1.
                // Null(M) ⊆ Null(M'): if no such v in Null(M'), none in Null(M) either.
                int chi_rows = non_zero_rows;
                int sketch_r = this_m + sketch_margin;
                bool sketch_rejected = false;
                total_pairs_tested++;
                
                if (use_random_sketch && chi_rows > sketch_r && sketch_r > 0) {
                    auto s_sketch = std::chrono::high_resolution_clock::now();
                    
                    // Fisher-Yates partial shuffle to pick r rows without replacement
                    std::vector<int> row_indices(chi_rows);
                    for (int i = 0; i < chi_rows; i++) row_indices[i] = i;
                    for (int i = 0; i < sketch_r; i++) {
                        std::uniform_int_distribution<int> dist(i, chi_rows - 1);
                        std::swap(row_indices[i], row_indices[dist(rng)]);
                    }
                    
                    // Create small sketch matrix by copying selected rows
                    mzd_t* sketch_mat = mzd_init(sketch_r, this_m);
                    for (int i = 0; i < sketch_r; i++) {
                        mzd_copy_row(sketch_mat, i, chi_A_experimental, row_indices[i]);
                    }
                    
                    // Compute nullspace of the sketch matrix
                    int sketch_d = 0;
                    bool** sketch_NS = M4RI_direct_nullspace(sketch_mat, sketch_d);
                    
                    if (sketch_NS) {
                        // Check: does ANY nullspace vector satisfy v[c1] XOR v[c2] = 1?
                        bool has_valid_vector = false;
                        for (int h = 0; h < sketch_d; h++) {
                            if ((sketch_NS[pair.c1][h] + sketch_NS[pair.c2][h]) % 2 == 1) {
                                has_valid_vector = true;
                                break;
                            }
                        }
                        LCL_Mat_GF2::destruct(sketch_NS, this_m, sketch_d);
                        
                        if (!has_valid_vector) {
                            // No vector satisfies pair condition => REJECT
                            sketch_rejected = true;
                            total_sketch_rejections++;
                        } else {
                            total_sketch_pass++;
                        }
                    } else {
                        // Nullspace empty (shouldn't happen for chi, but handle it)
                        sketch_rejected = true;
                        total_sketch_rejections++;
                    }
                    
                    mzd_free(sketch_mat);
                    total_sketch_duration += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - s_sketch);
                } else if (use_random_sketch) {
                    total_sketch_skipped++;
                }
                
                if (sketch_rejected) {
                    mzd_free(chi_A_experimental);
                    d_ns_cached = 0;
                    NS_cached = nullptr;
                    has_cached_ns = false;
                    continue;
                }
                // === End Random Sketch ===

                d_ns_cached = 0;
                auto s_ns = std::chrono::high_resolution_clock::now();
                NS_cached = M4RI_direct_nullspace(chi_A_experimental, d_ns_cached);
                total_ns_duration += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - s_ns);
                
                // Track whether full nullspace was empty (observability)
                if (d_ns_cached == 0) total_full_ns_empty++;
                
                mzd_free(chi_A_experimental);
                
                if (use_memoization) {
                    last_x_val = pair.x_val;
                    has_cached_ns = true;
                }
            }

            for (int h = 0; h < d_ns_cached; h++) {
                if ((NS_cached[pair.c1][h] + NS_cached[pair.c2][h]) % 2 == 1) {
                    for (int i = 0; i < n; i++) x_vec[i][0] = (A[i][pair.c1] + A[i][pair.c2]) % 2;

                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < this_m; j++) Anew[i][j] = (A[i][j] + x_vec[i][0] * NS_cached[j][h]) % 2;
                    }
                    int mp;
                    CleanupTrace cleanup_trace;
                    cleanup_with_trace(Anew, n, this_m, mp, cleanup_trace);
                    if (mp < this_m) {
                        std::cout << "  [HIT!] Pair(" << pair.c1 << "," << pair.c2 << ") Dist=" << pair.dist 
                                  << " | " << this_m << " -> " << mp << " columns" << std::endl;
                        LCL_Mat_GF2::copy((const bool**)Anew, n, mp, Abest);
                        m_best = mp;
                        best_cleanup_trace = cleanup_trace;
                        best_changed_rows.clear();
                        for (int i = 0; i < n; ++i) {
                            if (x_vec[i][0]) best_changed_rows.push_back(i);
                        }
                        found = true;
                        break;
                    }
                }
            }
        }
        
        if (NS_cached) {
            LCL_Mat_GF2::destruct(NS_cached, this_m, d_ns_cached);
            NS_cached = nullptr;
        }

        mzd_free_window(A_m4ri_win);

        if (found) {
            LCL_Mat_GF2::copy((const bool**)Abest, n, m_best, A);
            for (int r = 0; r < n; r++) {
                for (int c = 0; c < m_best; c++) mzd_write_bit(A_m4ri_full, r, c, A[r][c]);
                for (int c = m_best; c < m + 1; c++) mzd_write_bit(A_m4ri_full, r, c, 0);
            }
            if (use_aa_table) {
                auto s_aa = std::chrono::high_resolution_clock::now();
                update_aa_table_after_hit(AA_tab, A_m4ri_full, n, m_best, best_cleanup_trace, best_changed_rows);
                total_aa_tb_duration += std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::high_resolution_clock::now() - s_aa);
            }
            this_m = m_best;
        }
        auto round_end = std::chrono::high_resolution_clock::now();
        std::cout << "Round " << round << " Finished in " << std::chrono::duration_cast<std::chrono::milliseconds>(round_end - round_start).count() << " ms." << std::endl;
        round++;
    }

    omp = this_m;
    auto end_total = std::chrono::high_resolution_clock::now();
    auto total_dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);

    std::cout << "\n=== Optimization Summary ===" << std::endl;
    std::cout << "Algorithm       : LempelX2_M4RI_Experimental" << std::endl;
    std::cout << "Initial T-count : " << initial_m << std::endl;
    std::cout << "Final T-count   : " << omp << std::endl;
    std::cout << "Total Reduced   : " << (initial_m - omp) << " gates" << std::endl;
    std::cout << "Execution Time  : " << total_dur.count() << " ms" << std::endl;
    std::cout << "Chi calculation : " << total_chi_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "AA Table calc   : " << total_aa_tb_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Nullspace calc  : " << total_ns_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "--- Sketch Statistics ---" << std::endl;
    std::cout << "Sketch margin   : delta = " << sketch_margin << std::endl;
    std::cout << "Sketch time     : " << total_sketch_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Pairs tested    : " << total_pairs_tested << std::endl;
    std::cout << "Sketch fired    : " << (total_sketch_rejections + total_sketch_pass) << std::endl;
    std::cout << "Sketch skipped  : " << total_sketch_skipped << " (chi_rows <= m+delta)" << std::endl;
    std::cout << "Sketch rejected : " << total_sketch_rejections << " (nullspace empty)" << std::endl;
    std::cout << "Sketch passed   : " << total_sketch_pass << " (nullspace non-empty)" << std::endl;
    if (total_sketch_rejections + total_sketch_pass > 0) {
        double reject_rate = 100.0 * total_sketch_rejections / (total_sketch_rejections + total_sketch_pass);
        std::cout << "Reject rate     : " << std::fixed << std::setprecision(2) << reject_rate << "%" << std::endl;
    }
    std::cout << "Full NS empty   : " << total_full_ns_empty << " (of " << (total_pairs_tested - total_sketch_rejections) << " computed)" << std::endl;
    std::cout << "============================" << std::endl;

    if (AA_tab) mzd_free(AA_tab);
    LCL_Mat_GF2::destruct(x_vec, n, 1);
    LCL_Mat_GF2::destruct(Anew, n, m + 1);
    LCL_Mat_GF2::destruct(Abest, n, m + 1);
    mzd_free(A_m4ri_full);
}

void GateSynthesisMatrix::LempelX2_DynamicBasis(bool** A, int n, int m, int& omp) {
    auto start_total = std::chrono::high_resolution_clock::now();
    std::cout << "\n[Dynamic Basis TODD] Baseline target: original LempelX2" << std::endl;

    int this_m = m;
    int initial_m = m;
    bool** Anew = LCL_Mat_GF2::construct(n, m + 1);
    bool** Abest = LCL_Mat_GF2::construct(n, m + 1);
    LCL_Mat_GF2::copy((const bool**)A, n, m, Abest);
    int m_best = m;

    std::chrono::microseconds total_chi_update_duration(0);
    std::chrono::microseconds total_basis_rebuild_duration(0);
    std::chrono::microseconds total_ns_duration(0);

    long long total_pairs_tested = 0;
    long long total_filtered_miss = 0;
    long long total_ns_runs = 0;
    long long total_rebuilds = 0;
    long long total_affected_rows = 0;

    bool found = true;
    int round = 0;
    while (found && (round < m)) {
        found = false;
        std::cout << "--- Round " << round << " | Current Columns: " << this_m << " ---" << std::endl;

        std::vector<DynamicCandidate> candidates;
        candidates.reserve(this_m * (this_m - 1) / 2);
        for (int j1 = 0; j1 < this_m; ++j1) {
            for (int j2 = j1 + 1; j2 < this_m; ++j2) {
                std::vector<unsigned char> x_bits(n, 0);
                int x_weight = 0;
                for (int i = 0; i < n; ++i) {
                    x_bits[i] = (unsigned char)((A[i][j1] + A[i][j2]) % 2);
                    x_weight += x_bits[i];
                }
                DynamicCandidate cand;
                cand.c1 = j1;
                cand.c2 = j2;
                cand.x_weight = x_weight;
                cand.x_bits = x_bits;
                candidates.push_back(cand);
            }
        }

        std::sort(candidates.begin(), candidates.end(), [](const DynamicCandidate& a, const DynamicCandidate& b) {
            if (a.x_weight != b.x_weight) return a.x_weight < b.x_weight;
            if (a.x_bits != b.x_bits) return a.x_bits < b.x_bits;
            if (a.c1 != b.c1) return a.c1 < b.c1;
            return a.c2 < b.c2;
        });

        DynamicChiState state;
        std::vector<int> affected_rows;
        bool state_ready = false;

        for (int idx = 0; idx < (int)candidates.size() && !found; ++idx) {
            const DynamicCandidate& pair = candidates[idx];
            total_pairs_tested++;

            auto s_chi = std::chrono::high_resolution_clock::now();
            auto s_basis = std::chrono::high_resolution_clock::now();
            if (!state_ready) {
                dynamic_build_full_chi_state(state, A, pair.x_bits, n, this_m);
                state_ready = true;
                total_rebuilds++;
            } else {
                dynamic_update_chi_state(state, A, pair.x_bits, n, this_m, affected_rows);
                total_rebuilds++;
                total_affected_rows += (long long)affected_rows.size();
            }
            auto e_basis = std::chrono::high_resolution_clock::now();
            total_chi_update_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_basis - s_chi);
            total_basis_rebuild_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_basis - s_basis);

            std::vector<unsigned char> e_vec(this_m, 0);
            e_vec[pair.c1] = 1;
            e_vec[pair.c2] = 1;

            if (dynamic_membership_test(e_vec, state.basis)) {
                total_filtered_miss++;
                continue;
            }

            auto s_ns = std::chrono::high_resolution_clock::now();
            std::vector<std::vector<unsigned char>> ns_basis = dynamic_nullspace_basis(state.chi_rows, state.row_active, this_m);
            auto e_ns = std::chrono::high_resolution_clock::now();
            total_ns_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_ns - s_ns);
            total_ns_runs++;

            int good_idx = -1;
            for (int h = 0; h < (int)ns_basis.size(); ++h) {
                if ((ns_basis[h][pair.c1] ^ ns_basis[h][pair.c2]) == 1) {
                    good_idx = h;
                    break;
                }
            }
            if (good_idx < 0) continue;

            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < this_m; ++j) {
                    Anew[i][j] = (bool)((A[i][j] + pair.x_bits[i] * ns_basis[good_idx][j]) % 2);
                }
            }

            int mp = 0;
            GateSynthesisMatrix::cleanup(Anew, n, this_m, mp);
            if (mp < this_m) {
                std::cout << "  [HIT!] Pair(" << pair.c1 << "," << pair.c2 << ") x-weight=" << pair.x_weight
                          << " | " << this_m << " -> " << mp << " columns" << std::endl;
                LCL_Mat_GF2::copy((const bool**)Anew, n, mp, Abest);
                m_best = mp;
                found = true;
            }
        }

        if (found) {
            LCL_Mat_GF2::copy((const bool**)Abest, n, m_best, A);
            this_m = m_best;
        }
        round++;
    }

    omp = this_m;
    auto end_total = std::chrono::high_resolution_clock::now();
    auto total_dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);

    std::cout << "\n=== Dynamic Basis Summary ===" << std::endl;
    std::cout << "Algorithm       : LempelX2_DynamicBasis" << std::endl;
    std::cout << "Initial T-count : " << initial_m << std::endl;
    std::cout << "Final T-count   : " << omp << std::endl;
    std::cout << "Total Reduced   : " << (initial_m - omp) << " gates" << std::endl;
    std::cout << "Execution Time  : " << total_dur.count() << " ms" << std::endl;
    std::cout << "Pairs tested    : " << total_pairs_tested << std::endl;
    std::cout << "Filtered miss   : " << total_filtered_miss << std::endl;
    std::cout << "Nullspace runs  : " << total_ns_runs << std::endl;
    std::cout << "Chi update time : " << total_chi_update_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Basis rebuild   : " << total_basis_rebuild_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Nullspace time  : " << total_ns_duration.count() / 1000.0 << " ms" << std::endl;
    if (total_rebuilds > 1) {
        std::cout << "Avg affected    : " << (double)total_affected_rows / (double)(total_rebuilds - 1) << " rows" << std::endl;
    }
    std::cout << "=============================" << std::endl;

    LCL_Mat_GF2::destruct(Anew, n, m + 1);
    LCL_Mat_GF2::destruct(Abest, n, m + 1);
}

void GateSynthesisMatrix::LempelX2_DynamicBasisLocalRepair(bool** A, int n, int m, int& omp) {
    auto start_total = std::chrono::high_resolution_clock::now();
    std::cout << "\n[Dynamic Basis Local Repair TODD] Compare with lx2_dynamic" << std::endl;

    int this_m = m;
    int initial_m = m;
    bool** Anew = LCL_Mat_GF2::construct(n, m + 1);
    bool** Abest = LCL_Mat_GF2::construct(n, m + 1);
    LCL_Mat_GF2::copy((const bool**)A, n, m, Abest);
    int m_best = m;

    std::chrono::microseconds total_chi_update_duration(0);
    std::chrono::microseconds total_basis_rebuild_duration(0);
    std::chrono::microseconds total_ns_duration(0);

    long long total_pairs_tested = 0;
    long long total_filtered_miss = 0;
    long long total_ns_runs = 0;
    long long total_rebuilds = 0;
    long long total_local_add_rounds = 0;
    long long total_affected_rows = 0;

    bool found = true;
    int round = 0;
    while (found && (round < m)) {
        found = false;
        std::cout << "--- Round " << round << " | Current Columns: " << this_m << " ---" << std::endl;

        std::vector<DynamicCandidate> candidates;
        candidates.reserve(this_m * (this_m - 1) / 2);
        for (int j1 = 0; j1 < this_m; ++j1) {
            for (int j2 = j1 + 1; j2 < this_m; ++j2) {
                std::vector<unsigned char> x_bits(n, 0);
                int x_weight = 0;
                for (int i = 0; i < n; ++i) {
                    x_bits[i] = (unsigned char)((A[i][j1] + A[i][j2]) % 2);
                    x_weight += x_bits[i];
                }
                DynamicCandidate cand;
                cand.c1 = j1;
                cand.c2 = j2;
                cand.x_weight = x_weight;
                cand.x_bits = x_bits;
                candidates.push_back(cand);
            }
        }

        std::sort(candidates.begin(), candidates.end(), [](const DynamicCandidate& a, const DynamicCandidate& b) {
            if (a.x_weight != b.x_weight) return a.x_weight < b.x_weight;
            if (a.x_bits != b.x_bits) return a.x_bits < b.x_bits;
            if (a.c1 != b.c1) return a.c1 < b.c1;
            return a.c2 < b.c2;
        });

        DynamicChiState state;
        std::vector<int> affected_rows;
        bool state_ready = false;

        for (int idx = 0; idx < (int)candidates.size() && !found; ++idx) {
            const DynamicCandidate& pair = candidates[idx];
            total_pairs_tested++;

            auto s_all = std::chrono::high_resolution_clock::now();
            if (!state_ready) {
                dynamic_build_full_chi_state(state, A, pair.x_bits, n, this_m);
                state_ready = true;
                total_rebuilds++;
            } else {
                bool rebuilt_basis = false;
                dynamic_update_chi_state_local_repair(state, A, pair.x_bits, n, this_m, affected_rows, rebuilt_basis);
                total_affected_rows += (long long)affected_rows.size();
                if (rebuilt_basis) total_rebuilds++;
                else total_local_add_rounds++;
            }
            auto e_all = std::chrono::high_resolution_clock::now();
            total_chi_update_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_all - s_all);
            total_basis_rebuild_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_all - s_all);

            std::vector<unsigned char> e_vec(this_m, 0);
            e_vec[pair.c1] = 1;
            e_vec[pair.c2] = 1;

            if (dynamic_membership_test(e_vec, state.basis)) {
                total_filtered_miss++;
                continue;
            }

            auto s_ns = std::chrono::high_resolution_clock::now();
            std::vector<std::vector<unsigned char>> ns_basis = dynamic_nullspace_basis(state.chi_rows, state.row_active, this_m);
            auto e_ns = std::chrono::high_resolution_clock::now();
            total_ns_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_ns - s_ns);
            total_ns_runs++;

            int good_idx = -1;
            for (int h = 0; h < (int)ns_basis.size(); ++h) {
                if ((ns_basis[h][pair.c1] ^ ns_basis[h][pair.c2]) == 1) {
                    good_idx = h;
                    break;
                }
            }
            if (good_idx < 0) continue;

            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < this_m; ++j) {
                    Anew[i][j] = (bool)((A[i][j] + pair.x_bits[i] * ns_basis[good_idx][j]) % 2);
                }
            }

            int mp = 0;
            GateSynthesisMatrix::cleanup(Anew, n, this_m, mp);
            if (mp < this_m) {
                std::cout << "  [HIT!] Pair(" << pair.c1 << "," << pair.c2 << ") x-weight=" << pair.x_weight
                          << " | " << this_m << " -> " << mp << " columns" << std::endl;
                LCL_Mat_GF2::copy((const bool**)Anew, n, mp, Abest);
                m_best = mp;
                found = true;
            }
        }

        if (found) {
            LCL_Mat_GF2::copy((const bool**)Abest, n, m_best, A);
            this_m = m_best;
        }
        round++;
    }

    omp = this_m;
    auto end_total = std::chrono::high_resolution_clock::now();
    auto total_dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);

    std::cout << "\n=== Dynamic Basis Local Repair Summary ===" << std::endl;
    std::cout << "Algorithm       : LempelX2_DynamicBasisLocalRepair" << std::endl;
    std::cout << "Initial T-count : " << initial_m << std::endl;
    std::cout << "Final T-count   : " << omp << std::endl;
    std::cout << "Total Reduced   : " << (initial_m - omp) << " gates" << std::endl;
    std::cout << "Execution Time  : " << total_dur.count() << " ms" << std::endl;
    std::cout << "Pairs tested    : " << total_pairs_tested << std::endl;
    std::cout << "Filtered miss   : " << total_filtered_miss << std::endl;
    std::cout << "Nullspace runs  : " << total_ns_runs << std::endl;
    std::cout << "Rebuild count   : " << total_rebuilds << std::endl;
    std::cout << "Local add count : " << total_local_add_rounds << std::endl;
    std::cout << "Chi update time : " << total_chi_update_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Basis work time : " << total_basis_rebuild_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Nullspace time  : " << total_ns_duration.count() / 1000.0 << " ms" << std::endl;
    if ((total_rebuilds + total_local_add_rounds) > 1) {
        std::cout << "Avg affected    : " << (double)total_affected_rows / (double)(total_rebuilds + total_local_add_rounds - 1) << " rows" << std::endl;
    }
    std::cout << "=========================================" << std::endl;

    LCL_Mat_GF2::destruct(Anew, n, m + 1);
    LCL_Mat_GF2::destruct(Abest, n, m + 1);
}

void GateSynthesisMatrix::LempelX2_DynamicBasisLocalRepairK1(bool** A, int n, int m, int& omp) {
    auto start_total = std::chrono::high_resolution_clock::now();
    std::cout << "\n[Dynamic Basis Local Repair TODD k=1] Compare with lx2_dynamic_repair" << std::endl;

    int this_m = m;
    int initial_m = m;
    bool** Anew = LCL_Mat_GF2::construct(n, m + 1);
    bool** Abest = LCL_Mat_GF2::construct(n, m + 1);
    LCL_Mat_GF2::copy((const bool**)A, n, m, Abest);
    int m_best = m;

    std::chrono::microseconds total_chi_update_duration(0);
    std::chrono::microseconds total_basis_work_duration(0);
    std::chrono::microseconds total_ns_duration(0);

    long long total_pairs_tested = 0;
    long long total_filtered_miss = 0;
    long long total_ns_runs = 0;
    long long total_rebuilds = 0;
    long long total_local_add_rounds = 0;
    long long total_repair_attempts = 0;
    long long total_repair_success = 0;
    long long total_repair_fallback = 0;
    long long total_repaired_basis_rows = 0;
    long long total_affected_rows = 0;

    bool found = true;
    int round = 0;
    while (found && (round < m)) {
        found = false;
        std::cout << "--- Round " << round << " | Current Columns: " << this_m << " ---" << std::endl;

        std::vector<DynamicCandidate> candidates;
        candidates.reserve(this_m * (this_m - 1) / 2);
        for (int j1 = 0; j1 < this_m; ++j1) {
            for (int j2 = j1 + 1; j2 < this_m; ++j2) {
                std::vector<unsigned char> x_bits(n, 0);
                int x_weight = 0;
                for (int i = 0; i < n; ++i) {
                    x_bits[i] = (unsigned char)((A[i][j1] + A[i][j2]) % 2);
                    x_weight += x_bits[i];
                }
                DynamicCandidate cand;
                cand.c1 = j1;
                cand.c2 = j2;
                cand.x_weight = x_weight;
                cand.x_bits = x_bits;
                candidates.push_back(cand);
            }
        }

        std::sort(candidates.begin(), candidates.end(), [](const DynamicCandidate& a, const DynamicCandidate& b) {
            if (a.x_weight != b.x_weight) return a.x_weight < b.x_weight;
            if (a.x_bits != b.x_bits) return a.x_bits < b.x_bits;
            if (a.c1 != b.c1) return a.c1 < b.c1;
            return a.c2 < b.c2;
        });

        DynamicChiState state;
        std::vector<int> affected_rows;
        bool state_ready = false;

        for (int idx = 0; idx < (int)candidates.size() && !found; ++idx) {
            const DynamicCandidate& pair = candidates[idx];
            total_pairs_tested++;

            auto s_all = std::chrono::high_resolution_clock::now();
            if (!state_ready) {
                dynamic_build_full_chi_state(state, A, pair.x_bits, n, this_m);
                state_ready = true;
                total_rebuilds++;
            } else {
                bool rebuilt_basis = false;
                bool repaired_basis = false;
                int repaired_basis_rows = 0;
                dynamic_update_chi_state_local_repair_k1(
                    state, A, pair.x_bits, n, this_m, affected_rows,
                    rebuilt_basis, repaired_basis, repaired_basis_rows);
                total_affected_rows += (long long)affected_rows.size();
                if (repaired_basis_rows == 1) total_repair_attempts++;
                if (repaired_basis) {
                    total_repair_success++;
                    total_repaired_basis_rows += repaired_basis_rows;
                } else if (repaired_basis_rows == 1) {
                    total_repair_fallback++;
                }
                if (rebuilt_basis) total_rebuilds++;
                else if (repaired_basis) total_local_add_rounds++;
                else total_local_add_rounds++;
            }
            auto e_all = std::chrono::high_resolution_clock::now();
            total_chi_update_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_all - s_all);
            total_basis_work_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_all - s_all);

            std::vector<unsigned char> e_vec(this_m, 0);
            e_vec[pair.c1] = 1;
            e_vec[pair.c2] = 1;

            if (dynamic_membership_test(e_vec, state.basis)) {
                total_filtered_miss++;
                continue;
            }

            auto s_ns = std::chrono::high_resolution_clock::now();
            std::vector<std::vector<unsigned char>> ns_basis = dynamic_nullspace_basis(state.chi_rows, state.row_active, this_m);
            auto e_ns = std::chrono::high_resolution_clock::now();
            total_ns_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_ns - s_ns);
            total_ns_runs++;

            int good_idx = -1;
            for (int h = 0; h < (int)ns_basis.size(); ++h) {
                if ((ns_basis[h][pair.c1] ^ ns_basis[h][pair.c2]) == 1) {
                    good_idx = h;
                    break;
                }
            }
            if (good_idx < 0) continue;

            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < this_m; ++j) {
                    Anew[i][j] = (bool)((A[i][j] + pair.x_bits[i] * ns_basis[good_idx][j]) % 2);
                }
            }

            int mp = 0;
            GateSynthesisMatrix::cleanup(Anew, n, this_m, mp);
            if (mp < this_m) {
                std::cout << "  [HIT!] Pair(" << pair.c1 << "," << pair.c2 << ") x-weight=" << pair.x_weight
                          << " | " << this_m << " -> " << mp << " columns" << std::endl;
                LCL_Mat_GF2::copy((const bool**)Anew, n, mp, Abest);
                m_best = mp;
                found = true;
            }
        }

        if (found) {
            LCL_Mat_GF2::copy((const bool**)Abest, n, m_best, A);
            this_m = m_best;
        }
        round++;
    }

    omp = this_m;
    auto end_total = std::chrono::high_resolution_clock::now();
    auto total_dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);

    std::cout << "\n=== Dynamic Basis Local Repair k=1 Summary ===" << std::endl;
    std::cout << "Algorithm        : LempelX2_DynamicBasisLocalRepairK1" << std::endl;
    std::cout << "Initial T-count  : " << initial_m << std::endl;
    std::cout << "Final T-count    : " << omp << std::endl;
    std::cout << "Total Reduced    : " << (initial_m - omp) << " gates" << std::endl;
    std::cout << "Execution Time   : " << total_dur.count() << " ms" << std::endl;
    std::cout << "Pairs tested     : " << total_pairs_tested << std::endl;
    std::cout << "Filtered miss    : " << total_filtered_miss << std::endl;
    std::cout << "Nullspace runs   : " << total_ns_runs << std::endl;
    std::cout << "Rebuild count    : " << total_rebuilds << std::endl;
    std::cout << "Local add count  : " << total_local_add_rounds << std::endl;
    std::cout << "Repair attempts  : " << total_repair_attempts << std::endl;
    std::cout << "Repair success   : " << total_repair_success << std::endl;
    std::cout << "Repair fallback  : " << total_repair_fallback << std::endl;
    std::cout << "Repaired basis   : " << total_repaired_basis_rows << std::endl;
    std::cout << "Chi update time  : " << total_chi_update_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Basis work time  : " << total_basis_work_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Nullspace time   : " << total_ns_duration.count() / 1000.0 << " ms" << std::endl;
    if ((total_rebuilds + total_local_add_rounds) > 1) {
        std::cout << "Avg affected     : " << (double)total_affected_rows / (double)(total_rebuilds + total_local_add_rounds - 1) << " rows" << std::endl;
    }
    std::cout << "============================================" << std::endl;

    LCL_Mat_GF2::destruct(Anew, n, m + 1);
    LCL_Mat_GF2::destruct(Abest, n, m + 1);
}

void GateSynthesisMatrix::LempelX2_DynamicBasisLocalRepairChi(bool** A, int n, int m, int& omp) {
    auto start_total = std::chrono::high_resolution_clock::now();
    std::cout << "\n[Dynamic Basis Local Repair TODD + Chi Compression] Compare with lx2_dynamic_repair" << std::endl;

    int this_m = m;
    int initial_m = m;
    bool** Anew = LCL_Mat_GF2::construct(n, m + 1);
    bool** Abest = LCL_Mat_GF2::construct(n, m + 1);
    LCL_Mat_GF2::copy((const bool**)A, n, m, Abest);
    int m_best = m;

    std::chrono::microseconds total_chi_update_duration(0);
    std::chrono::microseconds total_basis_work_duration(0);
    std::chrono::microseconds total_ns_duration(0);

    long long total_pairs_tested = 0;
    long long total_filtered_miss = 0;
    long long total_ns_runs = 0;
    long long total_rebuilds = 0;
    long long total_local_add_rounds = 0;
    long long total_affected_rows = 0;
    long long total_active_rows = 0;
    long long total_active_snapshots = 0;

    bool found = true;
    int round = 0;
    while (found && (round < m)) {
        found = false;
        std::cout << "--- Round " << round << " | Current Columns: " << this_m << " ---" << std::endl;

        std::vector<DynamicCandidate> candidates;
        candidates.reserve(this_m * (this_m - 1) / 2);
        for (int j1 = 0; j1 < this_m; ++j1) {
            for (int j2 = j1 + 1; j2 < this_m; ++j2) {
                std::vector<unsigned char> x_bits(n, 0);
                int x_weight = 0;
                for (int i = 0; i < n; ++i) {
                    x_bits[i] = (unsigned char)((A[i][j1] + A[i][j2]) % 2);
                    x_weight += x_bits[i];
                }
                DynamicCandidate cand;
                cand.c1 = j1;
                cand.c2 = j2;
                cand.x_weight = x_weight;
                cand.x_bits = x_bits;
                candidates.push_back(cand);
            }
        }

        std::sort(candidates.begin(), candidates.end(), [](const DynamicCandidate& a, const DynamicCandidate& b) {
            if (a.x_weight != b.x_weight) return a.x_weight < b.x_weight;
            if (a.x_bits != b.x_bits) return a.x_bits < b.x_bits;
            if (a.c1 != b.c1) return a.c1 < b.c1;
            return a.c2 < b.c2;
        });

        DynamicChiState state;
        std::vector<int> affected_rows;
        bool state_ready = false;

        for (int idx = 0; idx < (int)candidates.size() && !found; ++idx) {
            const DynamicCandidate& pair = candidates[idx];
            total_pairs_tested++;

            auto s_all = std::chrono::high_resolution_clock::now();
            if (!state_ready) {
                dynamic_build_compressed_chi_state(state, A, pair.x_bits, n, this_m);
                state_ready = true;
                total_rebuilds++;
            } else {
                bool rebuilt_basis = false;
                dynamic_update_chi_state_local_repair_compressed(state, A, pair.x_bits, n, this_m, affected_rows, rebuilt_basis);
                total_affected_rows += (long long)affected_rows.size();
                if (rebuilt_basis) total_rebuilds++;
                else total_local_add_rounds++;
            }
            auto e_all = std::chrono::high_resolution_clock::now();
            total_chi_update_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_all - s_all);
            total_basis_work_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_all - s_all);

            long long active_now = 0;
            for (int r = 0; r < (int)state.row_active.size(); ++r) active_now += state.row_active[r] ? 1 : 0;
            total_active_rows += active_now;
            total_active_snapshots++;

            std::vector<unsigned char> e_vec(this_m, 0);
            e_vec[pair.c1] = 1;
            e_vec[pair.c2] = 1;

            if (dynamic_membership_test(e_vec, state.basis)) {
                total_filtered_miss++;
                continue;
            }

            auto s_ns = std::chrono::high_resolution_clock::now();
            std::vector<std::vector<unsigned char>> ns_basis = dynamic_nullspace_basis(state.chi_rows, state.row_active, this_m);
            auto e_ns = std::chrono::high_resolution_clock::now();
            total_ns_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_ns - s_ns);
            total_ns_runs++;

            int good_idx = -1;
            for (int h = 0; h < (int)ns_basis.size(); ++h) {
                if ((ns_basis[h][pair.c1] ^ ns_basis[h][pair.c2]) == 1) {
                    good_idx = h;
                    break;
                }
            }
            if (good_idx < 0) continue;

            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < this_m; ++j) {
                    Anew[i][j] = (bool)((A[i][j] + pair.x_bits[i] * ns_basis[good_idx][j]) % 2);
                }
            }

            int mp = 0;
            GateSynthesisMatrix::cleanup(Anew, n, this_m, mp);
            if (mp < this_m) {
                std::cout << "  [HIT!] Pair(" << pair.c1 << "," << pair.c2 << ") x-weight=" << pair.x_weight
                          << " | " << this_m << " -> " << mp << " columns" << std::endl;
                LCL_Mat_GF2::copy((const bool**)Anew, n, mp, Abest);
                m_best = mp;
                found = true;
            }
        }

        if (found) {
            LCL_Mat_GF2::copy((const bool**)Abest, n, m_best, A);
            this_m = m_best;
        }
        round++;
    }

    omp = this_m;
    auto end_total = std::chrono::high_resolution_clock::now();
    auto total_dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);

    std::cout << "\n=== Dynamic Basis Local Repair Chi Summary ===" << std::endl;
    std::cout << "Algorithm       : LempelX2_DynamicBasisLocalRepairChi" << std::endl;
    std::cout << "Initial T-count : " << initial_m << std::endl;
    std::cout << "Final T-count   : " << omp << std::endl;
    std::cout << "Total Reduced   : " << (initial_m - omp) << " gates" << std::endl;
    std::cout << "Execution Time  : " << total_dur.count() << " ms" << std::endl;
    std::cout << "Pairs tested    : " << total_pairs_tested << std::endl;
    std::cout << "Filtered miss   : " << total_filtered_miss << std::endl;
    std::cout << "Nullspace runs  : " << total_ns_runs << std::endl;
    std::cout << "Rebuild count   : " << total_rebuilds << std::endl;
    std::cout << "Local add count : " << total_local_add_rounds << std::endl;
    std::cout << "Chi update time : " << total_chi_update_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Basis work time : " << total_basis_work_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Nullspace time  : " << total_ns_duration.count() / 1000.0 << " ms" << std::endl;
    if ((total_rebuilds + total_local_add_rounds) > 1) {
        std::cout << "Avg affected    : " << (double)total_affected_rows / (double)(total_rebuilds + total_local_add_rounds - 1) << " rows" << std::endl;
    }
    if (total_active_snapshots > 0) {
        std::cout << "Avg active rows : " << (double)total_active_rows / (double)total_active_snapshots << std::endl;
    }
    std::cout << "============================================" << std::endl;

    LCL_Mat_GF2::destruct(Anew, n, m + 1);
    LCL_Mat_GF2::destruct(Abest, n, m + 1);
}

void GateSynthesisMatrix::LempelX2_DynamicBasisLocalRepairChiPacked(bool** A, int n, int m, int& omp) {
    auto start_total = std::chrono::high_resolution_clock::now();
    std::cout << "\n[Dynamic Basis Local Repair TODD + Chi Compression + Packed] Compare with lx2_dynamic_repair_chi" << std::endl;

    int this_m = m;
    int initial_m = m;
    bool** Anew = LCL_Mat_GF2::construct(n, m + 1);
    bool** Abest = LCL_Mat_GF2::construct(n, m + 1);
    LCL_Mat_GF2::copy((const bool**)A, n, m, Abest);
    int m_best = m;

    std::chrono::microseconds total_chi_update_duration(0);
    std::chrono::microseconds total_basis_work_duration(0);
    std::chrono::microseconds total_ns_duration(0);

    long long total_pairs_tested = 0;
    long long total_filtered_miss = 0;
    long long total_ns_runs = 0;
    long long total_rebuilds = 0;
    long long total_local_add_rounds = 0;
    long long total_affected_rows = 0;
    long long total_active_rows = 0;
    long long total_active_snapshots = 0;

    bool found = true;
    int round = 0;
    while (found && (round < m)) {
        found = false;
        std::cout << "--- Round " << round << " | Current Columns: " << this_m << " ---" << std::endl;

        std::vector<DynamicCandidate> candidates;
        candidates.reserve(this_m * (this_m - 1) / 2);
        for (int j1 = 0; j1 < this_m; ++j1) {
            for (int j2 = j1 + 1; j2 < this_m; ++j2) {
                std::vector<unsigned char> x_bits(n, 0);
                int x_weight = 0;
                for (int i = 0; i < n; ++i) {
                    x_bits[i] = (unsigned char)((A[i][j1] + A[i][j2]) % 2);
                    x_weight += x_bits[i];
                }
                DynamicCandidate cand;
                cand.c1 = j1;
                cand.c2 = j2;
                cand.x_weight = x_weight;
                cand.x_bits = x_bits;
                candidates.push_back(cand);
            }
        }

        std::sort(candidates.begin(), candidates.end(), [](const DynamicCandidate& a, const DynamicCandidate& b) {
            if (a.x_weight != b.x_weight) return a.x_weight < b.x_weight;
            if (a.x_bits != b.x_bits) return a.x_bits < b.x_bits;
            if (a.c1 != b.c1) return a.c1 < b.c1;
            return a.c2 < b.c2;
        });

        PackedChiState state;
        std::vector<int> affected_rows;
        bool state_ready = false;

        for (int idx = 0; idx < (int)candidates.size() && !found; ++idx) {
            const DynamicCandidate& pair = candidates[idx];
            total_pairs_tested++;

            auto s_all = std::chrono::high_resolution_clock::now();
            if (!state_ready) {
                packed_build_compressed_chi_state(state, A, pair.x_bits, n, this_m);
                state_ready = true;
                total_rebuilds++;
            } else {
                bool rebuilt_basis = false;
                packed_update_chi_state_local_repair_compressed(state, A, pair.x_bits, n, this_m, affected_rows, rebuilt_basis);
                total_affected_rows += (long long)affected_rows.size();
                if (rebuilt_basis) total_rebuilds++;
                else total_local_add_rounds++;
            }
            auto e_all = std::chrono::high_resolution_clock::now();
            total_chi_update_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_all - s_all);
            total_basis_work_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_all - s_all);

            long long active_now = 0;
            for (int r = 0; r < (int)state.row_active.size(); ++r) active_now += state.row_active[r] ? 1 : 0;
            total_active_rows += active_now;
            total_active_snapshots++;

            if (packed_membership_test_pair(pair.c1, pair.c2, this_m, state.basis)) {
                total_filtered_miss++;
                continue;
            }

            auto s_ns = std::chrono::high_resolution_clock::now();
            std::vector<std::vector<unsigned char>> ns_basis = packed_nullspace_basis(state.chi_rows, state.row_active, this_m);
            auto e_ns = std::chrono::high_resolution_clock::now();
            total_ns_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_ns - s_ns);
            total_ns_runs++;

            int good_idx = -1;
            for (int h = 0; h < (int)ns_basis.size(); ++h) {
                if ((ns_basis[h][pair.c1] ^ ns_basis[h][pair.c2]) == 1) {
                    good_idx = h;
                    break;
                }
            }
            if (good_idx < 0) continue;

            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < this_m; ++j) {
                    Anew[i][j] = (bool)((A[i][j] + pair.x_bits[i] * ns_basis[good_idx][j]) % 2);
                }
            }

            int mp = 0;
            GateSynthesisMatrix::cleanup(Anew, n, this_m, mp);
            if (mp < this_m) {
                std::cout << "  [HIT!] Pair(" << pair.c1 << "," << pair.c2 << ") x-weight=" << pair.x_weight
                          << " | " << this_m << " -> " << mp << " columns" << std::endl;
                LCL_Mat_GF2::copy((const bool**)Anew, n, mp, Abest);
                m_best = mp;
                found = true;
            }
        }

        if (found) {
            LCL_Mat_GF2::copy((const bool**)Abest, n, m_best, A);
            this_m = m_best;
        }
        round++;
    }

    omp = this_m;
    auto end_total = std::chrono::high_resolution_clock::now();
    auto total_dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);

    std::cout << "\n=== Dynamic Basis Local Repair Chi Packed Summary ===" << std::endl;
    std::cout << "Algorithm       : LempelX2_DynamicBasisLocalRepairChiPacked" << std::endl;
    std::cout << "Initial T-count : " << initial_m << std::endl;
    std::cout << "Final T-count   : " << omp << std::endl;
    std::cout << "Total Reduced   : " << (initial_m - omp) << " gates" << std::endl;
    std::cout << "Execution Time  : " << total_dur.count() << " ms" << std::endl;
    std::cout << "Pairs tested    : " << total_pairs_tested << std::endl;
    std::cout << "Filtered miss   : " << total_filtered_miss << std::endl;
    std::cout << "Nullspace runs  : " << total_ns_runs << std::endl;
    std::cout << "Rebuild count   : " << total_rebuilds << std::endl;
    std::cout << "Local add count : " << total_local_add_rounds << std::endl;
    std::cout << "Chi update time : " << total_chi_update_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Basis work time : " << total_basis_work_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Nullspace time  : " << total_ns_duration.count() / 1000.0 << " ms" << std::endl;
    if ((total_rebuilds + total_local_add_rounds) > 1) {
        std::cout << "Avg affected    : " << (double)total_affected_rows / (double)(total_rebuilds + total_local_add_rounds - 1) << " rows" << std::endl;
    }
    if (total_active_snapshots > 0) {
        std::cout << "Avg active rows : " << (double)total_active_rows / (double)total_active_snapshots << std::endl;
    }
    std::cout << "===================================================" << std::endl;

    LCL_Mat_GF2::destruct(Anew, n, m + 1);
    LCL_Mat_GF2::destruct(Abest, n, m + 1);
}

void GateSynthesisMatrix::LempelX2_DynamicBasisLocalRepairChiPackedMemo(bool** A, int n, int m, int& omp) {
    auto start_total = std::chrono::high_resolution_clock::now();
    std::cout << "\n[Dynamic Basis Local Repair TODD + Chi Compression + Packed + Memo] Compare with lx2_dynamic_repair_chi_packed" << std::endl;

    int this_m = m;
    int initial_m = m;
    bool** Anew = LCL_Mat_GF2::construct(n, m + 1);
    bool** Abest = LCL_Mat_GF2::construct(n, m + 1);
    LCL_Mat_GF2::copy((const bool**)A, n, m, Abest);
    int m_best = m;

    std::chrono::microseconds total_chi_update_duration(0);
    std::chrono::microseconds total_basis_work_duration(0);
    std::chrono::microseconds total_ns_duration(0);

    long long total_pairs_tested = 0;
    long long total_filtered_miss = 0;
    long long total_ns_runs = 0;
    long long total_rebuilds = 0;
    long long total_local_add_rounds = 0;
    long long total_affected_rows = 0;
    long long total_active_rows = 0;
    long long total_active_snapshots = 0;
    long long total_memo_x_hits = 0;
    long long total_memo_ns_hits = 0;

    bool found = true;
    int round = 0;
    while (found && (round < m)) {
        found = false;
        std::cout << "--- Round " << round << " | Current Columns: " << this_m << " ---" << std::endl;

        std::vector<DynamicCandidate> candidates;
        candidates.reserve(this_m * (this_m - 1) / 2);
        for (int j1 = 0; j1 < this_m; ++j1) {
            for (int j2 = j1 + 1; j2 < this_m; ++j2) {
                std::vector<unsigned char> x_bits(n, 0);
                int x_weight = 0;
                for (int i = 0; i < n; ++i) {
                    x_bits[i] = (unsigned char)((A[i][j1] + A[i][j2]) % 2);
                    x_weight += x_bits[i];
                }
                DynamicCandidate cand;
                cand.c1 = j1;
                cand.c2 = j2;
                cand.x_weight = x_weight;
                cand.x_bits = x_bits;
                candidates.push_back(cand);
            }
        }

        std::sort(candidates.begin(), candidates.end(), [](const DynamicCandidate& a, const DynamicCandidate& b) {
            if (a.x_weight != b.x_weight) return a.x_weight < b.x_weight;
            if (a.x_bits != b.x_bits) return a.x_bits < b.x_bits;
            if (a.c1 != b.c1) return a.c1 < b.c1;
            return a.c2 < b.c2;
        });

        PackedChiState state;
        std::vector<int> affected_rows;
        bool state_ready = false;
        bool has_cached_ns = false;
        std::vector<std::vector<unsigned char>> cached_ns_basis;

        for (int idx = 0; idx < (int)candidates.size() && !found; ++idx) {
            const DynamicCandidate& pair = candidates[idx];
            total_pairs_tested++;

            bool same_x = state_ready && (pair.x_bits == state.x_prev);
            auto s_all = std::chrono::high_resolution_clock::now();
            if (!state_ready) {
                packed_build_compressed_chi_state(state, A, pair.x_bits, n, this_m);
                state_ready = true;
                total_rebuilds++;
                has_cached_ns = false;
            } else if (!same_x) {
                bool rebuilt_basis = false;
                packed_update_chi_state_local_repair_compressed(state, A, pair.x_bits, n, this_m, affected_rows, rebuilt_basis);
                total_affected_rows += (long long)affected_rows.size();
                if (rebuilt_basis) total_rebuilds++;
                else total_local_add_rounds++;
                has_cached_ns = false;
            } else {
                total_memo_x_hits++;
            }
            auto e_all = std::chrono::high_resolution_clock::now();
            total_chi_update_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_all - s_all);
            total_basis_work_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_all - s_all);

            long long active_now = 0;
            for (int r = 0; r < (int)state.row_active.size(); ++r) active_now += state.row_active[r] ? 1 : 0;
            total_active_rows += active_now;
            total_active_snapshots++;

            if (packed_membership_test_pair(pair.c1, pair.c2, this_m, state.basis)) {
                total_filtered_miss++;
                continue;
            }

            std::vector<std::vector<unsigned char>> ns_basis;
            if (same_x && has_cached_ns) {
                ns_basis = cached_ns_basis;
                total_memo_ns_hits++;
            } else {
                auto s_ns = std::chrono::high_resolution_clock::now();
                ns_basis = packed_nullspace_basis(state.chi_rows, state.row_active, this_m);
                auto e_ns = std::chrono::high_resolution_clock::now();
                total_ns_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_ns - s_ns);
                total_ns_runs++;
                cached_ns_basis = ns_basis;
                has_cached_ns = true;
            }

            int good_idx = -1;
            for (int h = 0; h < (int)ns_basis.size(); ++h) {
                if ((ns_basis[h][pair.c1] ^ ns_basis[h][pair.c2]) == 1) {
                    good_idx = h;
                    break;
                }
            }
            if (good_idx < 0) continue;

            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < this_m; ++j) {
                    Anew[i][j] = (bool)((A[i][j] + pair.x_bits[i] * ns_basis[good_idx][j]) % 2);
                }
            }

            int mp = 0;
            GateSynthesisMatrix::cleanup(Anew, n, this_m, mp);
            if (mp < this_m) {
                std::cout << "  [HIT!] Pair(" << pair.c1 << "," << pair.c2 << ") x-weight=" << pair.x_weight
                          << " | " << this_m << " -> " << mp << " columns" << std::endl;
                LCL_Mat_GF2::copy((const bool**)Anew, n, mp, Abest);
                m_best = mp;
                found = true;
            }
        }

        if (found) {
            LCL_Mat_GF2::copy((const bool**)Abest, n, m_best, A);
            this_m = m_best;
        }
        round++;
    }

    omp = this_m;
    auto end_total = std::chrono::high_resolution_clock::now();
    auto total_dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);

    std::cout << "\n=== Dynamic Basis Local Repair Chi Packed Memo Summary ===" << std::endl;
    std::cout << "Algorithm       : LempelX2_DynamicBasisLocalRepairChiPackedMemo" << std::endl;
    std::cout << "Initial T-count : " << initial_m << std::endl;
    std::cout << "Final T-count   : " << omp << std::endl;
    std::cout << "Total Reduced   : " << (initial_m - omp) << " gates" << std::endl;
    std::cout << "Execution Time  : " << total_dur.count() << " ms" << std::endl;
    std::cout << "Pairs tested    : " << total_pairs_tested << std::endl;
    std::cout << "Filtered miss   : " << total_filtered_miss << std::endl;
    std::cout << "Nullspace runs  : " << total_ns_runs << std::endl;
    std::cout << "Rebuild count   : " << total_rebuilds << std::endl;
    std::cout << "Local add count : " << total_local_add_rounds << std::endl;
    std::cout << "Memo x hits     : " << total_memo_x_hits << std::endl;
    std::cout << "Memo ns hits    : " << total_memo_ns_hits << std::endl;
    std::cout << "Chi update time : " << total_chi_update_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Basis work time : " << total_basis_work_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Nullspace time  : " << total_ns_duration.count() / 1000.0 << " ms" << std::endl;
    if ((total_rebuilds + total_local_add_rounds) > 1) {
        std::cout << "Avg affected    : " << (double)total_affected_rows / (double)(total_rebuilds + total_local_add_rounds - 1) << " rows" << std::endl;
    }
    if (total_active_snapshots > 0) {
        std::cout << "Avg active rows : " << (double)total_active_rows / (double)total_active_snapshots << std::endl;
    }
    std::cout << "========================================================" << std::endl;

    LCL_Mat_GF2::destruct(Anew, n, m + 1);
    LCL_Mat_GF2::destruct(Abest, n, m + 1);
}

void GateSynthesisMatrix::LempelX2_DynamicBasisLocalRepairChiPackedM4RI(bool** A, int n, int m, int& omp) {
    auto start_total = std::chrono::high_resolution_clock::now();
    std::cout << "\n[Dynamic Basis Local Repair TODD + Chi Compression + Packed + M4RI] Compare with lx2_dynamic_repair_chi_packed" << std::endl;

    int this_m = m;
    int initial_m = m;
    bool** Anew = LCL_Mat_GF2::construct(n, m + 1);
    bool** Abest = LCL_Mat_GF2::construct(n, m + 1);
    LCL_Mat_GF2::copy((const bool**)A, n, m, Abest);
    int m_best = m;

    std::chrono::microseconds total_chi_update_duration(0);
    std::chrono::microseconds total_basis_work_duration(0);
    std::chrono::microseconds total_ns_duration(0);

    long long total_pairs_tested = 0;
    long long total_filtered_miss = 0;
    long long total_ns_runs = 0;
    long long total_rebuilds = 0;
    long long total_local_add_rounds = 0;
    long long total_affected_rows = 0;
    long long total_active_rows = 0;
    long long total_active_snapshots = 0;

    bool found = true;
    int round = 0;
    while (found && (round < m)) {
        found = false;
        std::cout << "--- Round " << round << " | Current Columns: " << this_m << " ---" << std::endl;

        std::vector<DynamicCandidate> candidates;
        candidates.reserve(this_m * (this_m - 1) / 2);
        for (int j1 = 0; j1 < this_m; ++j1) {
            for (int j2 = j1 + 1; j2 < this_m; ++j2) {
                std::vector<unsigned char> x_bits(n, 0);
                int x_weight = 0;
                for (int i = 0; i < n; ++i) {
                    x_bits[i] = (unsigned char)((A[i][j1] + A[i][j2]) % 2);
                    x_weight += x_bits[i];
                }
                DynamicCandidate cand;
                cand.c1 = j1;
                cand.c2 = j2;
                cand.x_weight = x_weight;
                cand.x_bits = x_bits;
                candidates.push_back(cand);
            }
        }

        std::sort(candidates.begin(), candidates.end(), [](const DynamicCandidate& a, const DynamicCandidate& b) {
            if (a.x_weight != b.x_weight) return a.x_weight < b.x_weight;
            if (a.x_bits != b.x_bits) return a.x_bits < b.x_bits;
            if (a.c1 != b.c1) return a.c1 < b.c1;
            return a.c2 < b.c2;
        });

        PackedChiState state;
        std::vector<int> affected_rows;
        bool state_ready = false;

        for (int idx = 0; idx < (int)candidates.size() && !found; ++idx) {
            const DynamicCandidate& pair = candidates[idx];
            total_pairs_tested++;

            auto s_all = std::chrono::high_resolution_clock::now();
            if (!state_ready) {
                packed_build_compressed_chi_state_m4ri(state, A, pair.x_bits, n, this_m);
                state_ready = true;
                total_rebuilds++;
            } else {
                bool rebuilt_basis = false;
                packed_update_chi_state_local_repair_compressed_m4ri(state, A, pair.x_bits, n, this_m, affected_rows, rebuilt_basis);
                total_affected_rows += (long long)affected_rows.size();
                if (rebuilt_basis) total_rebuilds++;
                else total_local_add_rounds++;
            }
            auto e_all = std::chrono::high_resolution_clock::now();
            total_chi_update_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_all - s_all);
            total_basis_work_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_all - s_all);

            long long active_now = 0;
            for (int r = 0; r < (int)state.row_active.size(); ++r) active_now += state.row_active[r] ? 1 : 0;
            total_active_rows += active_now;
            total_active_snapshots++;

            if (packed_membership_test_pair(pair.c1, pair.c2, this_m, state.basis)) {
                total_filtered_miss++;
                continue;
            }

            int d_ns = 0;
            auto s_ns = std::chrono::high_resolution_clock::now();
            bool** NS = packed_active_rows_to_m4ri_nullspace(state, this_m, d_ns);
            auto e_ns = std::chrono::high_resolution_clock::now();
            total_ns_duration += std::chrono::duration_cast<std::chrono::microseconds>(e_ns - s_ns);
            total_ns_runs++;

            int good_idx = -1;
            for (int h = 0; h < d_ns; ++h) {
                if ((NS[pair.c1][h] ^ NS[pair.c2][h]) == 1) {
                    good_idx = h;
                    break;
                }
            }
            if (good_idx >= 0) {
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < this_m; ++j) {
                        Anew[i][j] = (bool)((A[i][j] + pair.x_bits[i] * NS[j][good_idx]) % 2);
                    }
                }

                int mp = 0;
                GateSynthesisMatrix::cleanup(Anew, n, this_m, mp);
                if (mp < this_m) {
                    std::cout << "  [HIT!] Pair(" << pair.c1 << "," << pair.c2 << ") x-weight=" << pair.x_weight
                              << " | " << this_m << " -> " << mp << " columns" << std::endl;
                    LCL_Mat_GF2::copy((const bool**)Anew, n, mp, Abest);
                    m_best = mp;
                    found = true;
                }
            }
            if (NS) LCL_Mat_GF2::destruct(NS, this_m, d_ns);
        }

        if (found) {
            LCL_Mat_GF2::copy((const bool**)Abest, n, m_best, A);
            this_m = m_best;
        }
        round++;
    }

    omp = this_m;
    auto end_total = std::chrono::high_resolution_clock::now();
    auto total_dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);

    std::cout << "\n=== Dynamic Basis Local Repair Chi Packed M4RI Summary ===" << std::endl;
    std::cout << "Algorithm       : LempelX2_DynamicBasisLocalRepairChiPackedM4RI" << std::endl;
    std::cout << "Initial T-count : " << initial_m << std::endl;
    std::cout << "Final T-count   : " << omp << std::endl;
    std::cout << "Total Reduced   : " << (initial_m - omp) << " gates" << std::endl;
    std::cout << "Execution Time  : " << total_dur.count() << " ms" << std::endl;
    std::cout << "Pairs tested    : " << total_pairs_tested << std::endl;
    std::cout << "Filtered miss   : " << total_filtered_miss << std::endl;
    std::cout << "Nullspace runs  : " << total_ns_runs << std::endl;
    std::cout << "Rebuild count   : " << total_rebuilds << std::endl;
    std::cout << "Local add count : " << total_local_add_rounds << std::endl;
    std::cout << "Chi update time : " << total_chi_update_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Basis work time : " << total_basis_work_duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "Nullspace time  : " << total_ns_duration.count() / 1000.0 << " ms" << std::endl;
    if ((total_rebuilds + total_local_add_rounds) > 1) {
        std::cout << "Avg affected    : " << (double)total_affected_rows / (double)(total_rebuilds + total_local_add_rounds - 1) << " rows" << std::endl;
    }
    if (total_active_snapshots > 0) {
        std::cout << "Avg active rows : " << (double)total_active_rows / (double)total_active_snapshots << std::endl;
    }
    std::cout << "========================================================" << std::endl;

    LCL_Mat_GF2::destruct(Anew, n, m + 1);
    LCL_Mat_GF2::destruct(Abest, n, m + 1);
}
