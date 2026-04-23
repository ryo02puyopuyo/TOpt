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

#ifndef HEADER_GATESYNTHESISMATRIX
#define HEADER_GATESYNTHESISMATRIX

namespace GateSynthesisMatrix {
    // Basics
    bool** from_signature(bool*** S, int n, int& m);
    void cleanup(bool** A, int n, int m, int& mp); // Sets pairwise col duplicates to zero, moves all zero cols to rhs, returns number of non-zero cols

    // Advanced
    void LempelX(bool** A, int n, int m, int& mp); // Form x = c_i + c_j for all col pairs. Extend matrix with rows x_a r_b r_c + x_b r_c r_a + x_c r_a r_b. Find nullspace. If N_i + N_j = 1, add x to all cols in nullspace, and eliminate c_i and c_j. Repeat until all c_i and c_j pairs have been exhausted.
    void Chi(bool** A, bool** x, int n, int m, bool** Aext); // Constructs the chi matrix from GSM, A. Aext must be n^3 x m in dimension.
    void ChiPrime(bool** A, bool** x, int n, int m, bool** Aext); // Constructs the chi matrix from GSM, A. Aext must be n^3 x m in dimension.
    void LempelX2(bool** A, int n, int m, int& mp);
    void LempelX3(bool** A, int n, int m, int& mp);
    
    bool** M4RI_wrapper_for_nullspace(bool const** A, int n, int m, int& out_d);
    void LempelX2_M4RI(bool** A, int n, int m, int& omp);
    void LempelX2_M4RI_DetailedStats(bool** A, int n, int m, int& omp);
    //void Chi_M4RI(bool** A, bool** x, int n, int m, mzd_t* Aext); // これが必要
    void Chi_M4RI(mzd_t* A, bool** x, int n, int m, mzd_t* Aext);
    void LempelX2_M4RI_Hamming(bool** A, int n, int m, int& omp);

    void LempelX2_M4RI_BeamSearch(bool** A, int n, int m, int& omp);
    void LempelX2_M4RI_RandomBeamSearch(bool** A_init, int n, int m_init, int& omp);
    void LempelX2_M4RI_SequentialBeamSearch(bool** A_init, int n, int m_init, int& omp);
    void LempelX2_M4RI_Hamming_Preprocess(bool** A, int n, int m, int& omp);
    void SparsifyAndTrack_Bool(bool** A, int n, int m, std::vector<std::pair<int, int>>& cnot_history);
    void LempelX2_M4RI_GreedyPreprocess(bool** A, int n, int m, int& omp);
    void LempelX2_M4RI_GreedyPreprocess_Fast(bool** A, int n, int m, int& omp);
    void LempelX2_M4RI_Experimental(bool** A, int n, int m, int& omp, bool use_packed_chi, bool use_aa_table, bool use_memoization, bool use_random_sketch, int sketch_margin);
    void LempelX2_M4RI_Experimental_PackedLocal(bool** A, int n, int m, int& omp, bool use_packed_chi, bool use_aa_table, bool use_memoization, bool use_random_sketch, int sketch_margin);
    void LempelX2_DynamicBasis(bool** A, int n, int m, int& omp);
    void LempelX2_DynamicBasisLocalRepair(bool** A, int n, int m, int& omp);
    void LempelX2_DynamicBasisLocalRepairK1(bool** A, int n, int m, int& omp);
    void LempelX2_DynamicBasisLocalRepairChi(bool** A, int n, int m, int& omp);
    void LempelX2_DynamicBasisLocalRepairChiPacked(bool** A, int n, int m, int& omp);
    void LempelX2_DynamicBasisLocalRepairChiPackedAA(bool** A, int n, int m, int& omp);
    void LempelX2_DynamicBasisLocalRepairChiPackedMemo(bool** A, int n, int m, int& omp);
    void LempelX2_DynamicBasisLocalRepairChiPackedM4RI(bool** A, int n, int m, int& omp);
    void LempelX2_DynamicBasisLocalRepairChiPackedBasisM4RI(bool** A, int n, int m, int& omp);
    void GreedyWeightReduction_Bool(bool** A, int n, int m, std::vector<std::pair<int, int>>& cnot_history);
}

#endif // HEADER_GATESYNTHESISMATRIX
