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
#include <algorithm>
#include <map>
// (他にもあればここに追加)

// ステップ3: プロジェクト固有のヘッダ
#include "GateSynthesisMatrix.h"
#include "LCL/LCL_Mat_GF2.h"
#include "LCL/Core/LCL_ConsoleOut.h"
#include "LCL/LCL_Int.h"

// ステップ4: using namespace は全ての #include の「後」に置く
using namespace std;
using namespace LCL_ConsoleOut;
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
                    if(x_a && x_b) res ^= A->rows[gamma][w];
                    if(x_b && x_c) res ^= A->rows[alpha][w];
                    if(x_c && x_a) res ^= A->rows[beta][w];

                    // 3. 2次項 (x_a * A[beta] * A[gamma] など)
                    // ワード同士の AND をとってから XOR
                    if(x_a) res ^= (A->rows[beta][w] & A->rows[gamma][w]);
                    if(x_b) res ^= (A->rows[gamma][w] & A->rows[alpha][w]);
                    if(x_c) res ^= (A->rows[alpha][w] & A->rows[beta][w]);

                    // 結果を書き込み
                    Aext->rows[row_idx][w] = res;
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
        LempelX2_M4RI(A_copy, n, m, omp_m4ri);
        end_total = std::chrono::high_resolution_clock::now();
        duration_total = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
        std::cout << "Total LempelX2_M4RI time: " << duration_total.count() << " ms" << endl;
        std::cout << "M4RI result size: " << omp_m4ri << " columns" << endl;
        LCL_Mat_GF2::destruct(A_copy, n, m + 1);
        return;
    }
    //Hamming距離版
    if(true){
        std::cout << "\n--- Testing Hamming Distance Version ---" << endl;
        // 元の行列Aから新しいコピーを作成 (公平な比較のため)
        bool** A_copy = LCL_Mat_GF2::construct(n, m + 1);
        LCL_Mat_GF2::copy((const bool**)A, n, m, A_copy);
        
        int omp_hamming = m;

        // 時間計測開始
        start_total = std::chrono::high_resolution_clock::now();
        GateSynthesisMatrix::LempelX2_M4RI_Hamming(A_copy, n, m, omp_hamming);
        end_total = std::chrono::high_resolution_clock::now();
        
        duration_total = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
        std::cout << "Total LempelX2_M4RI_Hamming time: " << duration_total.count() << " ms" << endl;
        std::cout << "Hamming result size: " << omp_hamming << " columns" << endl;

        // メモリ解放
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
                      for(int c=0; c<chi_A->width; c++) chi_A->rows[r][c] = 0;
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
    std::cout << "\n====================================" << std::endl;
    std::cout << "nullspace time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(g_total_nullspace_duration).count()
              << " ms" 
              << ", chi time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(total_chi_duration).count()
              << " ms" << std::endl;
    std::cout << "Total attempts: " << total_attempts << ", successful: " << successful_attempts << std::endl;
    
    std::cout << "------------------------------------" << std::endl;
    std::cout << "Initial T-count: " << initial_total_m << std::endl;
    std::cout << "Final T-count  : " << this_m << std::endl;
    std::cout << "Total Reduced  : " << (initial_total_m - this_m) << " gates" << std::endl;
    std::cout << "====================================" << std::endl;

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
                    for(int c=0; c<chi_A->width; c++) chi_A->rows[r][c] = 0;
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