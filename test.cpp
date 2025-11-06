#include <iostream>

// C++からC言語の M4RI ヘッダを正しく読み込むために
extern "C" {
    #include <m4ri/m4ri.h>
}

int main() {
    // M4RIのデータ型 (mzd_t) を使えるか？
    mzd_t* A = mzd_init(10, 10);
    
    if (A == NULL) {
        std::cerr << "M4RI: mzd_init に失敗しました。" << std::endl;
        return 1;
    }

    // エラーが出ている関数 (mzd_get_rows) を使えるか？
    int rows = mzd_get_rows(A);

    std::cout << "M4RIテスト成功。" << std::endl;
    std::cout << "mzd_get_rows(A) = " << rows << std::endl;

    // メモリを解放
    mzd_free(A);
    
    return 0;
}