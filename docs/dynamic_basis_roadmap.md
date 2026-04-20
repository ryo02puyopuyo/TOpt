# TODD 動的基底更新 ロードマップ

## 目的

同一ラウンド内で `A` を固定したまま、TODD の候補ペア探索を高速化する新しい実装経路を追加する。

新しい実装では、

- ラウンド内では `A` を固定する
- ペアごとに変化する `x` に応じて `chi(A, x)` を差分更新する
- `chi` の行空間基底を自前で管理する
- その基底を使って `miss` 判定だけを高速化する
- `hit` 候補に対してのみ重い零空間計算を行う

ことを目指す。

これは既存の `LempelX2_M4RI_Experimental` を直接書き換えるのではなく、**新しい関数として切り出して実装する**。

---

## スコープ

### 今回の対象

- 同じラウンド内での最適化
- `A` が固定されている区間での工夫
- 新しい helper 関数・struct の追加
- GF(2) 上の自前の行空間基底管理
- 動的更新が難しい場合のフル再構築フォールバック

### 初版では対象外

- `A <- A xor x v^T` の後のラウンドをまたいだ再利用
- 零空間計算そのものの完全置換
- フォールバックなしの完全動的 RREF
- 初版での random sketch との統合

---

## 基本方針

ペア `(c1, c2)` に対して

- `x = A[:,c1] xor A[:,c2]`
- `chi = chi(A, x)`

を考える。

本当に知りたいのは、

- `Null(chi)` に `v[c1] xor v[c2] = 1` を満たすベクトルが存在するか

である。

これを毎回重い零空間計算で調べる代わりに、まず `chi` の**行空間基底**を使って

- `e_{c1,c2}` が `RowSpace(chi)` に入るか

を判定する。

このとき、

- `e_{c1,c2}` が行空間に入るなら、そのペアは **definite miss**
- 行空間に入らないなら、そのペアは **hit の可能性あり**

と判定できる。

したがって、

- `miss` は自前基底で高速に落とす
- 生き残ったペアだけ既存の零空間計算へ回す

という構成にする。

---

## フェーズ 1: 新しい実行経路を作る

まずは新しい関数を追加する。

例:

`GateSynthesisMatrix::LempelX2_M4RI_DynamicBasisExperimental(...)`

初版の挙動は次のようにする。

- ペアの列挙は既存 experimental をベースにする
- ラウンド先頭の最初のペアだけは `chi` をフル構築する
- その `chi` から自前の行空間基底を構築する
- 後続のペアでは `x` の変化に応じて `chi` を更新し、基底も更新する
- 基底で `miss` 判定を行う
- 通過したペアだけ既存の M4RI 零空間計算を実行する

---

## フェーズ 2: 主要データ構造

### 基底行

```cpp
struct DynamicBasisRow {
    std::vector<uint64_t> bits;
    int pivot;
    int source_row_id;
};
```

### ラウンド状態

```cpp
struct DynamicChiState {
    std::vector<std::vector<uint64_t>> chi_rows;
    std::vector<char> row_active;
    std::vector<DynamicBasisRow> basis;
    std::vector<int> row_to_basis_slot;
    std::vector<int> pivot_to_basis_slot;
    std::vector<uint64_t> x_prev;
};
```

意味は以下の通り。

- `chi_rows`: 現在の `chi(A, x)` の各行
- `row_active`: packed chi のとき、その行が有効かどうか
- `basis`: 行空間基底
- `row_to_basis_slot`: 各行が basis に入っているか
- `pivot_to_basis_slot`: pivot 列ごとの basis 行対応
- `x_prev`: 直前ペアの `x`

---

## フェーズ 3: 先に実装すべき helper 関数

### ビット演算系

- `find_leftmost_one(bits)`
- `xor_row(dst, src)`
- `is_zero_row(bits)`
- `read_bit(bits, col)`
- `write_bit(bits, col, val)`

### chi 行インデックス

- `row_id(alpha, beta, gamma, n)`

必要なら逆変換も用意する。

### 単一行生成

- `build_chi_row(A_m4ri, x_bits, alpha, beta, gamma, m_words)`

`chi` 全体ではなく、**1行だけ**生成する関数を作る。

### 基底操作

- `reduce_against_basis(row, basis)`
- `add_row_to_basis(row, source_row_id, basis, maps)`
- `membership_test(target, basis)`

---

## フェーズ 4: まずは正しいフル再構築経路を作る

動的削除をいきなり実装するのではなく、まずは安全な基準実装を作る。

関数:

`rebuild_basis_from_chi_rows(state)`

この関数は

- basis を空にする
- active な `chi_rows` を順に読む
- 各非零行を `add_row_to_basis` で追加する

ことで行空間基底を構築する。

これにより、

- 正しさ確認用の基準
- 動的更新が失敗したときのフォールバック

を確保できる。

---

## フェーズ 5: 影響行集合 `AffectedRows` の計算

`x_prev -> x_cur` の変化に対して

- `delta_x = x_prev xor x_cur`

を計算し、`delta_x[k] = 1` の各 `k` について

- `alpha = k`
- `beta = k`
- `gamma = k`

を含む行を `AffectedRows` として集める。

helper:

`collect_affected_rows(delta_x, n, out_rows)`

注意点:

- 重複除去する
- 順序は決定的にする

---

## フェーズ 6: 動的更新の基本戦略

最初から完全動的 RREF を狙わない。

初期戦略は次の通り。

1. `AffectedRows` を求める
2. `chi_rows` 内の該当行を `old_row -> new_row` に差し替える
3. その affected row の中に basis 行があったかを調べる
4. basis 行が触られていなければ、差分追加だけ試す
5. basis 行が触られていたら、局所修復を試す
6. 局所修復が難しければ `rebuild_basis_from_chi_rows` へフォールバック

これで実装リスクをかなり抑えられる。

---

## フェーズ 7: basis 行削除時の局所修復

フル再構築版が安定したあとで、局所修復を入れる。

basis 行が失われたら、

1. basis からその行を削除する
2. 対応する pivot 列を空席にする
3. まず changed rows から代替候補を探す
4. 見つからなければ非基底行から代替候補を探す
5. 候補がなければランク減少を許容する
6. 最後に `new_row` 群を basis に追加する

### 重要な設計判断

初版では **厳密な RREF ではなく row-echelon basis を管理する**。

理由:

- 動的削除が簡単になる
- membership 判定にはそれで十分
- strict な双方向掃き出しは不要

---

## フェーズ 8: `miss` フィルタの使い方

ペア `(c1, c2)` に対して

- `e_{c1,c2}` を作る
- 現在の行空間基底で reduce する

その結果、

- 0 になる -> definite miss
- 0 にならない -> hit の可能性あり

と判定する。

このフィルタは **true hit を落としてはいけない**。

---

## フェーズ 9: 探索順の方針

動的基底更新にとって重要なのは、「今のペアが良いか」だけでなく「前のペアからどれだけ `x` が変わるか」である。

したがって、探索順は列ハミング距離だけで決めるのではなく、

1. `popcount(x)` が小さい
2. `popcount(x xor x_prev)` が小さい
3. 列ハミング距離が小さい

の順で優先した方がよい。

例えばスコアは

`score = (popcount(x_cur), popcount(x_cur xor x_prev), column_dist)`

とする。

### 理由

- `popcount(x_cur)` が小さい
  - packed chi の非零行数が減る
- `popcount(x_cur xor x_prev)` が小さい
  - 影響行が減る
  - basis の壊れ方が小さくなる
- `column_dist` が小さい
  - 既存 TODD のヒューリスティックとして意味がある

### 結論

今の `0010` 的な「`x` の popcount 重視」は方向としては良い。  
ただし、**動的基底更新を活かすには、それだけでなく連続ペア間の `x` の近さも必要**。

つまり、

- `x` 自体が軽い
- 前の `x` からあまり飛ばない

の両方を満たす順番が望ましい。

---

## フェーズ 10: 計測項目

以下を必ず記録する。

- 試したペア数
- 動的更新を試した回数
- フル再構築にフォールバックした回数
- basis 行削除が起きた回数
- `AffectedRows` の平均サイズ
- `popcount(x xor x_prev)` の平均
- `miss` フィルタで落とせた回数
- フィルタを通過して零空間計算に進んだ回数
- 動的更新時間
- 基底再構築時間
- 零空間計算時間

これを取らないと、効いているのかどうかが判断できない。

---

## フェーズ 11: 検証方針

### 正しさ

同じラウンド・同じ候補順で、

- 動的基底フィルタの判定
- 従来のフル零空間計算結果

を比較し、**true hit を誤って弾いていないこと**を確認する。

### 性能

比較対象:

- 既存 experimental
- 動的基底 + 毎回再構築
- 動的基底 + 局所修復

記録するもの:

- 最終 T-count
- 総実行時間
- 零空間計算回数
- フォールバック頻度

---

## フェーズ 12: 実装順

推奨順序は以下。

1. ロードマップ作成
2. helper type 宣言
3. 単一 `chi` 行生成関数
4. 自前の行空間基底構築関数
5. `e_{c1,c2}` の membership 判定
6. 更新後 `chi_rows` からのフル再構築版を新関数として実装
7. `AffectedRows` 計算
8. 動的更新 + フォールバックを実装
9. basis 行削除の局所修復を実装
10. 探索順の改善を入れる

---

## 現在の実装済みアルゴリズムの流れ

ここでは、**現在 `lx2_dynamic` として実装済みの安全版**が、1 ラウンド内で実際にどう動いているかを説明する。

重要なのは、

- `chi` は差分更新している
- しかし basis は差分修復していない
- basis は毎回フル再構築している

という点である。

### 1. ラウンド開始時

- 現在の `A` から、すべての候補ペア `(c1, c2)` を列挙する
- 各候補について `x = A[:,c1] xor A[:,c2]` を計算する
- `x` の `popcount` が小さい順に候補を並べる

### 2. ラウンド内最初のペア

最初のペアだけは、

- `x_prev = x_cur`
- `chi(A, x_cur)` を全行フル構築
- `chi_rows` に保存
- `chi_rows` 全体から行空間基底 `basis` を構築

する。

### 3. 2 個目以降のペア

次のペアでは、

1. `x_cur = A[:,c1] xor A[:,c2]` を作る
2. `delta_x = x_prev xor x_cur` を計算する
3. `delta_x` に応じて `AffectedRows` を求める
4. `AffectedRows` に属する `chi` の行だけ再計算して `chi_rows` を上書きする
5. **更新後の `chi_rows` 全体から basis を毎回再構築する**
6. `x_prev = x_cur` に更新する

となる。

### 4. `miss` 判定

今回のペアに対して

- `e_{c1,c2}` を作り
- basis に対して reduce する

その結果、

- 0 になれば definite miss として棄却
- 0 にならなければ次へ進む

### 5. 零空間計算

`miss` 判定を通過したペアに対してのみ、

- 現在の `chi_rows` から自前の GF(2) 零空間基底を計算
- `v[c1] xor v[c2] = 1` を満たすベクトルを探す

### 6. HIT 判定

そのような `v` があれば

- `Anew = A xor x v^T`
- `cleanup(Anew)`

を行い、列数が減れば `HIT` として

- `Abest` に保存
- そのラウンドを終了
- 次ラウンドへ進む

となる。

### 7. 現在のボトルネック

現在の安全版では、

- `AffectedRows` のみ差分更新しているにもかかわらず
- basis は毎回全 rebuild

なので、主な時間は

- `Chi update`
- `Basis rebuild`

に使われる。

したがって、次の改良点は

- basis 行削除時の局所修復
- 失敗時のみ full rebuild

である。

---

## 現在の実装済みアルゴリズムの疑似コード

### ラウンド全体

```text
function LempelX2_DynamicBasis(A):
    this_m = m(A)
    while improved:
        candidates = all column pairs (c1, c2)
        for each candidate:
            x = A[:,c1] xor A[:,c2]
            score = popcount(x)
        sort candidates by score

        state_ready = false
        found = false

        for each pair in candidates:
            if not state_ready:
                build_full_chi_state(pair.x)
                state_ready = true
            else:
                update_chi_state(pair.x)

            e = unit_pair_vector(pair.c1, pair.c2)
            if e is in RowSpace(chi):
                continue   # definite miss

            NS = nullspace(chi)
            find v in NS with v[c1] xor v[c2] = 1
            if no such v:
                continue

            Anew = A xor x v^T
            cleanup(Anew)
            if column count reduced:
                A = Anew
                found = true
                break

        if not found:
            stop
```

### 最初の `chi` と basis の構築

```text
function build_full_chi_state(x):
    for each row id r = (alpha, beta, gamma):
        chi_rows[r] = build_chi_row(A, x, alpha, beta, gamma)
        row_active[r] = 1

    basis = rebuild_basis_from_all_active_rows(chi_rows)
    x_prev = x
```

### 2 個目以降の `chi` 更新

```text
function update_chi_state(x_cur):
    delta_x = x_prev xor x_cur
    affected_rows = collect_affected_rows(delta_x)

    for each r in affected_rows:
        (alpha, beta, gamma) = decode_row_id(r)
        chi_rows[r] = build_chi_row(A, x_cur, alpha, beta, gamma)

    basis = rebuild_basis_from_all_active_rows(chi_rows)
    x_prev = x_cur
```

### `miss` 判定

```text
function membership_test_for_pair(c1, c2, basis):
    e = zero vector of length m
    e[c1] = 1
    e[c2] = 1

    e_reduced = reduce(e, basis)
    if e_reduced == 0:
        return definite_miss
    else:
        return possible_hit
```

### 自前零空間基底

```text
function dynamic_nullspace_basis(chi_rows):
    mat = active chi rows only
    perform GF(2) elimination
    record pivot columns

    for each free column f:
        vec[f] = 1
        for each pivot row i:
            vec[pivot_i] = mat[i][f]
        append vec to nullspace basis

    return nullspace basis
```

### 現在の実装の要点

現在の `lx2_dynamic` は次のように整理できる。

- `chi` は差分更新している
- basis は差分更新していない
- basis は毎回フル再構築している
- `miss` だけ基底で高速判定している
- 通過ペアにだけ零空間計算をしている

すなわち、**完全な動的基底更新の前段階としての安全版**である。

---

## 最終方針

初手は「安全なハイブリッド構成」にする。

- 自前 basis
- 正しい miss フィルタ
- 影響行の差分更新
- 難しいときはフル再構築

そのうえで、安定したら完全な動的基底更新へ進む。
