# TODD 動的基底更新: 局所修復案メモ

## 目的

`lx2_dynamic` の現在の安全版では、

- `chi` の更新は `AffectedRows` のみ差分更新
- しかし行空間基底 `basis` は毎回フル再構築

となっている。

この方法でも `miss` フィルタ自体は非常に有効だが、実時間の大半が `basis` の再構築に使われている。

そこで次の段階として、

- 壊れた基底だけを局所的に修復し
- 難しい場合のみ full rebuild にフォールバックする

という方針を整理する。

---

## 現在の安全版の問題点

現在の `lx2_dynamic` は、

1. `x_prev -> x_cur` から `AffectedRows` を求める
2. `AffectedRows` に属する `chi` の行だけ更新する
3. 更新後の `chi_rows` 全体から basis を毎回再構築する
4. `e_{c1,c2}` の membership 判定で `miss` を落とす

という流れになっている。

このため、

- `chi` の差分更新による恩恵はある
- `nullspace` 計算回数も大きく減る

一方で、

- basis 再構築の計算量が依然として大きい

という状態になる。

特に `adder_8` のような大きい例では、現在の主ボトルネックが `Basis rebuild` である。

---

## 新しいアイデアの要点

基本発想は単純である。

- 変わった行だけが basis に影響を与える
- したがって、basis 全体を毎回作り直すのではなく
- 変わった行に関係する部分だけを修理できればよい

このとき鍵になるのは以下の 3 つである。

- `affected_rows`
  - 今回変わった `chi` の行番号集合
- `row_to_basis_slot`
  - ある行が現在 basis に入っているかどうか
- `basis`
  - 現在の行空間基底そのもの

これらが分かれば、

- 変わった行が basis 行だったか
- 非基底行だったか

をすぐ判定できる。

---

## 方針の全体像

### 基本思想

完全な動的 RREF をいきなり目指すのではなく、

- まず局所修復を試す
- うまくいかなければ full rebuild へ戻す

という安全なハイブリッド構成にする。

### 期待する効果

- basis 行が壊れていないケースでは、再構築を回避できる
- 壊れていても軽微な場合は局所修復で済む
- 難しいケースのみ full rebuild に落とせる

---

## 差分修復アルゴリズムの考え方

### ステップ A: 汚染チェック

まず `affected_rows` を走査し、各行について

- `row_to_basis_slot[r] != -1`

かどうかを確認する。

これにより、その行が

- basis 行だった
- 非基底行だった

のどちらかを判定する。

#### 解釈

- basis 行が 1 本も壊れていない
  - 大きな再構築は不要な可能性が高い
- basis 行が壊れている
  - pivot を支える行が失われるので修復が必要

---

## ステップ B: 基底の外科手術

### 現在の安全版との違い

現在の安全版では、少しでも行が変わると basis 全体をフル再構築する。

局所修復版ではそうしない。

やりたいことは、

- **affected な basis 行だけ落とす**
- 残り basis はなるべくそのまま維持する
- 新しい行を basis に再投入する
- 矛盾が出たところだけ直す

である。

### なぜ「全 basis を再正規化」してはいけないのか

もし毎回、

- 残った basis 全部を生データから再掃き出し
- きれいな形へ並べ直す

としてしまうと、結局 small rebuild に近づいてしまう。

それでは full rebuild との差が小さくなる。

したがって本当に必要なのは、

- **影響を受けた pivot 周辺だけを再正規化する**

ことである。

### 実装イメージ

1. affected basis 行を一旦 basis から外す
2. その pivot 列を「空席」にする
3. 残った basis はそのまま保持する
4. 更新後の新しい行 `chi_rows[r]` を `add_to_basis` し直す
5. 衝突や不整合が出たところだけ局所的に整理する
6. 整理できなければ full rebuild に落とす

---

## ステップ C: ランクの回復

basis 行が壊れて消えると、

- ある pivot 列を支える行がいなくなり
- ランクが一時的に落ちる

可能性がある。

このとき、候補としては

- 更新された `affected_rows` 側の新行
- 非基底行プール

を使う。

### 初版の現実的方針

最初から大規模な候補探索はしない方がよい。

まずは、

1. changed rows から回復を試す
2. それで埋まらなければ full rebuild

でも十分である。

必要なら次段階で、

- 列ごとの候補リスト
- 転置インデックス

を導入する。

---

## 非基底行の変化も重要

注意点として、変わった行が basis でなかった場合でも、

- その新しい行が basis に入れるようになる
- 新しい pivot 候補になる

ことがある。

したがって、

- 「basis 行が壊れていないから何もしない」

は危険である。

この場合は最低限、

- changed 非基底行を basis に `add_to_basis` してみる

必要がある。

---

## なぜ row-echelon basis を勧めるのか

局所修復を行うなら、厳密な RREF よりも

- row-echelon basis

の方が扱いやすい。

理由:

- basis 行 1 本の削除時に影響が広がりにくい
- 完全な双方向掃き出しを維持しなくてよい
- membership 判定には十分

つまり、

- 完全に美しい基底

よりも

- 動的更新しやすい基底

を持つ方が有利である。

---

## 推奨する安全な実装順

### 段階 1

- `affected_rows` に basis 行が含まれるか数える
- basis 行が 0 本なら full rebuild を避ける

### 段階 2

- changed 非基底行だけ `add_to_basis` する

### 段階 2.5: 比較用の保守的局所修復版

まずは完全な局所修復ではなく、比較しやすい単純版を別関数として実装する。

方針は以下の通り。

- `affected_rows` の中に basis 行が 1 本でも含まれる
  - full rebuild
- `affected_rows` の中に basis 行が 1 本も含まれない
  - rebuild はしない
  - 更新後の changed rows だけ `add_to_basis` を試す

この版の利点は、

- 実装が単純
- 安全版 `lx2_dynamic` と比較しやすい
- 「basis 行が壊れていないケースだけを差分処理する」効果を単独で測定できる

---

## 比較用局所修復版のアルゴリズム

### 実行タグ

- `lx2_dynamic_repair`
  - 従来の保守的版
  - basis 行が 1 本でも touched したら full rebuild
- `lx2_dynamic_repair_k1`
  - `affected_basis_count == 1` のときだけ局所修復を試す版
  - 失敗したら full rebuild にフォールバック
- `lx2_dynamic_repair_chi`
  - 保守的 repair 版に `χ` 行圧縮を加えた版
  - `x_alpha = x_beta = x_gamma = 0` の行を inactive として basis / nullspace から除外
- `lx2_dynamic_repair_chi_packed`
  - `lx2_dynamic_repair_chi` に bit-parallel を加えた版
  - basis / membership / nullspace の行演算を `uint64_t` 単位で処理
- `lx2_dynamic_repair_chi_packed_memo`
  - `lx2_dynamic_repair_chi_packed` に同一 `x` のメモ化を加えた版
  - 連続して同じ `x` が現れた場合に `chi/basis` 更新と nullspace を再利用
- `lx2_dynamic_repair_chi_packed_m4ri`
  - `lx2_dynamic_repair_chi_packed` の full rebuild と nullspace を M4RI に委譲した版
  - local add は packed 実装のまま維持
- `lx2_dynamic_repair_chi_packed_basis_m4ri`
  - `lx2_dynamic_repair_chi_packed_m4ri` の nullspace 入力を active 行全体ではなく `state.basis` に縮めた版
  - full rebuild は引き続き M4RI を使用

## Experimental 比較用タグ

- `todd_exp_1110`
  - 現在の Experimental 実装
  - packed chi + AA table + memoization + M4RI nullspace
- `todd_exp_packedlocal_1110`
  - 比較用
  - packed chi + AA table + memoization までは同じ
  - nullspace だけ自前 packed 実装

### 判定規則

`affected_rows` を走査し、

- `row_to_basis_slot[r] != -1`

な行があるかどうかを調べる。

#### ケース A: basis 行が含まれる

- その更新は危険とみなし、full rebuild する

#### ケース B: basis 行が含まれない

- basis は壊れていないとみなす
- 更新後の changed rows を順に `add_to_basis` する

このとき、

- 0 に reduce されれば従属行
- 非零が残れば新しい basis 行

として扱う。

### 疑似コード

```text
function dynamic_update_chi_state_local_repair(state, x_cur):
    affected_rows = collect_affected_rows(x_prev, x_cur)

    for r in affected_rows:
        chi_rows[r] = build_new_row(x_cur, r)

    basis_touched = false
    for r in affected_rows:
        if row_to_basis_slot[r] != -1:
            basis_touched = true
            break

    if basis_touched:
        rebuild_basis_from_chi_rows(state)
        return rebuild

    for r in affected_rows:
        add_to_basis(chi_rows[r])

    return local_add_only
```

### この版の意味

この版は完全な局所修復ではないが、

- basis 行が壊れないケースを安く処理する
- basis 行が壊れたケースは安全側に倒す

という意味で、非常に実験しやすい。

特に、

- `basis 行が壊れる頻度`
- `壊れないケースでどれだけ rebuild を回避できるか`

を測定する上で有効である。

### 段階 3

- basis 行が 1〜2 本壊れた場合だけ局所修復を試す

### 段階 4

- 壊れた basis 行が多い場合は即 full rebuild

### 段階 5

- すべての判定を安全版と比較して正しさを確認する

---

## 擬似コード

```text
function dynamic_repair_basis(state, affected_rows):
    affected_basis_rows = []
    changed_nonbasis_rows = []

    for r in affected_rows:
        if row_to_basis_slot[r] != -1:
            affected_basis_rows.push(r)
        else:
            changed_nonbasis_rows.push(r)

    if affected_basis_rows is empty:
        for r in changed_nonbasis_rows:
            try add chi_rows[r] into basis
        return success

    if size(affected_basis_rows) is too large:
        return fallback_rebuild

    remove affected_basis_rows from basis
    mark their pivot columns as vacant

    for r in affected_rows:
        try add chi_rows[r] into basis

    if basis consistency check fails:
        return fallback_rebuild

    return success
```

ここで `basis consistency check` は、例えば

- pivot の重複がないか
- basis 行が非零か
- membership 判定に必要な構造が壊れていないか

を確認する処理である。

---

## フォールバック戦略

局所修復は万能ではない。

したがって、

- basis 行が多く壊れた
- 局所修復後に整合性が保てない
- 回復候補が見つからない

といった場合は、潔く

- `rebuild_basis_from_chi_rows`

に戻るべきである。

この fallback があることで、

- 正しさを壊しにくい
- 実験もやりやすい
- どのケースで局所修復が効くか測定できる

という利点がある。

---

## 評価

この局所修復案は、

- 現在の安全版から自然に進化できる
- 既存の `affected_rows` / `row_to_basis_slot` / `basis` を活用できる
- full rebuild の回数を減らせる可能性が高い

という意味で有望である。

一方で、

- 実装は難しくなる
- RREF をそのまま持つと壊れやすい
- 候補探索を欲張ると重くなる

という注意もある。

したがって、最初は

- 小さな局所修復
- 早めの fallback

で入れるのがよい。

---

## まとめ

この案の本質は、

- `chi` を差分更新するだけでなく
- `basis` も差分的に修理する

ことである。

ただし、いきなり完全動的基底更新を狙うのではなく、

- 変わった basis 行だけ外す
- 新しい行を足してみる
- ダメなら full rebuild

という保守的な局所修復から始めるのが現実的である。

---

## 現在の実装: `lx2_dynamic_repair_chi_packed_aa`

現在の比較実装の中で、性能が最も良い系の一つが

- `lx2_dynamic_repair_chi_packed_aa`

である。

この実装は、以下を同時に使っている。

- `chi` のゼロ行圧縮
- packed bit-parallel 表現
- `AA table` による `A_i AND A_j` の前計算
- 行空間基底による membership filter

### 全体の流れ

各ラウンドでは、まず現在の `A` から全候補列ペア `(c1, c2)` を作り、

- `x = A[:,c1] XOR A[:,c2]`

を求める。

候補順は現在

- `x_weight`
- packed `x_words`
- `c1`
- `c2`

の順でソートしている。

### `AA table` の構築

ラウンド先頭で現在の `A` に対して

- packed された各行 `A_rows`
- その全行ペア AND を持つ `AA_rows`

を一度だけ構築する。

ここで `AA_rows[idx(i,j)]` は

- `A_rows[i] AND A_rows[j]`

を表す。

このテーブルは `x` ではなく `A` にだけ依存するため、
同一ラウンド内では全候補で使い回せる。

### `chi` 行列の差分更新

最初の候補に対しては、圧縮済み packed `chi_rows` を full build する。

2 個目以降は

- `x_prev XOR x_cur`

から `AffectedRows` を求め、その行だけを更新する。

すなわち、

- `chi` を毎回全再構築するのではなく
- 影響行だけ `build_chi_row_packed_aa(...)` で上書きする

という構成である。

### basis 更新

現在の basis 更新はまだ保守的である。

- `AffectedRows` の中に basis 行があれば rebuild
- そうでなければ changed rows を `add_to_basis`

とする。

したがって、現状の主ボトルネックは

- `chi` 差分更新そのもの

ではなく、

- touched basis 行があったときの rebuild

である。

### この実装の意味

この版はすでに

- `AA table`
- `chi` 差分更新
- packed
- membership filter

を同時に持っている。

そのため、次の主戦場は

- `basis` の局所修復によって rebuild をどれだけ減らせるか

に移っている。

---

## 次の局所修復ロードマップ

`lx2_dynamic_repair_chi_packed_aa` をさらに速くするための、
現実的な次段階を以下に示す。

### フェーズ 1: 計測の細分化

まず、現在の rebuild の性質を細かく測る。

少なくとも以下を記録したい。

- `affected_basis_count == 0` の回数
- `affected_basis_count == 1` の回数
- `affected_basis_count >= 2` の回数
- rebuild 回数
- local add 回数

これにより、`k=1` 局所修復の期待値が見える。

### フェーズ 2: `affected_basis_count == 1` の局所修復

最初の本命はこれである。

方針は、

1. touched した basis 行が 1 本だけなら、その basis 行だけを temporary basis から除く
2. `AffectedRows` の更新後行を temporary basis に再投入する
3. rank が元に戻れば採用する
4. 戻らなければ full rebuild に落とす

というものである。

この段階では、完全な動的 RREF は狙わず、

- 小さい成功ケースだけ拾う
- 失敗したら即 rebuild

でよい。

### フェーズ 3: 探索順の `x_prev` 近傍化

現在の候補順は

- `x_weight`
- `x_words`

を重視しているが、局所修復を効かせるには

- `popcount(x_cur XOR x_prev)`

も小さくする順番が有効である。

これにより

- `AffectedRows`
- touched basis 行数

が減り、local add と局所修復の成功率が上がる可能性が高い。

### フェーズ 4: `affected_basis_count <= 2` への拡張

`k=1` 修復が効くなら、次に

- touched basis 行が 2 本以下

まで局所修復対象を広げる。

ただしこの段階では

- temporary basis の構築コスト
- 修復成功率

を慎重に見なければならない。

### フェーズ 5: fallback 条件の最適化

最終的には、

- どの条件で即 rebuild
- どの条件で局所修復を試す

を最適化する。

たとえば

- touched basis 行数
- `AffectedRows` の大きさ
- 現在の rank

などで閾値を決めることが考えられる。

---

## 現時点でのまとめ

現在の `lx2_dynamic_repair_chi_packed_aa` は、

- `chi` 側の高速化はかなり進んでいる
- `AA table` も有効
- packed nullspace も十分軽い

という状態である。

したがって、次の本命は

- basis rebuild を減らす局所修復

である。

実装順としては、

1. `affected_basis_count` の分布を取る
2. `affected_basis_count == 1` だけ局所修復する
3. 効果が出たら `<= 2` に広げる

のが最も自然である。

---

## 現在の主要関数一覧

ここでは、現在このリポジトリで比較対象として使っている主要関数と、
その役割をまとめる。

### 関数の見方

- `include/GateSynthesisMatrix.h`
  - 実験用の主要エントリ関数の宣言
- `source/GateSynthesisMatrix.cpp`
  - 実際のアルゴリズム本体
- `include/TO_Decoder.h`
  - CLI から使うアルゴリズム文字列
- `source/TO_Decoder.cpp`
  - アルゴリズム文字列から各関数への dispatch

したがって、ある関数を追うときは

1. `TO_Decoder.h` で文字列タグを確認する
2. `TO_Decoder.cpp` でどの関数が呼ばれるかを見る
3. `GateSynthesisMatrix.h` と `GateSynthesisMatrix.cpp` で本体を確認する

という順で読むと分かりやすい。

### TODD / LempelX2 系の主要エントリ

- `GateSynthesisMatrix::LempelX2(...)`
  - 元の LempelX2 実装

- `GateSynthesisMatrix::LempelX2_M4RI(...)`
  - M4RI を使う TODD 系の基本実装

- `GateSynthesisMatrix::LempelX2_M4RI_Hamming(...)`
  - ハミング距離ベースの候補順を使う版

- `GateSynthesisMatrix::LempelX2_M4RI_Experimental(...)`
  - `todd_exp_****` 系の本体
  - `packed chi`
  - `AA table`
  - `memoization`
  - `random sketch`
  を切り替えられる

- `GateSynthesisMatrix::LempelX2_M4RI_Experimental_PackedLocal(...)`
  - `todd_exp_packedlocal_****` の本体
  - `Experimental` と同じ枠組みだが、nullspace を自前 packed 実装で計算する

- `GateSynthesisMatrix::LempelX2_M4RI_Experimental_PackedLocal_ChiDiff(...)`
  - `todd_exp_packedlocal_diff_****` の本体
  - `PackedLocal` のうち `chi` を差分更新する比較版

### 位置づけのまとめ

- `LempelX2(...)`
  - 元の比較基準
- `todd_exp_****`
  - TODD 側の既存高速化枠組み
- `todd_exp_packedlocal_****`
  - TODD 側の packed local nullspace 比較群
- `lx2_dynamic_****`
  - 行空間基底による membership filter を導入した比較群

### 動的基底系

- `GateSynthesisMatrix::LempelX2_DynamicBasis(...)`
  - 動的基底追跡の最初の安全版
  - `chi` は差分更新するが、basis は毎回 full rebuild

- `GateSynthesisMatrix::LempelX2_DynamicBasisLocalRepair(...)`
  - basis touched がないときだけ local add を使う保守版

- `GateSynthesisMatrix::LempelX2_DynamicBasisLocalRepairK1(...)`
  - `affected_basis_count == 1` の局所修復を試す比較版

- `GateSynthesisMatrix::LempelX2_DynamicBasisLocalRepairChi(...)`
  - `chi` の圧縮を入れた bool 版

- `GateSynthesisMatrix::LempelX2_DynamicBasisLocalRepairChiPacked(...)`
  - `chi` 圧縮 + packed 版

- `GateSynthesisMatrix::LempelX2_DynamicBasisLocalRepairChiPackedAA(...)`
  - `chi` 圧縮 + packed + `AA table` 版
  - 現在の主力比較対象

- `GateSynthesisMatrix::LempelX2_DynamicBasisLocalRepairChiPackedMemo(...)`
  - `same x` の再利用を試す比較版

- `GateSynthesisMatrix::LempelX2_DynamicBasisLocalRepairChiPackedM4RI(...)`
  - full rebuild と nullspace の一部を M4RI に寄せた版

- `GateSynthesisMatrix::LempelX2_DynamicBasisLocalRepairChiPackedBasisM4RI(...)`
  - nullspace 入力として active rows 全体ではなく basis だけを M4RI に渡す版

---

## 主要 helper 関数

### `chi` 行生成

- `build_chi_row_bool(...)`
  - bool 版の `chi` 行生成

- `build_chi_row_packed(...)`
  - packed 版の `chi` 行生成

- `build_chi_row_packed_aa(...)`
  - packed + `AA table` 版の `chi` 行生成

### AA table

- `build_packed_a_rows(...)`
  - 現在の `A` を packed row 群へ変換

- `build_packed_aa_rows(...)`
  - `A_i AND A_j` を全行ペアについて前計算

- `packed_get_aa_word(...)`
  - `AA table` 参照 helper
  - `i == j` のときは `A_i` を返す

### `chi` 全体構築 / 差分更新

- `packed_build_compressed_chi_state(...)`
  - packed 圧縮版 `chi` の full build

- `packed_build_compressed_chi_state_aa(...)`
  - packed + `AA table` 版 `chi` の full build

- `packed_update_chi_state_local_repair_compressed(...)`
  - packed 圧縮版の `chi` 差分更新 + basis 更新

- `packed_update_chi_state_local_repair_compressed_aa(...)`
  - packed + `AA table` 版の `chi` 差分更新 + basis 更新

- `packed_build_compressed_chi_rows_aa(...)`
  - `PackedLocal_ChiDiff` 用の、basis を持たない `chi` full build

- `packed_update_chi_rows_compressed_aa(...)`
  - `PackedLocal_ChiDiff` 用の、basis を持たない `chi` 差分更新

### basis / membership / nullspace

- `packed_add_row_to_basis(...)`
  - packed basis への 1 行追加

- `packed_rebuild_basis_from_state(...)`
  - packed basis の full rebuild

- `packed_membership_test_pair(...)`
  - `e_{c1,c2}` が行空間に入るかの判定

- `packed_nullspace_basis(...)`
  - packed `chi` から自前で nullspace 基底を作る

---

## 実行タグとの対応

主な CLI タグとの対応は次のとおり。

- `todd_exp_1110`
  - `LempelX2_M4RI_Experimental(...)`

- `todd_exp_packedlocal_1110`
  - `LempelX2_M4RI_Experimental_PackedLocal(...)`

- `todd_exp_packedlocal_diff_1110`
  - `LempelX2_M4RI_Experimental_PackedLocal_ChiDiff(...)`

- `lx2_dynamic_repair_chi_packed`
  - `LempelX2_DynamicBasisLocalRepairChiPacked(...)`

- `lx2_dynamic_repair_chi_packed_aa`
  - `LempelX2_DynamicBasisLocalRepairChiPackedAA(...)`

- `lx2_dynamic_repair_chi_packed_memo`
  - `LempelX2_DynamicBasisLocalRepairChiPackedMemo(...)`

- `lx2_dynamic_repair_chi_packed_m4ri`
  - `LempelX2_DynamicBasisLocalRepairChiPackedM4RI(...)`

- `lx2_dynamic_repair_chi_packed_basis_m4ri`
  - `LempelX2_DynamicBasisLocalRepairChiPackedBasisM4RI(...)`
