# TODD 小型 chi 行列実験の方針

## 概要

このメモは，`todd_exp_packedlocal_diff_1110` に対して，小型 chi 行列を用いた高速化を試す方針をまとめたものです．

対象の実装は以下です．

```text
todd_exp_packedlocal_diff_1111_...
```

4 桁目のフラグを `1` にすると，小型 chi 行列による sketch 機能が有効になります．

## 目的

TODD の探索が進むと，削減条件を満たす列ペアが少なくなります．その結果，削減に失敗するペアに対してもフルサイズの chi 行列の零空間計算を行うことになり，ここが計算時間の大きな割合を占めます．

そこで，フル chi 行列の零空間を計算する前に，一部の行だけを抜き出した小型 chi 行列を調べ，不要なフル零空間計算を省略することを狙います．

## 数学的な根拠

フル chi 行列を `M`，その一部の行だけを取り出した小型 chi 行列を `M'` とします．

`M` の零空間に含まれるベクトルは，当然 `M'` の全ての行も満たします．したがって，

```text
Null(M) subset Null(M')
```

が成り立ちます．

このため，もし小型 chi 行列について

```text
Null(M') = {0}
```

すなわち，零ベクトル以外の解が存在しないことが分かれば，

```text
Null(M) = {0}
```

も必ず成り立ちます．

したがって，小型 chi 行列が列フルランクであれば，フル chi 行列の零空間計算を安全にスキップできます．

ただし逆は成り立ちません．小型 chi 行列で非自明な解 `y` が見つかっても，サンプルされなかった行で条件を破る可能性があります．そのため，小型 chi 行列から得られた `y` をそのまま採用してはいけません．必ず full chi に対して検証します．

## 現在のアルゴリズム

各候補ペアに対して，以下の流れで処理します．

1. packed chi 行列を構築または差分更新する．
2. sketch が有効なら，active な chi 行から一部をサンプリングする．
3. 小型 chi 行列の rank を計算する．
4. rank が現在の列数 `this_m` と等しければ，そのペアを棄却し，フル零空間計算をスキップする．
5. rank が `this_m` 未満なら，小型 chi 行列の零空間ベクトルを計算する．
6. 小型 chi 行列の零空間に，現在のペア `(c1, c2)` を分離できる方向があるかを調べる．具体的には，零空間基底の中に

```text
y[c1] xor y[c2] = 1
```

を満たすベクトルがあるかを見る．これが存在しない場合，フル chi の零空間にもそのような `y` は存在しないため，そのペアを棄却してフル零空間計算をスキップする．

7. 小型 chi 行列から得られた分離可能な各 `y` について，

```text
full_chi * y = 0
```

を確認する．

8. full chi の検証を通った場合だけ，

```text
Anew = A xor (x * y^T)
```

を試す．

9. `cleanup` によって重複列を削除し，T-count が減った場合のみ採用する．
10. sketch で削減できなければ，従来通りフル chi 行列の零空間計算に進む．

この流れでは，小型 chi 行列だけで得た `y` を無検証で採用しないため，正しさを保てます．

## 行数指定の方法

小型 chi 行列の行数は，3 種類の方法で指定できます．

### margin 指定

```text
todd_exp_packedlocal_diff_1111_m32
```

この場合，

```text
sample_rows = this_m + 32
```

となります．

`this_m` は現在の A 行列の列数です．TODD の探索中では，ほぼ現在の T-count に対応します．探索が進んで T-count が減ると，`this_m` も減ります．

### 固定行数指定

```text
todd_exp_packedlocal_diff_1111_r256
```

この場合，

```text
sample_rows = 256
```

となります．

`this_m` が変化しても，基本的に 256 行を使います．条件比較をしやすい指定方法です．

### パーセント指定

```text
todd_exp_packedlocal_diff_1111_p10
```

この場合，

```text
sample_rows = active chi rows の 10%
```

となります．

ただし，実装では下限として `this_m` 行を保証しています．理由は，`sample_rows < this_m` の場合，小型 chi 行列の rank は最大でも `sample_rows` までしか上がらず，`rank == this_m` を示せないためです．つまり，早期棄却には使えません．

## 出力する統計

sketch 機能が有効な場合，以下の統計を出力します．

```text
Sketch mode
Sketch attempts
Sketch skipped
Sketch rejected
Sketch rank rej
Sketch pair rej
Sketch passed
Sketch rows avg
Sketch y vectors
Sketch pair y
Sketch y full ok
Sketch y full ng
Sketch y reduced
Sketch y no red
Full NS skipped
Full NS computed
Sketch time
```

重要な項目は以下です．

- `Sketch rejected`
  - 小型 chi 行列により，フル零空間計算なしでペアを棄却できた合計回数．
- `Sketch rank rej`
  - 小型 chi 行列が列フルランクとなり，`Null(small chi) = {0}` と分かって棄却できた回数．
- `Sketch pair rej`
  - 小型 chi 行列の零空間内に `y[c1] xor y[c2] = 1` を満たす方向が存在せず，ペア削減が不可能と分かって棄却できた回数．
- `Full NS skipped`
  - フル零空間計算を省略できた回数．
- `Sketch pair y`
  - 小型 chi 行列の零空間基底のうち，`y[c1] xor y[c2] = 1` を満たした候補数．
- `Sketch y full ok`
  - 小型 chi 行列から得た `y` が full chi 検証を通った回数．
- `Sketch y reduced`
  - その `y` により実際に T-count が削減できた回数．
- `Sketch time`
  - sketch 処理に追加でかかった時間．

## 実装上の注意点

`Sketch rank rej` は `Null(small chi) = {0}` を示すため，現在の `x` に対して pair に依存せず安全にキャッシュできます．

一方，`Sketch pair rej` は `y[c1] xor y[c2] = 1` という条件を使うため，pair `(c1, c2)` に依存します．したがって，同じ `x` の別 pair に対して，pair reject の結果をキャッシュして再利用してはいけません．実装では pair reject した候補だけをスキップし，`ns_cached` には入れないようにします．

高速化できるかどうかは，`Full NS skipped` による節約時間が `Sketch time` を上回るかで決まります．

## gf2^7_mult での初期観察

以下は，`Sketch pair rej` を導入する前の初期実験です．

条件は以下です．

```text
todd_exp_packedlocal_diff_1111_p10
```

結果は正しく，T-count は従来と同じ 164 まで削減されました．

```text
Initial T-count : 217
Final T-count   : 164
Fail count      : 0
```

しかし，高速化には失敗しました．

```text
1110 baseline : 約 30.5 s
1111_p10      : 約 51.5 s
```

主な理由は以下です．

```text
Sketch rejected : 0
Full NS skipped : 0
Sketch time     : 10835 ms
```

つまり，`p10` では早期棄却が一度も発生せず，結局ほぼ全ての試行でフル零空間計算を行っています．そのため，sketch の計算時間がほぼ純粋な追加コストになりました．

一方で，小型 chi 行列から有効な `y` を見つけた例は 1 回ありました．

```text
Sketch y full ok : 1
Sketch y reduced : 1
```

このことから，早期棄却よりも「小型 chi 行列から有効な `y` を早く見つける」方向には可能性があります．

この結果を受けて，棄却条件を `Null(small chi) = {0}` だけでなく，

```text
Null(small chi) 内に y[c1] xor y[c2] = 1 を満たす方向が存在しない
```

場合にも拡張しました．この条件は TODD のペア削減条件に直接対応しており，`Full NS empty = 0` の状況でもフル零空間計算を省略できる可能性があります．

## 今後の実験

以下のバッチファイルで，`gf2^7_mult` に対する複数条件の比較実験を行います．

```text
forme/batch_gf27_sketch.sh
```

実行する条件は以下です．

```text
todd_exp_packedlocal_diff_1110
todd_exp_packedlocal_diff_1111_m16
todd_exp_packedlocal_diff_1111_m32
todd_exp_packedlocal_diff_1111_m64
todd_exp_packedlocal_diff_1111_r128
todd_exp_packedlocal_diff_1111_r256
todd_exp_packedlocal_diff_1111_r512
todd_exp_packedlocal_diff_1111_p5
todd_exp_packedlocal_diff_1111_p10
todd_exp_packedlocal_diff_1111_p20
```

確認すべき点は以下です．

1. `Sketch rejected` と `Full NS skipped` が増える条件があるか．
2. 小型 chi 行列から有効な `y` を早期発見できるか．
3. `Sketch time` よりも，フル零空間計算の節約時間が大きいか．
4. 最終 T-count が baseline と一致するか．

この機能は，正しさを保ったまま `todd_exp_packedlocal_diff_1110` より実行時間が短くなった場合に有効と判断します．
