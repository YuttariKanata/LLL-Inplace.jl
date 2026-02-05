# LLL-Inplace.jl

Juliaの標準的な `BigFloat` 演算は計算のたびに新しいオブジェクトを生成するため、LLLアルゴリズムのような反復計算ではGC（ガベージコレクション）がボトルネックになります。

この実装は、`libmpfr` と `libgmp` の関数を `ccall` で直接叩き、メモリを破壊的に再利用（In-place mutation）することで、GCのオーバーヘッドを極限まで排除した爆速のLLL実装です。

## 特徴
- **Zero-allocation loop**: メインループ内での新規メモリ確保を回避。
- **Direct MPFR/GMP calls**: Juliaのオーバーヘッドを通さず、Cレイヤーで計算。
- **Performance**: 標準的な実装に比べ、大幅な高速化を実現。

## 使い方
`LLL4.jl` を include して `LLL4` 関数を呼び出してください。

```julia
include("LLL4.jl")

# x: 入力ベクトル, C: 拡大係数
# b_int: 縮小基底（整数行列）, b_real: 実数残差, iter: 反復回数
b_int, b_real, iter = LLL4(x, C)
```
