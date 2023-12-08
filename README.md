# dsu_equation

UnionFind に 一次方程式 をのせます。
Weighted Union Find の上位互換のようなものです。
分数ライブラリが必要です。

## ドキュメント

### `dsu_equation(int N)`

長さ $N$ の数列 $A$ として初期化します。 $O(N)$ 時間

### `int leader(int x)`

$x$ が含まれるグループのリーダーを返します。$O(\log N)$ 時間

### `int size(int x)`

$x$ が含まれているグループのサイズを返します。$O(\log N)$ 時間

### `bool same(int x,int y)`

$x,y$ が含まれているグループが同じかどうかを返します。 $O(\log N)$ 時間

### `tuple<int,frac,frac> relationship(int x,int p=1,int q=0)`

$x$ 含まれているグループのリーダーを $root_x$ としたときに

$A[root_x]=p\cdot A[x]+q$ と表せられる $p,q$ を求め $root_x,p,q$ を返します。 $O(\log N)$ 時間

### `bool merge(int x,int y,int a,int b,int c)`

$a\cdot A[x]+b\cdot A[y]=c$ として $x,y$ をマージします。

矛盾しないなら `true` , しているなら `false` を返します。 $O(\log N)$ 時間

### `bool set_val(int x,frac f)`

$A[x] = f$ としてセットします。

これまでの情報に矛盾しないなら`true` , しているなら `false` を返します。 $O(\log N)$ 時間 

### `pair<bool,frac> specify(int x)`

$A[x]$ の値を求めます。一意に定まらないなら `{false,0}` を返します。

そうでないなら `{true,A[x]}` を返します。 $O(\log N)$ 時間

### `tuple<bool,frac,frac> solve_coef(int x,int y)`

$A[y]=p\cdot A[x]+q$ を満たす $(p,q)$ を返します。

一意に定まらないなら `{false,0,0}` を返します。

そうでないなら `{true,p,q}` を返します。 $O(\log N)$ 時間

## verify
- https://atcoder.jp/contests/abc328/submissions/48274354
- https://atcoder.jp/contests/abc087/submissions/48274365
- https://atcoder.jp/contests/typical90/submissions/48274381
- https://atcoder.jp/contests/abc320/submissions/48274314
