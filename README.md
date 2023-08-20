# dsu_equation

WeightedUnionFind木の拡張のようなものです。

$N$ 元 $1$ 次方程式を高速に扱います

分数ライブラリが必要です。別で用意してください。

なお分数のデータ構造をここでは frac と定義します

fracは四則演算と整数型での値の取得する関数 (ここではgetNumber)が使える必要があります

## ドキュメント

### `dsu_equation d(N)`

$N$ を初期化します $O(N)$ 時間です

### `int root(x)`

$x$ の root を返します $O(logN)$ 時間です


### `bool same(x,y)`

$A[x],A[y]$ お互いの値が関係しているか

つまり、ある$a(\not ={0}),b$を用いて $A[x]=a A[y] +b$ と表すことができるかを返します

$O(logN)$ 時間です


### `int size(x)`

$x$ の連結成分のサイズを返します $O(logN)$ 時間です

### `int d.merge(int x,int y,int a,int b,int c)`

$A[x]$ と $A[y]$ を $a A[x] + b A[y] =c$としてmergeします

それまでの merge に矛盾していれば $0$を返しmergeをやめます

矛盾がなければ　$1$を返します

$O(log N)$ 時間です

### `pair<int,int> calc(int x,int y,int v)`

$1.$ $A[x]$ , $A[y]$ が元々一意に定まっていれば $1$ $A[y]$

$2.$ $A[y]$ のみが元々一意に定まっていれば $2$ $A[y]$

$3.$ $A[x]$ の値によってx_qが一意に定まるならば $v$ を $A[x]$ に代入して $3$ $A[y]$

$4.$ $A[x]$ の値によって $A[y]$ が一意に定まらない場合は $4$ $0$

を返します

$O(logN)$ 時間です

### `pair<frac,frac> relationship(x)`

$x$ の root を $rx$ としたときに $A[rx]=a A[x]+b$

を満たす $a,b$ を返します 

$O(logN)$ 時間です


### `pair<bool,int> getval(x)`

$A[x]$ の値が一意に定まっているならば true $A[x]$

そうでなければ false $0$

を返します

$O(logN)$ 時間です

