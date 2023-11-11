template<class T>
class frac{
    T bunsi,bunbo;
    constexpr void reduct()noexcept{
        T g=gcd(bunsi,bunbo);
        bunsi/=g;bunbo/=g;
        if(bunbo<0){
            bunsi=-bunsi;bunbo=-bunbo;
        }
    }
public:
    constexpr frac(T Bunsi=0,T Bunbo=1)noexcept{
        bunsi=Bunsi;bunbo=Bunbo;
        reduct();
    }
    constexpr T &Bunsi(){return bunsi;}
    constexpr T &Bunbo(){return bunbo;}
    constexpr frac<T>&operator+=(const frac<T>&rhs)noexcept{
        bunsi=bunsi*rhs.bunbo+bunbo*rhs.bunsi;
        bunbo*=rhs.bunbo;
        reduct();
        return *this;
    }
    constexpr frac<T>&operator-=(const frac<T>&rhs)noexcept{
        bunsi=bunsi*rhs.bunbo-bunbo*rhs.bunsi;
        bunbo*=rhs.bunbo;
        reduct();
        return *this;
    }
    constexpr frac<T>&operator*=(const frac<T>&rhs)noexcept{
        bunbo*=rhs.bunbo;
        bunsi*=rhs.bunsi;
        reduct();
        return *this;
    }
    constexpr frac<T>&operator/=(const frac<T>&rhs)noexcept{
        bunbo*=rhs.bunsi;
        bunsi*=rhs.bunbo;
        reduct();
        return *this;
    }
    constexpr frac<T>operator+(const frac<T>&rhs)const noexcept{return frac(*this)+=rhs;}
    constexpr frac<T>operator-(const frac<T>&rhs)const noexcept{return frac(*this)-=rhs;}
    constexpr frac<T>operator*(const frac<T>&rhs)const noexcept{return frac(*this)*=rhs;}
    constexpr frac<T>operator/(const frac<T>&rhs)const noexcept{return frac(*this)/=rhs;}
    constexpr bool operator<(const frac<T>&rhs)const noexcept{return bunsi*rhs.bunbo<bunbo*rhs.bunsi;}
    constexpr bool operator>(const frac<T>&rhs)const noexcept{return bunsi*rhs.bunbo>bunbo*rhs.bunsi;}
    constexpr bool operator>=(const frac<T>&rhs)const noexcept{return bunsi*rhs.bunbo>=bunbo*rhs.bunsi;}
    constexpr bool operator<=(const frac<T>&rhs)const noexcept{return bunsi*rhs.bunbo<=bunbo*rhs.bunsi;}
    constexpr bool operator==(const frac<T>&rhs)const noexcept{return bunsi*rhs.bunbo==bunbo*rhs.bunsi;}
    constexpr bool operator!=(const frac<T>&rhs)const noexcept{return bunsi*rhs.bunbo!=bunbo*rhs.bunsi;}
};

#define frac frac<int>

struct dsu_equation{
    private:
    vector<int> par;
    vector<int> siz;
    vector<frac> coef_x,coef_c;
    vector<bool> isvaild;
    vector<frac> val;
    public:
    dsu_equation(int n=1){
        par.resize(n);
        siz.resize(n);
        coef_x.resize(n);
        coef_c.resize(n);
        isvaild.resize(n);
        val.resize(n);
        for(int i=0;i<n;i++)par[i]=i;
        for(int i=0;i<n;i++){
            siz[i]=1;
            coef_x[i]=1;
            coef_c[i]=0;
        }
    }
    int root(int x){
        if(par[x]==x)return x;
        return root(par[x]);
    }
    
    int size(int x){
        return siz[root(x)];
    }

    // A_rx=pA_x+q となるp,qを返す
    pair<frac,frac> relationship(int x,frac cx=1,frac cc=0){
        if(par[x]==x)return make_pair(cx,cc);
        return relationship(par[x],cx*coef_x[x],cc*coef_x[x]+coef_c[x]);
    }

    bool same(int x,int y){
        return root(x)==root(y);
    }

    //ax+by=cとして結合
    int merge(int x,int y,int a,int b,int c){
        if(a==0||b==0){
            if(a==0&&b==0){
                if(c==0)return 1;
                else return 0;
            }else if(a!=0){
                frac res=c/a;
                int rx=root(x);
                pair<frac,frac> lx=relationship(x);
                if(isvaild[rx]&&lx.first*res+lx.second!=val[rx])return 0;
                else{
                    isvaild[rx]=true;
                    val[rx]=lx.first*res+lx.second;
                    return 1;
                }
            }else if(b!=0){
                frac res=c/b;
                int ry=root(y);
                pair<frac,frac> ly=relationship(y);
                if(isvaild[ry]&&ly.first*res+ly.second!=val[ry])return 0;
                else{
                    isvaild[ry]=true;
                    val[ry]=ly.first*res+ly.second;
                    return 1;
                }
            }
        }
        //xを根に持ってくる
        if(same(x,y)){
            pair<frac,frac> lx=relationship(x);
            pair<frac,frac> ly=relationship(y);
            int rx=root(x);
            if(isvaild[rx]){
                frac cx=frac(a)/lx.first+frac(b)/ly.first;
                frac cc=frac(c)+frac(a)*lx.second/lx.first+frac(b)*ly.second/ly.first;
                if(cx==0){
                    if(cc==0)return 1;
                    else return 0;
                }
                else if(val[rx]==cc/cx)return 1;
                else return 0;
            }else{
                frac cx=frac(a)/lx.first+frac(b)/ly.first;
                frac cc=frac(c)+frac(a)*lx.second/lx.first+frac(b)*ly.second/ly.first;
                if(cx==0){
                    if(cc==0)return 1;
                    else return 0;
                }
                else{
                    isvaild[rx]=true;
                    val[rx]=cc/cx;
                    return 1;
                }
            }
        }
        if(size(x)<size(y)){
            swap(x,y);swap(a,b);
        }
        pair<frac,frac> lx=relationship(x);
        pair<frac,frac> ly=relationship(y);
        int rx=root(x);
        int ry=root(y);
        if(isvaild[ry]){
            frac res=(frac(c)-frac(b)*val[ry]/ly.first+frac(b)*ly.second/ly.first+frac(a)*lx.second/lx.first)*lx.first/a;
            if(isvaild[rx]&&res!=val[rx])return 0;
            isvaild[rx]=true;
            val[rx]=res;
        }
        par[ry]=rx;
        siz[rx]+=siz[ry];
        coef_x[ry]=(frac(-b)*lx.first)/(frac(a)*ly.first);
        coef_c[ry]=(frac(c)+frac(b)*ly.second/ly.first+frac(a)*lx.second/lx.first)*lx.first/a;
        return 1;
    }

    pair<bool,int> getval(int x){
        pair<frac,frac> lx=relationship(x);
        int rx=root(x);
        if(!isvaild[rx])return make_pair(false,0);
        else{
            frac res=(val[rx]-lx.second)/lx.first;
            return make_pair(true,res.Bunsi());
        }
    }
    
    pair<int,int> calc(int x,int y,int xval){
        auto xv=getval(x);
        auto yv=getval(y);
        if(xv.first&&yv.first)return make_pair(1,yv.second);
        if(yv.first)return make_pair(2,yv.second);
        if(same(x,y)){
            pair<frac,frac> lx=relationship(x),ly=relationship(y);
            frac res=(lx.first*xval+lx.second-ly.second)/ly.first;
            return make_pair(3,res.Bunsi());
        }else{
            return make_pair(4,0);
        }
    }

};
