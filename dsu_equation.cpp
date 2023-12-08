template<class T>
struct Fraction{
    T Bunsi,Bunbo;
    Fraction(T Bunsi=0,T Bunbo=1)noexcept:Bunsi(Bunsi),Bunbo(Bunbo){reduce();}
    void reduce(){
        if(Bunbo<0){Bunbo=-Bunbo;Bunsi=-Bunsi;}
        if(Bunsi==0){Bunbo=1;return;}
        T GCD=gcd(abs(Bunsi),abs(Bunbo));
        Bunsi/=GCD;Bunbo/=GCD;
    }
    const Fraction<T>& operator++(){
        Bunsi+=Bunbo;reduce();return *this;
    }
    const Fraction<T>& operator--(){
        Bunsi-=Bunbo;reduce();return *this;
    }
    const Fraction<T>& operator+=(const Fraction<T>& rhs)noexcept{
        Bunsi=Bunsi*rhs.Bunbo+Bunbo*rhs.Bunsi;
        Bunbo=Bunbo*rhs.Bunbo;reduce();
        return *this;
    }
    const Fraction<T>& operator-=(const Fraction<T>& rhs)noexcept{
        Bunsi=Bunsi*rhs.Bunbo-Bunbo*rhs.Bunsi;
        Bunbo=Bunbo*rhs.Bunbo;reduce();
        return *this;
    }
    const Fraction<T>& operator*=(const Fraction<T>& rhs)noexcept{
        Bunsi=Bunsi*rhs.Bunsi;
        Bunbo=Bunbo*rhs.Bunbo;reduce();
        return *this;
    }
    const Fraction<T>& operator/=(const Fraction<T>& rhs)noexcept{
        Bunsi*=rhs.Bunbo;
        Bunbo*=rhs.Bunsi;reduce();
        return *this;
    }

    friend Fraction<T> operator+(const Fraction<T> &lhs,const Fraction<T> &rhs)noexcept{
        return Fraction(lhs)+=rhs;
    }
    friend Fraction<T> operator-(const Fraction<T> &lhs,const Fraction<T> &rhs)noexcept{
        return Fraction(lhs)-=rhs;
    }
    friend Fraction<T> operator*(const Fraction<T> &lhs,const Fraction<T> &rhs)noexcept{
        return Fraction(lhs)*=rhs;
    }
    friend Fraction<T> operator/(const Fraction<T> &lhs,const Fraction<T> &rhs)noexcept{
        return Fraction(lhs)/=rhs;
    }
    constexpr bool operator<(const Fraction<T> &rhs)const noexcept{
        return Bunsi*rhs.Bunbo<Bunbo*rhs.Bunsi;
    }
    constexpr bool operator<=(const Fraction<T> &rhs)const noexcept{
        return Bunsi*rhs.Bunbo<=Bunbo*rhs.Bunsi;
    }
    constexpr bool operator>(const Fraction<T> &rhs)const noexcept{
        return Bunsi*rhs.Bunbo>Bunbo*rhs.Bunsi;
    }
    constexpr bool operator>=(const Fraction<T> &rhs)const noexcept{
        return Bunsi*rhs.Bunbo>=Bunbo*rhs.Bunsi;
    }
    constexpr bool operator==(const Fraction<T> &rhs)const noexcept{
        return Bunsi*rhs.Bunbo==Bunbo*rhs.Bunsi;
    }
    constexpr bool operator!=(const Fraction<T> &rhs)const noexcept{
        return Bunsi*rhs.Bunbo!=Bunbo*rhs.Bunsi;
    }
    
};



struct dsu_equation{
private:
    using frac=Fraction<long long>;
    vector<int> par; //parent
    vector<int> siz; //size
    vector<frac> coef_p,coef_q; //coefficient
    vector<bool> specify_; //is value uniquely determined
    vector<frac> val; //value
public:
    dsu_equation(int n=0){
        par.resize(n);
        siz.resize(n);
        coef_p.resize(n);
        coef_q.resize(n);
        specify_.resize(n);
        val.resize(n);
        for(int i=0;i<n;i++){
            par[i]=i;
            siz[i]=1;
        }
    }

    //Find a leader of the group contains x
    int leader(int x){
        if(par[x]==x)return x;
        return leader(par[x]);
    }

    //Find the size of the group contains x
    int size(int x){
        return siz[leader(x)];
    }

    //Check if x y are in the same group
    bool same(int x,int y){
        return leader(x)==leader(y);
    }
    
    //return rootx(leader of x) , p , q satisfying A[rootx] = p * A[x] + q 
    tuple<int,frac,frac> relationship(int x,frac p=1,frac q=0){
        if(par[x]==x)return {x,p,q};
        return relationship(par[x],p*coef_p[x],q*coef_p[x]+coef_q[x]);
    }

    //unite as a * A[x] + b * A[y] = c
    //return false if contradiction else true
    //must satisfy a != 0 and b != 0
    bool merge(int x,int y,int a,int b,int c){
        assert(a!=0&&b!=0);
        auto [rx,px,qx]=relationship(x);
        auto [ry,py,qy]=relationship(y);
        if(specify_[rx]&&specify_[ry]){
            return (a*(val[rx]/px-qx/px)+b*(val[ry]/py-qy/py)==c);
        }
        if(rx==ry){
            //value is determined uniquely or other
            frac ca=a/px+b/py,cb=c+a*qx/px+b*qy/py;
            if(ca==0)return cb==0;
            specify_[rx]=true;
            val[rx]=cb/ca;
            return true;
        }else{
            if(siz[rx]<siz[ry]){
                swap(rx,ry),swap(a,b),swap(px,py),swap(qx,qy);
            }
            par[ry]=rx;
            siz[rx]+=siz[ry];
            coef_p[ry]=(-b*px)/(a*py);
            coef_q[ry]=(qx+b*px*qy/(a*py)+(c*px)/a);
            if(specify_[ry]){
                val[rx]=coef_p[ry]*val[ry]+coef_q[ry];
            }
            return true;
        }
    }

    //set as A[x]=f , return false if contradiction else true
    bool set_val(int x,frac f){
        auto [rx,px,qx]=relationship(x);
        if(specify_[rx])return px*f+qx==val[rx];
        else {val[rx]=px*f+qx;return specify_[rx]=true;}
    }
    
    //return {is A[x] uniquely determined ?, its value}
    pair<bool,frac> specify(int x){
        auto [rx,px,qx]=relationship(x);
        if(!specify_[rx])return {false,0};
        return {true,val[rx]/px-qx/px};
    }

    //return {p,q} satisfying A[y] = p * A[x] + q 
    tuple<bool,frac,frac> solve_coef(int x,int y){
        auto [rx,px,qx]=relationship(x);
        auto [ry,py,qy]=relationship(y);
        if(rx!=ry)return {false,0,0};
        return {true,px/py,(qx-qy)/py};
    }
};
