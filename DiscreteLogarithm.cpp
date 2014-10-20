#include <gmp.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>

using namespace std;

class ZZ
{
public:
  ZZ(ZZ const & a) { mpz_init_set(data, a.data); }
  ZZ(long a = 0) { mpz_init_set_si(data, a); }

  ZZ(string const & str, int base = 10) { mpz_init_set_str(data, str.c_str(), base); }

  ~ZZ() noexcept { mpz_clear(data); }

  ZZ & operator=(ZZ const & a) { mpz_set(data, a.data); return *this; } // self-assignment ok

  ZZ & operator+=(ZZ const & a) { mpz_add(data, data, a.data); return *this; }
  ZZ & operator-=(ZZ const & a) { mpz_sub(data, data, a.data); return *this; }
  ZZ & operator*=(ZZ const & a) { mpz_mul(data, data, a.data); return *this; }

  ZZ operator+(ZZ const & a) const { return ZZ(*this) += a; }
  ZZ operator-(ZZ const & a) const { return ZZ(*this) -= a; }
  ZZ operator*(ZZ const & a) const { return ZZ(*this) *= a; }

  ZZ operator-() const { ZZ res; mpz_neg(res.data, data); return res; }

  ZZ & operator/=(ZZ const & a) { mpz_fdiv_q(data, data, a.data);     return *this; }
  ZZ & operator%=(ZZ const & a) { mpz_fdiv_r(data, data, a.data);     return *this; }

  ZZ operator/(ZZ const & a) const { return ZZ(*this) /= a; }
  ZZ operator%(ZZ const & a) const { return ZZ(*this) %= a; }

  ZZ operator^(unsigned long int e) const { ZZ res; mpz_pow_ui(res.data, data, e); return res; }

  bool operator==(ZZ const & a) const { return mpz_cmp(data, a.data) == 0; }
  bool operator!=(ZZ const & a) const { return mpz_cmp(data, a.data) != 0; }

  bool operator>(ZZ const& a) const { return mpz_cmp(data, a.data) > 0; }
  bool operator<(ZZ const& a) const { return mpz_cmp(data, a.data) < 0; }
  bool operator>=(ZZ const& a) const { return mpz_cmp(data, a.data) >= 0; }
  bool operator<=(ZZ const& a) const { return mpz_cmp(data, a.data) <= 0; }

  // only prefix
  ZZ & operator++() { mpz_add_ui(data, data, 1); return *this; }
  ZZ & operator--() { mpz_sub_ui(data, data, 1); return *this; }

  // output
  string ToString(int base = 10) const
  {
    char * cstr = new char[mpz_sizeinbase(data, base) + 10];
    mpz_get_str(cstr, base, data);
    string str(cstr);
    delete [] cstr;
    return str;
  }

  friend ostream & operator<<(ostream & out, ZZ const & a) { out << a.ToString(); return out; }

  // friends
  friend ZZ Gcd(ZZ const & a, ZZ const & b) { ZZ res; mpz_gcd(res.data, a.data, b.data); return res; }
  friend ZZ PowMod(ZZ const & a, ZZ const & e, ZZ const & n) { ZZ res; mpz_powm(res.data, a.data, e.data, n.data); return res; }
  friend ZZ InvMod(ZZ const & a, ZZ const & n) { ZZ res; if (!mpz_invert(res.data, a.data, n.data)) return 0; return res; }
  friend bool IsPrime(ZZ const& n, int reps = 10) { return (mpz_probab_prime_p(n.data, reps) > 0); }
  friend ZZ Sqrt(ZZ const & a) { ZZ res; mpz_sqrt(res.data, a.data); return res; }

  // random
  static void Seed(unsigned long int seed) { gmp_randseed_ui(state.s, seed); }
  friend ZZ Random(ZZ const & n) { ZZ res; mpz_urandomm(res.data, ZZ::state.s, n.data); return res; }

  // other
  bool operator[](size_t index) const { return mpz_tstbit(data, (mp_bitcnt_t)index); }
  size_t BitLength() const { return mpz_sizeinbase(data, 2); }

private:
  mpz_t data;

  // random generation
    class RandState
    {
    public:
      RandState() { gmp_randinit_default(s); gmp_randseed_ui(s, 0xf0123456789abcdeul); }
      ~RandState() { gmp_randclear(s); }

      gmp_randstate_t s;
    };

  static RandState state;
};

ZZ::RandState ZZ::state;

ZZ pollard(ZZ n)
{
  if (n<2) return n;
  if (n%2==0) return 2;
  ZZ x=Random(n),y(x),d=1;
  do{
    x=(x*x+1)%n;
    y=(y*y+1)%n;
    y=(y*y+1)%n;
    d=Gcd(x-y,n);
  }while(d==1);
  return d;
}

ZZ babyGiant(ZZ b,ZZ a,ZZ n,ZZ p)
{
  ZZ m=Sqrt(n)+1,j(0),am(1),i(0),g=b;
  map<ZZ,ZZ> hsTb;
  while(j<m)
    {
      hsTb[am]=j;
      am*=a;
      am%=p;
      j+=1;
    }
  am=InvMod(PowMod(a,m,p),p);
  while(i<=m)
    {
      if (hsTb.find(g)!=hsTb.end())
      return i*m+hsTb[g];
      g*=am;
      g%=p;
      i+=1;
    }
  while(1);
}

ZZ discLog(ZZ b,ZZ a,ZZ ord,ZZ p)
{
  ZZ m=pollard(ord);
  if (ord==m)
    return babyGiant(b,a,ord,p);
  ZZ n=ord/m,x2=discLog(PowMod(b,n,p),PowMod(a,n,p),m,p);
  ZZ x1=discLog((b*InvMod(PowMod(a,x2,p),p))%p,PowMod(a,m,p),n,p);
  return m*x1+x2;
}

int main()
{
  ios_base::sync_with_stdio(0);
  string s;
  cin>>s;
  ZZ p(s);
  cin>>s;
  ZZ g(s);
  int x;
  cin>>x;
  while(x--)
    {
      cin>>s;
      ZZ b(s);
      cout<<discLog(b,g,p-1,p)<<endl;
    }
  return 0;
}
