#include <bits/stdc++.h>
using namespace std;
constexpr int MOD=998244353,ROOT=3;
int powmod(int a,int b){long long r=1,x=a;while(b){if(b&1)r=r*x%MOD;x=x*x%MOD;b>>=1;}return (int)r;}
void ntt(vector<int>& a,bool inv){int n=a.size();for(int i=1,j=0;i<n;++i){int b=n>>1;for(;j&b;b>>=1)j^=b;j^=b;if(i<j)swap(a[i],a[j]);}for(int len=2;len<=n;len<<=1){int wlen=powmod(ROOT,(MOD-1)/len);if(inv)wlen=powmod(wlen,MOD-2);for(int i=0;i<n;i+=len){long long w=1;for(int j=0;j<len/2;j++){int u=a[i+j],v=(long long)a[i+j+len/2]*w%MOD;int x=u+v;if(x>=MOD)x-=MOD;int y=u-v;if(y<0)y+=MOD;a[i+j]=x;a[i+j+len/2]=y;w=w*wlen%MOD;}}}if(inv){int w=powmod(n,MOD-2);for(int &x:a)x=(long long)x*w%MOD;}}
int main(){const int N=1000000, L=2*N+20;vector<unsigned char>p(L+1,1);p[0]=p[1]=0;for(int i=2;i*i<=L;i++)if(p[i])for(int j=i*i;j<=L;j+=i)p[j]=0;int len=1;while(len<2*(L+1))len<<=1;vector<int> a(len);for(int q=2;q<=L-6;q++)if(p[q]&&p[q+6])a[q]=1;ntt(a,0);for(int &v:a)v=(long long)v*v%MOD;ntt(a,1);vector<int>failure;int maxfailure=0;for(int n=12;n<=N;n++){int m=2*n-6;int count=a[m];int mid=m/2;if(2*mid==m && mid>=0 && mid<=L-6 && p[mid]&&p[mid+6])count--;int lo=(m-6)/2; if(2*lo==m-6 && lo>=2 && lo+12<=L && p[lo]&&p[lo+6]&&p[lo+12])count-=2;
 if(count<=0){failure.push_back(n);maxfailure=n;}}
cout<<"quad fail count "<<failure.size()<<" last "<<maxfailure<<"\n";for(int n:failure)cout<<n<<",";cout<<"\n";}
