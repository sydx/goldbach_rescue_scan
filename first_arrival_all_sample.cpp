#include <bits/stdc++.h>
using namespace std;
int main(int argc,char**argv){
 int maxN=argc>1?atoi(argv[1]):30000000;
 vector<uint8_t> isp(maxN+1,1);isp[0]=isp[1]=0;
 for(long long i=2;i*i<=maxN;i++)if(isp[i])for(long long j=i*i;j<=maxN;j+=i)isp[j]=0;
 vector<int> primes; for(int i=2;i*i<=maxN;i++)if(isp[i])primes.push_back(i);
 cout<<"N,t,d,p,q,H,div5,div7,div11,div13,midpoint_even\n";
 long long observed=0,missed=0;
 for(int N=16;N<=maxN;N+=6){
  int t=(N-4)/6;
  if(isp[N-5])continue;
  observed++;
  if(t%101>=20)continue;
  int m=N/2;int d=m%2==0?3:0;
  int p=0,q=0;
  for(;m-d>=5;d+=6)if(isp[m-d] && isp[m+d]){p=m-d;q=m+d;break;}
  if(p==0){if(isp[N-3]){p=3;q=N-3;d=m-3;}else {d=-1;missed++;}}
  int left=N;double H=1.;
  for(int ell:primes){if((long long)ell*ell>left)break;if(left%ell==0){if(ell>2)H*=(double)(ell-1)/(ell-2);do{left/=ell;}while(left%ell==0);}}
  if(left>2)H*=(double)(left-1)/(left-2);
  cout<<N<<','<<t<<','<<d<<','<<p<<','<<q<<','<<setprecision(12)<<H<<','<<(N%5==0)<<','<<(N%7==0)<<','<<(N%11==0)<<','<<(N%13==0)<<','<<(m%2==0)<<'\n';
 }
 cerr<<"reset-target rows examined="<<observed<<", missed="<<missed<<"\n";
}
