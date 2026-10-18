#include <bits/stdc++.h>
using namespace std;
int main(int argc,char**argv){
 int maxN = argc>1 ? atoi(argv[1]):30000000;
 vector<uint8_t> isprime(maxN+1,1);isprime[0]=isprime[1]=0;
 for(long long i=2;i*i<=maxN;i++) if(isprime[i])for(long long j=i*i;j<=maxN;j+=i)isprime[j]=0;
 vector<int> primes;for(int i=2;i*i<=maxN;i++)if(isprime[i])primes.push_back(i);
 cout<<"N,u,half_gap,j,p,q,H,div7,div11,div13,div17,div19\n";
 for(int N=70;N<=maxN;N+=60){
   int m=N/2;
   int j=1;
   for(; m-6*j>=5; j++) if(isprime[m-6*j]&&isprime[m+6*j])break;
   int p,q,d;
   if(m-6*j<5){
     if(isprime[N-3]) {p=3;q=N-3;d=m-3;j=-1;}
     else {p=0;q=0;d=-1;j=-1;}
   } else {d=6*j;p=m-d;q=m+d;}
   int left=N;double H=1.0;
   for(int ell:primes){if((long long)ell*ell>left)break;if(left%ell==0){if(ell>2)H*= (double)(ell-1)/(ell-2);do{left/=ell;}while(left%ell==0);}}
   if(left>2)H*= (double)(left-1)/(left-2);
   cout<<N<<','<<(N+2)/12<<','<<d<<','<<j<<','<<p<<','<<q<<','<<setprecision(12)<<H<<','<<(N%7==0)<<','<<(N%11==0)<<','<<(N%13==0)<<','<<(N%17==0)<<','<<(N%19==0)<<'\n';
 }
}
