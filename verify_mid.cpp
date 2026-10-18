// Exact direct integer-loop verification of the restricted convolution
// for 2 <= t <= 100000; sample counts also printed at t=1000000.
#include <bits/stdc++.h>
using namespace std;
int main(){
 const int T=1000000, L=6*(T+1)+5;
 vector<uint8_t> p(L+1,1),f(T+2),g(T+2);
 p[0]=p[1]=0;
 for(int i=2;i*i<=L;i++)if(p[i])for(int j=i*i;j<=L;j+=i)p[j]=0;
 for(int i=1;i<=T+1;i++)f[i]=p[6*i-1];
 for(int i=1;i<=T;i++)g[i]=f[i]&&f[i+1];
 int extinctionFailures=0,allFailures=0;
 for(int t=2;t<=100000;t++){
   int r=0;
   for(int a=(t+3)/4;a<=t/2;a++)r+=g[a]*f[t-a];
   if(r==0){allFailures++;if(!f[t])extinctionFailures++;}
 }
 cout<<"Restricted t<=100000: all failures="<<allFailures
     <<", extinction failures="<<extinctionFailures<<"\n";
 for(int t: {92,1000,10000,100000,1000000}){
   int r=0,k=0;
   for(int a=(t+3)/4;a<=t/2;a++){ k+=g[a];r+=g[a]*f[t-a]; }
   cout<<"t="<<t<<" K="<<k<<" Rmid="<<r<<"\n";
 }
}
