#include <bits/stdc++.h>
using namespace std;
int main(){const int X=30000000, K=600;vector<uint8_t> prime(X+1,1);prime[0]=prime[1]=0;for(int p=2;p*p<=X;p++)if(prime[p])for(int q=p*p;q<=X;q+=p)prime[q]=0;
vector<int> widths={10,20,40,60,100,150,250,400,600};long long cnt=0;vector<long double> means(widths.size()),sq(widths.size());
for(int N=1000030;N<=X;N+=60){cnt++;int m=N/2,c=0,id=0;for(int j=1;j<=K;j++){if(j%5 && prime[m-6*j] && prime[m+6*j])c++; if(j==widths[id]){means[id]+=c;sq[id]+=(long double)c*c;id++;if(id==(int)widths.size()) break;}}}
for(int z=0;z<(int)widths.size();z++){double mu=(double)(means[z]/cnt),var=(double)(sq[z]/cnt-means[z]*means[z]/cnt/cnt);cout<<widths[z]<<","<<mu<<","<<var<<","<<var/mu<<"\n";}}
