#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
using namespace std;
int main(){
 const int M=1000000;
 const int MAXN=2*M+7;
 vector<unsigned char> prime(MAXN+1,1); prime[0]=prime[1]=0;
 for(long long p=2;p*p<=MAXN;p++) if(prime[p])for(long long q=p*p;q<=MAXN;q+=p)prime[q]=0;
 vector<int> ps;for(int x=3;x<=MAXN;x+=2)if(prime[x])ps.push_back(x);
 vector<int> earliest(MAXN+1,0);earliest[0]=3;
 for(int d=2;d<=2*M-6;d+=2){
   for(int p:ps){if(p+d>MAXN)break;if(prime[p+d]){earliest[d]=p;break;}}
 }
 long long resets=0,no_ancestor=0,without_goldbach=0,all_without_ancestor=0; long long distinct_m=0;int max_min_gap=0;
 vector<int> exc;
 for(int m=4;m<=M;m++){
   int n=m-1;bool reset=(n>=7 && n%3==1 && !prime[2*n-3]);
   if(reset)resets++;
   bool any=false,ancestral=false;
   for(int p:ps){if(p>m)break;int q=2*m-p;if(!prime[q])continue;
       any=true; if(earliest[q-p]<p){ancestral=true;break;}
   }
   if(!any) {if(reset)without_goldbach++; continue;}
   if(!ancestral){all_without_ancestor++; if(reset){no_ancestor++;exc.push_back(m);} if(m<100)cerr<<"ancestry-free m="<<m<<"\n";}
 }
 cout<<"m_max="<<M<<" reset_rows="<<resets<<" reset_goldbach_empty="<<without_goldbach<<" reset_nonempty_ancestryfree="<<no_ancestor<<" all_nonempty_ancestryfree="<<all_without_ancestor<<"\n";
 cout<<"reset_exceptions:";for(int m:exc)cout<<" "<<m;cout<<"\n";
}
