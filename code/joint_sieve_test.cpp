#include <bits/stdc++.h>
using namespace std;
int main(int argc, char**argv){
 int X=30000000,J=120,L=26, minN=1000000;
 vector<uint8_t> prime(X+1,1); prime[0]=prime[1]=0;
 for(int p=2;p*p<=X;p++)if(prime[p])for(int q=p*p;q<=X;q+=p)prime[q]=0;
 long long K=0;
 vector<long long> hits(J+1,0),both(L+1,0), opp(L+1,0), binned(8,0), eligible(8,0);
 vector<long long> paircount(8,0), moments(8,0), any(8,0), groupcount(8,0);
 vector<int> widths={5,10,15,20,30,40,60,100};
 vector<long long> residue_hits(5,0),residue_trials(5,0);
 for(int N=70;N<=X;N+=60) {
  if(N<minN) continue;
  ++K;
  int m=N/2;int a[121]={};
  for(int j=1;j<=J;j++){
   bool e=(j%5!=0);
   if(e) {
      a[j]=prime[m-6*j]&&prime[m+6*j];
      hits[j]+=a[j]; residue_hits[j%5]+=a[j]; residue_trials[j%5]++;
   }
  }
  for(int d=1;d<=L;d++) for(int j=1;j+d<=J;j++) if(j%5 && (j+d)%5) {
    opp[d]++;
    both[d]+=a[j]*a[j+d];
  }
  for(int h=0;h<(int)widths.size();h++) {
     int mcount=0;
     for(int j=1;j<=widths[h];j++)mcount+=a[j];
     groupcount[h]++; paircount[h]+=mcount;
     moments[h]+=(long long)mcount*mcount; any[h]+=mcount>0;
  }
 }
 cout<<setprecision(10)<<"K="<<K<<"\n";
 cout<<"width,mean,variance,p_zero,poisson_zero,fano\n";
 for(int h=0;h<(int)widths.size();h++) {
   double av=(double)paircount[h]/groupcount[h], var=(double)moments[h]/groupcount[h]-av*av;
   cout<<widths[h]<<","<<av<<","<<var<<","<<1.-(double)any[h]/groupcount[h]<<","<<exp(-av)<<","<<var/av<<"\n";
 }
 cout<<"index residue hit rates: \n";
 for(int j=0;j<5;j++) cout<<j<<" "<<residue_hits[j]<<" "<<residue_trials[j]<<" "<<(residue_trials[j]?(double)residue_hits[j]/residue_trials[j]:0)<<"\n";
 cout<<"lag,eligible_trials,both_hits,rate_joint,p_independent,ratio\n";
 for(int d=1;d<=L;d++){
   double prediction=0;int nc=0;
   for(int j=1;j+d<=J;j++) if(j%5 && (j+d)%5) {prediction+=(double)hits[j]/K*((double)hits[j+d]/K);nc++;}
   prediction/=nc;
   double obs=(double)both[d]/opp[d];
   cout<<d<<","<<opp[d]<<","<<both[d]<<","<<obs<<","<<prediction<<","<<obs/prediction<<"\n";
 }
}
