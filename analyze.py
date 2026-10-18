import pandas as pd,numpy as np
from scipy.stats import pearsonr, spearmanr
f='/mnt/data/imbalance_first_arrival/reset_first_arrivals_30m.csv'
df=pd.read_csv(f)
C2=0.6601618158468696
df['logm2']=np.log(df.N/2)**2
df['rawnorm']=df.half_gap/df.logm2
df['adj']=df.half_gap*2*C2*df.H/df.logm2
df['normN']=df.half_gap/np.log(df.N)**2
print('count',len(df),'missing', (df.half_gap<0).sum())
for low,hi in [(1e3,1e4),(1e4,1e5),(1e5,1e6),(1e6,1e7),(1e7,3e7+1),(1e6,3e7+1)]:
 s=df[(df.N>=low)&(df.N<hi)];
 print('bin',int(low),int(hi),'N',len(s),'mean d',round(s.half_gap.mean(),2),'median d',s.half_gap.median(),'q90 d',s.half_gap.quantile(.9),'mean d/log^2(m)',round(s.rawnorm.mean(),3),'median d/log^2m',round(s.rawnorm.median(),3),'mean adj',round(s.adj.mean(),3),'median adj',round(s.adj.median(),3),'p90 adj',round(s.adj.quantile(.9),3),'mean H',round(s.H.mean(),3),'j1',(s.j==1).mean())
print('max',df.loc[df.half_gap.idxmax()].to_dict());print('max norm',df.loc[df.adj.idxmax()].to_dict())
s=df[df.N>=1e6]
for ell in [7,11,13,17,19]:
 print('prime',ell)
 for val in (0,1):
  z=s[s[f'div{ell}']==val]
  print(val,'N',len(z),'mean H',round(z.H.mean(),4),'mean j',round(z.j.mean(),3),'median j',z.j.median(),'mean d/logm2',round(z.rawnorm.mean(),4),'median',round(z.rawnorm.median(),4),'mean adj',round(z.adj.mean(),4),'median adj',round(z.adj.median(),4))
print('correlations n>=1e6 H with d/logm² Pearson Spearman',pearsonr(s.H,s.rawnorm),spearmanr(s.H,s.rawnorm))
print('adj distribution percentiles',s.adj.quantile([.05,.1,.25,.5,.75,.9,.95,.99]).to_dict())
for zz in [.25,.5,1,1.5,2,3]:
 print('survival',zz,'observed',np.mean(s.adj>zz),'exp',np.exp(-zz))
print('earliest j distribution N>=1e6',s.j.value_counts().sort_index().head(25).to_dict())
print('proportion j%5==0',np.mean(s.j%5==0))
