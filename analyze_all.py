import pandas as pd,numpy as np
from scipy.stats import ks_1samp
from scipy.stats import spearmanr
D='/mnt/data/imbalance_first_arrival/reset_first_arrivals_all_sample_30m.csv.gz'
df=pd.read_csv(D)
C2=0.6601618158468696
m=df.N/2
df['rate_z']=df.d*2*C2*df.H/np.log(m)**2
df['raw_z']=df.d/np.log(m)**2
print('total',len(df),'missing',sum(df.d<0),'diagonal',sum(df.d==0))
for l,h in [(1e4,1e5),(1e5,1e6),(1e6,1e7),(1e7,3e7+1)]:
 s=df[(df.N>=l)&(df.N<h)]
 print('bin',int(l),int(h),'rows',len(s),'mean d',round(s.d.mean(),1),'median d',s.d.median(),'mean d/log²',round(s.raw_z.mean(),3),'median d/log²',round(s.raw_z.median(),3),'mean norm',round(s.rate_z.mean(),3),'med norm',round(s.rate_z.median(),3),'mean H',round(s.H.mean(),3))
s=df[df.N>=1e6]
for factor in ('div5','div7','div11','div13'):
 print('\n',factor)
 for b in (0,1):
  z=s[s[factor]==b]
  print(b, len(z), 'meanH',round(z.H.mean(),3),'mean_d',round(z.d.mean(),2),'raw',round(z.raw_z.mean(),3),'median raw',round(z.raw_z.median(),3),'z',round(z.rate_z.mean(),3),'med z',round(z.rate_z.median(),3),'d=0%',round(100*np.mean(z.d==0),2))
print('spearman H versus raw',spearmanr(s.H,s.raw_z))
print('all q of normalized',s.rate_z.quantile([.1,.25,.5,.75,.9,.99]).to_dict())
for z in (.25,.5,1,2,3):
 print('tail',z,'obs',round(np.mean(s.rate_z>z),4),'exp',round(np.exp(-z),4))
print('largest gap',df.loc[df.d.idxmax()].to_dict())
