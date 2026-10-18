import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
P=Path(__file__).resolve().parent
D=pd.read_csv(P/'reset_first_arrivals_30m.csv')
C2=0.6601618158468696
D['dlog']=D.half_gap/np.log(D.N/2)**2
D['z']=2*C2*D.H*D.dlog
D['imbalance']=2*D.half_gap/D.N
# Tail distribution.
D2=D[D.N>=1000000]
grid=np.linspace(0,5,101)
survival=np.array([(D2.z>x).mean() for x in grid])
fig,ax=plt.subplots(figsize=(8.8,5.3));ax.plot(grid,survival,label='Observed (N ≥ 1 million)',linewidth=2)
ax.plot(grid,np.exp(-grid),'--',label='Exponential benchmark',linewidth=1.8)
ax.set(xlabel=r'Singular-series-scaled distance  $z=2C_2 H(N)d/\log^2(N/2)$',ylabel=r'$\Pr(Z>z)$',title='First-arrival distance: observed survival versus exponential heuristic')
ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(P/'survival_vs_exponential.png',dpi=180);plt.close(fig)
# Geometric bins for normalized distance mean and median
edges=np.geomspace(1e3,3e7,26)
mid=[];avg=[];med=[]
for a,b in zip(edges[:-1],edges[1:]):
 s=D[(D.N>=a)&(D.N<b)]
 if len(s)>10:mid.append(np.sqrt(a*b));avg.append(s.dlog.mean());med.append(s.dlog.median())
fig,ax=plt.subplots(figsize=(8.8,5.3));ax.plot(mid,avg,marker='o',markersize=3,label='Mean')
ax.plot(mid,med,marker='.',markersize=3,label='Median')
ax.set(xscale='log',xlabel='Even target N (guaranteed reset, N ≡ 10 mod 60)',ylabel=r'Nearest distance / $\log^2(N/2)$',title='First-arrival scale across midpoint sizes');ax.grid(alpha=.2);ax.legend();fig.tight_layout();fig.savefig(P/'scaling_with_N.png',dpi=180);plt.close(fig)
# Sample of the dynamics, consecutive reset targets around 10,000,000
s=D[(D.N>=10_000_030)&(D.N<=10_060_000)]
fig,ax=plt.subplots(figsize=(9,4.8));ax.plot(s.N,s.half_gap,linewidth=.8)
ax.set(xlabel='Successive reset targets N (spacing 60)',ylabel='Nearest Goldbach half-gap d(N)',title='Irregular first-arrival dynamics over 1,000 successive resets');ax.grid(alpha=.2);fig.tight_layout();fig.savefig(P/'successive_reset_dynamics.png',dpi=180);plt.close(fig)
