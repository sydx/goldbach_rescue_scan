import numpy as np, json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":8.5,"axes.edgecolor":"#52514e","axes.labelcolor":"#0b0b0b",
    "xtick.color":"#52514e","ytick.color":"#52514e","axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":0.6,"mathtext.fontset":"cm"})
C=["#2a78d6","#eb6834","#1baf7a","#eda100","#e87ba4","#008300"]
N=8000
s=np.ones(N+40,bool);s[:2]=False
for i in range(2,int((N+40)**.5)+1):
    if s[i]: s[i*i::i]=False
def Pd(d):
    up=s[d:d+N+1];dn=np.zeros(N+1,bool);dn[d:]=s[:N+1-d];return s[:N+1]&(up|dn)
P2,P4=Pd(2),Pd(4)
r22=np.convolve(P2.astype(int),P2.astype(int))[:N+1]; ev=np.arange(4,N+1,2); r22u=(r22[ev]+P2[ev//2])/2
r24=np.convolve(P2.astype(int),P4.astype(int))[:N+1][ev]
fig,ax=plt.subplots(1,2,figsize=(6.6,2.85),gridspec_kw=dict(width_ratios=[1.25,1],wspace=0.32))
a=ax[0]
a.scatter(ev,r24,s=2.2,c=C[1],lw=0,label=r"$r_{2,4}(n)$ (twin + cousin)",rasterized=True)
a.scatter(ev,r22u,s=2.2,c=C[0],lw=0,label=r"$r_{2,2}(n)$ (twin + twin)",rasterized=True)
z=ev[r22u==0]
a.scatter(z,np.zeros(len(z)),s=9,marker="x",c="#0b0b0b",lw=0.7,label=r"$r_{2,2}(n)=0$",zorder=5)
a.set_xlim(0,N);a.set_ylim(-1.5,None);a.set_xlabel(r"$n$");a.set_ylabel("representations")
a.grid(axis="y",color="#e4e3df",lw=0.5);a.set_axisbelow(True)
a.legend(frameon=False,loc="upper left",markerscale=2.5,handletextpad=0.3,fontsize=7.5)
a.set_title("(a) Counts for small $n$",fontsize=8.5,loc="left")
b=ax[1]
d=json.load(open("../data/hl2.json"))
order=[("2-2","twin+twin","o"),("4-4","cousin+cousin","s"),("2-4","twin+cousin","^"),
       ("2-6","twin+sexy","D"),("4-6","cousin+sexy","v"),("6-6","sexy+sexy","P")]
xs=[1e3,1e4,1e5,1e6,1e7]; keys=["1000","10000","100000","1000000",str(10**7-20000)]
for i,(k,lab,m) in enumerate(order):
    b.plot(xs,[d[k][q]["ratio"] for q in keys],marker=m,ms=4,lw=1.4,c=C[i],label=lab,mec="#fcfcfb",mew=0.6)
b.plot(xs,[d["6-6"][q]["ratio_lead"] for q in keys],ls=(0,(3,2)),lw=1.1,c=C[5],marker="P",ms=4,mfc="none",label="sexy+sexy, leading terms only")
b.axhline(1,c="#52514e",lw=0.6)
b.set_xscale("log");b.set_xlabel(r"$n$ (window start)");b.set_ylabel("actual / predicted")
b.set_ylim(0.6,1.05);b.grid(axis="y",color="#e4e3df",lw=0.5);b.set_axisbelow(True)
b.legend(frameon=False,fontsize=6.6,loc="upper center",bbox_to_anchor=(0.45,-0.19),ncol=2,handlelength=2.2,labelspacing=0.25,columnspacing=1.0)
b.set_title("(b) Fit of the conjectured asymptotic",fontsize=8.5,loc="left")
fig.savefig("../fig_counts.pdf",bbox_inches="tight",dpi=400); pass
