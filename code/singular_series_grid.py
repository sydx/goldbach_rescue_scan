import numpy as np, json
G=[2,4,6,8,10,12,14,16,18,20,22,24,26,28,30]
d=json.load(open('../data/cert_grid.json'))
s=np.ones(1001,bool);s[:2]=False
for i in range(2,32):
    if s[i]: s[i*i::i]=False
P=[l for l in np.nonzero(s)[0] if l<=200]
n=np.arange(2_000_000,2_000_000+2*30030*2,2)   # covers all residues mod 2*3*5*7*11*13
rows=[]
for i,a in enumerate(G):
    for b in G[i:]:
        S=np.zeros(len(n))
        for ea in (a,-a):
            for db in (b,-b):
                Sx=np.ones(len(n))
                for l in P:
                    tab=np.array([((l-len({0,(-ea)%l,r,(r+db)%l}))/l)/(1-1/l)**4 for r in range(l)])
                    Sx*=tab[n%l]
                S+=Sx
        if a==b: S/=2
        big=d[f'{a}-{b}']['threshold']-2
        rows.append((a,b,float(S.min()),float(S.mean()),big))
rows.sort(key=lambda r:r[2])
for r in rows:
    if r[0]==r[1] or (r[0],r[1]) in [(2,4),(2,16)]: print("%2d-%-2d minS %5.2f meanS %6.2f n0 %d"%(r[0],r[1],r[2],r[3],r[4]+2))
x=np.log([r[2] for r in rows]); y=np.log([r[4]+2 for r in rows]); print("corr", np.corrcoef(x,y)[0,1])
