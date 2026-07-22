#!/usr/bin/env python3
"""Derive a wrap-correct mixed-radix butterfly factorization of the TRUE
HElib Step2 dim-0 block (m=50731, p=65537, D=96, slot ring F_p[X]/G, deg G=18).

Pipeline: y = P_rader( INTTchain( diag( NTTchain(v) ) ) ) with DC fixup.
NTT chain: radices [6,4,4] DIF, offsets per stage:
  stage1: {0,16,32,48,64,80}          (6 diagonals, wrap-safe subgroup)
  stage2: {0,+-4,+-8,+-12}            (7 diagonals)
  stage3: {0,+-1,+-2,+-3}             (7 diagonals)
All stage matrices built as explicit sparse matrices, verified vs DFT, and the
full pipeline verified against dump-matrix matvec on random vectors.
Outputs butterfly_stages.json for the C++ encrypted prototype.
"""
import json, random
from pathlib import Path
ROOT = Path("/home/luck/xzy/0424project")
OUT = ROOT/"order4bgv/tower_linear_transform/encrypted_prototype"
D, P = 96, 65537
cfg = json.load(open(ROOT/"order4bgv/tower_linear_transform/sage_verification/tower_stage_coefficients.json"))
G = cfg["G"]; d = len(G)-1

def pmul(a,b):
    res=[0]*(2*d-1)
    for i,ai in enumerate(a):
        if ai:
            for j,bj in enumerate(b):
                if bj: res[i+j]=(res[i+j]+ai*bj)%P
    for i in range(len(res)-1,d-1,-1):
        c=res[i]
        if c:
            for j in range(d+1): res[i-d+j]=(res[i-d+j]-c*G[j])%P
    return tuple(x%P for x in res[:d])
def padd(a,b): return tuple((x+y)%P for x,y in zip(a,b))
def psub(a,b): return tuple((x-y)%P for x,y in zip(a,b))
ZERO=tuple([0]*d); ONE=tuple([1]+[0]*(d-1))
def ppow(a,e):
    r=ONE; b=tuple(a)
    while e>0:
        if e&1: r=pmul(r,b)
        b=pmul(b,b); e>>=1
    return r
def pinv(a):
    def deg(x):
        for i in range(len(x)-1,-1,-1):
            if x[i]%P: return i
        return -1
    r0,r1=list(G),list(a)+[0]*(d+1-len(a))
    s0,s1=[0],[1]
    while deg(r1)>0:
        db=deg(r1); binv=pow(r1[db],P-2,P)
        rr=list(r0); q=[0]*(deg(r0)-db+1)
        while deg(rr)>=db:
            da=deg(rr); c=rr[da]*binv%P
            q[da-db]=c
            for i2 in range(db+1): rr[da-db+i2]=(rr[da-db+i2]-c*r1[i2])%P
        r0,r1=r1,rr
        qs=[0]*(len(q)+len(s1))
        for i2,qi in enumerate(q):
            if qi:
                for j2,sj in enumerate(s1):
                    if sj: qs[i2+j2]=(qs[i2+j2]+qi*sj)%P
        L=max(len(s0),len(qs))
        ns=[((s0[i2] if i2<len(s0) else 0)-(qs[i2] if i2<len(qs) else 0))%P for i2 in range(L)]
        s0,s1=s1,ns
    ci=pow(r1[0],P-2,P)
    inv=[x*ci%P for x in s1]
    return tuple((inv+[0]*d)[:d])

# ---- load dump ----
lines=(ROOT/"github_order4_cleaner/step2_matrix_dim0_sz96.txt").read_text().strip().split("\n")
S=[]
for i in range(1,D+1):
    row=[]
    for e in lines[i].split("|"):
        cs=[int(c)%P for c in e.strip()[1:-1].split()]
        row.append(tuple((cs+[0]*d)[:d]))
    S.append(row)

# ---- find omega: primitive 96th root of unity in F_p^18 ----
random.seed(42)
pd_minus1 = P**d - 1
assert pd_minus1 % 96 == 0
expo = pd_minus1 // 96
omega=None
for trial in range(50):
    a=tuple(random.randrange(P) for _ in range(d))
    w=ppow(a,expo)
    if w!=ONE and ppow(w,48)!=ONE and ppow(w,32)!=ONE:
        assert ppow(w,96)==ONE
        omega=w; break
assert omega is not None
print("omega found, order-96 verified")
wp=[ONE]
for k in range(1,96): wp.append(pmul(wp[-1],omega))
def wpow(e): return wp[e%96]

# ---- DIF stages, radices [6,4,4] ----
# stage: block length L, radix r, m=L/r. positions x = b + j + i*m
# DIF butterfly with twiddle:
#   out[b+j+i*m] = sum_i2 in[b+j+i2*m] * w_r^(i*i2) * w_L^(i*j)
# where w_r = omega^(96/r), w_L = omega^(96/L)
def stage_matrix(L,r):
    m=L//r
    M={}
    for b in range(0,D,L):
        for j in range(m):
            for i in range(r):
                row=b+j+i*m
                tw=wpow((96//L)*i*j)
                for i2 in range(r):
                    col=b+j+i2*m
                    val=pmul(wpow((96//r)*i*i2), tw)
                    M[(row,col)]=val
    return M
stages=[stage_matrix(96,6), stage_matrix(16,4), stage_matrix(4,4)]
def offsets_of(M):
    return sorted(set((c-r)%D for (r,c) in M))
for si,M in enumerate(stages):
    print("stage",si+1,"offsets:",offsets_of(M))

def apply_sparse(M,v):
    out=[ZERO]*D
    for (r,c),val in M.items():
        out[r]=padd(out[r],pmul(val,v[c]))
    return out
def chain(v,Ms):
    for M in Ms: v=apply_sparse(M,v)
    return v

# ---- verify chain == R * DFT ----
Fcols=[]
for j in range(D):
    e=[ZERO]*D; e[j]=ONE
    Fcols.append(chain(e,stages))
Fmat=[[Fcols[j][i] for j in range(D)] for i in range(D)]
W=[[wpow(i*j) for j in range(D)] for i in range(D)]
Wrows={tuple(W[t]):t for t in range(D)}
perm=[Wrows.get(tuple(Fmat[x])) for x in range(D)]
assert all(t is not None for t in perm), "chain is not a permuted DFT"
print("NTT chain == digit-reversed DFT, perm sample:",perm[:8])

# ---- inverse stages ----
def invert_stage(M):
    from collections import defaultdict
    rows=defaultdict(dict)
    for (r,c),v in M.items(): rows[r][c]=v
    groups=defaultdict(list)
    for r,cols in rows.items():
        groups[tuple(sorted(cols))].append(r)
    Minv={}
    for colset,rws in groups.items():
        rws=sorted(rws); cols=list(colset)
        n=len(rws)
        A=[[rows[r][c] for c in cols] for r in rws]
        I=[[ONE if i==j else ZERO for j in range(n)] for i in range(n)]
        for col in range(n):
            piv=None
            for rr in range(col,n):
                if A[rr][col]!=ZERO: piv=rr;break
            A[col],A[piv]=A[piv],A[col]; I[col],I[piv]=I[piv],I[col]
            pi=pinv(A[col][col])
            A[col]=[pmul(x,pi) for x in A[col]]; I[col]=[pmul(x,pi) for x in I[col]]
            for rr in range(n):
                if rr!=col and A[rr][col]!=ZERO:
                    f=A[rr][col]
                    A[rr]=[psub(x,pmul(f,y)) for x,y in zip(A[rr],A[col])]
                    I[rr]=[psub(x,pmul(f,y)) for x,y in zip(I[rr],I[col])]
        for a2,c in enumerate(cols):
            for b2,r in enumerate(rws):
                if I[a2][b2]!=ZERO: Minv[(c,r)]=I[a2][b2]
    return Minv
inv_stages=[invert_stage(M) for M in reversed(stages)]
for si,M in enumerate(inv_stages):
    print("inv stage",si+1,"offsets:",offsets_of(M))
v=[tuple(random.randrange(P) for _ in range(d)) for _ in range(D)]
vv=chain(chain(v,stages),inv_stages)
print("inverse chain roundtrip:", vv==v)

# ---- circulant middle diag ----
# correlation matrix C[t][j] = c[(j-t)%96], kernel c = S[1] (row 1 of dump)
c=[S[1][u] for u in range(D)]
cols=[]
ok_diag=True
for j in range(D):
    e=[ZERO]*D; e[j]=ONE
    x=chain(e,inv_stages)
    y=[ZERO]*D
    for t in range(D):
        acc=ZERO
        for jj in range(D):
            acc=padd(acc,pmul(c[(jj-t)%D],x[jj]))
        y[t]=acc
    z=chain(y,stages)
    cols.append(z)
for i in range(D):
    for j in range(D):
        if i!=j and cols[j][i]!=ZERO: ok_diag=False
print("F*C*Finv diagonal:", ok_diag)
Dmid=[cols[i][i] for i in range(D)]

# ---- full pipeline vs dump matvec ----
g=5
dlog={}
gp=1
for t in range(96):
    dlog[gp]=t; gp=(gp*g)%97
def full_pipeline(v):
    x=chain(v,stages)
    x=[pmul(Dmid[i],xi) for i,xi in enumerate(x)]
    w=chain(x,inv_stages)
    y=[ZERO]*D
    dc=ZERO
    for j in range(D): dc=padd(dc,v[j])
    y[0]=dc
    for i in range(1,D):
        y[i]=w[dlog[i]]
    return y,w
ok_all=True
for trial in range(3):
    v=[tuple(random.randrange(P) for _ in range(d)) for _ in range(D)]
    y,w=full_pipeline(v)
    ref=[ZERO]*D
    for i in range(D):
        acc=ZERO
        for j in range(D):
            acc=padd(acc,pmul(S[i][j],v[j]))
        ref[i]=acc
    if y!=ref:
        ok_all=False
        bad=[i for i in range(D) if y[i]!=ref[i]]
        print("trial",trial,"MISMATCH slots:",len(bad), bad[:10])
    else:
        print("trial",trial,"full pipeline == dump matvec: OK")
print("ALL OK" if ok_all else "FAILED")

if ok_all:
    def exp_stage(M):
        offs=offsets_of(M)
        diags={}
        for o in offs:
            dv=[list(ZERO) for _ in range(D)]
            for (r,cc),val in M.items():
                if (cc-r)%D==o: dv[r]=list(val)
            diags[str(o)]=dv
        return {"offsets":offs,"diags":diags}
    stage3f={k:pmul(Dmid[k[0]],v2) for k,v2 in stages[2].items()}
    exp={
      "p":P,"D":D,"G":G,"d":d,
      "ntt_stages":[exp_stage(stages[0]),exp_stage(stages[1]),exp_stage(stage3f)],
      "intt_stages":[exp_stage(M) for M in inv_stages],
      "rader_perm_w_index_for_output":[(-1 if i==0 else dlog[i]) for i in range(D)],
      "notes":"y[i]=w[dlog5(i)] for i>=1; y[0]=sum_j v[j]; Dmid folded into ntt stage 3"
    }
    json.dump(exp,open(OUT/"butterfly_stages.json","w"))
    print("exported butterfly_stages.json")
