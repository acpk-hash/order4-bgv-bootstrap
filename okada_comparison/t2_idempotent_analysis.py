import sys
# Check whether the Lagrange-style idempotent e0(X^d) is a valid CRT idempotent
# modulo X^N+1 over F_p (i.e. u^2 == u mod X^N+1), and sparsity in X^(d/2).
def polmulmod(a,b,p,N):
    # multiply a,b mod (X^N + 1) over F_p ; a,b are dicts idx->coeff (sparse)
    res={}
    for i,ai in a.items():
        for j,bj in b.items():
            k=i+j
            c=ai*bj
            if k>=N:
                k-=N; c=-c  # X^N = -1
            res[k]=(res.get(k,0)+c)%p
    return {k:v for k,v in res.items() if v%p!=0}

def build_eY(p,n):
    roots=[y for y in range(p) if pow(y,n,p)==(p-1)]
    assert len(roots)==n
    z0=roots[0]; others=roots[1:]
    num=[1]
    for zj in others:
        new=[0]*(len(num)+1)
        for i,c in enumerate(num):
            new[i]=(new[i]-c*zj)%p
            new[i+1]=(new[i+1]+c)%p
        num=new
    denom=1
    for zj in others: denom=denom*((z0-zj)%p)%p
    dinv=pow(denom,p-2,p)
    return [(c*dinv)%p for c in num], roots

def test(p,N,d,n):
    eY,roots=build_eY(p,n)
    u={k*d:c for k,c in enumerate(eY) if c%p!=0}  # sparse dict in X
    # idempotent test: u*u == u mod X^N+1
    u2=polmulmod(u,u,p,N)
    ok = (u2==u)
    # sparsity in X^(d/2): all indices multiple of d/2
    sparse=all(k%(d//2)==0 for k in u)
    print("p=%d N=%d d=%d n=%d: idempotent u^2==u? %s ; sparse in X^(d/2)? %s ; #nz=%d"%(p,N,d,n,ok,sparse,len(u)))
    return u,roots

test(257,32768,256,128)
u1297,roots1297=test(1297,32768,4096,8)
print("1297 roots:", roots1297)
