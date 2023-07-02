from cmath import sqrt
import scipy
from scipy.optimize import minimize
class ring:
    rList=[]
    def __init__(self,S,m,p):
        # x=(R1,R2,g,omega)
        self.sa=lambda x,*args:S*(x[p]+x[1+p])/x[p+1]+2*3.14*(x[1+p]*x[1+p]-x[p]*x[p])
        
        # limit
        self.lim1=lambda x:-sqrt(x[2+p]/x[2+p])+0.03
        self.con1={"type":"ineq","fun":self.lim1}
        self.lim4=lambda x:x[1+p]-1000*x[2+p]/9
        self.con4={"type":"ineq","fun":self.lim4}
        self.lim5=lambda x:-x[1+p]+S/314
        self.con5={"type":"ineq","fun":self.lim5}
        self.lim6=lambda x:x[1+p]-x[p]-70
        self.con6={"type":"ineq","fun":self.lim6}
        self.cons=(self.con1,self.con4,self.con5,self.con6)
        self.bounds=((0,None),(0,None),(4.9,7.84),(None,None))
        self.args=(500,600,6,0.02)
        self.m=m
        self.p=p
        self.S=S
        self.rList.append(self)
def s(rings):
    # return lambda x:rings[0].sa(x)+rings[1].sa(x)+rings[2].sa(x)+rings[3].sa(x)
    def func(x):
        a=0
        for i in range(0,4):
          a+=rings[i].sa(x)
        return a
    return func
def eq(rs):
    def eq1(x):
        a=0
        for i in range(0,4):
            r1=x[rs[i].p]
            r2=x[rs[i].p+1]
            omega=x[rs[i].p+3]
            a+=rs[i].m/2*(r2*r2+r1*r1)*omega
        return a
    return eq1

def gArgs(rs):
    a=()
    for i in range(0,4):
        a+=rs[i].args
    return a
con_eq=({"type":"eq","fun":eq(ring.rList)},)

def gCons(rs):
    a=()
    for i in range(0,4):
        a+=rs[i].cons
    return a
def gBound(rs):
    a=()
    for i in range(0,4):
        a+=rs[i].bounds
    return a

residential = ring(1369240,10000,0)
tra_in=ring(1620240,90000,4)
arg=ring(235206,10000,8)
ent=ring(235206,10000,12)
solution = minimize(s(ring.rList),gArgs(ring.rList),bounds=gBound(ring.rList),method='SLSQP',constraints=gCons(ring.rList)+con_eq)
print(solution.success)
print(solution.x)
print(solution.message)