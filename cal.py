from cmath import sqrt
import cmath
import traceback
#import scipy
#import cyipopt
from scipy.optimize import minimize
def nzero(a):
    return 0.01 if (a==0) else a
class ring:
    rList=[]
    def __init__(self,S,m,p,init_r1,init_r2,init_g,scalar):
        # x=(R1,R2,g,omega)
        self.sa=lambda x:S*(x[p]+x[1+p])/nzero(x[p+1])+4*3.14*(x[1+p]*x[1+p]-x[p]*x[p])
        #self.sa=wraper(S,m,p)
        # limit

        self.lim6=lambda x:x[1+p]-x[p]-70
        self.con6={"type":"ineq","fun":self.lim6}

        self.lim7=lambda x:x[p+3]**2-(x[p+2]/nzero(x[p+1]))
        self.con7={"type":"eq","fun":self.lim7}

        self.cons=(self.con6,self.con7)
        self.bounds=((50,S/314),(50,S/314),(4.9,7.84),(None,None))

        self.args=(init_r1,init_r2,init_g,scalar*sqrt(init_g/init_r2).real)# the value should not be on the edge 
        self.m=m
        self.p=p
        self.S=S
        self.rList.append(self)
        self.id=len(self.rList)

def s(rings):
    # return lambda x:rings[0].sa(x)+rings[1].sa(x)+rings[2].sa(x)+rings[3].sa(x)
    def func(x):
        a=0
        for i in range(0,4):
          a+=rings[i].sa(x)
        #print(a)
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
            #print("R1:",r1)
            #print("R2:",r2)
        #print("L:",a)
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

def gInterval(rs):
    def interval(x):
        a={}
        for i in rs:
            a[x[i.p]]=i.id# R1
            a[x[i.p+1]]=i.id# R2
        lst=sorted(a)
        activeID=[]
        length=0
        for element in lst:
            if(a[element] in activeID):
                del activeID[activeID.index(a[element])]
            else:
                activeID.append(a[element])
            current_r1=x[rs[a[element]-1].p]
            current_r2=x[rs[a[element]-1].p+1]
            for i in activeID:
                if(i!=a[element]):
                    r1=x[rs[i-1].p]
                    r2=x[rs[i-1].p+1]
                    r1=current_r1 if(r1<current_r1) else r1
                    r2=current_r2 if(r2>current_r2) else r2
                    length+=r2-r1
            #length-=current_r2-current_r1
        #print(length)
        return nzero(length/2)
    return interval

con_eq+=({"type":"eq","fun":gInterval(ring.rList)},)
residential=ring(1369240,100,0,200,300,6,1)
tra_in=ring(1620240,100,4,300,400,6,-1)
arg=ring(235206,100,8,400,500,6,1)
ent=ring(235206,26.7,12,500,600,6,-1)
solution = minimize(s(ring.rList),gArgs(ring.rList),bounds=gBound(ring.rList),constraints=gCons(ring.rList)+con_eq,options={'maxiter': 1000},tol=1e-02)

print(solution)
print(solution.jac)
print(solution.x)