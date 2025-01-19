#%% Importing Dependencies
#%% Importing Dependencies
import pyomo.environ as pyo
import pandas as pd
import numpy as np
from pyomo.opt import SolverFactory
#%% Model Definition
m = pyo.ConcreteModel()
#%% RangeSet Definition
m.I = pyo.RangeSet(10)   #suppliers
m.J = pyo.RangeSet(10)  #buyer(hospitals)
m.P = pyo.RangeSet(6)   #product(medicine)
m.S = pyo.RangeSet(6)   #scenarios

#%% Read Data from excel file

d_jps_ = pd.read_excel(r"problem3.xlsx", sheet_name="d_jps",header=None)
as_ip_ = pd.read_excel(r"problem3.xlsx", sheet_name="As_ip",header=None)
ar_jp_ = pd.read_excel(r"problem3.xlsx", sheet_name="Ar_jp",header=None)
hr_jp_ = pd.read_excel(r"problem3.xlsx", sheet_name="hr_jp",header=None)
ff_p_ = pd.read_excel(r"problem3.xlsx", sheet_name="ff_p",header=None)
Ocu_is_= pd.read_excel(r"problem3.xlsx", sheet_name="ocu2",header=None)
Em_is_=  pd.read_excel(r"problem3.xlsx", sheet_name="em2",header=None)
Eg_is_= pd.read_excel(r"problem3.xlsx", sheet_name="eg2",header=None)
pp_s_ = pd.read_excel(r"problem3.xlsx", sheet_name="pp_s",header=None)
#%% convert to dict
ffp = list(ff_p_[0].array)
pps = list(pp_s_[0].array)
d_jps = {(j,p,s):d_jps_.iloc[((j-1)*len(m.P) + p)-1][s-1] for j in m.J for p in m.P for s in m.S}
as_ip = {(i,j):as_ip_.iloc[i-1][j-1] for i in m.I for j in m.P} 
ar_jp = {(i,j):ar_jp_.iloc[i-1][j-1] for i in m.J for j in m.P}
hr_jp = {(i,j):hr_jp_.iloc[i-1][j-1] for i in m.J for j in m.P}
ff_p = {i:ffp[i-1] for i in m.P}
pp_s = {i:pps[i-1] for i in m.S} 
Ocu_is ={(i,j):Ocu_is_.iloc[i-1][j-1] for i in m.I for j in m.S}
Em_is={(i,j):Em_is_.iloc[i-1][j-1] for i in m.I for j in m.S}
Eg_is={(i,j):Eg_is_.iloc[i-1][j-1] for i in m.I for j in m.S}


#%% parameters
m.d = pyo.Param(m.J,m.P,m.S, initialize=d_jps)
m.As = pyo.Param(m.I,m.P, initialize=as_ip) 
m.Ar = pyo.Param(m.J,m.P, initialize=ar_jp)
m.hr = pyo.Param(m.J,m.P, initialize=hr_jp)
m.ff = pyo.Param(m.P, initialize=ff_p)
m.ps = pyo.Param(m.S, initialize=pp_s)
m.ocu = pyo.Param(m.I,m.S,initialize=Ocu_is)
m.em = pyo.Param(m.I,m.S,initialize=Em_is)
m.eg = pyo.Param(m.I,m.S,initialize=Eg_is)
m.phi = pyo.Param(m.I,m.P,initialize=0)
m.phi_hat = pyo.Param(m.I,m.P,initialize=3)
m.e = pyo.Param(m.I,initialize=158)
m.en = pyo.Param(m.I,initialize=1)
m.o = pyo.Param(m.I,initialize=20)

z_alfa = 4
a = 100    
F     = 210
M     = 50
landa = 1

ee    = 0.03


en    = 0.01

oo    = 0.05

bb    = 0.1


#%% decision variables
m.q = pyo.Var(m.I,m.J,m.P,m.S,domain=pyo.NonNegativeReals)

m.b = pyo.Var(m.I,m.J,m.P,m.S,domain=pyo.NonNegativeReals)

m.dd = pyo.Var(m.I,m.J,m.P,m.S,domain=pyo.Reals)

m.Cs = pyo.Var(m.S,domain=pyo.Reals)

m.z = pyo.Var(domain=pyo.NonNegativeReals)

m.va = pyo.Var(m.S,domain=pyo.NonNegativeReals)

m.vb = pyo.Var(m.S,domain=pyo.NonNegativeReals)

#%% Objective
average_costs = sum(m.ps[s]*m.Cs[s] for s in m.S)
coefficient_of_mean_deviation = z_alfa* sum(m.ps[s]*m.va[s]+m.ps[s]*m.vb[s] for s in m.S)
expression = average_costs + coefficient_of_mean_deviation
#def objective_rule(model):
   # return sum(m.ps[s]*m.Cs[s] for s in m.S) + z_alfa * sum(m.ps[s]*m.va[s]+m.ps[s]*m.vb[s] for s in m.S)
m.obj = pyo.Objective(expr= expression , sense=pyo.minimize)

#%% Constraints
m.Constraints = pyo.ConstraintList()

#1
for s in m.S :
    m.Constraints.add(m.Cs[s] == sum((m.As[i,p] + m.Ar[j,p]) * m.dd[i,j,p,s] / m.q[i,j,p,s] + m.hr[j,p] * (m.q[i,j,p,s] - m.b[i,j,p,s])**2 / (2 * m.q[i,j,p,s]) + m.phi_hat[i,p] * m.b[i,j,p,s] / (2 * m.q[i,j,p,s]) + (m.phi[i,p] * m.b[i,j,p,s] * m.dd[i,j,p,s]) / m.q[i,j,p,s] for i in m.I for j in m.J for p in m.P))

#2 Constraints of total order volumes:
for i in m.I:
    for j in m.J:
        for p in m.P:
            for s in m.S:
                m.Constraints.add(m.ff[p] * m.q[i,j,p,s] <= F)
                
#3 Constraints of resiliency in orders:
for i in m.I:
    for j in m.J:
        for p in m.P:
            for s in m.S:
                m.Constraints.add((m.dd[i,j,p,s]/m.q[i,j,p,s]) <= landa*M)
                
#4 Constraints of shortage:
for i in m.I:
    for j in m.J:
        for p in m.P:
            for s in m.S:
                m.Constraints.add(m.b[i,j,p,s]<=bb*m.q[i,j,p,s])
                
#5 Constraints of demand assignment:
for i in m.I:
    for j in m.J:
        for p in m.P:
            for s in m.S:
                m.Constraints.add(m.dd[i,j,p,s] == m.d[j,p,s]/len(m.I))





#6 Constraints of sustainability (environments and energy consumption):
for i in m.I:
    for s in m.S:
        m.Constraints.add(m.e[i] + ee*sum(m.q[i,j,p,s] for j in m.J for p in m.P)<=m.em[i,s])
        
#7
for i in m.I:
    for s in m.S:
        m.Constraints.add(m.en[i] + en*sum(m.q[i,j,p,s] for j in m.J for p in m.P)<=m.eg[i,s])

#8
for i in m.I:
    for s in m.S:
        m.Constraints.add(m.o[i] + oo*sum(m.q[i,j,p,s] for j in m.J for p in m.P)>=m.ocu[i,s])
        
#9
for s in m.S:
    m.Constraints.add(m.va[s] - m.vb[s] == m.Cs[s] - sum(m.ps[s] * m.Cs[s] for s in m.S))
    
#%%        
opt = SolverFactory('ipopt')
results = opt.solve(m,tee=True)
print(results)
print(pyo.value(m.obj))
     
                
                




