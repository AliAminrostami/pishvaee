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


#%% Read Data from excel file

d_jps_ = pd.read_excel(r"problem8.xlsx", sheet_name="d_jps",header=None)
as_ip_ = pd.read_excel(r"problem8.xlsx", sheet_name="As_ip",header=None)
ar_jp_ = pd.read_excel(r"problem8.xlsx", sheet_name="Ar_jp",header=None)
hr_jp_ = pd.read_excel(r"problem8.xlsx", sheet_name="hr_jp",header=None)
ff_p_ = pd.read_excel(r"problem8.xlsx", sheet_name="ff_p",header=None)

#%% convert to dict
ffp = list(ff_p_[0].array)

d_jps = {(j,p):d_jps_.iloc[j-1][p-1] for j in m.J for p in m.P}
as_ip = {(i,j):as_ip_.iloc[i-1][j-1] for i in m.I for j in m.P} 
ar_jp = {(i,j):ar_jp_.iloc[i-1][j-1] for i in m.J for j in m.P}
hr_jp = {(i,j):hr_jp_.iloc[i-1][j-1] for i in m.J for j in m.P}
ff_p = {i:ffp[i-1] for i in m.P}


#%% parameters
m.d = pyo.Param(m.J,m.P, initialize=d_jps)
m.As = pyo.Param(m.I,m.P, initialize=as_ip) 
m.Ar = pyo.Param(m.J,m.P, initialize=ar_jp)
m.hr = pyo.Param(m.J,m.P, initialize=hr_jp)
m.ff = pyo.Param(m.P, initialize=ff_p)

m.ocu = pyo.Param(m.I,initialize={1:60,2:60,3:60,4:60,5:60,6:60,7:60,8:60,9:60,10:60})
m.em = pyo.Param(m.I,initialize={1:522,2:522,3:522,4:522,5:522,6:522,7:522,8:522,9:522,10:522})
m.eg = pyo.Param(m.I,initialize={1:70,2:70,3:70,4:70,5:70,6:70,7:70,8:70,9:70,10:70})
m.phi = pyo.Param(m.I,m.P,initialize=0)
m.phi_hat = pyo.Param(m.I,m.P,initialize=3)
m.e = pyo.Param(m.I,initialize=158)
m.en = pyo.Param(m.I,initialize=1)
m.o = pyo.Param(m.I,initialize=20)
#m.z_alf = pyo.Param(initialize=confidence_level_)


F     = 210
M     = 50
landa = 1

ee    = 0.03


en    = 0.01

oo    = 0.05

bb    = 0.1


#%% decision variables
m.q = pyo.Var(m.I,m.J,m.P,domain=pyo.NonNegativeReals)

m.b = pyo.Var(m.I,m.J,m.P,domain=pyo.NonNegativeReals)

m.dd = pyo.Var(m.I,m.J,m.P,domain=pyo.NonNegativeReals)

m.Cs = pyo.Var(domain=pyo.NonNegativeReals)

#m.z = pyo.Var(domain=pyo.Reals) 

#m.va = pyo.Var(m.S,domain=pyo.NonNegativeReals)

#m.vb = pyo.Var(m.S,domain=pyo.NonNegativeReals)

#%% Objective
#average_costs = sum(m.ps[s]*m.Cs[s] for s in m.S)
#coefficient_of_mean_deviation = z_alfa*a * sum(m.ps[s]*m.va[s]+m.ps[s]*m.vb[s] for s in m.S)
expression = m.Cs
#def objective_rule(model):
   # return sum(m.ps[s]*m.Cs[s] for s in m.S) + z_alfa * sum(m.ps[s]*m.va[s]+m.ps[s]*m.vb[s] for s in m.S)
m.obj = pyo.Objective(expr= expression , sense=pyo.minimize)

#%% Constraints
m.Constraints = pyo.ConstraintList()

#1

m.Constraints.add(m.Cs == sum((m.As[i,p] + m.Ar[j,p]) * m.dd[i,j,p] / m.q[i,j,p] + m.hr[j,p] * (m.q[i,j,p] - m.b[i,j,p])**2 / (2 * m.q[i,j,p]) + m.phi_hat[i,p] * m.b[i,j,p] / (2 * m.q[i,j,p]) + (m.phi[i,p] * m.b[i,j,p] * m.dd[i,j,p]) / m.q[i,j,p] for i in m.I for j in m.J for p in m.P))

#2 Constraints of total order volumes:
for i in m.I:
    for j in m.J:
        for p in m.P:
            m.Constraints.add(m.ff[p] * m.q[i,j,p] <= F)
                
#3 Constraints of resiliency in orders:
for i in m.I:
    for j in m.J:
        for p in m.P:
            m.Constraints.add((m.dd[i,j,p]/m.q[i,j,p]) <= landa*M)
                
#4 Constraints of shortage:
for i in m.I:
    for j in m.J:
        for p in m.P:
            m.Constraints.add(m.b[i,j,p]<=bb*m.q[i,j,p])
                
#5 Constraints of demand assignment:
for i in m.I:
    for j in m.J:
        for p in m.P:
            m.Constraints.add(m.dd[i,j,p] == m.d[j,p]/len(m.I))





#6 Constraints of sustainability (environments and energy consumption):
for i in m.I:
    m.Constraints.add(m.e[i] + ee*sum(m.q[i,j,p] for j in m.J for p in m.P)<=m.em[i])
        
#7
for i in m.I:
    m.Constraints.add(m.en[i] + en*sum(m.q[i,j,p] for j in m.J for p in m.P)<=m.eg[i])

#8
for i in m.I:
    m.Constraints.add(m.o[i] + oo*sum(m.q[i,j,p] for j in m.J for p in m.P)>=m.ocu[i])
        
#9
#for s in m.S:
 #   m.Constraints.add(m.va[s] - m.vb[s] == m.Cs[s] - sum(m.ps[s] * m.Cs[s] for s in m.S))
    
#%%        
opt = SolverFactory('ipopt')
results = opt.solve(m,tee=True)
print(results)
print(pyo.value(m.obj))
     
                
                




