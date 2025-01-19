#%% Importing Dependencies
import pyomo.environ as pyo
import pandas as pd
import numpy as np
from pyomo.opt import SolverFactory
#%% Model Definition
m = pyo.ConcreteModel()
#%% RangeSet Definition
#%% RangeSet Definition
m.I = pyo.RangeSet(10)   #suppliers
m.J = pyo.RangeSet(10)  #buyer(hospitals)
m.P = pyo.RangeSet(6)   #product(medicine)
m.S = pyo.RangeSet(6)   #scenarios
#%% Read Data from excel file
d_jps_ = pd.read_excel(r"problem3.xlsx", sheet_name="d_jps",header=None)
as_ip_ = pd.read_excel(r"problem3.xlsx", sheet_name="As_ip",header=None)
as_hat_= 0.1*as_ip_
ar_jp_ = pd.read_excel(r"problem3.xlsx", sheet_name="Ar_jp",header=None)
ar_hat_= 0.1*ar_jp_
hr_jp_ = pd.read_excel(r"problem3.xlsx", sheet_name="hr_jp",header=None)
hr_hat_ = 0.1*hr_jp_
ff_p_ = pd.read_excel(r"problem3.xlsx", sheet_name="ff_p",header=None)
ff_hat_ = 0.1 * ff_p_
pp_s_ = pd.read_excel(r"problem3.xlsx", sheet_name="pp_s",header=None)
Ocu_1_is_= pd.read_excel(r"problem3.xlsx", sheet_name="ocu1",header=None)
Ocu_2_is_= pd.read_excel(r"problem3.xlsx", sheet_name="ocu2",header=None)
Ocu_3_is_= pd.read_excel(r"problem3.xlsx", sheet_name="ocu3",header=None)
Em_1_is_= pd.read_excel(r"problem3.xlsx", sheet_name="em1",header=None)
Em_2_is_= pd.read_excel(r"problem3.xlsx", sheet_name="em2",header=None)
Em_3_is_= pd.read_excel(r"problem3.xlsx", sheet_name="em3",header=None)
Eg_1_is_= pd.read_excel(r"problem3.xlsx", sheet_name="eg1",header=None)
Eg_2_is_= pd.read_excel(r"problem3.xlsx", sheet_name="eg2",header=None)
Eg_3_is_= pd.read_excel(r"problem3.xlsx", sheet_name="eg3",header=None)
#%% convert to dict
ffp = list(ff_p_[0].array)
ffhat = list(ff_hat_[0].array)

pps = list(pp_s_[0].array)
d_jps = {(j,p,s):d_jps_.iloc[((j-1)*len(m.P) + p)-1][s-1] for j in m.J for p in m.P for s in m.S}

as_ip = {(i,j):as_ip_.iloc[i-1][j-1] for i in m.I for j in m.P} 
as_hat={(i,j):as_hat_.iloc[i-1][j-1] for i in m.I for j in m.P}

ar_jp = {(i,j):ar_jp_.iloc[i-1][j-1] for i in m.J for j in m.P}
ar_hat = {(i,j):ar_hat_.iloc[i-1][j-1] for i in m.J for j in m.P}

hr_jp = {(i,j):hr_jp_.iloc[i-1][j-1] for i in m.J for j in m.P}
hr_hatt = {(i,j):hr_hat_.iloc[i-1][j-1] for i in m.J for j in m.P}

ff_p = {i:ffp[i-1] for i in m.P}
ff_hatt = {i:ffhat[i-1] for i in m.P}

pp_s = {i:pps[i-1] for i in m.S} 
Ocu_1_is ={(i,j):Ocu_1_is_.iloc[i-1][j-1] for i in m.I for j in m.S}
Ocu_2_is ={(i,j):Ocu_2_is_.iloc[i-1][j-1] for i in m.I for j in m.S}
Ocu_3_is ={(i,j):Ocu_3_is_.iloc[i-1][j-1] for i in m.I for j in m.S}
Em_1_is={(i,j):Em_1_is_.iloc[i-1][j-1] for i in m.I for j in m.S}
Em_2_is={(i,j):Em_2_is_.iloc[i-1][j-1] for i in m.I for j in m.S}
Em_3_is={(i,j):Em_3_is_.iloc[i-1][j-1] for i in m.I for j in m.S}
Eg_1_is={(i,j):Eg_1_is_.iloc[i-1][j-1] for i in m.I for j in m.S}
Eg_2_is={(i,j):Eg_2_is_.iloc[i-1][j-1] for i in m.I for j in m.S}
Eg_3_is={(i,j):Eg_3_is_.iloc[i-1][j-1] for i in m.I for j in m.S}
#%% parameters
m.d = pyo.Param(m.J,m.P,m.S, initialize=d_jps)

m.As = pyo.Param(m.I,m.P, initialize=as_ip) 
m.As_hat = pyo.Param(m.I,m.P, initialize=as_hat)

m.Ar = pyo.Param(m.J,m.P, initialize=ar_jp)
m.Ar_hat = pyo.Param(m.J,m.P, initialize=ar_hat)

m.hr = pyo.Param(m.J,m.P, initialize=hr_jp)
m.hr_hat = pyo.Param(m.J,m.P, initialize=hr_hatt)

m.ff = pyo.Param(m.P, initialize=ff_p)
m.ff_hat = pyo.Param(m.P, initialize=ff_hatt)

m.ps = pyo.Param(m.S, initialize=pp_s)
m.ocu_1 = pyo.Param(m.I,m.S,initialize=Ocu_1_is)
m.ocu_2 = pyo.Param(m.I,m.S,initialize=Ocu_2_is)
m.ocu_3 = pyo.Param(m.I,m.S,initialize=Ocu_3_is)
m.em_1 = pyo.Param(m.I,m.S,initialize=Em_1_is)
m.em_2 = pyo.Param(m.I,m.S,initialize=Em_2_is)
m.em_3 = pyo.Param(m.I,m.S,initialize=Em_3_is)
m.eg_1 = pyo.Param(m.I,m.S,initialize=Eg_1_is)
m.eg_2 = pyo.Param(m.I,m.S,initialize=Eg_2_is)
m.eg_3 = pyo.Param(m.I,m.S,initialize=Eg_3_is)
m.c_ip = pyo.Param(m.I,m.P,initialize=0.05)
m.c_jp = pyo.Param(m.J,m.P,initialize=0.05)
m.c_p = pyo.Param(m.P,initialize=0.05)
m.c_i = pyo.Param(m.I,initialize=0.05) 
m.c = pyo.Param(initialize=0.05)

m.phi = pyo.Param(m.I,m.P,initialize=0)
m.phi_hatt = pyo.Param(m.I,m.P,initialize=0)

m.phi_hat = pyo.Param(m.I,m.P,initialize=3)
m.phi_hat_hat = pyo.Param(m.I,m.P,initialize=0.3)

m.e = pyo.Param(m.I,initialize=158)
m.e_hat = pyo.Param(m.I,initialize=15.8)

m.en = pyo.Param(m.I,initialize=1)
m.en_hat = pyo.Param(m.I,initialize=0.1)

m.o = pyo.Param(m.I,initialize=20)
m.o_hat = pyo.Param(m.I,initialize=2)

F1 = 189
F2 = 210
F3 = 231


M     = 50
landa = 1  

ee    = 0.03
ee_hat= 0.003

en    = 0.01
en_hat= 0.001

oo    = 0.05
oo_hat = 0.0005

bb    = 0.1
#%%
alfa_cut_ = [0,0.2,0.4,0.6,0.8,1]
objective_values = []
check_optimalities_condition = []
#%% decision variables
m.q = pyo.Var(m.I,m.J,m.P,m.S,domain=pyo.NonNegativeReals)

m.b = pyo.Var(m.I,m.J,m.P,m.S,domain=pyo.NonNegativeReals)

m.dd = pyo.Var(m.I,m.J,m.P,m.S,domain=pyo.Reals)
 
m.Cs = pyo.Var(m.S,domain=pyo.Reals)

m.z = pyo.Var(domain=pyo.Reals)

m.va = pyo.Var(m.S,domain=pyo.NonNegativeReals)

m.vb = pyo.Var(m.S,domain=pyo.NonNegativeReals)

#%% Objective
for alfa_cut in alfa_cut_:

    expression = m.z
    
    m.obj = pyo.Objective(expr= expression , sense=pyo.minimize)
    
    
    m.Constraints = pyo.ConstraintList()
    
    #0
    for s in m.S:
        m.Constraints.add(m.z >= m.Cs[s])
        
    
    #1
    for s in m.S :
         m.Constraints.add(m.Cs[s] == sum((m.As[i,p] + m.c_ip[i,p]*m.As_hat[i,p] + m.Ar[j,p] + m.c_jp[j,p]*m.Ar_hat[j,p]) * m.dd[i,j,p,s] / m.q[i,j,p,s] + (m.hr[j,p]+ m.c_jp[j,p]*m.hr_hat[j,p]) * (m.q[i,j,p,s] - m.b[i,j,p,s])**2 / (2 * m.q[i,j,p,s]) + (m.phi_hat[i,p]+ m.c_ip[i,p]*m.phi_hat_hat[i,p]) * m.b[i,j,p,s] / (2 * m.q[i,j,p,s]) + ((m.phi[i,p]+m.c_ip[i,p]*m.phi_hatt[i,p]) * m.b[i,j,p,s] * m.dd[i,j,p,s]) / m.q[i,j,p,s] for i in m.I for j in m.J for p in m.P))
    
    #2 Constraints of total order volumes:
    for i in m.I:
        for j in m.J:
            for p in m.P:
                for s in m.S:
                    m.Constraints.add((m.ff[p] + m.c_p[p]*m.ff_hat[p]) * m.q[i,j,p,s]<= alfa_cut*(F1+F2)/2 + (1-alfa_cut)*(F2+F3)/2)
                    
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
            m.Constraints.add((m.e[i]+m.c_i[i]*m.e_hat[i]) + (ee + m.c*ee_hat)*sum(m.q[i,j,p,s] for j in m.J for p in m.P)<= alfa_cut*(m.em_1[i,s]+ m.em_2[i,s])/2+(1-alfa_cut)*(m.em_2[i,s] + m.em_3[i,s])/2)
    
    #7
    for i in m.I:
        for s in m.S:
            m.Constraints.add((m.en[i] + m.c_i[i]*m.en_hat[i]) + (en + m.c*en_hat)*sum(m.q[i,j,p,s] for j in m.J for p in m.P)<= alfa_cut*(m.eg_1[i,s]+ m.eg_2[i,s])/2+(1-alfa_cut)*(m.eg_2[i,s]+m.eg_3[i,s])/2)
            
    #8
    for i in m.I:
        for s in m.S:
            m.Constraints.add((m.o[i] + m.c_i[i]*m.o_hat[i]) + (oo+ m.c * oo_hat)*sum(m.q[i,j,p,s] for j in m.J for p in m.P) >=  (1-alfa_cut)*(m.ocu_1[i,s]+ m.ocu_2[i,s])/2 + alfa_cut*(m.ocu_2[i,s]+m.ocu_3[i,s])/2)
            
    #9
    for s in m.S:
        m.Constraints.add(m.va[s] - m.vb[s] == m.Cs[s] - sum(m.ps[s] * m.Cs[s] for s in m.S))
        
        
    opt = SolverFactory('ipopt')
    opt.options['tol']=1e-6
    results = opt.solve(m,tee=True ) 
    objective_values.append(pyo.value(m.obj))
    check_optimalities_condition.append(pyo.check_optimal_termination(results))

print(f'objective_values : {objective_values}')
print(f'check_optimalities_condition : {check_optimalities_condition}')




















