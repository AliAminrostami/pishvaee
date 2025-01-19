#%%
import numpy as np
import pandas as pd
#%% d_jps 
def scenrio(s,k_s):
    return  0.1/(((k_s + 1)/2)- 1)*(s - 1)+ 0.9

x = np.array([scenrio(s,6) for s in range(1,7)]) # جای عدد 3 ، کاردینالیتی مجموعه ی اس رو قرار بدع . ودر قسمت رنج از عدد یک شرو کن تا کاردینالیتی اس به اضافه یک
rand_d = np.random.uniform(300,600, 60) # جای عدد 12 حاصلضرب تعداد مشتری در تعداد محصول رو قرار بدع
d_jps= np.array([i*j for i in rand_d for j in x]).reshape((60 ,6)) # در قسمت ریشیپ میشع تعداد سطر در تعداد سناریو

#%%pps
def p_ps(s,k_s):
    return s / (k_s * (k_s + 1) / 2)
pp_s = np.array([p_ps(s,6) for s in range(1,7)])
#%%eg1
su1 = np.array([63 for i in range(10)])# با این رنج مشخص میکنی که چنتا تامین کننده داری
eg1 = np.array([i*j for i in su1 for j in x]).reshape((10,6))
#%%eg2
su2 = np.array([70 for i in range(10)])# با این رنج مشخص میکنی که چنتا تامین کننده داری
eg2 = np.array([i*j for i in su2 for j in x]).reshape((10,6))
#%%eg3
su3 = np.array([77 for i in range(10)])# با این رنج مشخص میکنی که چنتا تامین کننده داری
eg3 = np.array([i*j for i in su3 for j in x]).reshape((10,6))
#%%em1
tu1= np.array([469.8 for i in range(10)])
em1= np.array([i*j for i in tu1 for j in x]).reshape((10,6))
#%%em2
tu2= np.array([522 for i in range(10)])
em2= np.array([i*j for i in tu2 for j in x]).reshape((10,6))
#%%em3
tu3= np.array([574.2 for i in range(10)])
em3= np.array([i*j for i in tu3 for j in x]).reshape((10,6))
#%%ocu1
xu1= np.array([54 for i in range(10)])
ocu1= np.array([i*j for i in xu1 for j in x]).reshape((10,6))
#%%ocu1
xu2= np.array([60 for i in range(10)])
ocu2= np.array([i*j for i in xu2 for j in x]).reshape((10,6))
#%%ocu1
xu3= np.array([66 for i in range(10)])
ocu3= np.array([i*j for i in xu3 for j in x]).reshape((10,6))
#%%
# Set the random seed for reproducibility
np.random.seed(0)

# Define the distributions based on Table 2
#d_jps = np.random.uniform(300, 600, 12)

A_sip = np.random.uniform(1, 5, size=(10, 6)) # Assuming 3 suppliers and 5 products
A_rjp = np.random.uniform(1, 4, size=(10, 6)) # Assuming 10 buyers and 5 products
h_rjp = np.random.uniform(2, 10, size=(10, 6)) # Assuming 10 buyers and 5 products
f_fp = np.random.uniform(1, 5, 6) # Assuming 5 products


#def scenrio(s,k_s):
 #   return  0.1/(((k_s + 1)/2)- 1)*(s - 1)+ 0.9


# p_ps calculation example
#def p_ps(s, s_j):
   # return s / (s_j * (s_j - 1) / 2)

# Create a dictionary to hold data frames
data_dict = {
    "d_jps": pd.DataFrame(d_jps),#, columns=["d_jps"]),
    "As_ip": pd.DataFrame(A_sip),#, columns=[f"A_sip_{i+1}" for i in range(A_sip.shape[1])]),
    "Ar_jp": pd.DataFrame(A_rjp),#, columns=[f"A_rjp_{i+1}" for i in range(A_rjp.shape[1])]),
    "hr_jp": pd.DataFrame(h_rjp),#, columns=[f"h_rjp_{i+1}" for i in range(h_rjp.shape[1])]),
    "ff_p": pd.DataFrame(f_fp),#, columns=["f_fp"]),
    "pp_s":pd.DataFrame(pp_s),
    "em1":pd.DataFrame(em1),
    "em2":pd.DataFrame(em2),
    "em3":pd.DataFrame(em3),
    "eg1":pd.DataFrame(eg1),
    "eg2":pd.DataFrame(eg2),
    "eg3":pd.DataFrame(eg3),
    "ocu1":pd.DataFrame(ocu1),
    "ocu2":pd.DataFrame(ocu2),
    "ocu3":pd.DataFrame(ocu3),
}

# Write data to Excel file
with pd.ExcelWriter("problem8.xlsx", engine="openpyxl") as writer:
    for sheet_name, df in data_dict.items():
        df.to_excel(writer, sheet_name=sheet_name, index=False,header=None)

print("Data successfully written to data_driven_model.xlsx")

