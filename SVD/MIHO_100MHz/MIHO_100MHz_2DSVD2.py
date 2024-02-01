import numpy as np
import matplotlib.pyplot as plt
import math

#PATH
input_path="/home/csun/MODULE/SVD/MIHO_100MHz/BSCAN_GPR_rmbgr.txt"

#READ DATASET
data=np.loadtxt(input_path)

rows = data.shape[0]
dis = data.shape[1]

print(data.shape[0])
print(data.shape[1])

#RESAHPE THE INPUT DATA
print("input_shape (z,x)=",data.shape)    
#100MHz
print("+++++++++++++++++++++")
print("\n")

plt.rcParams['figure.figsize'] = (20,13)
plt.rcParams['font.size'] = 15
#=========================================================
#Full matrix
#U, s, Vt = np.linalg.svd(data_1[:,line,:].T, full_matrices=True)

#Economy-sized decomposition 
U = np.zeros( shape = (rows, dis) )
s = np.zeros( shape = (dis ) )
Vt = np.zeros( shape = (dis,  dis) )

print("U_shape=",U.shape)

U2 = U
s2 = s
Vt2 = Vt

print(U.shape, s.shape, Vt.shape)
print(U2.shape, s2.shape, Vt2.shape)
print(data.shape)
U2[:,:], s2[:], Vt2[:,:] = np.linalg.svd(data[:,:], full_matrices=False)

#U, s, Vt = np.linalg.svd(bscan_flat2, full_matrices=False)
print("U.shape=", U.shape,",", "s.shape=",",", s.shape,",", "Vt.shape=", Vt.shape)

egi = np.zeros(shape=(rows,dis))
egi2 = np.zeros(shape=(rows,dis))
egi2_lp = np.zeros(shape=(rows,dis))
egi2_bp = np.zeros(shape = (rows,dis))
egi2_hp = np.zeros(shape = (rows,dis))

#Band-pass filter
b = 10
#p = dis
p = 100


print("U[:,",b,",",p,"]=", U[:,b:p].shape )
print("diag(s[",b,":",p,"]",np.diag(s[b:p]).shape )
print("Vt[",b,":",p,":]=", Vt[b:p,:].shape )

egi2[:,:] =  np.dot(np.dot(U2[:,:], np.diag(s2[:])), Vt2[:,:]) 


#bad-pass filter
egi2_bp[:,:] =  np.dot(np.dot(U2[:,b:p], np.diag(s2[b:p])), Vt2[b:p,:]) 

#low-pass filter
egi2_lp[:,:] =  np.dot(np.dot(U2[:,0:b], np.diag(s2[0:b])), Vt2[0:b,:]) 

#high-pass filter
egi2_hp[:,:] =  np.dot(np.dot(U2[:,b:], np.diag(s2[b:])), Vt2[b:,:]) 


cbar_min = -3*10**4
cbar_max = +3*10**4

#row_sample
rs = 500
re = 3000


plt.subplot(1,4,1)
plt.imshow(data[rs:re,:], cmap = 'gray', vmin=cbar_min, vmax=cbar_max)

plt.subplot(1,4,2)
plt.imshow(egi2_bp[rs:re,:], cmap = 'gray', vmin=cbar_min, vmax=cbar_max)

plt.subplot(1,4,3)
plt.imshow(egi2_lp[rs:re,:], cmap = 'gray', vmin=cbar_min, vmax=cbar_max)

plt.subplot(1,4,4)
plt.imshow(egi2_hp[rs:re,:], cmap = 'gray', vmin=cbar_min, vmax=cbar_max)


plt.show()

