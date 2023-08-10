import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from scipy.interpolate import interp1d

#import ipyvolume as ipv
#import cv2
import math

#PATH
#input_path="/home/changwan/GPR_DATA/KOREA/MIHO_ri/40MHz/2023/EW/attribute/3D_HILBERT_powerdB.txt"
input_path="/home/changwan/GPR_DATA/KOREA/MIHO_ri/3D_trench/40MHz/CSCAN3/EW/3D_CUBE_GPR.txt"


#READ DATASET
inputdata = np.loadtxt(input_path)

#extract the values and their coordinates
values = inputdata[:, 0]
coords = inputdata[:, 1:].astype(int) - 1 #subtract 1 due to indices starting at 0

#create your matrix
#(dis, tra, rows)
dis = 41
tra = 21
rows = 4096

data_1 = np.zeros( shape = (dis, tra, rows) )

#fill in your matrix at the appropriate coordinates
data_1.flat[np.ravel_multi_index(coords.T, data_1.shape)] = values


#CALCULATE THE DISTANCE INTERVAL:
#c_sl = 3*10**8   # speed of light
#t_window = 10**(-7)
#permit = 16.0
depth_range = 50
#depth_range = t_window * (c_sl/math.sqrt(permit)) * 0.5 
 
dis_int = 0.5

depth_int = (depth_range/ rows) 
depth_int2 = round(depth_int, 2)

print("\n")
print("++++++DATA INFO++++++")
print("depth_interval=",depth_int)

#C_SCAN_IMAGE

#RESAHPE THE INPUT DATA
print("input_shape=",data_1.shape)    
#40MHz

print("3D_shape (x,y,z) =",data_1.shape)
print("+++++++++++++++++++++")
print("\n")

print(data_1.shape[2])

ax1_min=0
#ax1_max=data_1.shape[0]*0.01  #This is Northing.
ax1_max=data_1.shape[0] * dis_int   #This is Northing.

ay1_min=0
ay1_max=data_1.shape[1] * dis_int   #This is Easting.

ay2_max=0
ay2_min=data_1.shape[2] * depth_int #This is Easting.


#cmap=colors.ListedColormap(["white","blue","red"])

#     ++++++++++++++++++++++
#++++++Before interploation++++++
#     ++++++++++++++++++++++
start = 1037
end   = 2000


rows = list(range(start,end,1))

#samples=list(range(0,sample,1,))

dis_s = 20
dis_e = data_1.shape[0]
dis=list(range(dis_s,dis_e,1))
#dis = dis_e *dis_int

line = 1

#print(data_1[:,line,:].shape)
#print(data_1[:,line,:].T.shape)


#Full matrix
#U, s, Vt = np.linalg.svd(data_1[:,line,:].T, full_matrices=True)

#Economy-sized decomposition 
U, s, Vt = np.linalg.svd(data_1[:,line,:].T, full_matrices=False)


 
print("U.shape=", U.shape,",", "s.shape=",",", s.shape,",", "Vt.shape=", Vt.shape)
#print( np.allclose(data_1[:,line,:].T, np.dot(U[:, :41] * s, Vt)) ) 


t = np.dot(U[:, :41] * s, Vt) 
#fig, axes = plt.subplots(nrows=1, ncols=4)
#axes[0].imshow(U)
#axes[1].imshow(np.diag(s))
#axes[2].imshow(Vt)
#axes[3].imshow(eigenimage)
#plt.tight_layout()
#plt.show()

#print(s)
#plt.plot(s)
#plt.show()


bscan_extent = (ax1_min,ax1_max,ay2_min,ay2_max)
egi = np.zeros(shape=(4096,41,41))
egi2 = np.zeros(shape=(4096,41))


egi_ftd = np.zeros(shape = (4096,41))

# Define the desired color limits for the colorbar
cbar_min = -10**6
cbar_max = +10**6


#compress image by truncating SVD
for i in np.arange(0,41,1):

#High-pass filter
#    egi[:,:,i] =  np.dot(np.dot(U[:,:i], np.diag(s[:i])), Vt[:i,:]) 

#Low-pass filter 
    egi[:,:,i] =  np.dot(np.dot(U[:,i:], np.diag(s[i:])), Vt[i:,:]) 
#total eigenimage check

#for i in np.arange(0,41,1):
#    plt.subplot(2,21,i+1)
    #plt.rcParams['figure.figsize'] = [20, 10]
#    plt.imshow(egi[:,:,i],extent = bscan_extent
#           ,cmap="Greys_r", vmin=cbar_min, vmax=cbar_max 
#           ,aspect= 3
#)

#    plt.clim(-10**(6),10**(6))
#    plt.clim(-10**(5),10**(5))


#Band-pass filter
b = 0
p = 5

print("U[:,",b,",",p,"]=", U[:,b:p].shape )
print("diag(s[",b,":",p,"]",np.diag(s[b:p]).shape )
print("Vt[",b,":",p,":]=", Vt[b:p,:].shape )



for i in np.arange(b, p, 1):
    egi2[:,:] =  np.dot(np.dot(U[:,b:p], np.diag(s[b:p])), Vt[b:p,:]) 


#total eigenimage check
plt.subplot(1,2,1)
plt.imshow(data_1[:,line,:].T,extent = bscan_extent
           ,cmap="Greys_r", vmin=cbar_min, vmax=cbar_max 
           ,aspect= 3
)

plt.clim(-10**(5),10**(5))


plt.subplot(1,2,2)
 #plt.rcParams['figure.figsize'] = [20, 10]
plt.imshow(egi2[:,:],extent = bscan_extent
           ,cmap="Greys_r", vmin=cbar_min, vmax=cbar_max 
           ,aspect= 3
)

#plt.clim(-10**(6),10**(6))
plt.clim(-10**(5),10**(5))


    #plt.ylim(21.5,11.5)
   
plt.show()






#plt.imshow(egi[:,:,3],extent = bscan_extent,cmap="Greys_r",
#                        vmin=cbar_min, vmax=cbar_max)
#plt.show()


#"""
#band-pass eigenimage
#p = 0
#q = 16
p = 10
q = 20


for i in np.arange(p,q,1):
    egi_ftd = egi_ftd + egi[:,:,i]  

bscan_extent = (ax1_min,ax1_max,ay2_min,ay2_max)

fig, axes = plt.subplots(nrows=1, ncols=5, figsize = (20,15) )
im0 = axes[0].imshow(data_1[:,line,:].T
           ,extent=bscan_extent
           #,cmap='gist_rainbow'
           ,cmap="Greys_r", vmin=cbar_min, vmax=cbar_max 
           ,aspect= 3
)


im1 = axes[1].imshow(egi[:,:,0] ,extent=bscan_extent
           ,cmap="Greys_r", vmin=cbar_min, vmax=cbar_max
           ,aspect= 3
)
  
im2 = axes[2].imshow(egi[:,:,1] ,extent=bscan_extent
         ,cmap="Greys_r", vmin=cbar_min, vmax=cbar_max
           ,aspect= 3
)

egi_ftd = axes[3].imshow(egi_ftd ,extent=bscan_extent
         ,cmap="Greys_r", vmin=cbar_min, vmax=cbar_max
           ,aspect= 3
)


t = axes[4].imshow(t, extent=bscan_extent 
         ,cmap="Greys_r", vmin=cbar_min, vmax=cbar_max
           ,aspect= 3
)





#Amp
fig.colorbar(im0, ax=axes[0])
fig.colorbar(im1, ax=axes[1])
fig.colorbar(im2, ax=axes[2])
fig.colorbar(egi_ftd, ax=axes[3])
fig.colorbar(t, ax=axes[4])

im0.axes.set_ylim(21.5, 11.5)
im1.axes.set_ylim(21.5, 11.5)
im2.axes.set_ylim(21.5, 11.5)
egi_ftd.axes.set_ylim(21.5, 11.5)
t.axes.set_ylim(21.5,11.5)






plt.show()
#"""
