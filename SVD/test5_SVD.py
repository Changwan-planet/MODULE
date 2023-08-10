import numpy as np
import matplotlib.pyplot as plt
import math
#import dtw

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



line = 0


bscan = data_1[:,line,:].T 
print( np.allclose(bscan, data_1[:,line,:].T) )
print(bscan.shape)


bscan_flat = np.zeros( shape = (5096,41) )
bscan_flat2 = np.zeros( shape = (5096,41) )
bscan_flat3 = np.zeros( shape = (5096,41) )



bscan_flat[1000:,:] = bscan 
print(bscan_flat.shape)


k_x = 0
k_y = 0 

l_p   = 5

hp_p  = 6.0
hp_xd = 4   #hyperpoblic functiom x-direction movment

base  = 1.15

for j in np.arange(0,41,1):
#linear
#   k_y = k_y + l_p
 
#hyperbolic
  k_y = 10 +  int( (1.0/hp_p)**2 * ((k_x-hp_xd)**2) )   #for line 1

#logarimthic
#    k_y = int (np.log10(k_x))

#exponential
#    k_y = int (np.exp((k_x-10)) )                      
#    k_y = int (np.power(base,k_x))

#print("k_y=",k_y)
  for i in np.arange(1000,5096,1):
       bscan_flat2[i-k_y,j] = bscan_flat[i,j]
  
  k_x = k_x + 1



bscan_extent = (ax1_min,ax1_max,ay2_min,ay2_max)


plt.rcParams['figure.figsize'] = (10,13)
plt.rcParams['font.size'] = 15
#=========================================================
#Full matrix
#U, s, Vt = np.linalg.svd(data_1[:,line,:].T, full_matrices=True)

#Economy-sized decomposition 
#U, s, Vt = np.linalg.svd(data_1[:,line,:].T, full_matrices=False)
U, s, Vt = np.linalg.svd(bscan_flat2, full_matrices=False)

print("U.shape=", U.shape,",", "s.shape=",",", s.shape,",", "Vt.shape=", Vt.shape)
#print( np.allclose(data_1[:,line,:].T, np.dot(U[:, :41] * s, Vt)) ) 


t = np.dot(U[:, :41] * s, Vt) 

bscan_extent = (ax1_min,ax1_max,ay2_min,ay2_max)
egi = np.zeros(shape=(5096,41,41))
egi2 = np.zeros(shape=(5096,41))
egi_ftd = np.zeros(shape = (5096,41))



#=========================================================
#compress image by truncating SVD
for i in np.arange(0,41,1):

#High-pass filter
#    egi[:,:,i] =  np.dot(np.dot(U[:,:i], np.diag(s[:i])), Vt[:i,:]) 

#Low-pass filter 
    egi[:,:,i] =  np.dot(np.dot(U[:,i:], np.diag(s[i:])), Vt[i:,:]) 

#Band-pass filter
b = 0
p = 1

print("U[:,",b,",",p,"]=", U[:,b:p].shape )
print("diag(s[",b,":",p,"]",np.diag(s[b:p]).shape )
print("Vt[",b,":",p,":]=", Vt[b:p,:].shape )


for i in np.arange(b, p, 1):
    egi2[:,:] =  np.dot(np.dot(U[:,b:p], np.diag(s[b:p])), Vt[b:p,:]) 

#=========================================================



#define the exponent of color limit
cl = 5.5

#=========================================================
plt.subplot(1,4,1)
plt.imshow(bscan_flat[1000:,:],  extent = bscan_extent
          ,cmap="Greys_r",aspect= 3
)
plt.clim(-10**(cl),10**(cl))


#---------------------------------------------------------
plt.subplot(1,4,2)
#plt.imshow(data_1[:,line,:].T,extent = bscan_extent
plt.imshow(bscan_flat2[1000:,:],extent = bscan_extent
           ,cmap="Greys_r",aspect= 3
)

plt.clim(-10**(cl),10**(cl))


#---------------------------------------------------------
plt.subplot(1,4,3)
plt.imshow(egi2[1000:,:],extent = bscan_extent
           ,cmap="Greys_r",aspect= 3
)
plt.clim(-10**(cl),10**(cl))

#---------------------------------------------------------
k_x = 0
k_y = 0

for j in np.arange(0,41,1):
#linear
#   k_y = k_y + l_p

#hyperbolic 
   k_y =  int( (1.0/hp_p)**2 * ((k_x-hp_xd)**2) )

#logarimthic
#    k_y = int (np.log10(k_x))

#exponential
#    k_y = int (np.exp(k_x))                       
#     k_y = int (np.power(base,k_x))

   for i in np.arange(1000,5096,1):
       bscan_flat3[i,j] = egi2[i-k_y,j]

   k_x = k_x + 1

plt.subplot(1,4,4)
plt.imshow(bscan_flat3[1000:,:],extent = bscan_extent
           ,cmap="Greys_r",aspect= 3
) 
plt.clim(-10**(cl),10**(cl))

#=========================================================



plt.show()


