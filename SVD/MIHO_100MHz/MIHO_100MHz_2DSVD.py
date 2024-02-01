import numpy as np
import matplotlib.pyplot as plt
import math
#from colormap_creator import create_od_seismic1_cmap
#import gprpyTools as gT

#create_od_seismic1_cmap()
#import dtw

#PATH
#input_path="/home/changwan/GPR_DATA/KOREA/MIHO_ri/3D_trench/40MHz/CSCAN3/EW/3D_CUBE_GPR.txt"

input_path="/home/csun/MODULE/SVD/MIHO_100MHz/BSCAN_GPR_rmbgr.txt"
#input_path="/home/changwan/GPR_DATA/KOREA/MIHO_ri/2022_BOMIN/Edit_data/Miho11/CUBE/3D_CUBE_GPR_EDIT.txt"

#input_path = "3D_HILBERT_powerdB_MaxdB.txt"

#READ DATASET
#inputdata = np.loadtxt(input_path)

#READ DATASET
data=np.loadtxt(input_path)
#inputdata=np.loadtxt(input_path)

#extract the coordinates and the value
#coords = inputdata[:,0:3].astype(int)
#values = inputdata[:,3]

rows = data.shape[0]
dis = data.shape[1]

print(data.shape[0])
print(data.shape[1])


#print(coords)

#data_1 = np.zeros( shape = (dis, tra, rows) )
#data_1.flat[np.ravel_multi_index(coords.T, data_1.shape)] = values

#CALCULATE THE DISTANCE INTERVAL:
c_sl = 3*10**8   # speed of light
#t_window = 10**(-7)
permit = 16.0
depth_range = 50
#depth_range = t_window * (c_sl/math.sqrt(permit)) * 0.5 
 
dis_int = 0.05

depth_int = (depth_range/np.sqrt(permit)) /rows 
depth_int2 = round(depth_int, 2)

print("\n")
print("++++++DATA INFO++++++")
print("depth_interval=",depth_int)

#C_SCAN_IMAGE

#RESAHPE THE INPUT DATA
print("input_shape (z,x)=",data.shape)    
#100MHz
print("+++++++++++++++++++++")
print("\n")

#plt.imshow(data)
#plt.show()


ax1_min=0
#ax1_max=data_1.shape[0]*0.01  #This is Northing.
ax1_max=data.shape[0] * dis_int   #This is Northing.

ay1_min=0
ay1_max=data.shape[1] * dis_int   #This is Easting.

#ay2_max=0
#ay2_min=data_1.shape[2] * depth_int #This is Easting.

#bscan_extent = (ax1_min,ax1_max,ay2_min,ay2_max)


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
#U2[:,:], s2[:], Vt2[:,:] = np.linalg.svd(data[:,:], full_matrices=True)


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

#for line in np.arange(0,tra,1):
#    for i in np.arange(b, p, 1):
#        egi2[:,line,:] =  np.dot(np.dot(U2[:,line,b:p], np.diag(s2[b:p,line])), Vt2[b:p,line,:]) 

#    for i in np.arange(b, p, 1):
#        egi2_lp[:,line,:] =  np.dot(np.dot(U2[:,line,0:b], np.diag(s2[0:b,line])), Vt2[0:b,line,:]) 

egi2[:,:] =  np.dot(np.dot(U2[:,:], np.diag(s2[:])), Vt2[:,:]) 


#bad-pass filter
egi2_bp[:,:] =  np.dot(np.dot(U2[:,b:p], np.diag(s2[b:p])), Vt2[b:p,:]) 


#low-pass filter
egi2_lp[:,:] =  np.dot(np.dot(U2[:,0:b], np.diag(s2[0:b])), Vt2[0:b,:]) 

#high-pass filter
egi2_hp[:,:] =  np.dot(np.dot(U2[:,b:], np.diag(s2[b:])), Vt2[b:,:]) 



#automatic gain control
#agc_w = 3
#data_agc = gT.agcGain(data_1[:,line,:].T,agc_w)
#egi2_agc = gT.agcGain(egi2,agc_w)
#egi2_lp_agc = gT.agcGain(egi2_lp,agc_w)

#x, y, z = np.indices(egi2.shape)

#x_flat = x.flatten(order='F').astype(int)
#y_flat = y.flatten(order='F').astype(int)
#z_flat = z.flatten(order='F').astype(int)
#data_flat = egi2.flatten(order='F')

#data_out = np.column_stack((x_flat, y_flat, z_flat, data_flat))


#transpose from (rows, tra, dis) into (dis, tra, rows)
#egi3= egi2.T

shape = egi2.shape
print("export shape=",shape)
data_out = []

#for k in range(shape[2]):
#    for j in range(shape[1]):
#        for i in range(shape[0]):
#                value = egi2[i,j,k]              
#                data_out.append([i, j, k, value])

#data_out = np.array(data_out)
#np.savetxt('Miho11_SVD.txt', data_out, fmt=('%d', '%d', '%d', '%.6f') )


#np.savetxt('test.txt', data_out, fmt=('%d', '%d', '%d', '%.6f') )

print("complete to export SVD-filterd data")


cbar_min = -3*10**4
cbar_max = +3*10**4
#bscan_extent = (ax1_min,ax1_max,ay2_min,ay2_max)

#row_sample
rs = 500
re = 1500


plt.subplot(1,4,1)
plt.imshow(data[rs:re,:], cmap = 'gray', vmin=cbar_min, vmax=cbar_max)

plt.subplot(1,4,2)
plt.imshow(egi2_bp[rs:re,:], cmap = 'gray', vmin=cbar_min, vmax=cbar_max)

plt.subplot(1,4,3)
plt.imshow(egi2_lp[rs:re,:], cmap = 'gray', vmin=cbar_min, vmax=cbar_max)

plt.subplot(1,4,4)
plt.imshow(egi2_hp[rs:re,:], cmap = 'gray', vmin=cbar_min, vmax=cbar_max)


plt.show()

"""
#plt.rcParams['figure.figsize'] = [20, 10]
plt.subplot(2,3,1)
plt.imshow(data_1[:,line,:].T, extent = bscan_extent
#           ,cmap="od_seismic1", vmin=cbar_min, vmax=cbar_max, aspect=7) 
#           ,cmap="Greys_r", vmin=cbar_min, vmax=cbar_max, aspect= 7)
           ,cmap="Greys_r", aspect= 7)

plt.subplot(2,3,2)
#plt.imshow(egi2_agc, extent = bscan_extent
plt.imshow(egi2, extent = bscan_extent
#           ,cmap="od_seismic1", vmin=cbar_min, vmax=cbar_max, aspect=7) 
#           ,cmap="Greys_r", vmin=cbar_min, vmax=cbar_max, aspect= 7)
            ,cmap="Greys_r", aspect= 7)
   
plt.subplot(2,3,3)
#plt.imshow(egi2_lp_agc, extent = bscan_extent
plt.imshow(egi2_lp, extent = bscan_extent
#           ,cmap="od_seismic1", vmin=cbar_min, vmax=cbar_max, aspect=7) 
#           ,cmap="Greys_r", vmin=cbar_min, vmax=cbar_max, aspect= 7)
            ,cmap="Greys_r", aspect= 7)

plt.subplot(2,3,4)
plt.plot(data_1[10,line,:])

plt.subplot(2,3,5)
plt.plot(egi2[:,0])
plt.plot(egi2_agc[:,0])



plt.subplot(2,3,6)
plt.plot(egi2_lp[:,0])


#define the exponent of color limit
#cl = 5.5
cl = 6


#=========================================================
plt.subplot(1,4,1)
plt.imshow(bscan_flat[:,:],  extent = bscan_extent
          ,cmap="Greys_r",aspect= 3
)
plt.clim(-10**(cl),10**(cl))


#---------------------------------------------------------
plt.subplot(1,4,2)
plt.imshow(data_1[:,line,:].T,extent = bscan_extent
#plt.imshow(bscan_flat2[100:,:],extent = bscan_extent
           ,cmap="Greys_r",aspect= 3
)

plt.clim(-10**(cl),10**(cl))


#---------------------------------------------------------
plt.subplot(1,4,3)
plt.imshow(egi2[:,:],extent = bscan_extent
           ,cmap="Greys_r",aspect= 3

)
plt.clim(-10**(cl),10**(cl))

#---------------------------------------------------------

k_x = 0
k_y = 0

for j in np.arange(0,dis-k,1):
#linear
   k_y = k_y + l_p

#hyperbolic 
#   k_y =  int( (1.0/hp_p)**2 * ((k_x-hp_xd)**2) )

#logarimthic
#    k_y = int (np.log10(k_x))

#exponential
#    k_y = int (np.exp(k_x))                       
#     k_y = int (np.power(base,k_x))

   for i in np.arange(100,612,1):
       bscan_flat3[i,j] = egi2[i-k_y,j]

   k_x = k_x + 1

plt.subplot(1,4,4)
plt.imshow(bscan_flat3[100:,:],extent = bscan_extent
           ,cmap="Greys_r",aspect= 3
) 
plt.clim(-10**(cl),10**(cl))

#=========================================================


plt.show()

"""
