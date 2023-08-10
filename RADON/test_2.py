import matplotlib.pyplot as plt
import numpy as np

import pylops
from pylops.utils.wavelets import ricker

plt.close("all")
#np.random.seed(0)


#PATH
#Miho
#input_path="/home/changwan/GPR_DATA/KOREA/MIHO_ri/40MHz/2023/NS/3D_CUBE_GPR.txt"
input_path="/home/changwan/GPR_DATA/KOREA/MIHO_ri/40MHz/2023/NS/attribute/3D_HILBERT_powerdB.txt"

#READ DATASET
inputdata = np.loadtxt(input_path)

#extract the values and their coordinates
values = inputdata[:, 0]
coords = inputdata[:, 1:].astype(int) - 1 #subtract 1 due to indices starting at 0

#create your matrix
#(dis, tra, rows)
dis = 41
tra = 18
rows = 4096

data_1 = np.zeros( shape = (dis, tra, rows) )

#fill in your matrix at the appropriate coordinates
data_1.flat[np.ravel_multi_index(coords.T, data_1.shape)] = values


#CALCULATE THE DISTANCE INTERVAL:
c_sl = 3*10**8   # speed of light
depth_range = 50 # 90 m 
permit = 16
sample = 4096
depth_int = depth_range/ sample 
depth_int = round(depth_int, 2)

print("\n")
print("++++++DATA INFO++++++")
print("depth_interval=",depth_int)

#C_SCAN_IMAGE

#RESAHPE THE INPUT DATA
print("input_shape=",data_1.shape)    

#data2_2=data2.reshape(1900,17,4096)

#Miho-ri
#202303

print("3D_shape (x,y,z) =",data_1.shape)
print("+++++++++++++++++++++")
print("\n")

fig,host =plt.subplots()

print(data_1.shape[2])

ax1_min=0
#ax1_max=data_1.shape[0]*0.01  #This is Northing.
ax1_max=data_1.shape[0]*0.5  #This is Northing.


ay1_min=0
ay1_max=data_1.shape[1]*0.5 #This is Easting.

ay2_max=0
ay2_min=data_1.shape[2]*depth_int #This is Easting.


par = {"ox": 0, "dx": 2, "nx": 121, "ot": 0, "dt": 0.004, "nt": 100, "f0": 30}
# create axis
taxis, taxis2, xaxis, yaxis = pylops.utils.seismicevents.makeaxis(par)

# radon operator
npx = 61
pxmax = 5e-4  # s/m
px = np.linspace(-pxmax, pxmax, npx)

#plt.plot(px)
#plt.show()

Rop = pylops.signalprocessing.Radon2D(
    taxis, xaxis, px, kind="linear", interp=False, centeredh=False, dtype="float64"
)

print(Rop)

#plt.imshow(Rop)
#plt.show()

# adjoint Radon transform
xadj = Rop.H * y

# sparse Radon transform
xinv, niter, cost = pylops.optimization.sparsity.fista(
    Rop, y.ravel(), niter=15, eps=1e1
)
xinv = xinv.reshape(Rop.dims)

# filtering
xfilt = np.zeros_like(xadj)
xfilt[npx // 2 - 3 : npx // 2 + 4] = xadj[npx // 2 - 3 : npx // 2 + 4]

yfilt = Rop * xfilt

# filtering on sparse transform
xinvfilt = np.zeros_like(xinv)
xinvfilt[npx // 2 - 3 : npx // 2 + 4] = xinv[npx // 2 - 3 : npx // 2 + 4]

yinvfilt = Rop * xinvfilt




pclip = 0.7
fig, axs = plt.subplots(1, 5, sharey=True, figsize=(12, 5))
axs[0].imshow(
    y.T,
    cmap="gray",
    vmin=-pclip * np.abs(y).max(),
    vmax=pclip * np.abs(y).max(),
    extent=(xaxis[0], xaxis[-1], taxis[-1], taxis[0]),
)
axs[0].set(xlabel="$x$ [m]", ylabel="$t$ [s]", title="Data")
axs[0].axis("tight")
axs[1].imshow(
    xadj.T,
    cmap="gray",
    vmin=-pclip * np.abs(xadj).max(),
    vmax=pclip * np.abs(xadj).max(),
    extent=(px[0], px[-1], taxis[-1], taxis[0]),
)
axs[1].axvline(px[npx // 2 - 3], color="r", linestyle="--")
axs[1].axvline(px[npx // 2 + 3], color="r", linestyle="--")
axs[1].set(xlabel="$p$ [s/m]", title="Radon")
axs[1].axis("tight")
axs[2].imshow(
    yfilt.T,
    cmap="gray",
    vmin=-pclip * np.abs(yfilt).max(),
    vmax=pclip * np.abs(yfilt).max(),
    extent=(xaxis[0], xaxis[-1], taxis[-1], taxis[0]),
)
axs[2].set(xlabel="$x$ [m]", title="Filtered data")
axs[2].axis("tight")
axs[3].imshow(
    xinv.T,
    cmap="gray",
    vmin=-pclip * np.abs(xinv).max(),
    vmax=pclip * np.abs(xinv).max(),
    extent=(px[0], px[-1], taxis[-1], taxis[0]),
)
axs[3].axvline(px[npx // 2 - 3], color="r", linestyle="--")
axs[3].axvline(px[npx // 2 + 3], color="r", linestyle="--")
axs[3].set(xlabel="$p$ [s/m]", title="Sparse Radon")
axs[3].axis("tight")
axs[4].imshow(
    yinvfilt.T,
    cmap="gray",
    vmin=-pclip * np.abs(y).max(),
    vmax=pclip * np.abs(y).max(),
    extent=(xaxis[0], xaxis[-1], taxis[-1], taxis[0]),
)

axs[4].set(xlabel="$x$ [m]", title="Sparse filtered data")
axs[4].axis("tight")
#plt.tight_layout()

plt.show()

