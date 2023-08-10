import matplotlib.pyplot as plt
import numpy as np

import pylops
from pylops.utils.wavelets import ricker
from icecream import ic

plt.close("all")
np.random.seed(0)


par = {"ox": 0, "dx": 2, "nx": 121, "ot": 0, "dt": 0.004, "nt": 100, "f0": 30}

# linear events
v = 1500  # m/s
t0 = [0.1, 0.2, 0.3]  # s
theta = [0, 0, 0]
amp = [1.0, -2, 0.5]

# parabolic event
tp0 = [0.13]  # s
px = [0]  # s/m
pxx = [5e-7]  # s²/m²
ampp = [0.7]

# create axis
taxis, taxis2, xaxis, yaxis = pylops.utils.seismicevents.makeaxis(par)

#print(taxis)
#print(taxis2)
#print(xaxis)
#print(yaxis)


#plt.plot(taxis)
#plt.plot(taxis2)
#plt.plot(xaxis)
#plt.plot(yaxis)
#plt.show()
ic(taxis)
ic(taxis2)
ic(xaxis)
ic(yaxis)
print("taxis.shape=",taxis.shape)
print("taxis2.shape=",taxis2.shape)
print("xaxis.shape=",xaxis.shape)






# create wavelet : source
wav = ricker(taxis[:41], f0=par["f0"])[0]


# generate model :synthetic b_scan image
y = (
    pylops.utils.seismicevents.linear2d(xaxis, taxis, v, t0, theta, amp, wav)[1]
    + pylops.utils.seismicevents.parabolic2d(xaxis, taxis, tp0, px, pxx, ampp, wav)[1]
)


#print(y)
#plt.imshow(y.T)
#plt.imshow(y)
#plt.show()
#print(y.shape)

# radon operator
npx = 61
pxmax = 5e-4  # s/m
px = np.linspace(-pxmax, pxmax, npx)

ic(px)
print(px.shape)

Rop = pylops.signalprocessing.Radon2D(
    taxis, xaxis, px, kind="linear", interp=False, centeredh=False, dtype="float64"
)

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





'''
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
'''
