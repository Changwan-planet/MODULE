# Re-importing necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import chatgpt_gpr_algorithm as cga


# Loading the data from the provided file
file_path = "/home/changwan/GPR_DATA/KOREA/MIHO_ri/3D_trench/40MHz/CSCAN3/EW/3D_BSCAN_GPR.txt"

data = np.loadtxt(file_path).T

# Apply AGC
agc_data = np.apply_along_axis(cga.agc, axis=1, arr=data)

# Compute structure tensor
st_xx, st_xy, st_yy = cga.structure_tensor(agc_data)

# Compute eigenvectors and ratio
_, eigenvectors = cga.eigenvectors(st_xx, st_xy, st_yy)
ratio = eigenvectors[:, :, 0, 1] / (eigenvectors[:, :, 0, 0] + 1e-9)

# Plot all graphs in a single figure
cga.plot([data, agc_data, ratio], [
    'Original Transposed Data',
    'AGC Transposed Data (Window Length = 11)',
    'Ratio of the Components of the Smaller Eigenvector'
])


