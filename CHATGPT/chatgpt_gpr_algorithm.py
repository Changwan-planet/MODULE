import numpy as np
import scipy.signal
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt

def agc(trace, window_length=11):
    if window_length % 2 == 0: window_length += 1
    envelope = np.abs(scipy.signal.hilbert(trace))
    gain = np.convolve(envelope, np.ones(window_length) / window_length, mode='same')
    return trace / (gain + 1e-9)

def structure_tensor(data, sigma_v=2.5, sigma_i=1):
    grad_v = scipy.ndimage.gaussian_filter1d(data, sigma_v, axis=0, order=1)
    grad_i = scipy.ndimage.gaussian_filter1d(data, sigma_i, axis=1, order=1)
    st_xx, st_xy, st_yy = grad_v * grad_v, grad_v * grad_i, grad_i * grad_i
    return [scipy.ndimage.gaussian_filter(x, (sigma_v, sigma_i)) for x in [st_xx, st_xy, st_yy]]

def eigenvectors(st_xx, st_xy, st_yy):
    tensor_shape, tensor_matrix = st_xx.shape, np.zeros(st_xx.shape + (2, 2))
    tensor_matrix[:, :, 0, 0], tensor_matrix[:, :, 0, 1], tensor_matrix[:, :, 1, 0], tensor_matrix[:, :, 1, 1] = st_xx, st_xy, st_xy, st_yy
    return np.linalg.eigh(tensor_matrix)

def plot(data_list, titles):
    fig, axes = plt.subplots(1, 3, figsize=[18, 4])
    for i, ax in enumerate(axes):
        im = ax.imshow(data_list[i], cmap='RdGy', extent=(0, 41 * 0.5, 50.0, 0))
        ax.set_title(titles[i])
        ax.set_xlabel('Dimension 2 (meters)')
        ax.set_ylabel('Dimension 1 (meters)')
        plt.colorbar(im, ax=ax, label='Amplitude')
    plt.tight_layout()
    plt.show()
