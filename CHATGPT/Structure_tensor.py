from scipy.ndimage.filters import gaussian_filter

# Function to compute the structure tensor for a 2D array
def compute_structure_tensor(data, sigma=1):
    # Compute the gradients in the x (horizontal) and y (vertical) directions
    Ix, Iy = np.gradient(data)
    
    # Compute the components of the structure tensor
    Ixx = gaussian_filter(Ix * Ix, sigma=sigma)
    Ixy = gaussian_filter(Ix * Iy, sigma=sigma)
    Iyy = gaussian_filter(Iy * Iy, sigma=sigma)
    
    return Ixx, Ixy, Iyy

# Compute the structure tensor for the original transposed data
Ixx, Ixy, Iyy = compute_structure_tensor(transposed_data, sigma=2)

# Visualize the components of the structure tensor
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))
ax1.imshow(Ixx, cmap='RdGy', aspect='auto')
ax1.set_title('Ixx (Horizontal Gradient Squared)')
ax1.set_xlabel('Trace Number')
ax1.set_ylabel('Depth / Time')

ax2.imshow(Ixy, cmap='RdGy', aspect='auto')
ax2.set_title('Ixy (Horizontal-Vertical Gradient Product)')
ax2.set_xlabel('Trace Number')
ax2.set_ylabel('Depth / Time')

ax3.imshow(Iyy, cmap='RdGy', aspect='auto')
ax3.set_title('Iyy (Vertical Gradient Squared)')
ax3.set_xlabel('Trace Number')
ax3.set_ylabel('Depth / Time')

plt.show()
