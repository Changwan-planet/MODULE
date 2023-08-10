# Function to perform eigen-decomposition on the structure tensor at each point
# and extract the dominant eigenvectors
def dominant_eigenvectors_decomposition(Ixx, Ixy, Iyy):
    dominant_eigenvectors_x = np.zeros_like(Ixx)
    dominant_eigenvectors_y = np.zeros_like(Ixx)
    for i in range(Ixx.shape[0]):
        for j in range(Ixx.shape[1]):
            J = np.array([[Ixx[i, j], Ixy[i, j]], [Ixy[i, j], Iyy[i, j]]])
            _, eigvecs = np.linalg.eigh(J)
            dominant_eigenvectors_x[i, j] = eigvecs[0, 1]
            dominant_eigenvectors_y[i, j] = eigvecs[1, 1]
    return dominant_eigenvectors_x, dominant_eigenvectors_y

# Compute the dominant eigenvectors of the structure tensor
dominant_eigenvectors_x, dominant_eigenvectors_y = dominant_eigenvectors_decomposition(Ixx, Ixy, Iyy)

# Calculate the local slopes using the negative ratio of the elements of the dominant eigenvector
print(dominant_eigenvectors_y.shape)
local_slopes_dominant = -dominant_eigenvectors_y / dominant_eigenvectors_x

# Visualize the local slopes with the specified extent
plt.figure(figsize=(8, 5))
plt.imshow(local_slopes_dominant, cmap='RdGy', extent=plot_extent_eigenvectors)
plt.title('Local Slopes (Dominant Eigenvector)')
plt.xlabel('Trace Number')
plt.ylabel('Depth / Time')
plt.colorbar(label='Slope')
plt.tight_layout()
plt.show()

# Return the local slopes for further analysis if needed
local_slopes_dominant

