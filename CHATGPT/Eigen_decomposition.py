# Function to perform eigen-decomposition on the structure tensor at each point
def eigen_decomposition(Ixx, Ixy, Iyy):
    eigenvalues1 = np.zeros_like(Ixx)
    eigenvalues2 = np.zeros_like(Ixx)
    orientations = np.zeros_like(Ixx)
    for i in range(Ixx.shape[0]):
        for j in range(Ixx.shape[1]):
            J = np.array([[Ixx[i, j], Ixy[i, j]], [Ixy[i, j], Iyy[i, j]]])
            eigvals, eigvecs = np.linalg.eigh(J)
            eigenvalues1[i, j] = eigvals[0]
            eigenvalues2[i, j] = eigvals[1]
            orientations[i, j] = np.arctan2(eigvecs[1, 1], eigvecs[0, 1])
    return eigenvalues1, eigenvalues2, orientations

