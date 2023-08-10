import numpy as np


a = np.random.randn(9, 6) + 1j*np.random.randn(9, 6)
b = np.random.randn(2, 7, 8, 3) + 1j*np.random.randn(2, 7, 8, 3)


U, S, Vh = np.linalg.svd(a, full_matrices=True)


print( U.shape, S.shape, Vh.shape )
print(U[:,:6].shape)
print( np.allclose(a, np.dot(U[:, :6] * S, Vh)) )


smat = np.zeros((9, 6), dtype=complex)
smat[:6, :6] = np.diag(S)


print(smat[:6, :6])
print(np.allclose(a, np.dot(U, np.dot(smat, Vh))))

