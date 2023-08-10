import numpy as np


 # 행렬 A는 이전 포스트에 있는 예제에서 사용된 것이다.
A = [[i**2, i, 1] for i in range(1930, 2020, 10)]


 # 행렬 A를 numpy 함수들이 다룰 수 있는 numpy array으로 변환한다.
matA = np.array(A).astype(np.float64)


# svd함수를 사용하여 3개의 반환값(U,s,V)를 저장한다.
U, s, V = np.linalg.svd(matA, full_matrices = True)


 # s는 matA의 고유값(eigenvalue) 리스트이다.
 # svd를 이용하여 근사한(approximated) 결과를 원본과 비교하기 위해
 # s를 유사대각행렬로 변환한다.
 # 유사행렬에 대한 내용은 이전 포스트 참조.

S = np.zeros(matA.shape)

for i in range(len(s)):
     S[i][i] = s[i]


 # 근사한 결과를 계산한다.
appA = np.dot(U, np.dot(S, V))


 # 원래 행렬인 matA와 근사한 행렬인 appA가 서로 비슷하다면
 # 두 행렬의 차이는 영행렬(zero matrix)에 가까울 것이다.
 # 즉, 다음의 결과가 대부분 0으로 채워져있다면 성공적인 svd가 이루어진 것이다.
 # 파이썬 버전에 따라 출력방식이 다르므로 주의한다.

print(matA - appA)
