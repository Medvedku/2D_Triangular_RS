from lib.functions import D_mat, B_mat_RM, np

maT = np.zeros((3,3))
matB = np.ones((3,4))
g=D_mat(.2,3, .2)
matC = np.matmul(maT, matB)
pass