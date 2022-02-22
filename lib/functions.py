import numpy as np


def B_mat_RM(x: list, y: list, xg:float, yg:float):
    """Computes strain-displacement matrix (B matrix) for bending moments
    and shear forces.

        Arg:
            x:     X coordinates of the element  [list]
            y:     Y coordinates of the element  [list]
            xg:    Gauss point - local X coordinate [float]
            yg:    Gauss point - local X coordinate [float]
        Returns:
            b_mat_b:    strain-displacement matrix for bending moments
            b_mat_s:    strain-displacement matrix for shear forces
            area:       area of the element
    """

    N = np.zeros((3, 1))
    N[0] = 1.0 - xg - yg
    N[1] = xg
    N[2] = yg

    loc_dxN = np.zeros((3, 1))
    loc_dxN[0] = -1.0
    loc_dxN[1] = 1.0
    loc_dxN[2] = 0.0

    loc_dyN = np.zeros((3, 1))
    loc_dyN[0] = -1.0
    loc_dyN[1] = 0.0
    loc_dyN[2] = 1.0

    jac_m = np.zeros((2, 2))
    jac_m[0, 0] = x[0]*loc_dxN[0] + x[1]*loc_dxN[1] + x[2]*loc_dxN[2]
    jac_m[0, 1] = y[0]*loc_dxN[0] + y[1]*loc_dxN[1] + y[2]*loc_dxN[2]
    jac_m[1, 0] = x[0]*loc_dyN[0] + x[1]*loc_dyN[1] + x[2]*loc_dyN[2]
    jac_m[1, 1] = y[0]*loc_dyN[0] + y[1]*loc_dyN[1] + y[2]*loc_dyN[2]

    jac_im = np.linalg.inv(jac_m)

    area = abs(jac_m[0, 0]*jac_m[1, 1] - jac_m[1, 0]*jac_m[0, 1]) / 2

    Ndx = np.zeros((3, 1))
    Ndx[0] = jac_im[0, 0]*loc_dxN[0] + jac_im[0, 1]*loc_dyN[0]
    Ndx[1] = jac_im[0, 0]*loc_dxN[1] + jac_im[0, 1]*loc_dyN[1]
    Ndx[2] = jac_im[0, 0]*loc_dxN[2] + jac_im[0, 1]*loc_dyN[2]

    Ndy = np.zeros((3, 1))
    Ndy[0] = jac_im[1, 0]*loc_dxN[0] + jac_im[1, 1]*loc_dyN[0]
    Ndy[1] = jac_im[1, 0]*loc_dxN[1] + jac_im[1, 1]*loc_dyN[1]
    Ndy[2] = jac_im[1, 0]*loc_dxN[2] + jac_im[1, 1]*loc_dyN[2]

    bmat_b1 = np.array(
        [[0, -Ndx[0, 0],         0],
         [0,         0,  -Ndy[0, 0]],
         [0, -Ndy[0, 0], -Ndx[0, 0]]])
    bmat_b2 = np.array(
        [[0, -Ndx[1, 0],         0],
         [0,         0,  -Ndy[1, 0]],
         [0, -Ndy[1, 0], -Ndx[1, 0]]])
    bmat_b3 = np.array(
        [[0, -Ndx[2, 0],         0],
         [0,         0,  -Ndy[2, 0]],
         [0, -Ndy[2, 0], -Ndx[2, 0]]])
    bmat_b = np.concatenate((bmat_b1, bmat_b2, bmat_b3), axis=1)

    bmat_s1 = np.array(
        [[Ndx[0, 0], -N[0, 0],        0],
         [Ndy[0, 0],        0, -N[0, 0]]])
    bmat_s2 = np.array(
        [[Ndx[1, 0], -N[1, 0],        0],
         [Ndy[1, 0],        0, -N[1, 0]]])
    bmat_s3 = np.array(
        [[Ndx[2, 0], -N[2, 0],        0],
         [Ndy[2, 0],        0, -N[2, 0]]])
    bmat_s = np.concatenate((bmat_s1, bmat_s2, bmat_s3), axis=1)

    return bmat_b, bmat_s, area

def NodeStress(D_mat_b, D_mat_s, xg, yg, u):
    """Computes stresses at the Gauss points and evaluates values to the nodes.

        Arg:
            D_mat_b:    X coordinates of the element  [list]
            D_mat_s:    Y coordinates of the element  [list]
            xg:         Gauss point - local X coordinate [float]
            yg:         Gauss point - local X coordinate [float]
            u:          Nodal displacements
        Returns:
            ndStress:   Nodal stress matrix
    """

    pass


def D_mat(h, E, mi):
    """Computes constitutive matrices for bending
    moments and shear forces of element.

        Arg:
            h:          Thickness of element  [m, float]
            E:          Young modulus of elasticity  [Pa, int]
            mi:         Poissons number [float]
        Returns:
            D_mat_b:    Constitutive matrix for bending moments
            D_mat_s:    Constitutive matrix for shear forces
    """

    D_mat_b = (h**3 / 12) * (E / (1 - mi**2)) * np.array(
        [[1, mi,   0],
         [mi, 1,   0],
         [0,  0,   (1-mi)/2]])

    D_mat_s = (h*5 / 6) * (E / (2*(1 + mi))) * np.array(
        [[1,  0],
         [0,  1]])

    return D_mat_b, D_mat_s
