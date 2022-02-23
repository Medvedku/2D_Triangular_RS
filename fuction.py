import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as plt3d
from matplotlib import cm

import surf2stl #import from local

def monkey(x,y,a):
    z = (x *(x**2 - 3*y**2))/a**2
    return z


def mobius(u,v,a,c=1):
    x = np.cos(v)*(a + c*u*np.cos(v/2))
    y = np.sin(v)*(a + c*u*np.cos(v/2))
    z = u*np.sin(v/2)
    return x, y, z


def astroidal(u, v, a=1, b=1, c=1):
    x = a*np.cos(u)**3 *np.cos(v)**3
    y = b*np.sin(u)**3 *np.cos(v)**3
    z = c*np.sin(v)**3 
    return x, y, z


def sea_shell(u, v, a=2, b=6*np.pi, c=3*np.pi):
    x = a * (np.e**(u/b) - 1) * np.cos(u) * np.cos(v/2)**2
    y = a * (1- np.e**(u/b)) * np.sin(u) * np.cos(v/2)**2
    z = np.e**(u/b) *np.sin(v) - np.e**(u/c) - np.sin(v) + 1
    return x, y, z

def breather(u, v, a=2/5):
    x = -u + ((2*(1-a**2)*np.cosh(a*u)*np.sinh(a*u)) / (a*((1-a**2)*np.cosh(a*u)**2+a**2*np.sin((1-a**2)**(1/2)*v)**2)))
    y = (2*(1-a**2)**(1/2)*np.cosh(a*u) * (-(1-a**2)**(1/2)*np.cos(v) *
         np.cos((1-a**2)**(1/2)*v) - np.sin(v)*np.sin((1-a**2)**(1/2)*v))) / (a*((1-a**2)*np.cosh(a*u)**2+a**2*np.sin((1-a**2)**(1/2)*v)**2))
    z = (2*(1-a**2)**(1/2)*np.cosh(a*u) * (-(1-a**2)**(1/2)*np.sin(v) *
         np.cos((1-a**2)**(1/2)*v) + np.cos(v)*np.sin((1-a**2)**(1/2)*v))) / (a*((1-a**2)*np.cosh(a*u)**2+a**2*np.sin((1-a**2)**(1/2)*v)**2))
    return x, y, z

mesh = 250

# <<<<<<<<<<<<<  MONKEY SLUT >>>>>>>>>>>>>>>>>>
# x = np.arange(-0.5, 0.5+1/mesh, 1/mesh)
# y = np.arange(-0.5, 0.5+1/mesh, 1/mesh)
# X, Y = np.meshgrid(x, y)
# Z = np.zeros(np.shape(X))

# for i in range(len(X)):
#     for j in range(len(X[0])):
#         Z[i,j] = monkey(X[i,j], Y[i,j], 1)

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# ax.set_xlabel("x")
# ax.set_ylabel("y")
# ax.set_zlabel("z")

# surf_def   = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
#                 linewidth=0, antialiased=True, alpha = 0.5 )

# plt.show()


# # <<<<<<<<<<<<< MOBIUM STRIPTEASE >>>>>>>>>>>>>>>>>>
# u = np.arange(-0.5, 0.5+1/mesh, 1/mesh)
# v = np.arange(0, 2*np.pi+2*np.pi/mesh, 2*np.pi/mesh)
# U, V = np.meshgrid(u, v)

# X = np.zeros(np.shape(U))
# Y = np.zeros(np.shape(U))
# Z = np.zeros(np.shape(U))
# acko = 1
# for i in range(len(U)):
#     for j in range(len(U[0])):
#         X[i,j], Y[i,j], Z[i,j] = mobius(U[i,j], V[i,j], acko)

# # <<<<<<<<<<<<< aSTEROIDS >>>>>>>>>>>>>>>>>>
# u = np.arange(-np.pi/2, np.pi/2+np.pi/mesh, np.pi/mesh)
# v = np.arange(-np.pi, np.pi+2*np.pi/mesh, 2*np.pi/mesh)
# U, V = np.meshgrid(u, v)

# X = np.zeros(np.shape(U))
# Y = np.zeros(np.shape(U))
# Z = np.zeros(np.shape(U))
# acko = 1
# for i in range(len(U)):
#     for j in range(len(U[0])):
#         X[i,j], Y[i,j], Z[i,j] = astroidal(U[i,j], V[i,j])

# <<<<<<<<<<<<< Sea Shell >>>>>>>>>>>>>>>>>>
# u = np.arange(0, 6*np.pi+6*np.pi/mesh, 6*np.pi/mesh)
# v = np.arange(0, 2*np.pi+2*np.pi/mesh, 2*np.pi/mesh)
# U, V = np.meshgrid(u, v)

# X = np.zeros(np.shape(U))
# Y = np.zeros(np.shape(U))
# Z = np.zeros(np.shape(U))
# acko = 1
# for i in range(len(U)):
#     for j in range(len(U[0])):
#         X[i,j], Y[i,j], Z[i,j] = sea_shell(U[i,j], V[i,j])

# <<<<<<<<<<<<< Breather >>>>>>>>>>>>>>>>>>
u = np.arange(-14, 14+28/mesh, 28/mesh)
v = np.arange(-38, 38+(2*38)/mesh, (2*38)/mesh)
U, V = np.meshgrid(u, v)

X = np.zeros(np.shape(U))
Y = np.zeros(np.shape(U))
Z = np.zeros(np.shape(U))
acko = 1
for i in range(len(U)):
    for j in range(len(U[0])):
        X[i, j], Y[i, j], Z[i, j] = breather(U[i, j], V[i, j])


fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

surf_def   = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                linewidth=0, antialiased=True, alpha = 0.5 )

# plt.show()

surf2stl.write('breather.stl', X, Y, Z)

pass