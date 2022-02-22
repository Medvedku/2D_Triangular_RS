import numpy as np
import pandas as pd

np.set_printoptions(precision=1)

spiel = 0

d_con = pd.read_csv('Concrete_data.csv')
d_reb = pd.read_csv('Rebar_data.csv')

gamma_s = 1.15
gamma_c = 1.5

## USED MATERIALS:
#  Concrete
concrete = "C20/25"
filt_c = (d_con["Class"]==concrete)
f_ck   = d_con.loc[filt_c, "f_ck"].values[0]
f_cd   = f_ck / gamma_c
f_ctm  = d_con.loc[filt_c, "f_ctm"].values[0]
E_cm   = d_con.loc[filt_c, "E_cm"].values[0]

#  Steel
steel  = "B500B"
filt_s = (d_reb["Class"]==steel)
f_yk   = d_reb.loc[filt_s, "f_yk"].values[0]
f_yd   = f_yk / gamma_s

def check_symmetric(a, rtol=1e-05, atol=1e-05):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)

class Node:
    num_of_nodes = 0
    def __init__(self, x, y):
        self.nd_id    = Node.num_of_nodes
        self.w        = Node.num_of_nodes*3 + 0
        self.fi_x     = Node.num_of_nodes*3 + 1
        self.fi_y     = Node.num_of_nodes*3 + 2
        self.co_x     = x
        self.co_y     = y
        self.def_z    = 0

        Node.num_of_nodes += 1

    def apply_deformation(self, defm):
        self.def_z    -= defm


class Element:
    num_of_elements = 0
    a  = None
    b  = None
    E  = None
    h  = None
    mi = None
    k  = None
    d_mat  = None
    b_mat  = None

    aa = 1 + 3**(1/2)
    bb = 1 - 3**(1/2)
    mstress = (1/4) * np.array([[aa*aa, aa*bb, bb*bb, aa*bb],
                                [aa*bb, aa*aa, aa*bb, bb*bb],
                                [bb*bb, aa*bb, aa*aa, aa*bb],
                                [aa*bb, bb*bb, aa*bb, aa*aa] ])

    def __init__(self, n0, n1, n2, n3, E = 200*1e9, mi = 0.3, h = 0.01):
        self.Flag   = False

        self.el_id  = Element.num_of_elements
        self.n_0    = n0.nd_id
        self.n_1    = n1.nd_id
        self.n_2    = n2.nd_id
        self.n_3    = n3.nd_id
        self.a      = abs(n0.co_x - n1.co_x)/2
        self.b      = abs(n1.co_y - n2.co_y)/2
        self.E      = E
        self.h      = h
        self.mi     = mi
        self.vec    = [ n0.w, n0.fi_x, n0.fi_y,
                        n1.w, n1.fi_x, n1.fi_y,
                        n2.w, n2.fi_x, n2.fi_y,
                        n3.w, n3.fi_x, n3.fi_y, ]

        if Element.num_of_elements == 0:
            Element.a  = self.a
            Element.b  = self.b
            Element.E  = self.E
            Element.h  = self.h
            Element.mi = self.mi
            if spiel:
                print("First element created", self.el_id)
            self.Flag   = True

        if abs((self.a / Element.a)-1)*100 > 1:
            if spiel:
                print("a changed, from: {} to: {}, element: {}".format( Element.a,self.a, self.el_id  ))
            Element.a  = self.a
            self.Flag   = True

        if abs((self.b / Element.b)-1)*100 > 1:
            if spiel:
                print("b changed, from: {} to: {}, element: {}".format( Element.b,self.b, self.el_id  ))
            Element.b  = self.b
            self.Flag   = True

        if abs((self.E / Element.E)-1)*100 > 1:
            if spiel:
                print("E changed, from: {} to: {}, element: {}".format( Element.E,self.E, self.el_id  ))
            Element.E  = self.E
            self.Flag   = True

        if abs((self.h / Element.h)-1)*100 > 1:
            if spiel:
                print("h changed, from: {} to: {}, element: {}".format( Element.h,self.h, self.el_id  ))
            Element.h  = self.h
            self.Flag   = True

        if abs((self.mi / Element.mi)-1)*100 > 1:
            if spiel:
                print("mi changed, from: {} to: {}, element: {}".format( Element.mi,self.mi, self.el_id  ))
            Element.mi  = self.mi
            self.Flag   = True

        if self.Flag:
            Element.k, Element.b_mat, Element.d_mat = self.k_elem()
            self.k_el  = Element.k
            self.b_mat = Element.b_mat
            self.d_mat = Element.d_mat
        else:
            self.k_el  = Element.k
            self.b_mat = Element.b_mat
            self.d_mat = Element.d_mat

        self.e_loads = []

        Element.num_of_elements += 1

    def k_elem(self):
        a = self.a
        b = self.b
        E = self.E
        h = self.h
        mi = self.mi

        d_mat = np.zeros( (12,12) )
        D_con = E*h**3/(12*(1-mi**2))
        d_mat[0][0], d_mat[1][1], d_mat[3][3], d_mat[4][4]   = D_con * 1, D_con * 1, D_con * 1, D_con * 1
        d_mat[6][6], d_mat[7][7], d_mat[9][9], d_mat[10][10] = D_con * 1, D_con * 1, D_con * 1, D_con * 1

        d_mat[0][1], d_mat[1][0], d_mat[3][4], d_mat[4][3]   = D_con * mi, D_con * mi, D_con * mi, D_con * mi
        d_mat[6][7], d_mat[7][6], d_mat[9][10], d_mat[10][9] = D_con * mi, D_con * mi, D_con * mi, D_con * mi

        d_mat[2][2], d_mat[5][5], d_mat[8][8], d_mat[11][11] = D_con * (1-mi)/2, D_con * (1-mi)/2, D_con * (1-mi)/2, D_con * (1-mi)/2

        gauss_x=[-1/3**0.5, +1/3**0.5, +1/3**0.5, -1/3**0.5]
        gauss_y=[-1/3**0.5, -1/3**0.5, +1/3**0.5, +1/3**0.5]

        b_mat = None
        b_matrix = []
        for i in range(4):
            x_g = gauss_x[i]
            y_g = gauss_y[i]
            b_mat_b = self.B_matrix(a, b, x_g, y_g)
            if i == 0:
                b_mat = b_mat_b
                b_matrix.append(b_mat_b)
            else:
                b_mat = np.concatenate( ( b_mat, b_mat_b) )
                b_matrix.append(b_mat_b)

        K_elem = (np.transpose(b_mat)@d_mat@b_mat)*a*b

        return K_elem, b_matrix, d_mat

    def B_matrix(self, a, b, x, y):
        b_matrix = np.zeros( (3, 12) )
        b_matrix[0][0]  = -3*(x-x*y)/(4*a**2)
        b_matrix[0][1]  = -((3*a*x-3*a*x*y-a+a*y)/4)/a**2
        b_matrix[0][2]  = 0
        b_matrix[0][3]  = -3*(-x+x*y)/(4*a**2)
        b_matrix[0][4]  = -((3*a*x-3*a*x*y+a-a*y)/4)/a**2
        b_matrix[0][5]  = 0
        b_matrix[0][6]  = -3*(-x-x*y)/(4*a**2)
        b_matrix[0][7]  = -((3*a*x+3*a*x*y+a+a*y)/4)/a**2
        b_matrix[0][8]  = 0
        b_matrix[0][9]  = -3*(x+x*y)/(4*a**2)
        b_matrix[0][10] = -((3*a*x+3*a*x*y-a-a*y)/4)/a**2
        b_matrix[0][11] = 0
        b_matrix[1][0]  = -3*(y-x*y)/(4*b**2)
        b_matrix[1][1]  = 0
        b_matrix[1][2]  = -((3*b*y-3*b*x*y-b+b*x)/4)/b**2
        b_matrix[1][3]  = -3*(y+x*y)/(4*b**2)
        b_matrix[1][4]  = 0
        b_matrix[1][5]  = -((3*b*y+3*b*x*y-b-b*x)/4)/b**2
        b_matrix[1][6]  = -3*(-y-x*y)/(4*b**2)
        b_matrix[1][7]  = 0
        b_matrix[1][8]  = -((3*b*y+3*b*x*y+b+b*x)/4)/b**2
        b_matrix[1][9]  = -3*(-y+x*y)/(4*b**2)
        b_matrix[1][10] = 0
        b_matrix[1][11] = -((3*b*y-3*b*x*y+b-b*x)/4)/b**2
        b_matrix[2][0]  = -2*(1/2-3*x**2/8-3*y**2/8)/(a*b)
        b_matrix[2][1]  = -2*(-3/8*a*x**2+a*x/4+a/8)/(a*b)
        b_matrix[2][2]  = -2*(-3/8*b*y**2+b*y/4+b/8)/(a*b)
        b_matrix[2][3]  = -2*(-1/2+3*x**2/8+3*y**2/8)/(a*b)
        b_matrix[2][4]  = -2*(-3/8*a*x**2-a*x/4+a/8)/(a*b)
        b_matrix[2][5]  = -2*(3/8*b*y**2-b*y/4-b/8)/(a*b)
        b_matrix[2][6]  = -2*(1/2-3*x**2/8-3*y**2/8)/(a*b)
        b_matrix[2][7]  = -2*(3/8*a*x**2+a*x/4-a/8)/(a*b)
        b_matrix[2][8]  = -2*(3/8*b*y**2+b*y/4-b/8)/(a*b)
        b_matrix[2][9]  = -2*(-1/2+3*x**2/8+3*y**2/8)/(a*b)
        b_matrix[2][10] = -2*(3/8*a*x**2-a*x/4-a/8)/(a*b)
        b_matrix[2][11] = -2*(-3/8*b*y**2-b*y/4+b/8)/(a*b)
        return b_matrix

    def get_load_vector(self, q):
        """Computes load vector from surface load for each element. Can be called multiple times
            Arg:
                q:       Value of surface load in [N/m^2].    [constant]
        """
        a = self.a
        b = self.b
        q = 4*q*a*b * np.array([1/4, a/12, b/12, 1/4, -a/12, b/12, 1/4, -a/12, -b/12, 1/4, a/12, -b/12])
        self.e_loads.append(q)

    def get_internal_forces(self, c_n, r_t):
        """Computes internal forces (bending moments) at the centre of element.
        Also computes design internal forces.

            Arg:
                c_n:       Reduced list of code numbers.    [list]
                r_t:       Reduced vector of deformations.  [list]
        """
        a = []
        for i in range(12):
            if (self.vec[i] in c_n):
                indx = list(c_n).index(self.vec[i])
                a.append(r_t[indx])
            else:
                a.append(0)
        a = np.array(a)

        d_mtrx =    self.E*self.h**3/(12*(1-self.mi**2)) * np.array([
                    [1,       self.mi, 0             ],
                    [self.mi, 1,       0             ],
                    [0,       0,       (1-self.mi)/2 ] ] )

        stress = [ 0, 0, 0 ]
        for i in range(4):
            stress += d_mtrx @ self.b_mat[i] @ a

        self.moments = stress/4

        self.m_x  = self.moments[0]
        self.m_y  = self.moments[1]
        self.m_xy = -self.moments[2]

        m_x = self.m_x
        m_y = self.m_y
        m_xy = self.m_xy

        if m_x >= -abs(m_xy):
            m_x_bot = m_x + abs(m_xy)
            m_y_bot = m_y + abs(m_xy)
            self.m_c_bot = - 2 * abs(m_xy)
        else:
            m_x_bot = 0
            # m_y_bot = m_y + (m_xy**2)/(abs(m_x))
            m_y_bot = 0
            self.m_c_bot = -abs(m_x)*(1+(m_xy/m_x)**2)

        if m_y <= abs(m_xy):
            m_x_top = m_x - abs(m_xy)
            m_y_top = m_y - abs(m_xy)
            self.m_c_top = -2 * abs(m_xy)
        else:
            # m_x_top = -m_x + (m_xy**2)/(abs(m_y))
            m_x_top = 0
            m_y_top = 0
            self.m_c_top = -abs(m_y)*(1+((m_xy/m_y)**2))

        if m_x_bot > 0:
            self.M_x_bot = m_x_bot
        else:
            self.M_x_bot = 0
        if m_y_bot > 0:
            self.M_y_bot = m_y_bot
        else:
            self.M_y_bot = 0
        if m_x_top < 0:
            self.M_x_top = m_x_top
        else:
            self.M_x_top = 0
        if m_y_top < 0:
            self.M_y_top = m_y_top
        else:
            self.M_y_top = 0

    def get_bot_rebar(self, layer_b0="X", d_b0=10, d_b1=12, c_n=None):
        """Takes argument layer_b0 as a first layer of rebar at the bottom of slab.

            Arg:
                layer_b0:   Direction of the first layer of rebars placed at the bottom of the slab.
                            Input must be string X or Y
                d_b0:       Diameter of the rebar - first layer
                            Input in milimeters [mm]
                d_b1:       Diameter of the rebar - second layer
                            Input in milimeters [mm]

            Optional Arg:
                c_n:      Concrete reinforcement cover layer
                            Direct input in milimeters or cumputed based on EN if left void.
        """
        b_0 = 1
        b_1 = 1
        if layer_b0 == "X":
            # b_0 = 2*self.b
            # b_1 = 2*self.a
            self.bot_layer_0 = layer_b0
            self.bot_layer_1 = "Y"
            self.M_Ed_0_b = self.M_x_bot
            self.M_Ed_1_b = self.M_y_bot
        elif layer_b0 == "Y":
            # b_0 = 2*self.a
            # b_1 = 2*self.b
            self.bot_layer_0 = layer_b0
            self.bot_layer_1 = "X"
            self.M_Ed_0_b = self.M_y_bot
            self.M_Ed_1_b = self.M_x_bot
        else:
            print("Non valid direction of botom first layer -> Aborting computation")
            quit()
        self.d_b_0 = d_b0
        self.d_b_1 = d_b1
        if c_n:
            self.c_b = c_n / 1000
        else:
            if self.d_b_0%5 == 0:
                self.c_b = ( ((self.d_b_0//5+0)*5) + 10) / 1000
            else:
                self.c_b = ( ((self.d_b_0//5+1)*5) + 10) / 1000
        h_ = self.h
        d_0 = h_ - (self.c_b + 0.5*self.d_b_0*0.001)
        As_min_0 = 0.26 * (f_ctm/f_yk) * b_0 * d_0
        As_min_1 = 0.0013 * b_0 * d_0
        self.As_min_b_0 = max(As_min_0, As_min_1)

        # First layer
        x_lim_0 = (700*d_0)/(700+f_yd)
        if self.M_Ed_0_b > 0:
            x_b_0 = d_0 - (d_0**2 - ( self.M_Ed_0_b / (0.5*b_0*f_cd*1e6) ))**(1/2)
            x_0 = x_b_0/0.8
            if x_0 >= x_lim_0:
                print("x is greater then x_lim")
                self.As_b_0 = 0.0099
            else:
                self.As_req_b_0 = (x_b_0 * b_0 * f_cd) / f_yd
                self.As_b_0 = max( self.As_min_b_0 , self.As_req_b_0 )
        else:
            self.As_b_0 = self.As_min_b_0

        # Second layer
        d_1 = h_ - (self.c_b + (0.5*self.d_b_0+self.d_b_1)*0.001)
        As_min_0 = 0.26 * (f_ctm/f_yk) * b_1 * d_1
        As_min_1 = 0.0013 * b_1 * d_1
        self.As_min_b_1 = max(As_min_0, As_min_1)
        x_lim_1 = (700*d_1)/(700+f_yd)
        if self.M_Ed_1_b > 0:
            x_b_1 = d_1 - (d_1**2 - ( self.M_Ed_1_b / (0.5*b_1*f_cd*1e6) ))**(1/2)
            x_1 = x_b_1/0.8
            if x_1 >= x_lim_1:
                print("x is greater then x_lim")
                self.As_b_1 = 0.0099
            else:
                self.As_req_b_1 = (x_b_1 * b_1 * f_cd) / f_yd
                self.As_b_1 = max( self.As_min_b_1 , self.As_req_b_1 )
        else:
            self.As_b_1 = self.As_min_b_1


    def get_top_rebar(self, layer_t0="X", d_t0=10, d_t1=12, c_n=None):
        """Takes argument layer_b0 as a first layer of rebar at the top of slab.

            Arg:
                layer_t0:   Direction of the first layer of rebars placed at the top of the slab.
                            Input must be string X or Y
                d_t0:       Diameter of the rebar - first layer
                            Input in milimeters [mm]
                d_t1:       Diameter of the rebar - second layer
                            Input in milimeters [mm]

            Optional Arg:
                c_n:      Concrete reinforcement cover layer
                            Direct input in milimeters or cumputed based on EN if left void.
        """
        b_0 = 1
        b_1 = 1
        if layer_t0 == "X":
            # b_0 = 2*self.b
            # b_1 = 2*self.a
            self.top_layer_0 = layer_t0
            self.top_layer_1 = "Y"
            self.M_Ed_0_t = self.M_x_top
            self.M_Ed_1_t = self.M_y_top
        elif layer_t0 == "Y":
            # b_0 = 2*self.a
            # b_1 = 2*self.b
            self.top_layer_0 = layer_t0
            self.top_layer_1 = "X"
            self.M_Ed_0_t = self.M_y_top
            self.M_Ed_1_t = self.M_x_top
        else:
            print("Non valid direction of botom first layer -> Aborting computation")
            quit()
        self.d_t_0 = d_t0
        self.d_t_1 = d_t1
        if c_n:
            self.c_t = c_n / 1000
        else:
            if self.d_t_0%5 == 0:
                self.c_t = ( ((self.d_t_0//5+0)*5) + 10) / 1000
            else:
                self.c_t = ( ((self.d_t_0//5+1)*5) + 10) / 1000
        h_ = self.h
        d_0 = h_ - (self.c_t + 0.5*self.d_t_0*0.001)
        As_min_0 = 0.26 * (f_ctm/f_yk) * b_0 * d_0
        As_min_1 = 0.0013 * b_0 * d_0
        self.As_min_t_0 = max(As_min_0, As_min_1)

        # First layer
        x_lim_0 = (700*d_0)/(700+f_yd)
        if self.M_Ed_0_t < 0:
            x_t_0 = d_0 - (d_0**2 - ( abs(self.M_Ed_0_t) / (0.5*b_0*f_cd*1e6) ))**(1/2)
            x_0 = x_t_0/0.8
            if x_0 >= x_lim_0:
                print("x is greater then x_lim")
                self.As_t_1 = 0.0099
            else:
                self.As_req_t_0 = (x_t_0 * b_0 * f_cd) / f_yd
                self.As_t_0 = max( self.As_min_t_0 , self.As_req_t_0 )
        else:
            self.As_t_0 = 0

        # Second layer
        d_1 = h_ - (self.c_t + (0.5*self.d_t_0+self.d_t_1)*0.001)
        As_min_0 = 0.26 * (f_ctm/f_yk) * b_1 * d_1
        As_min_1 = 0.0013 * b_1 * d_1
        self.As_min_t_1 = max(As_min_0, As_min_1)
        x_lim_1 = (700*d_1)/(700+f_yd)
        if self.M_Ed_1_t < 0:
            x_t_1 = d_1 - (d_1**2 - ( abs(self.M_Ed_1_t) / (0.5*b_1*f_cd*1e6) ))**(1/2)
            x_1 = x_t_1/0.8
            if x_1 >= x_lim_1:
                print("x is greater then x_lim")
                self.As_t_1 = 0.0099
            else:
                self.As_req_t_1 = (x_t_1 * b_1 * f_cd) / f_yd
                self.As_t_1 = max( self.As_min_t_1 , self.As_req_t_1 )
        else:
            self.As_t_1 = 0


# Definition of construction

LX = 5
LY = 9


mesh = .5

nx = int(LX/mesh)
ny = int(LY/mesh)


coor_x = [0+i*(LX/nx) for i in range(nx+1)]
coor_y = [0+i*(LY/ny) for i in range(ny+1)]

Nodes_ = np.empty( (ny+1, nx+1), dtype = "int" )
cunt = 0
for i in range(ny+1):
    for j in range(nx+1):
        Nodes_[i][j] = cunt
        cunt += 1

l_nodes = []
for i in range(ny+1):
    for j in range(nx+1):
        x_coor = coor_x[j]
        y_coor = coor_y[i]
        l_nodes.append(Node(x_coor,y_coor))

l_elems = []
for i in range(ny):
    for j in range(nx):
        n0 = Nodes_[i+1][j]
        n1 = Nodes_[i+1][j+1]
        n2 = Nodes_[i][j+1]
        n3 = Nodes_[i][j]
        l_elems.append( Element(l_nodes[n0], l_nodes[n1], l_nodes[n2], l_nodes[n3], h = 0.200, E = E_cm*1e9, mi=0.2) )

## LOAD CASES

# for i in l_elems:
#     i.get_load_vector(2000)

for i in l_elems:
    i.get_load_vector(20*1e3)

# Assembly of Global stiffness B_matrix
n_e       = len(l_elems)
n_n       = len(l_nodes)
dofs      = n_n*3
code_nums = list(range(dofs))

load = [0 for i in range(dofs)]
for e in l_elems:
    for l in e.e_loads:
        for i in range(12):             # because element matrix is 12x12 shape
            load[e.vec[i]] += l[i]

# # 4 points around
# boundary = []
# for i in l_nodes:
#     if i.co_x == 0 and i.co_y == LY:
#         boundary.append(i.w)
#     if i.co_x == LX and i.co_y == 0:
#         boundary.append(i.w)
#     if i.co_x == LX and i.co_y == LY:
#         boundary.append(i.w)
#     if i.co_x == 0 and i.co_y == 0:
#         boundary.append(i.w)
# deleto = boundary

# Sample
boundary = []
for i in l_nodes:
    if i.co_x == 0 and i.co_y >= 0.5*LY:
        boundary.append(i.w)
    if i.co_x == LX and i.co_y == 0:
        boundary.append(i.w)
    if i.co_x == LX and i.co_y == LY:
        boundary.append(i.w)
    if i.co_x == 0 and i.co_y == 0:
        boundary.append(i.w)
    if i.co_y == LY/2 and i.co_x <= LX/2:
        boundary.append(i.w)
    if i.co_x == LX and i.co_y <= 0.5*LY:
        boundary.append(i.w)
deleto = boundary

# # 4 Edges around
# boundary = []
# for i in l_nodes:
#     if i.co_x == 0:
#         boundary.append(i.w)
#     # if i.co_x == LX:
#     #     boundary.append(i.w)
#     if i.co_y == 0:
#         boundary.append(i.w)
#     if i.co_y == LY:
#         boundary.append(i.w)
#     if i.co_y == (5/9)*LY:
#         boundary.append(i.w)
# deleto = boundary


K_gl = np.zeros( (dofs, dofs) )
for i in l_elems:
    for j in range(12):                 # because element matrix is 12x12 shape
        for k in range(12):             # because element matrix is 12x12 shape
            K_gl[i.vec[j],i.vec[k]] += i.k_el[j,k]

# np.savetxt("2elems_full.csv",K_gl, fmt = "%.2e")

K_gl = np.delete(K_gl, deleto, axis = 0)
K_gl = np.delete(K_gl, deleto, axis = 1)
load = np.delete(load, deleto, axis = 0)

code_nums = np.delete(code_nums, deleto, axis = 0)

delta = np.linalg.inv(K_gl)
r_tot = np.matmul(delta, load)

## APPLICATION OF DEFORMATIONS
for i in range(len(r_tot)):
    if code_nums[i]%3 == 0:
        l_nodes[code_nums[i]//3].apply_deformation(r_tot[i])

## GET INTERNAL FORCES ON ELEMENTS
for i in l_elems:
    i.get_internal_forces(code_nums, r_tot)

for i in l_elems:
    i.get_bot_rebar(layer_b0="X", d_b0=8, d_b1=8, c_n=None)

for i in l_elems:
    i.get_top_rebar(layer_t0="X", d_t0=8, d_t1=8, c_n=None)


#
# print(l_elems[0].moments)

x = np.arange(0+0.5*mesh, LX, mesh)
y = np.arange(0+0.5*mesh, LY, mesh)
X, Y = np.meshgrid(x, y)
Z = np.zeros(np.shape(X))
A_b_0 = np.zeros(np.shape(X))
A_b_1 = np.zeros(np.shape(X))
A_t_0 = np.zeros(np.shape(X))
A_t_1 = np.zeros(np.shape(X))

M_b_0 = np.zeros(np.shape(X))
M_b_1 = np.zeros(np.shape(X))
M_t_0 = np.zeros(np.shape(X))
M_t_1 = np.zeros(np.shape(X))

mx  = np.zeros(np.shape(X))
my  = np.zeros(np.shape(X))
mxy = np.zeros(np.shape(X))

cunt = 0
for i in range(len(Z)):
    for j in range(len(Z[0])):
        Z[i,j] = int(l_elems[cunt].moments[0])
        A_t_0[i,j] = l_elems[cunt].As_t_0 * 1e4
        A_t_1[i,j] = l_elems[cunt].As_t_1 * 1e4
        A_b_0[i,j] = l_elems[cunt].As_b_0 * 1e4
        A_b_1[i,j] = l_elems[cunt].As_b_1 * 1e4
        M_t_0[i,j] = l_elems[cunt].M_x_top / 1000
        M_t_1[i,j] = l_elems[cunt].M_y_top / 1000
        M_b_0[i,j] = l_elems[cunt].M_x_bot / 1000
        M_b_1[i,j] = l_elems[cunt].M_y_bot / 1000

        mx[i,j]  = l_elems[cunt].m_x
        my[i,j]  = l_elems[cunt].m_y
        mxy[i,j] = l_elems[cunt].m_xy
        cunt +=1

bot_0 = l_elems[0].bot_layer_0
bot_1 = l_elems[0].bot_layer_1
top_0 = l_elems[0].top_layer_0
top_1 = l_elems[0].top_layer_1

# print("Top_0", top_0)
# print(A_t_0)
# print("Top_1", top_1)
# print(A_t_1)
#
# print("Bot_0", bot_0)
# print(A_b_0)
# print("Bot_1", bot_1)
# print(A_b_1)
#
# print("M_Ed Top_0", top_0)
# print(M_t_0)
# print("M_Ed Top_1", top_1)
# print(M_t_1)
#
# print("M_Ed Bot_0", top_0)
# print(M_b_0)
# print("M_Ed Bot_1", top_1)
# print(M_b_1)
#
# print("M_x")
# print(mx*0.001)
# print("M_y")
# print(my*0.001)
# print("M_xy")
# print(mxy*0.001)

As_min_b_0 = l_elems[0].As_min_b_0
As_min_b_1 = l_elems[0].As_min_b_1
As_min_t_0 = l_elems[0].As_min_t_0
As_min_t_1 = l_elems[0].As_min_t_1

s_b_0_req = 1 / ((As_min_b_0*1e4) / (((l_elems[0].d_b_0**2)*np.pi/4)*1e-2))
s_b_1_req = 1 / ((As_min_b_1*1e4) / (((l_elems[0].d_b_1**2)*np.pi/4)*1e-2))
s_t_0_req = 1 / ((As_min_t_0*1e4) / (((l_elems[0].d_t_0**2)*np.pi/4)*1e-2))
s_t_1_req = 1 / ((As_min_t_1*1e4) / (((l_elems[0].d_t_1**2)*np.pi/4)*1e-2))

s_b_0_prov = int(s_b_0_req*100)/100
s_b_1_prov = int(s_b_1_req*100)/100
s_t_0_prov = int(s_t_0_req*100)/100
s_t_1_prov = int(s_t_1_req*100)/100

As_0_b_0 = (1/s_b_0_prov) * (((l_elems[0].d_b_0**2)*np.pi/4)*1e-2)
As_0_b_1 = (1/s_b_1_prov) * (((l_elems[0].d_b_1**2)*np.pi/4)*1e-2)
As_0_t_0 = (1/s_t_0_prov) * (((l_elems[0].d_t_0**2)*np.pi/4)*1e-2)
As_0_t_1 = (1/s_t_1_prov) * (((l_elems[0].d_t_1**2)*np.pi/4)*1e-2)


l_A_b_0 = [0, As_0_b_0]
l_A_b_1 = [0, As_0_b_1]
l_A_t_0 = [0, As_0_t_0]
l_A_t_1 = [0, As_0_t_1]

cunt = 0
while l_A_b_0[-1] < np.max(A_b_0):
    l_A_b_0.append( l_A_b_0[1] + (cunt+1) * (((l_elems[0].d_b_0**2)*np.pi/4)*1e-2) )
    cunt += 1

cunt = 0
while l_A_b_1[-1] < np.max(A_b_1):
    l_A_b_1.append( l_A_b_1[1] + (cunt+1) * (((l_elems[0].d_b_1**2)*np.pi/4)*1e-2) )
    cunt += 1

cunt = 0
while l_A_t_0[-1] < np.max(A_t_0):
    l_A_t_0.append( l_A_t_0[1] + (cunt+1) * (((l_elems[0].d_t_0**2)*np.pi/4)*1e-2) )
    cunt += 1

cunt = 0
while l_A_t_1[-1] < np.max(A_t_1):
    l_A_t_1.append( l_A_t_1[1] + (cunt+1) * (((l_elems[0].d_t_1**2)*np.pi/4)*1e-2) )
    cunt += 1


# print(l_A_b_0)
# print(l_A_b_1)
# print(l_A_t_0)
# print(l_A_t_1)


red_r_tot  = []
red_c_nums = []
for i in range(len(code_nums)):
    if code_nums[i]%3 == 0:
        red_r_tot.append(r_tot[i])
        red_c_nums.append(code_nums[i])

r_max = max(red_r_tot)
r_min = min(red_r_tot)

c_n_max = red_c_nums[red_r_tot.index(r_max)]
c_n_min = red_c_nums[red_r_tot.index(r_min)]


_3D = 0
if _3D:
    scale = 1
    # PLOT
    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d as plt3d
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    red_r_tot  = []
    red_c_nums = []
    for i in range(len(code_nums)):
        if code_nums[i]%3 == 0:
            red_r_tot.append(r_tot[i])
            red_c_nums.append(code_nums[i])

    r_max = max(red_r_tot)
    r_min = min(red_r_tot)
    c_n_max = red_c_nums[red_r_tot.index(r_max)]
    c_n_min = red_c_nums[red_r_tot.index(r_min)]

    if LX > LY:
        ax.set_xlim( (0, LX) )
        ax.set_ylim( (-abs(LX-LY)/2, LY+(abs(LX-LY)/2)) )
        ax.set_zlim( (min(-scale*r_max,0), max(-scale*r_min,0)) )
    elif LX < LY:
        ax.set_xlim( (-abs(LX-LY)/2, LX+(abs(LX-LY)/2)) )
        ax.set_ylim( (0,LY) )
        ax.set_zlim( (min(-scale*r_max,0), max(-scale*r_min,0)) )
    else:
        ax.set_xlim( (0,LX) )
        ax.set_ylim( (0,LY) )
        ax.set_zlim( (min(-scale*r_max,0), max(-scale*r_min,0)) )

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    # Remove gray panes and axis grid
    ax.xaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('white')
    ax.yaxis.pane.fill = False
    ax.yaxis.pane.set_edgecolor('white')
    ax.zaxis.pane.fill = False
    ax.zaxis.pane.set_edgecolor('white')
    ax.grid(False)
    # Remove z-axis
    ax.w_zaxis.line.set_lw(0.)
    ax.set_zticks([])

    X = np.zeros( (ny+1, nx+1) )
    Y = np.zeros( (ny+1, nx+1) )
    Z = np.zeros( (ny+1, nx+1) )
    Z_= np.zeros( (ny+1, nx+1) )

    for i in range(ny+1):
        for j in range(nx+1):
            n = Nodes_[i,j]
            X[i,j] = l_nodes[n].co_x
            Y[i,j] = l_nodes[n].co_y
            Z[i,j] = 0
            Z_[i,j]= l_nodes[n].def_z

    surf_undef = ax.plot_surface(X, Y, Z,
                 linewidth=1, antialiased=True, alpha = 0.1 )

    surf_def   = ax.plot_surface(X, Y, Z_, cmap=cm.coolwarm,
                 linewidth=0, antialiased=True, alpha = 0.5 )

    for i in deleto:
        ax.plot([l_nodes[i//3].co_x], [l_nodes[i//3].co_y], [0.], markerfacecolor='k', markeredgecolor='k', marker='o', markersize=1, alpha=1.)

    ax.plot([l_nodes[c_n_max//3].co_x], [l_nodes[c_n_max//3].co_y], [l_nodes[c_n_max//3].def_z], markeredgecolor='b', marker='x', markersize=7, alpha=0.7)
    ax.plot([l_nodes[c_n_min//3].co_x], [l_nodes[c_n_min//3].co_y], [l_nodes[c_n_min//3].def_z], markeredgecolor='r', marker='x', markersize=7, alpha=0.7)

    plt.show()

_2D = 1
if _2D:
    import matplotlib.pyplot as plt
    from matplotlib import gridspec as grd
    # plt.style.use("seaborn")
    scale = 11.5/max(LX,LY)/2
    plt.figure(figsize=(2*LX*scale,2*LY*scale))
    G = grd.GridSpec(2,2)

    X = [0, LX, LX, 0 , 0]
    Y = [0, 0,  LY, LY, 0]
    # fig0, ax0 = plt.subplots(figsize=(LX, LY))
    ax0 = plt.subplot(G[0, 0])
    ax0.plot(X, Y, linewidth=0.5)

    # fig1, ax1 = plt.subplots(figsize=(LX, LY))
    ax1 = plt.subplot(G[0, 1])
    ax1.plot(X, Y, linewidth=0.5)

    # fig2, ax2 = plt.subplots(figsize=(LX, LY))
    ax2 = plt.subplot(G[1, 0])
    ax2.plot(X, Y, linewidth=0.5)

    # fig3, ax3 = plt.subplots(figsize=(LX, LY))
    ax3 = plt.subplot(G[1, 1])
    ax3.plot(X, Y, linewidth=0.5)

    for i in deleto:
        ax0.plot([l_nodes[i//3].co_x], [l_nodes[i//3].co_y], markerfacecolor='k', markeredgecolor='k', marker='o', markersize=1, alpha=1.)
        ax1.plot([l_nodes[i//3].co_x], [l_nodes[i//3].co_y], markerfacecolor='k', markeredgecolor='k', marker='o', markersize=1, alpha=1.)
        ax2.plot([l_nodes[i//3].co_x], [l_nodes[i//3].co_y], markerfacecolor='k', markeredgecolor='k', marker='o', markersize=1, alpha=1.)
        ax3.plot([l_nodes[i//3].co_x], [l_nodes[i//3].co_y], markerfacecolor='k', markeredgecolor='k', marker='o', markersize=1, alpha=1.)

    x = np.arange(0+0.5*mesh, LX, mesh)
    y = np.arange(0+0.5*mesh, LY, mesh)
    X, Y = np.meshgrid(x, y)
    Z0 = np.zeros(np.shape(X))
    Z1 = np.zeros(np.shape(X))
    Z2 = np.zeros(np.shape(X))
    Z3 = np.zeros(np.shape(X))

    cunt = 0
    for i in range(len(Z)):
        for j in range(len(Z[0])):
            # Z[i,j] = l_elems[cunt].moments[2]
            Z0[i,j] = l_elems[cunt].As_b_0 * 1e4
            Z1[i,j] = l_elems[cunt].As_b_1 * 1e4
            Z2[i,j] = l_elems[cunt].As_t_0 * 1e4
            Z3[i,j] = l_elems[cunt].As_t_1 * 1e4
            cunt +=1

    ax0.contour(X, Y, Z0, levels=l_A_b_0)
    ax1.contour(X, Y, Z1, levels=l_A_b_1)
    ax2.contour(X, Y, Z2, levels=l_A_t_0)
    ax3.contour(X, Y, Z3, levels=l_A_t_1)
    # ax.clabel(CS, inline=True, fontsize=5)

    bot_0 = l_elems[0].bot_layer_0
    bot_1 = l_elems[0].bot_layer_1
    top_0 = l_elems[0].top_layer_0
    top_1 = l_elems[0].top_layer_1

    ax0.set_title('As bottom - first layer '+ l_elems[0].bot_layer_0)
    ax1.set_title('As bottom - second layer '+ l_elems[0].bot_layer_1)
    ax2.set_title('As top - first layer '+ l_elems[0].top_layer_0)
    ax3.set_title('As top - second layer '+ l_elems[0].top_layer_1)
    # ax0.set_xlabel('X-axis')
    # ax0.set_ylabel('Y-axis')

    # plt.savefig("As2.eps")
    plt.show()


_2D_M = 0
if _2D_M:
    import matplotlib.pyplot as plt
    from matplotlib import gridspec as grd
    # plt.style.use("seaborn")
    scale = 11.5/max(LX,LY)/2
    plt.figure(figsize=(2*LX*scale,2*LY*scale))
    G = grd.GridSpec(2,2)

    X = [0, LX, LX, 0 , 0]
    Y = [0, 0,  LY, LY, 0]
    # fig0, ax0 = plt.subplots(figsize=(LX, LY))
    ax0 = plt.subplot(G[0, 0])
    ax0.plot(X, Y, linewidth=0.5)

    # fig1, ax1 = plt.subplots(figsize=(LX, LY))
    ax1 = plt.subplot(G[0, 1])
    ax1.plot(X, Y, linewidth=0.5)

    # fig2, ax2 = plt.subplots(figsize=(LX, LY))
    ax2 = plt.subplot(G[1, 0])
    ax2.plot(X, Y, linewidth=0.5)

    # fig3, ax3 = plt.subplots(figsize=(LX, LY))
    ax3 = plt.subplot(G[1, 1])
    ax3.plot(X, Y, linewidth=0.5)

    for i in deleto:
        ax0.plot([l_nodes[i//3].co_x], [l_nodes[i//3].co_y], markerfacecolor='k', markeredgecolor='k', marker='o', markersize=1, alpha=1.)
        ax1.plot([l_nodes[i//3].co_x], [l_nodes[i//3].co_y], markerfacecolor='k', markeredgecolor='k', marker='o', markersize=1, alpha=1.)
        ax2.plot([l_nodes[i//3].co_x], [l_nodes[i//3].co_y], markerfacecolor='k', markeredgecolor='k', marker='o', markersize=1, alpha=1.)
        ax3.plot([l_nodes[i//3].co_x], [l_nodes[i//3].co_y], markerfacecolor='k', markeredgecolor='k', marker='o', markersize=1, alpha=1.)

    x = np.arange(0+0.5*mesh, LX, mesh)
    y = np.arange(0+0.5*mesh, LY, mesh)
    X, Y = np.meshgrid(x, y)
    Z0 = np.zeros(np.shape(X))
    Z1 = np.zeros(np.shape(X))
    Z2 = np.zeros(np.shape(X))
    Z3 = np.zeros(np.shape(X))

    cunt = 0
    for i in range(len(Z)):
        for j in range(len(Z[0])):
            # Z[i,j] = l_elems[cunt].moments[2]
            Z0[i,j] = l_elems[cunt].M_x_bot
            Z1[i,j] = l_elems[cunt].M_y_bot
            Z2[i,j] = l_elems[cunt].M_x_top
            Z3[i,j] = l_elems[cunt].M_y_top
            cunt +=1

    ax0.contour(X, Y, Z0)
    ax1.contour(X, Y, Z1)
    ax2.contour(X, Y, Z2)
    ax3.contour(X, Y, Z3)
    # ax.clabel(CS, inline=True, fontsize=5)

    bot_0 = l_elems[0].bot_layer_0
    bot_1 = l_elems[0].bot_layer_1
    top_0 = l_elems[0].top_layer_0
    top_1 = l_elems[0].top_layer_1

    ax0.set_title("M_x_bot")
    ax1.set_title("M_y_bot")
    ax2.set_title("M_x_top")
    ax3.set_title("M_y_top")
    # ax0.set_xlabel('X-axis')
    # ax0.set_ylabel('Y-axis')

    plt.show()
