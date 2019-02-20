from matplotlib import pyplot as plt
import numpy as np



a = np.array([[3,1], [1,2]])
b = np.array([9,8])

x = np.linalg.solve(a, b)

print np.allclose(np.dot(a, x), b)

print x

x1 = 0.174
x2 = 1.051
x3 = 2.512
d1 = 1.034
d3 = 2.066
xa = 0.3
theta = 25.0*(np.pi)/180
E = 1.0*(10**9)
Iyy = 1.0*(10**(-5))
Izz = 1.0*(10**(-5))
q = 1.0*(10**3)
L = 2.691
P = 20.6*(10**3)
ch = 0.515
h = 0.248

def reactionforces(E, Izz, Iyy, q, P, theta, xa, x1, x2, x3, d1, d3, ch, h, L):

    ## H2x == 0

    ## Rotated forces and deflections

    qy = q*np.cos(theta)
    qz = q*np.sin(theta)
    Py = P*np.sin(theta)
    Pz = P*np.cos(theta)
    dy1 = d1*np.cos(theta)
    dz1 = d1*np.sin(theta)
    dy2 = 0.0
    dz2 = 0.0
    dy3 = d3*np.cos(theta)
    dz3 = d3*np.sin(theta)

    ## Moment about the hinge line
    ## Actuator force

    A = (qy*((0.25*ch) - (h/2)) + Pz*(h/2) - Py*(h/2))/((np.sin(theta)*(h/2)) - (np.cos(theta)*(h/2)))

    Ay = A*np.sin(theta)
    Az = A*np.cos(theta)

## -------------------------------------------------------------------------------------------------------------------- ##

    ## [H1y,H2y,H3y,cf1,cf2]

    ## Maccauly at Hinge 1

    C11 = 0.0
    C12 = 0.0
    C13 = 0.0
    C14 = -x1
    C15 = -1.0

    ## Maccauly at Hinge 2

    C21 = ((x2 - x1)**3)/6
    C22 = 0.0
    C23 = 0.0
    C24 = -x2
    C25 = -1.0

    ## Maccauly at Hinge 3

    C31 = ((x3 - x1)**3)/6
    C32 = ((x3 - x2)**3)/6
    C33 = 0.0
    C34 = -x3
    C35 = -1.0

    ## Force equilibrium

    C41 = 1.0
    C42 = 1.0
    C43 = 1.0
    C44 = 0.0
    C45 = 0.0

    ## Moment equilibrium
    
    C51 = x1
    C52 = x2
    C53 = x3
    C54 = 0.0
    C55 = 0.0

    ## Matric of coefficients

    C = np.array([[C11, C12, C13, C14, C15],
                  [C21, C22, C23, C24, C25],
                  [C31, C32, C33, C34, C35],
                  [C41, C42, C43, C44, C45],
                  [C51, C52, C53, C54, C55]])

    ## Matrix of knowns K

    K1 = (qy/24)*(x1**4) - (E*Iyy*dy1)
    K2 = (qy/24)*(x2**4) - (Py/6)*((x2+(xa/2)-x2)**3) - (E*Iyy*dy2)
    K3 = (qy/24)*(x3**4) - (Py/6)*((x3+(xa/2)-x2)**3) - (Ay/6)*((x3-(xa/2)-x2)**3) - (E*Iyy*dy3)
    K4 = qy*L - Py - Ay
    K5 = (qy/2)*(L**2) - Ay*(x2-(xa/2)) - Py*(x2+(xa/2))

    K = np.array([K1, K2, K3, K4, K5])

## -------------------------------------------------------------------------------------------------------------------- ##

    ## [H1z,H2z,H3z,cf1,cf2]

    ## Maccauly at Hinge 1

    D11 = 0.0
    D12 = 0.0
    D13 = 0.0
    D14 = x1
    D15 = 1.0

    ## Maccauly at Hinge 2
    
    D21 = ((x2 - x1)**3)/6
    D22 = 0.0
    D23 = 0.0
    D24 = x2
    D25 = 1.0

    ## Maccauly at Hinge 3

    D31 = ((x3 - x1)**3)/6
    D32 = ((x3 - x2)**3)/6
    D33 = 0.0
    D34 = x3
    D35 = 1.0

    ## Force equilibrium

    D41 = 1.0
    D42 = 1.0
    D43 = 1.0
    D44 = 0.0
    D45 = 0.0

    ## Moment equilibrium
    
    D51 = x1
    D52 = x2
    D53 = x3
    D54 = 0.0
    D55 = 0.0

    ## Matric of coefficients

    D = np.array([[D11, D12, D13, D14, D15],
                  [D21, D22, D23, D24, D25],
                  [D31, D32, D33, D34, D35],
                  [D41, D42, D43, D44, D45],
                  [D51, D52, D53, D54, D55]])

    ## Matrix of knowns K

    M1 = (qz/24)*(x1**4) - (E*Izz*dz1)
    M2 = (qz/24)*(x2**4) - (Pz/6)*((x2+(xa/2)-x2)**3) - (E*Izz*dz2)
    M3 = (qz/24)*(x3**4) - (Pz/6)*((x3+(xa/2)-x2)**3) - (Az/6)*((x3-(xa/2)-x2)**3) - (E*Izz*dz3)
    M4 = qz*L - Pz - Az
    M5 = (qz/2)*(L**2) - Az*(x2-(xa/2)) - Pz*(x2+(xa/2))

    M = np.array([M1, M2, M3, M4, M5])

    return (np.linalg.solve(C, K), np.linalg.solve(D, M), A)

print reactionforces(E, Izz, Iyy, q, P, theta, xa, x1, x2, x3, d1, d3, ch, h, L)
