#--------SHEAR FLOW AND TWIST FINDER----------

def twistrate(Sz,Sy,Mz,My,P):
    from math import sqrt, pi, radians, cos, sin
    import numpy as np
    from cs import crosssec
    from baseshear import qbase
    from torsion import torque

    #IMPORT INTERNAL SHEAR FORCES

    c = 0.515  # aileron chord
    h = 0.248  # aileron height
    r = h / 2  # leading edge section radius
    le = sqrt((c - r) ** 2 + r ** 2)  # length of linear section
    circ = pi * r + 2 * le  # circumference of cross section
    phi = circ / 11  # stiffener spacing
    Astiff = 5.4 * 10 ** (-5)  # Stiffener point area (m^2)
    tsk = 0.0011  # skin thickness
    tsp = 0.0022  # spar thickness
    G = 28*10**9 #Shear modulus Al-2024 T3
    A1 = 0.5*h*(c-r) #Area of section 1
    A2 = 0.5*pi*r**2 #Area of section 2
    theta = radians(25) #Aileron deflection angle (rad)


    T = 4000
    Mb = []
    Mbs = 1500 #REPLACE WITH Szt*r IN REAL PROGRAM!
    qb1 = []
    qb2 = []
    qtot = []
    Mz = 20*10**3
    My = 20*10**3
    Sy = 10*10**3
    Sz = 10*10**3
    P = 20.6*10**3

    ##Szt = Sz * cos(theta) - Sy * sin(theta)  # shear force rotated in our coordinate system

    nstiff, zpos, ypos = crosssec(11)
    qt,qr,ql = qbase(Sz,Sy,Mz,My)
    T = torque(Sz,Sy,Mz,My,P)


    for i in range(len(qr)):
        if i==4:
            qbval = qr[i]*h/(2*tsp*A1)
            qb1.append(qbval)
        elif i==3 or i==5:
            qbval = qr[i]*(sqrt((ypos[4]-ypos[3])**2+(zpos[4]-zpos[3])**2))/(2*tsk*A1)
            qb1.append(qbval)
        else:
            qbval = qr[i]*phi/(2*tsk*A1)
            qb1.append(qbval)

    for j in range(len(ql)):
        if j == 0 or j == 3:
            qbval = ql[j]*(phi-sqrt((ypos[4]-ypos[3])**2+(zpos[4]-zpos[3])**2))/(2*tsk*A2)
            qb2.append(qbval)
        elif j == 4:
            qbval = ql[j]*h/(2*tsp*A2)
            qb2.append(qbval)
        else:
            qbval = ql[j]*phi/(2*tsk*A2)
            qb2.append(qbval)

    for k in range(len(nstiff)-1):
        if k == 7 or k == 8:
            Mb.append(0)
        else:
            qbm = qt[k]*(ypos[k+1]-ypos[k])*((c-r)-zpos[k])+qt[k]*(zpos[k+1]-zpos[k])*(ypos[k]+r)
            Mb.append(qbm)




    c1a = (7*phi+2*sqrt((ypos[4]-ypos[3])**2+(zpos[4]-zpos[3])**2))/(2*tsk*A1)
    c1b = h/(2*tsp*A1)
    c2a = (2*phi+2*(phi-sqrt((ypos[4]-ypos[3])**2+(zpos[4]-zpos[3])**2)))/(2*tsk*A2)
    c2b = h/(2*tsp*A2)

    torco = np.array([[c1a+c1b,-c1b,-G],[c2a+c2b,-c2b,-G],[2*A1,2*A2,0]])
    torde = np.array([0,0,T])

    sol = list(np.linalg.solve(torco,torde))



    qco = np.array([[c1a+c1b,-c1b,-G],[c2a+c2b,-c2b,-G],[2*A1,2*A2,0]])
    qde = np.array([(-sum(qb1)-(c1a+c1b)*sol[0]+c1b*sol[1]),(-sum(qb2)-(c2a-c2b)*sol[1]+c2b*sol[0]),Mbs-sum(Mb)])

    solfin = list(np.linalg.solve(qco,qde))


    for x in range(len(qr)):
        if x == 4:
            qt1 = solfin[0]+sol[0]+qr[x]-solfin[1]-sol[1]
            qtot.append(qt)

        else:
            qt1 = solfin[0]+sol[0]+qr[x]
            qtot.append(qt)

    for y in range(len(ql)):
        if y == 4:
            qt1 = solfin[1]+sol[1]+ql[y]-solfin[0]-sol[0]
            qtot.append(qt)
        else:
            qt1 = solfin[1]+sol[1]+ql[y]
            qtot.append(qt)

    return solfin[2] #rate of twist at the section