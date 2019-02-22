#------SHEAR CENTER FINDER-----------


def shearcent(Mz, My):
    from math import sqrt, pi
    from cs import crosssec
    from simmoi import moisim
    from simcent import centsim
    from Boomarea import aboom
    import numpy as np
    #------Input functions: moments of inertia, aileron cross section geometry

    c = 0.515  # aileron chord
    h = 0.248  # aileron height
    r = h / 2  # leading edge section radius
    le = sqrt((c - r) ** 2 + r ** 2)  # length of linear section
    circ = pi * r + 2 * le  # circumference of cross section
    phi = circ / 11  # stiffener spacing
    Astiff = 5.4 * 10 ** (-5)  # Stiffener point area (m^2)
    tsk = 0.0011  # skin thickness
    tsp = 0.0022  # spar thickness
    A1 = 0.5*h*(c-r) #Area of section 1
    A2 = 0.5*pi*r**2 #Area of section 2

    Izz, Iyy = moisim(Mz, My)
    B = aboom(Mz,My)
    zc, yc = centsim(Mz,My)
    nstiff, zpos, ypos = crosssec(11)

    c1a = (7*phi+2*sqrt((ypos[4]-ypos[3])**2+(zpos[4]-zpos[3])**2))/(2*tsk*A1)
    c1b = h/(2*tsp*A1)
    c2a = (2*phi+2*(phi-sqrt((ypos[4]-ypos[3])**2+(zpos[4]-zpos[3])**2)))/(2*tsk*A2)
    c2b = h/(2*tsp*A2)

    siy = 1/Izz
    siz = 1/Iyy

    lstqbright = []
    lstqbleft = []
    lstqbright1 = []
    lstqbleft1 = []
    total1 = []
    total2 = []
    rightcell = B[0:5] + B[8:13]
    leftcell = B[4:9]
    qb1 = []
    qb2 = []
    qb11 = []
    qb22 = []
    Mb1 = []
    Mb2 = []

    #--------- y'-location--------
    for n in range(len(rightcell)):
        if n == 4 or n == 9:
            lstqbright.append(0)
        elif n<5:
            lstqbright.append(siy*B[n]*ypos[n])
        else:
            lstqbright.append(siy*B[n]*ypos[n+3])
    for n in range(len(leftcell)):
        if n == 4:
            lstqbleft.append(0)
        else:
            lstqbleft.append(siy*B[n]*ypos[n+5])

    for n in range(len(B)):
        qb = siy*B[n]*ypos[n]
        total1.append(qb)

    for i in range(len(rightcell)):
        if i==4:
            qbval = lstqbright[i]*h/(2*tsp*A1)
            qb1.append(qbval)
        elif i==3 or i==5:
            qbval = lstqbright[i]*(sqrt((ypos[4]-ypos[3])**2+(zpos[4]-zpos[3])**2))/(2*tsk*A1)
            qb1.append(qbval)
        else:
            qbval = lstqbright[i]*phi/(2*tsk*A1)
            qb1.append(qbval)

    for j in range(len(leftcell)):
        if j == 0 or j == 3:
            qbval = lstqbleft[j]*(phi-sqrt((ypos[4]-ypos[3])**2+(zpos[4]-zpos[3])**2))/(2*tsk*A2)
            qb2.append(qbval)
        elif j == 4:
            qbval = lstqbleft[j]*h/(2*tsp*A2)
            qb2.append(qbval)
        else:
            qbval = lstqbleft[j]*phi/(2*tsk*A2)
            qb2.append(qbval)

    qar1 = np.array([[c1a+c1b,-c1b],[c2a+c2b,-c2b]])
    qde1 = np.array([-sum(qb1),-sum(qb2)])

    sol1 = np.linalg.solve(qar1,qde1)


    for x in range(len(total1)-1):
        val = total1[x]*(-(ypos[x+1]-ypos[x])*zpos[x]+(zpos[x+1]-zpos[x])*ypos[x])
        Mb1.append(val)

    ys = (sum(Mb1)+sol1[0]*A1+sol1[1]*A2)-r


    #--------z' location-------
    for n in range(len(rightcell)):
        if n == 4 or n == 9:
            lstqbright1.append(0)
        elif n<5:
            lstqbright1.append(siz*B[n]*zpos[n])
        else:
            lstqbright1.append(siz*B[n]*zpos[n+3])
    for n in range(len(leftcell)):
        if n == 4:
            lstqbleft1.append(0)
        else:
            lstqbleft1.append(siz*B[n]*zpos[n+5])

    for n in range(len(B)):
        qb = siz*B[n]*zpos[n]
        total2.append(qb)

    for i in range(len(rightcell)):
        if i==4:
            qbval = lstqbright1[i]*h/(2*tsp*A1)
            qb11 .append(qbval)
        elif i==3 or i==5:
            qbval = lstqbright1[i]*(sqrt((ypos[4]-ypos[3])**2+(zpos[4]-zpos[3])**2))/(2*tsk*A1)
            qb11.append(qbval)
        else:
            qbval = lstqbright1[i]*phi/(2*tsk*A1)
            qb11.append(qbval)

    for j in range(len(leftcell)):
        if j == 0 or j == 3:
            qbval = lstqbleft1[j]*(phi-sqrt((ypos[4]-ypos[3])**2+(zpos[4]-zpos[3])**2))/(2*tsk*A2)
            qb22.append(qbval)
        elif j == 4:
            qbval = lstqbleft[j]*h/(2*tsp*A2)
            qb22.append(qbval)
        else:
            qbval = lstqbleft[j]*phi/(2*tsk*A2)
            qb22.append(qbval)

    qar2 = np.array([[c1a+c1b,-c1b],[c2a+c2b,-c2b]])
    qde2 = np.array([-sum(qb11),-sum(qb22)])

    sol2 = np.linalg.solve(qar2,qde2)


    for x in range(len(total2)-1):
        val = total2[x]*(-(ypos[x+1]-ypos[x])*zpos[x]+(zpos[x+1]-zpos[x])*ypos[x])
        Mb2.append(val)

    zs = (sum(Mb2)+sol2[0]*A1+sol2[1]*A2)
    return zs, ys
