#---------------CROSS SECTION DEFINITION------------#
from math import*
from matplotlib import pyplot as plt

##NOTE: coordinate system is defined with origin at the trailing edge and
##positive y' pointing up and positive z' pointing towards the leading edge.
##Stiffeners are numbered starting at trailing edge, moving over positive y'
##towards LE, then back to TE over negative y'. Stiffeners closest to TE are
##both half a stiffener pitch away from the TE. Stiffeners are equally spaced
##over cross section.

#------STARTING VALUES-----------
c = 0.515 #aileron chord
h = 0.248 #aileron height
r = h/2 #leading edge section radius
le = sqrt((c-r)**2+r**2) #length of linear section
circ = pi*r + 2*le #circumference of cross section
phi = circ/11 #stiffener spacing

ypos = [] #y' position of stiffener
zpos = [] #z' position of stiffener
nstiff = [] #Stiffener numbering
n = 0

while n<=10: #Collecting coordinate data of stiffeners
    if n==0:
        dely = (phi/2)*(r/le)
        delz = (phi/2)*((c-r)/le)
        ypos.append(dely)
        zpos.append(delz)
        nstiff.append(n+1)
        

    elif n>0 and ypos[n-1]>0:
        if (ypos[n-1]+(phi*r)/le)<r:
            dely = (phi*r)/le
            delz = (phi*(c-r))/le
            ynew = dely + ypos[n-1]
            znew = delz + zpos[n-1]
            ypos.append(ynew)
            zpos.append(znew)
            nstiff.append(n+1)
            
        elif (ypos[n-1]+(phi*r)/le)>r and zpos[n-1]<(c-r):
            eta = phi - (le*(r-ypos[n-1]))/r
            ynew = r*cos(eta/r)
            znew= ((c-r)+r*sin(eta/r))
            ypos.append(ynew)
            zpos.append(znew)
            nstiff.append(n+1)

        elif zpos[n-1]>(c-r):
            ypos.append(0)
            zpos.append(c)
            nstiff.append(n+1)

    elif n>0 and ypos[n-1]<=0:
        if ypos[n-1]==0:
            dely = -(phi/2)*(r/le)
            delz = (phi/2)*((c-r)/le)
            ypos.append(dely)
            zpos.append(delz)
            nstiff.append(n+1)
        
        if (ypos[n-1]-(phi*r)/le)>-r and ypos[n-1]!=0:
            dely = -(phi*r)/le
            delz = (phi*(c-r))/le
            ynew = dely + ypos[n-1]
            znew = delz + zpos[n-1]
            ypos.append(ynew)
            zpos.append(znew)
            nstiff.append(n+1)
            
        elif (ypos[n-1]-(phi*r)/le)<-r and zpos[n-1]<(c-r):
            eta = phi - (le*(r+ypos[n-1]))/r
            ynew = -r*cos(eta/r)
            znew= ((c-r)+r*sin(eta/r))
            ypos.append(ynew)
            zpos.append(znew)
            nstiff.append(n+1)
    n+=1

print(zpos,ypos,nstiff)
plt.plot(zpos, ypos, "ro")
plt.gca().invert_xaxis()
plt.axis("equal")
plt.grid(linestyle='-', linewidth='0.5')
plt.show()
