from pydefaults import *
import rayABCD as ABCD
import pydefaults
import math as m

global zL, Elem, focL

#Spaces = [0.30,0.42,0.45]
Spaces = [0.36, 0.27, 0.45]
zL = np.cumsum(Spaces)
print ('positions:', zL)

# specify elements and focal lengths
Elem   = ["TL","TL","D"]
#focL   = [0.05,0.20,0]
focL = [0.04, 0.15, 0]
# specify thick lens parameters (thickness, index of refraction, radius of curvature)

'''
from ThorLabs:

System 1
f = 0.05: d = 0.0058, n = 1.460, R = 0.023
f = 0.20: d = 0.0066, n = 1.460, R = 0.092
Li = 0.30 m, L1 = 0.42 m, LC = 0.45 m

System 2
f = 0.04: d = 0.0071, n = 1.460, R = 0.0184
f = 0.15: d = 0.0078, n = 1.460, R = 0.069
Li = 0.36 m, L1 = 0.27 m, Lc = 0.45 m

'''

#d = [0.0058, 0.0066, 0]
#R = [0.023, 0.092, 0]

d = [0.0071, 0.0078, 0]
R = [0.0184, 0.069, 0]
n = [1.460, 1.460, 0]

zScreen= zL[0]

###################################################
###################### INPUT BEAM PARAMETERS ######
###################################################

# number of rays
N     = 100000
# beam size in m
sigx  = [20e-6,50e-6,75e-6,100e-6,125e-6,150e-6,175e-6,200e-6]

# m^2 values; choose which to use
m2 = [1,1.5,2,2.5,3,3.5,4,4.5]
m2_ind = 0

lamb = 262e-9 # wavelength of laser in m

# beam divergence in m
def divergence(sigx):
	return m2[m2_ind] * lamb / (np.pi * sigx)

sigxp = np.zeros(len(sigx))
for i in range(len(sigxp)):
	sigxp[i] = divergence(sigx[i])

# choose which divergence and beam size to use
ind = 3

###################################################
###################### beamline ###################
###################################################
# number of tracking step along beamline
NStep        =  151


# 
# ### lattice definition
# 


def Lattice(zinit, dz,X):
# propagate the beam at location zinit by an incremental length dz
# lens position
   global Elem, zL, focL
   print("zin=",zin)


# go up to lens entrance
# I need to bread the loop
   for i in range(len(Elem)):
      if (zinit<=zL[i]) & (zinit+dz>zL[i]):
        print("At element #",i)
        X1=ABCD.Drift(zL[i]-zinit, X)
        
        if Elem[i]=="L": # thin lens
           print("------------------------ lens")
           print("position [m]=",zL[i])
           print("f [m]=       ",focL[i])
           X2=ABCD.ThinLens(focL[i],X1,aper=0.01)
           Xfin=ABCD.Drift(zinit+dz-zL[i], X2)           
        if Elem[i]=="D": # drift
           print("------------------------ drift")
           print("end of Drift [m]=",zL[i])
           print("L [m]=       ",dz)
           Xfin=ABCD.Drift(zinit+dz-zL[i], X1)
        if Elem[i]=="DMD": # DMD
           print("------------------------ DMD")
           print("position [m]=",zL[i])
           print("aperture [m]=       ",focL[i])
           Xfin=ABCD.DMD(focL[i],X1, angle=10)
        if Elem[i]=="TL": # thick plano-convex lens
           print("------------------------ thick lens")
           print("position [m]=",zL[i])
           print("f [m]=       ",focL[i])
           X2=ABCD.ThickLensPC(d[i],n[i],R[i],X1,aper=0.01)
           Xfin=ABCD.Drift(zinit+dz-zL[i], X2)   
        return(Xfin)
      
   Xfin=ABCD.Drift(dz, X)
   return(Xfin)



###################################################
##################### diagnostics #################
###################################################
# profile savings
Nbins        =  501
xmin         = -2e-2
xmax         =  2e-2
# length of beamline in meter
TotalLength  =  zL[-1]


######################################################################################################
######################################################################################################
######################################################################################################

# initialize array
Profile   = np.zeros((Nbins-1, NStep))
xProfile  = np.linspace(xmin, xmax, Nbins)


      
# generate rays

Xinit= ABCD.MkBuddle_G(sigx[ind], sigxp[ind], N)


XOld = Xinit
Xnew = Xinit
zin  = 0.

if NStep==1:
	dz=TotalLength
else:
	dz=TotalLength/(NStep*1.)

print("step size:s", dz)

distance=np.zeros((NStep))
sx=np.zeros((NStep))
XOld=Xinit.copy()

	
scr=0
for i in range(NStep):
	Xnew=Lattice(zin, dz, XOld)
	if (zin>=zScreen and scr==0):
		Xdmd=Xnew
		scr=1
	xtmp=Xnew[0,:]
	xtmp=xtmp[(np.abs(xtmp)<xmax)]
	cc=np.histogram(xtmp, bins=xProfile)
	Profile[:,i]=cc[0]
	distance[i] =zin
	sx[i]=np.std(xtmp)
	zin+=dz
	XOld=Xnew.copy()



############################
###-----### plotting
############################
#print(len(distance))
#print(len(sx))
sx_rounded = np.round(sx[-1]*10**6,3)

plt.figure()
#ax=plt.subplot (3,2,1)
Profile=Profile/np.max(np.max(Profile))
plt.imshow (Profile, aspect='auto', \
            extent=[0, TotalLength, xmin*1e3, xmax*1e3],cmap = 'RdPu',norm=LogNorm())
plt.plot (distance,  sx*1000, "C1")
plt.plot (distance, -sx*1000, "C1")
for i in range(len(Elem)):
    if Elem[i]=="L":
        plt.plot (zL[i]*np.ones(10), 0.75*1e3*np.linspace(xmin, xmax,10), \
	         "C7-", linewidth=3)
    if Elem[i]=="TL":
        plt.plot (zL[i]*np.ones(10), 0.75*1e3*np.linspace(xmin, xmax,10), \
	         "C6-", linewidth=3)
    if Elem[i]=="D":
        plt.plot (zL[i]*np.ones(10), 0.5*1e3*np.linspace(xmin, xmax,10), \
	          "C8-", linewidth=3)
    if Elem[i]=="DMD":
        plt.plot (zL[i]*np.ones(10), 0.5*1e3*np.linspace(xmin, xmax,10), \
	          "C9-", linewidth=3)

plt.colorbar()
plt.xlabel('distance from cathode (m)')
plt.ylabel('$x$, $y$ (mm)')
#plt.ylim(-2,2)
plt.title("$M^2$ = "+ str(m2[m2_ind]) + " ; final beam size = "+ str(sx_rounded) + "($\mu$m)")
plt.tight_layout()
print(sx_rounded)

plt.show()

