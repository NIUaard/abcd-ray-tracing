from pydefaults import *
import rayABCD as ABCD
import pydefaults

global zL, Elem, focL

Spaces = [0.01,0.005,0.010,0.03]
zL = np.cumsum(Spaces)
#zL     = np.cumsum(Spaces)
print ('positions:', zL)

# our test MLA is MLA_example in Downloads
Elem   = ["MLA","D","MLA","D"]
focL   = 0.0467

zScreen= zL[0]  # location of DMD It is the 3rd element

# number of rays
N     = 100000
# beam size in m
sigx  = [3e-3,9e-4,5e-3,2.5e-3,5e-3]
# beam divergence in m
m2 = 1
lamb = 262e-9

def divergence(sigx):
	return m2 * lamb / (np.pi * sigx)

sigxp = np.zeros(len(sigx))
for i in range(len(sigxp)):
	sigxp[i] = divergence(sigx[i])

# choose beam size
ind = 4

###################################################
###################### beamline ###################
###################################################
# number of tracking step along beamline
NStep        =  151

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
        if Elem[i]=="D":
           print("------------------------ drift")
           print("end of Drift [m]=",zL[i])
           print("L [m]=       ",dz)
           Xfin=ABCD.Drift(zinit+dz-zL[i], X1)
        if Elem[i]=="MLA":
           print("------------------------ MLA")
           print("position [m]=",zL[i])
           print("aperture [m]=       ",focL)
           Xfin=ABCD.MLA(focL,0.5e-3, 0.25e-3, 0, 20, X1)
        return(Xfin)
      
   Xfin=ABCD.Drift(dz, X)
   #print(Xfin[0])
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
	
print("shape is ",np.shape(XOld))
############################
###-----### plotting
############################
#print(len(distance))
#print(len(sx))


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
    if Elem[i]=="D":
        plt.plot (zL[i]*np.ones(10), 0.5*1e3*np.linspace(xmin, xmax,10), \
	          "C8-", linewidth=3)
    if Elem[i]=="MLA":
        plt.plot (zL[i]*np.ones(10), 0.5*1e3*np.linspace(xmin, xmax,10), \
	          "C9-", linewidth=3)

plt.colorbar()
plt.xlabel('distance from cathode (m)')
plt.ylabel('$x$, $y$ (mm)')
#plt.ylim(-2,2)
plt.title("beam size = "+ str(sigx[ind]*10**6) + " ($\mu$m)")
plt.tight_layout()


plt.show()
