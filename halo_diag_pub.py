from pydefaults import *
import rayABCD as ABCD
import pydefaults

global zL, Elem, focL

Spaces = [0.25, 0.1, 0.125, 0.150, 0.075]
zL     = np.cumsum(Spaces)
print ('positions:', zL)
# [0.25, 0.25+0.1, 0.25+0.1+0.125, 0.25+0.1+0.125+0.150, 0.25+0.1+0.125+0.150+0.075]
# DMD off###########
#Elem   = ["L", "L", "D", "L", "D"]
#focL   = [0.25, 0.125, 0, 0.05, 0 ]
###################

# DMD ON###########
Elem   = ["L", "L", "DMD", "L", "D"]
focL   = [0.25, 0.125, 0.001, 0.05, 0 ]

zScreen= zL[2]  # location of DMD It is the 3rd element

###################################################
###################### INPUT BEAM PARAMETERS ######
###################################################

# number of rays
N     = 100000
# beam size in m
sigx  = 5e-3
# beam divergence in m
sigxp = 1.25e-2

###################################################
###################### beamline ###################
###################################################
# number of tracking step along beamline
NStep        =  151


# 
# ### lattice definion (3 example)
# 


def Lattice(zinit, dz,X):
# propagate the beam at location zinit by an increemental length dz
# lens position
   global Elem, zL, focL
   print("zin=",zin)

#### lattice definion (3 lenses)
#   Elem   = ["L" , "L",  "L"]                # element type "L"=lens
############next 3 line is a working example
#   Elem   = ["L", "D"]
#   zL     = [0.211277, 0.5]
#   focL   = [0.211277, 0.0]
#
# this 3 lines go up to the DMD
##   Driftto= [0.1, 0.13715


# go up to lens entrance
# I need to bread the loop
   for i in range(len(Elem)):
      if (zinit<=zL[i]) & (zinit+dz>zL[i]):
        print("At element #",i)
        X1=ABCD.Drift(zL[i]-zinit, X)
        
        if Elem[i]=="L":
           print("------------------------ lens")
           print("position [m]=",zL[i])
           print("f [m]=       ",focL[i])
           X2=ABCD.ThinLens(focL[i],X1,aper=0.01)
           Xfin=ABCD.Drift(zinit+dz-zL[i], X2)           
        if Elem[i]=="D":
           print("------------------------ drift")
           print("end of Drift [m]=",zL[i])
           print("L [m]=       ",dz)
           Xfin=ABCD.Drift(zinit+dz-zL[i], X1)
        if Elem[i]=="DMD":
           print("------------------------ DMD")
           print("position [m]=",zL[i])
           print("aperture [m]=       ",focL[i])
           Xfin=ABCD.DMD(focL[i],X1, angle=10)
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

Xinit= ABCD.MkBuddle_G(sigx, sigxp, N)
#Xinit= ABCD.MkBuddle_T(sigx, sigxp, 5, N)

#plt.figure()
#plt.plot(Xinit[0,:], Xinit[2,:],'.')

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
print(len(distance))
print(len(sx))


plt.figure()
ax=plt.subplot (2,1,1)
Profile=Profile/np.max(np.max(Profile))
plt.imshow (Profile, aspect='auto', \
            extent=[0, TotalLength, xmin*1e3, xmax*1e3],\
	    cmap = 'RdPu',\
            norm=LogNorm())
plt.plot (distance,  sx*1000, "C1")
plt.plot (distance, -sx*1000, "C1")
for i in range(len(Elem)):
    if Elem[i]=="L":
        plt.plot (zL[i]*np.ones(10), 0.75*1e3*np.linspace(xmin, xmax,10), \
	         "C7-", linewidth=3)
    if Elem[i]=="D":
        plt.plot (zL[i]*np.ones(10), 0.5*1e3*np.linspace(xmin, xmax,10), \
	          "C8-", linewidth=3)
    if Elem[i]=="DMD":
        plt.plot (zL[i]*np.ones(10), 0.5*1e3*np.linspace(xmin, xmax,10), \
	          "C9-", linewidth=3)

plt.colorbar()
plt.xlabel('distance from YAG screen (m)')
plt.ylabel('$x$, $y$ (mm)')
plt.text(0.075, 0.85,r'$\bf{(a)}$', ha='center', va='center', \
         transform=ax.transAxes, fontsize=18)


ax=plt.subplot (2,3,4)
plt.hist2d(Xinit[0,:]*1e3, Xinit[2,:]*1e3, bins=151, cmap = 'RdPu')
plt.xlabel('$x$ (mm)')
plt.ylabel('$y$ (mm)')
plt.text(0.2, 0.85,r'$\bf{(b)}$', ha='center', va='center', \
         transform=ax.transAxes, fontsize=18)

ax=plt.subplot (2,3,5)
plt.hist2d(Xdmd[0,:]*1e3,  Xdmd[2,:]*1e3,  bins=151, cmap = 'RdPu')
plt.xlabel('$x$ (mm)')
plt.ylabel('$y$ (mm)')
plt.text(0.2, 0.85,r'$\bf{(c)}$', ha='center', va='center', \
         transform=ax.transAxes, fontsize=18)

ax=plt.subplot (2,3,6)
plt.hist2d(Xnew[0,:]*1e3,  Xnew[2,:]*1e3,  bins=151, cmap = 'RdPu')
plt.xlabel('$x$ (mm)')
plt.ylabel('$y$ (mm)')
plt.text(0.2, 0.85,r'$\bf{(d)}$', ha='center', va='center', \
         transform=ax.transAxes, fontsize=18)
plt.tight_layout()
#plt.figure()
#plt.plot (distance, sx)

bbb=plt.figure()

ax=plt.subplot (2,3,1)
plt.hist2d(Xinit[0,:]*1e3, Xinit[2,:]*1e3, bins=151, cmap = 'RdPu')
plt.xlabel('$x$ (mm)')
plt.ylabel('$y$ (mm)')
plt.text(0.2, 0.85,r'$\bf{(a)}$', ha='center', va='center', \
         transform=ax.transAxes, fontsize=18)

ax=plt.subplot (2,3,2)
plt.hist2d(Xdmd[0,:]*1e3,  Xdmd[2,:]*1e3,  bins=151, cmap = 'RdPu')
plt.xlabel('$x$ (mm)')
plt.ylabel('$y$ (mm)')
plt.text(0.2, 0.85,r'$\bf{(b)}$', ha='center', va='center', \
         transform=ax.transAxes, fontsize=18)

ax=plt.subplot (2,3,3)
plt.hist2d(Xnew[0,:]*1e3,  Xnew[2,:]*1e3,  bins=151, cmap = 'RdPu')
plt.xlabel('$x$ (mm)')
plt.ylabel('$y$ (mm)')
plt.text(0.2, 0.85,r'$\bf{(c)}$', ha='center', va='center', \
         transform=ax.transAxes, fontsize=18)
plt.tight_layout()

plt.savefig('target.pdf', bbox_inches='tight')
#plt.figure()
#plt.plot (distance, sx)

plt.show()


#X=Drift(L,X)
#X=ThinLens(f, X)
