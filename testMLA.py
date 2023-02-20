from pydefaults import *
import rayABCD as ABCD

N     = 100000
sigx  = 3.e-3
sigxp = 0.0001e-3

# beamline
TotalLength  =  1
NStep        =  100

# profile savings
Nbins        =  501
xmin         = -1e-2
xmax         =  1e-2

# initialize array
Profile   = np.zeros((Nbins-1, NStep))
xProfile  = np.linspace(xmin, xmax, Nbins)

# -------------------------------------------------------------
# ---------------------  define lattice  ----------------------
# -------------------------------------------------------------

def Lattice(z,X):
# lens position
   zMLA1  = 0.5
   zMLA2  = zMLA1+0.5
   focMLA = 0.1  
   zL     = 0.6
   focL   = 0.5  
# go up to MLA entrance
   if (z<zMLA1) & (z<zMLA1):
      print ("before MLA1")
      X=ABCD.Drift(z, X)
# after MLA1 go to L entrance
   if (z>=zMLA1) & (z<zMLA2):
      print ("between MLA1 & MLA2")
      X=ABCD.Drift(zMLA1, X)
      X=ABCD.MLA(focL, 500e-6, 90e-6, 0, 10, X)
      X=ABCD.Drift(z-zMLA1, X)
# after MLA2 go to L entrance
#   if (z>=zMLA2) & (z<zL):
#      X0=ABCD.Drift(zMLA2, X)
#      X1=ABCD.ThinLens(focL,X0)
#      X1=ABCD.MLA(focL, 500e-6, 90e-6, 0, 10, X0)
#      X2=ABCD.Drift(z-zMLA2, X1)
# after L just drift
   if (z>=zL):
      print ("after length")
      X=ABCD.ThinLens(focL,X)
      X=ABCD.Drift(z-zL, X)

   return(X)

def LatticeOld(z,X):
# lens position
   zL     = 0.5
   focL   = 0.25  
# go up to MLA entrance
   if (z<zL):
      return(ABCD.Drift(z, X))
# after MLA go to L entrance
   if (z>=zL):
      X0=ABCD.Drift(zL, X)
      X1=ABCD.ThinLens(focL,X0)
      X2=ABCD.Drift(z-zL, X1)
      return(X2)

# -------------------------------------------------------------
# ------------------ end lattice definition -------------------
# -------------------------------------------------------------

      
# generate rays

Xinit=ABCD.MkBuddle(sigx, sigxp, N)

#XOld=Xinit

for i in range(NStep):
    L=TotalLength/(NStep*1.)*i
    print (L)
    Xnew=Lattice(L, Xinit)
    cc=np.histogram(Xnew[0,:], bins=xProfile)
    Profile[:,i]=cc[0]
#    XOld=Xnew
#norm_x, mean_x, rms_x =ABCD.coord_stat(X[0,:])
#print (rms_x)

plt.figure()
plt.hist2d(Xinit[0,:], Xinit[2,:], bins=501)
plt.colorbar()
plt.figure()
plt.hist2d(Xnew[0,:],  Xnew[2,:],  bins=501)
plt.colorbar()


plt.figure()
plt.imshow (Profile, aspect='auto',cmap='Spectral', extent=[0, TotalLength, xmin*1e3, xmax*1e3])
plt.colorbar()
plt.xlabel('distance (m)')
plt.ylabel('x (mm')

plt.show()


#X=Drift(L,X)
#X=ThinLens(f, X)
