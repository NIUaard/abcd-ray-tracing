'''
this script defines a set of function to ray trace optical ray
P.P., 11-17-2015

a buddle of ray is arranged as the 4xn matrix

/                             \
| x_0    x_1  ........ x_n    |
| x'_0   x'_1......... x'_n   |
| y_0    y_1...........y_n    |
| y'_0   y'_1..........y'_n   |
\                             /

'''


import numpy as np
import numpy.random as random
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap






##################### end cmap definition ####################################

def MkBuddle(sigx, sigxp, N):
   '''
     generate rays
   '''
   X=np.zeros((4,N))
   X[0,:]=random.normal(0., sigx,  N)
   X[1,:]=random.normal(0., sigxp, N)
   X[2,:]=random.normal(0., sigx,  N)
   X[3,:]=random.normal(0., sigxp, N)

   return(X)

def MkBuddle_G(sigx, sigxp, N):
   '''
     generate rays distributed to follow a Gaussian distribution
   '''
   X=np.zeros((4,N))
   X[0,:]=random.normal(0., sigx,  N)
   X[1,:]=random.normal(0., sigxp, N)
   X[2,:]=random.normal(0., sigx,  N)
   X[3,:]=random.normal(0., sigxp, N)

   return(X)

def MkBuddle_T(sigx, sigxp, NCirc, N):
   '''
     generate rays distributed to follow a grid pattern
   '''
   sigx = 2*sigx
   sigxp= 2*sigxp
   X=np.zeros((4,N))
   partperline=int(N/NCirc)
   print(partperline)
   for i in range(NCirc):
      tmp=random.uniform(-1, 1,  partperline)*(i+1)*sigx/NCirc
      X[0,i*partperline:(i+1)*partperline]=\
                    np.hstack([tmp[0::2], tmp[1::2]])
      tmpp= np.sqrt(((i+1)*sigx/NCirc)**2-tmp[0::2]**2)
      tmpm=-np.sqrt(((i+1)*sigx/NCirc)**2-tmp[1::2]**2)

      X[2,i*partperline:(i+1)*partperline]= \
                    np.hstack([tmpp, tmpm])
      
   X[1,:]=random.normal(0., sigxp, N)
   X[3,:]=random.normal(0., sigxp, N)

   return(X)

def Drift(L,X):
   '''
      drift matrix and ray tranformation
   '''
   M=np.zeros((4,4))
   M[0,0]=M[2,2]=1.; M[0,1]=M[2,3]=L; 
   M[1,0]=M[3,2]=0.; M[1,1]=M[3,3]=1.;
   
   Y=np.dot(M,X)

   return(Y)

def heaviside(x):
   if x==0:
      return (0.5)
   if x<0:
      return (1.0)
   if x>0:
      return (0.0) 

def ThinLens(f,X, aper=0, cen=[0,0]):
   '''
     thin lens matrix and ray tranformation
     ray with radius>aper are left untouched
   '''
   M=np.zeros((4,4))
   M[0,0]=M[2,2]=1.;     M[0,1]=M[2,3]=0.; 
   M[1,0]=M[3,2]=-1./f ; M[1,1]=M[3,3]=1.;


   if aper==0:
     Y=np.dot(M,X)
   else:
     R=np.sqrt((X[0,:]-cen[0])**2+(X[2,:]-cen[1])**2)
     X2=X[:,(R<aper)]
     Y=np.dot(M,X2)

   return(Y)

def DMD(R ,X, angle=10 , cen=[0,0]):
   '''
     DMD simulation
     ray with radius>aper are left untouched
   '''   
   if np.abs(angle)>0:
     for i in range(np.shape(X)[1]):
        if ((X[0,i]-cen[0])**2+(X[2,i]-cen[1])**2< R**2):
           X[1,i]= X[1,i]+angle*np.pi/180
#d         X[3,i]= X[3,i]+angle*np.pi/180
       
   return(X)

'''   
   if aper==0:
     Y=np.dot(M,X)
   else:
     R=(1.+np.sign(aper**2-((X[0,:]-cen[0])**2+(X[2,:]-cen[1])**2)))/2.
     T=np.diag(R)
     I=np.diag(np.ones(len(R)))
     Y=np.dot(M,X)*R+X*(1-R)
'''

   
def MLA(fmic, pitch, radius, thick, num, X):
   '''

      Simulate a MLA via piecewise matrices

        fmic        focal lens of microlens
        pitch        pitch of MLA
        radius        radius of the MLA (radius< pitch)
        thick        thickness of substrate // not used
        num     makes a numXnum array of micro lenses
        X        input rays

   '''
   print((np.shape(X)))
   Y=X
   Xc=Yc=np.linspace(-num/2.*pitch,num/2.*pitch, num)
   for j in range(num):
      for i in range(num):
          X=ThinLens(fmic,X, radius, cen=[Xc[i], Yc[j]])
   return(X)          


#
def MLAold(fmic, pitch, radius, thick, num, X):
   '''

      Simulate a MLA via piecewise matrices

        fmic        focal lens of microlens
        pitch        pitch of MLA
        radius        radius of the MLA (radius< pitch)
        pitch        thickness of substrate // not used
        num     makes a numXnum array of micro lenses
        X        input rays

   '''
   print((np.shape(X)))
   Y=X
   Xc=Yc=np.linspace(-num/2.*pitch,num/2.*pitch, num)
   for j in range(num):
      for i in range(num):
#          print (Xc[i], Yc[j])
          Rho=np.sqrt((X[0,:]-Xc[i])**2+(X[2,:]-Yc[j])**2)
#          print "here"
#          print (np.shape(Rho))
          index=np.where(Rho<radius)
#          print index[:]
          if len(index)>0:
             X=X
# how to vectorize this loop ???
#          for k in range(len(index)):
#             if len(index)>0:
#                 print (len(index))
#               Y[:,index[k]]=ThinLens(fmic, X[:,index[k]])
#              Y[:,index]=ThinLens(fmic, X[:,index])

#             print (np.shape(zip(*zip(*X[:,index]))))
#          Y[:,index]=ThinLens(fmic, X[:,index])
#          print (X[:,index])
   return(X)          

#
def coord_stat(x):
#
# comptute stat parameters of a 1D array (x)
# cal is a muliplier factor (pixel to real length)
#
# returns norm, mean, rms, fwhm
#
# 
   Nbin=21
   profile, pix_p  =np.histogram(x)
   profile = np.hstack((profile, 0.0))
   norm_p  = sum(profile)
   mean_p  = 1/norm_p*sum(pix_p*profile)
   var_p   = 1/norm_p*sum(pix_p**2*profile)
   rms_p   = np.sqrt(var_p-mean_p**2)

   statistics=[norm_p, mean_p, rms_p]

   return (statistics)
#
def dump_prof(x, xmin, xmax, Nbin=21):
#
# comptute stat parameters of a 1D array (x)
# cal is a muliplier factor (pixel to real length)
#
# returns norm, mean, rms, fwhm
#
# 
   profile, pix_p  =np.histogram(x)
   profile = np.hstack((profile, 0.0))
   norm_p  = sum(profile)
   mean_p  = 1/norm_p*sum(pix_p*profile)
   var_p   = 1/norm_p*sum(pix_p**2*profile)
   rms_p   = np.sqrt(var_p-mean_p**2)

   statistics=[norm_p, mean_p, rms_p]

   return (statistics)
