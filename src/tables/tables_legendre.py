#=============================================================================
## computes the Legendre roots and weights for Gaussian quadrature
#=============================================================================

from typing import List
import mpmath as mp

roots_precision=1e-32

mp.prec= mp.inf 
#=============================================================================
## we provide an initial guess (due to Tricomi) for the roots
## and then root polish with a Newton-Raphson method
#=============================================================================
def roots_weights_Legendre(n:int)->(List[mp.mpf],List[mp.mpf]):
   assert(n%2==0)

   roots=   [mp.mpf(0) for i in range(n)]
   weights= [mp.mpf(0) for i in range(n)]

   for i in range(1,int(n/2)+1):
      z= (
      (1-mp.fdiv(1,8)*mp.power(n,-2)+mp.fdiv(1,8)*mp.power(n,-3))
      *   mp.cos(mp.pi*mp.fdiv(4*i-1,4*n+2))
      )   
      z_old= mp.mpf(z) 
      f= mp.mpf(0)
      der_f= mp.mpf(0)
      while True:
         f=     mp.legendre(n,z) 
         der_f= mp.diff(lambda z:mp.legendre(n,z),z) 
         z_old= z
         z= mp.fsub(z,mp.fdiv(f,der_f)) 
         if mp.fabs(z-z_old)<roots_precision:
            break
      roots[n-i]= z
      roots[i-1]= -roots[n-i] 
      weights[n-i]= 2/((1-mp.power(z,2))*mp.power(der_f,2));
      weights[i-1]= weights[n-i];
   return (roots, weights)
#=============================================================================
def save_roots_weights_Legendre(dir_name:str,gauss_pts:int)->None:
   (roots,weights)= roots_weights_Legendre(gauss_pts)

   with open(dir_name+'/roots_legendre.txt','w') as f:
      for i in range(gauss_pts):
         f.write(mp.nstr(roots[i],32)+'\n')

   with open(dir_name+'/weights_legendre.txt','w') as f:
      for i in range(gauss_pts):
         f.write(mp.nstr(weights[i],32)+'\n')
