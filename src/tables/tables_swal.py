#=============================================================================
## computes values for normalized spin-weighted associated Legendre (swaL)
## functions over the Legendre nodes for Gaussian quadrature
#=============================================================================
import mpmath as mp
from typing import List
from tables_legendre import roots_weights_Legendre

mp.prec= +mp.inf 
#=============================================================================
def swaL(spin:int,m_ang:int,l_ang:int,x:float)->mp.mpf:
   al= mp.mpf(abs(m_ang-spin))
   be= mp.mpf(abs(m_ang+spin))
   assert((al+be)%2==0)
   n= mp.mpf(l_ang - (al+be)/2)

   if n<0:
      return mp.mpf(0)

   norm= mp.sqrt(
       (2*n+al+be+1)*mp.power(2,-al-be-1)
   *	mp.fdiv(mp.fac(n+al+be),mp.fac(n+al))
   *	mp.fdiv(mp.fac(n)      ,mp.fac(n+be))
   )
   norm*= mp.power(-1,max(m_ang,-spin))

   return norm*mp.power(1-x,al/2.)*mp.power(1+x,be/2.)*mp.jacobi(n,al,be,x)
#=============================================================================
def diff_swaL(spin:int,m_ang:int,l_ang:int,x:mp.mpf)->mp.mpf:
   return mp.diff(lambda x: swaL(spin,m_ang,l_ang,x),x)
#=============================================================================
def write_to_file(name:str,arr:List[List[float]])->None:
   with open(name,'w') as f:
      for line in arr:
         for val in line:
            f.write(mp.nstr(val,32)+' ')
         f.write('\n')
#=============================================================================
def save_Gauss_quad_vals_swaL(
   dir_name:str,spin:int,m_ang:int,nl:int,gauss_pts:int)->None:

   roots, weights= roots_weights_Legendre(gauss_pts)

   cy = [-pt                                   for pt in roots]
   sy = [mp.sqrt(mp.fadd(1,-pt)*mp.fadd(1,pt)) for pt in roots]

   with open(dir_name+'/cos.txt','w') as f:
      for i in range(len(cy)):
         f.write(mp.nstr(cy[i],32)+'\n')

   with open(dir_name+'/sin.txt','w') as f:
      for i in range(len(sy)):
         f.write(mp.nstr(sy[i],32)+'\n')

   swaL_vals= [ 
      [swaL(spin,m_ang,l_ang,root) for l_ang in range(0,nl)]
      for root in roots 
   ]
   write_to_file(
      "{}/s_{}_m_{}.txt".format(dir_name,spin,m_ang),
      swaL_vals
   )
