#=============================================================================
## computes values for normalized spin-weighted associated Legendre (swaL)
## functions over the Legendre nodes for Gaussian quadrature
#=============================================================================
import mpmath as mp
from typing import List
from tables_Legendre import roots_weights_Legendre

mp.prec= +mp.inf 
#=============================================================================
def swaL(spin:int,m_ang:int,l_ang:int,x:float)->mp.mpf:
   al= mp.mpf(abs(m_ang-spin))
   be= mp.mpf(abs(m_ang+spin))
   assert((al+be)%2==0)
   n= mp.mpf(l_ang - (al+be)/2)

   norm= mp.sqrt(
       (2*n+al+be+1)*mp.power(2,-al-be-1)
   *	mp.fdiv(mp.fac(n+al+be),mp.fac(n+al))
   *	mp.fdiv(mp.fac(n)      ,mp.fac(n+be))
   )
   norm*= pow(-1,max(m_ang,-spin))

   return norm*mp.power(1-x,al/2.)*mp.power(1+x,be/2.)*mp.jacobi(n,al,be,x)
#=============================================================================
def diff_swaL(spin:int,m_ang:int,l_ang:int,x:mp.mpf)->mp.mpf:
   return mp.diff(lambda x: swaL(spin,m_ang,l_ang,x),x)
#=============================================================================
def write_to_file(name:str,arr:List[List[float]])->None:
   with open(name,'w') as f:
      for line in arr:
         for val in line:
            f.write(mp.nstr(val,16)+' ')
            f.write('\n')
#=============================================================================
def save_Gauss_quad_vals_swaL(
   dir_name:str,spin:int,m_ang:int,nl:int,gauss_pts:int)->None:

   roots, weights= roots_weights_Legendre(gauss_pts)

   lmin= max(abs(m_ang),abs(spin))

   swaL_vals= [ 
      [swaL(spin,m_ang,l_ang,root) for root in roots]
      for l_ang in range(lmin,lmin+nl+1)
   ]
   write_to_file(
      "{}/s_{}_m_{}.txt".format(dir_name,spin,m_ang),
      swaL_vals
   )
#=============================================================================
def write_module(
   dir_name:str,spin:int,m_ang:int,nl:int,gauss_pts:int)->None:

   roots, weights= roots_weights_Legendre(gauss_pts)

   lmin= max(abs(m_ang),abs(spin))

   swaL_vals= [ 
      [swaL(spin,m_ang,l_ang,root) for root in roots]
      for l_ang in range(lmin,lmin+nl+1)
    ]
   write_to_file(
      "{}/s_{}_m_{}.txt".format(dir_name,spin,m_ang),
      swaL_vals
   )
   name= 'src/mod_coords.f90'
   with open(name,'w') as f:
      for y in roots:
         for l in range(0,nl):
            for m in range(min_m,max_m+1):
               for s in [-3,-2,-1,0,1,2]:
                  f.write(mp.str(swal(s,m,l,y),16)+'_rp &\n')
