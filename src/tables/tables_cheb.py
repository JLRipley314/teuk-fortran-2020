#=============================================================================
## computes derivative matrices for Chebyshev polynomials T_n 
## at Chebyshev points  
#=============================================================================

import mpmath as mp
from typing import List

mp.prec= +mp.inf 

#=============================================================================
def cheb_pts(n:int)->List[mp.mpf]:
   return [mp.cos(mp.pi*i/(n-1)) for i in range(n)]
#=============================================================================
def diff_cheb(n:int,x:float)->mp.mpf:
   return mp.diff(lambda x: mp.chebyt(n,x), x)
#=============================================================================
def cheb_D(n:int)->List[mp.mpf]:
   pts= cheb_pts(n)
   return [
      [diff_cheb(n,x) for x in pts]
      for n in range(n)
   ] 
#=============================================================================
def save_cheb(dir_name:str,n:int)->None:
   pts= cheb_pts(n)

   cheb_D_matrix= [
      [diff_cheb(n,x) for x in pts]
      for n in range(n)
   ] 

   with open(dir_name+"/cheb_pts.txt","w") as f:
      for val in pts:
         f.write(mp.nstr(val,32)+'\n')

   with open(dir_name+"/cheb_D.txt","w") as f:
      for line in cheb_D_matrix:
         for val in line:
            f.write(mp.nstr(val,32)+' ')
         f.write('\n')
