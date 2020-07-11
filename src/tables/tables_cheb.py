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
   return mp.diff(lambda x: mp.chebyt(n,x) x)
#=============================================================================
def write_to_file(name:str,arr:List[List[mp.mpf]])->None:
   with open(name,'w') as f:
      for line in arr:
         for val in line:
            f.write(mp.nstr(val,16)+' ')
            f.write('\n')
#=============================================================================
def save_cheb_D_matrix(dir_name:str,n:int)->None:
   pts= cheb_pts(n)

   cheb_D_matrix= [
      [diff_cheb(n,x) for x in pts]
      for n in range(n)
   ] 

   write_to_file(
      "{}/cheb_D_matrix_{}.txt".format(dir_name,n),
      cheb_D_matrix
   )
