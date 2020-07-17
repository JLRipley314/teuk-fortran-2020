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
   Dmat = [[0 for i in range(n)] for j in range(n)] 

   N = n-1

   Dmat[0][0] =   mp.fdiv(mp.mpf(2)*mp.power(N,2)+mp.mpf(1),mp.mpf(6))
   Dmat[N][N] = - mp.fdiv(mp.mpf(2)*mp.power(N,2)+mp.mpf(1),mp.mpf(6))
   Dmat[0][N] =   mp.mpf(0.5) * mp.power(-1,N)
   Dmat[N][0] = - mp.mpf(0.5) * mp.power(-1,N)

   for i in range(1,N):
      Dmat[0][i] =   mp.mpf(2.0) * mp.power(-1,i  ) * mp.fdiv(1,mp.mpf(1)-pts[i])
      Dmat[N][i] = - mp.mpf(2.0) * mp.power(-1,i+N) * mp.fdiv(1,mp.mpf(1)+pts[i])
      Dmat[i][0] = - mp.mpf(0.5) * mp.power(-1,i  ) * mp.fdiv(1,mp.mpf(1)-pts[i])
      Dmat[i][N] =   mp.mpf(0.5) * mp.power(-1,i+N) * mp.fdiv(1,mp.mpf(1)+pts[i])

      Dmat[i][i] = - 0.5 * pts[i] * mp.fdiv(1,mp.mpf(1)-mp.power(pts[i],2))

      for j in range(1,N):
         if i!=j:
            Dmat[i][j] = mp.power(-1,i+j) * mp.fdiv(1,pts[i]-pts[j])

   return Dmat
#=============================================================================
def save_cheb(dir_name:str,n:int)->None:
   pts= cheb_pts(n)

   cheb_D_matrix= cheb_D(n)

   with open(dir_name+"/cheb_pts.txt","w") as f:
      for val in pts:
         f.write(mp.nstr(val,8)+'\n')

   with open(dir_name+"/cheb_D.txt","w") as f:
      for line in cheb_D_matrix:
         for val in line:
            f.write(mp.nstr(val,8)+' ')
         f.write('\n')
