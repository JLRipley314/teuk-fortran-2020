#=============================================================================
import subprocess, os, sys, time
from typing import List 

sys.path.insert(1,os.getcwd()+'/src/tables/')
#=============================================================================
from tables_Legendre import save_roots_weights_Legendre
from tables_swaL import save_Gauss_quad_vals_swaL 
#=============================================================================
class Sim:
#=============================================================================
   def __init__(self,args:List[str])->None:
      self.home_dir= str(os.getcwd())

      assert len(args) > 1, (
         'argv[1] is empty-meed a run_type to run!'
      )	
      self.run_type = args[1]
   if (len(args)>2 and args[2]=='debug'):
      self.debug=True
   else:
      self.debug=False
#=============================================================================
   def set_derived_params(self)->None:
#-----------------------------------------------------------------------------
## put the R_max slightly inside the outer horizon.
## if near extremal limit then put R_max a m, which is in between
## the outer and inner horizons
      sqrt_term= pow(
         pow(self.black_hole_mass,2)
      -   pow(self.black_hole_spin,2)
      ,0.5)
      self.horizon=  self.black_hole_mass+sqrt_term

      self.R_max= float(
         pow(self.compactification_length,2)
         /self.horizon
      )
#-----------------------------------------------------------------------------
## Gauss points for integration
## want to exactly integrate polynomials of order
## 2l + 2m(i.e. lmin) + alpha + beta (so being a bit conservative here) 
      self.lmin= max(abs(self.pm_ang),abs(self.spin))
      self.ny= (self.nl
      +	int(abs(2*self.lmin))
      + int(abs(2*self.pm_ang+self.spin))
      +	int(abs(2*self.pm_ang-self.spin))
      )
      if (self.ny%2!=0):
         self.ny+= 1
#-----------------------------------------------------------------------------
## consider characteristic speeds: it looks like max is ~2/3 so can multiply
## the factor 6/N^2 by up to (3/2) (need to experiment)
#-----------------------------------------------------------------------------
      self.dt= float(
         8.*pow(max(self.nx,self.ny),-2)
      )
#-----------------------------------------------------------------------------
      self.nt= int(
         self.evolve_time*self.black_hole_mass/self.dt
      )
#-----------------------------------------------------------------------------
      self.t_step_save= int(
         self.nt/float(self.num_saved_times)
      )	
      if (self.t_step_save==0):
         self.t_step_save= 1
#=============================================================================
   def make_output_dir(self)->None:
      self.output_dir= str(
         'output/'
      +	'_'.join('_'.join(time.asctime().split(' ')).split(':'))
      +	'_a'+str(self.black_hole_spin/self.black_hole_mass)
      +	'_nx'+str(self.nx)
      +	'_ny'+str(self.ny)
      +	'_nl'+str(self.nl)
      +	'_s'+str(self.spin)
      +	'_lpm'+str(self.l_ang_pm)
      +	'_lnm'+str(self.l_ang_nm)
      +	'_pm'+str( self.pm_ang)
      )
      os.makedirs(self.output_dir)
#=============================================================================
   def make_swaL_dir(self)->None:
      self.swaL_dir= self.output_dir+"/swaL_tables"
      os.makedirs(self.swaL_dir)
#=============================================================================
   def write_sim_params(self)->None:
      with open(self.output_dir+'/sim_params.txt','w') as f:
         attrs= vars(self)
         for param in attrs:
            f.write('{} {}\n'.format(param,attrs[param]))	
#=============================================================================
   def compile(self)->None:
      subprocess.call('make '+self.bin,shell='True')
#=============================================================================
   def launch_run(self)->None:
      if (self.computer=='home'):
         output_file= self.output_dir+'/output.txt'
         run_str= (
            './bin/ '+self.bin+self.home_dir+'/'+self.output_dir
         +  '>'+output_file+' 2>&1 &'
         )
         if (self.debug):
            run_str= 'valgrind -v --track-origins=yes --leak-check=full '+run_str
            subprocess.call(run_str,shell='True') 
      else:
         raise ValueError('computer= '+self.computer+' not yet supported')
#=============================================================================
   def make_Legendre_pts(self)->None:
      save_roots_weights_Legendre(
         self.swaL_dir,
         self.ny
      )
#=============================================================================
   def make_Gauss_pts(self)->None:
      for spin in [-3,-2,-1,0,1,2]:
         for m_ang in [-2*self.pm_ang,-self.pm_ang,0,self.pm_ang,2*self.pm_ang]:
            save_Gauss_quad_vals_swaL(self.swaL_dir,spin,m_ang,self.nl,self.ny) 
