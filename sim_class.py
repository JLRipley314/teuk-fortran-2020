#=============================================================================
import subprocess, os, sys, time, shutil
from typing import List 
from math import log

sys.path.insert(1,os.getcwd()+'/src/tables/')
#=============================================================================
from tables_cheb     import save_cheb
from tables_legendre import save_roots_weights_Legendre
from tables_swal     import save_Gauss_quad_vals_swaL 
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
   def make_output_dir(self)->None:
      time_of_day= time.asctime().split(' ')
      day = time_of_day[0]
      hr_min_sec = (
         time_of_day[4] 
         if time_of_day[2]=='' else
         time_of_day[3]
      )
      self.output_stem= str(
        day 
      + '_'+hr_min_sec.split(':')[0] 
      + '_'+hr_min_sec.split(':')[1] 
      +	'_bhm'+str(self.black_hole_mass)
      +	'_bhs'+str(self.black_hole_spin)
      +	'_nx'+str(self.nx)
      +	'_nl'+str(self.nl)
      +	'_m_'+str(self.pm1_ang)+'_'+str(self.pm2_ang)
      )
      if (self.computer=="home"):
         self.output_dir= "output/"+self.output_stem
      elif (self.computer=="feynman"):
         self.output_dir=self.feyn_out_stem+self.output_stem
      else:
         raise ValueError("self.computer="+self.computer+" not supported")
      os.makedirs(self.output_dir)
#=============================================================================
   def set_derived_params(self)->None:
      self.max_l = int(self.nl - 1)
#-----------------------------------------------------------------------------
## put the R_max at the location of the outer horizon.
## if near extremal limit then put R_max a m, which is in between
## the outer and inner horizons
      sqrt_term = pow(
         pow(self.black_hole_mass,2)
      -  pow(self.black_hole_spin,2)
      ,0.5)

      self.horizon= self.black_hole_mass+(0.999*sqrt_term)

      self.R_max= float(
         pow(self.compactification_length,2)
         /self.horizon
      )

      self.rl_pm1= self.horizon*self.rl_pm1_0
      self.ru_pm1= self.horizon*self.ru_pm1_0

      self.rl_nm1= self.horizon*self.rl_nm1_0
      self.ru_nm1= self.horizon*self.ru_nm1_0

      self.rl_pm2= self.horizon*self.rl_pm2_0
      self.ru_pm2= self.horizon*self.ru_pm2_0

      self.rl_nm2= self.horizon*self.rl_nm2_0
      self.ru_nm2= self.horizon*self.ru_nm2_0 
#-----------------------------------------------------------------------------
      absa = abs(self.black_hole_spin/self.black_hole_mass)
      self.constraint_damping = abs(
         (10.0/self.black_hole_mass)*pow(abs(1.000000001-absa),-0.5)
      )
#-----------------------------------------------------------------------------
## when to begin metric reconstruction
      self.scd_order_start_time = max(
            self.start_multiple*(
               (2.0/self.black_hole_mass)*(self.ru_nm1 - self.horizon)
            +  4.0*log(self.ru_nm1/self.horizon)
            ),
            self.start_multiple*(
               (2.0/self.black_hole_mass)*(self.ru_pm1 - self.horizon)
            +  4.0*log(self.ru_pm1/self.horizon)
            ),
            self.start_multiple*(
               (2.0/self.black_hole_mass)*(self.ru_nm2 - self.horizon)
            +  4.0*log(self.ru_nm2/self.horizon)
            ),
            self.start_multiple*(
               (2.0/self.black_hole_mass)*(self.ru_pm2 - self.horizon)
            +  4.0*log(self.ru_pm2/self.horizon)
            )
      )
#-----------------------------------------------------------------------------
## Gauss points for integration
## want to exactly integrate polynomials of order
## 2l + 2m(i.e. lmin) + alpha + beta (so being a bit conservative here) 
      lmin= max(abs(self.pm1_ang),abs(self.pm2_ang),abs(self.psi_spin))
      self.ny= (self.nl
      +	int(abs(2*lmin))
      + int(abs(2*self.pm1_ang+self.psi_spin))
      +	int(abs(2*self.pm1_ang-self.psi_spin))
      )
      if (self.ny%2!=0):
         self.ny+= 1
#-----------------------------------------------------------------------------
## consider characteristic speeds: it looks like max is ~2/3 so can multiply
## the factor 6/N^2 by up to (3/2) (need to experiment)
#-----------------------------------------------------------------------------
      self.dt= float(
         9.*pow(max(self.nx,self.ny),-2)
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
#-----------------------------------------------------------------------------
      self.max_m=  max(abs(2*self.pm1_ang),abs(2*self.pm2_ang),1)
      self.min_m= -max(abs(2*self.pm1_ang),abs(2*self.pm2_ang),1)
#-----------------------------------------------------------------------------
      self.max_s=  3
      self.min_s= -3
#-----------------------------------------------------------------------------
      assert(self.l_ang_nm1>=max(abs(self.pm1_ang),abs(self.psi_spin)))
      assert(self.l_ang_pm1>=max(abs(self.pm1_ang),abs(self.psi_spin)))
      assert(self.l_ang_nm2>=max(abs(self.pm2_ang),abs(self.psi_spin)))
      assert(self.l_ang_pm2>=max(abs(self.pm2_ang),abs(self.psi_spin)))
#=============================================================================
   def make_tables_dir(self)->None:
      self.tables_dir= self.output_dir+"/tables"
      os.makedirs(self.tables_dir)
#=============================================================================
   def write_sim_params(self)->None:
      with open(self.output_dir+'/sim_params.txt','w') as f:
         attrs= vars(self)
         for param in attrs:
            f.write('{} {}\n'.format(param,attrs[param]))	
#=============================================================================
   def write_slurm_script(self):
      with open('{}/run.slurm'.format(self.home_dir), 'w') as f:
         f.write('#!/bin/sh\n')
         f.write('#SBATCH -J fteuk\t\t# job name\n')
         f.write('#SBATCH -t {}\t\t# walltime (dd:hh:mm:ss)\n'.format(self.walltime))
         f.write('#SBATCH -p dept\t\t# partition/queue name\n')
         f.write('#SBATCH --mem={}MB\t\t# memory in MB\n'.format(self.memory))
         f.write('#SBATCH --output={}\t\t# file for STDOUT\n'.format(self.output_file))
         f.write('#SBATCH --mail-user=jripley@princeton.edu\t\t# Mail  id of the user\n')
         #------------
         ## for openmp
         #------------
         f.write('#SBATCH -N 1\t\t# nodes= 1\n')
         f.write('#SBATCH -c {}\n'.format(self.num_threads))

         f.write('if [ -n "$SLURM_CPUS_PER_TASK" ]; then\n')
         f.write('  omp_threads=$SLURM_CPUS_PER_TASK\n')
         f.write('else\n')
         f.write('  omp_threads=1\n')
         f.write('fi\n')

         f.write('export OMP_NUM_THREADS=$omp_threads')
         #------------
         ## executable
         #------------
         run_str= './bin/{} {}\n\n'.format(self.bin_name, self.output_dir)
         if (self.debug):
            run_str= 'valgrind -v --track-origins=yes --leak-check=full '+run_str
         f.write('\n'+run_str)

         shutil.copyfile(
         '{}/run.slurm'.format(self.home_dir),
         '{}/run.slurm'.format(self.output_dir)
         )
#=============================================================================
   def compile(self)->None:
      subprocess.call('make '+self.bin_name,shell=True)
#=============================================================================
   def launch_run(self)->None:
      self.set_derived_params()

      self.make_output_dir()
      self.bin_name= self.output_stem+'.run'

      self.make_tables_dir()
      save_cheb(self.tables_dir,self.nx)
      save_roots_weights_Legendre(self.tables_dir,self.ny)
      self.make_Gauss_pts()

      self.write_sim_params()
      self.write_mod_params()

      self.output_file= self.output_dir+'/output.txt'
      if (self.computer=='home'):
         run_str= (
            'time ./bin/'+self.bin_name+' > '+self.output_file+' 2>&1 &'
         )
         if (self.debug):
            run_str= 'valgrind -v --track-origins=yes --leak-check=full '+run_str
         os.environ['OMP_NUM_THREADS']= str(self.num_threads)
         subprocess.call(run_str,shell=True) 
      elif (self.computer=='feynman'):
         self.write_slurm_script()
         subprocess.call('sbatch run.slurm', shell='True')		
      else:
         raise ValueError('computer= '+self.computer+' not yet supported')
#=============================================================================
   def make_Gauss_pts(self)->None:
      for spin in [-3,-2,-1,0,1,2,3]:
         for m_ang in range(self.min_m,self.max_m+1):
            save_Gauss_quad_vals_swaL(self.tables_dir,spin,m_ang,self.nl,self.ny) 
#=============================================================================
## write mod_params.f90 and recompile so everything is a module
## will probably rewrite at some point...
   def write_mod_params(self)->None:
      attrs= vars(self)
      with open('src/mod_params.f90','w') as f:
         f.write('!\n')
         f.write('! automatically generated from sim_class.py\n')
         f.write('!\n')
         f.write('module mod_params\n')
         f.write('use mod_prec\n')
         f.write('implicit none\n')
         attrs= vars(self)
         for param in attrs:
            if (type(attrs[param])==str):
               f.write("   character(*), parameter :: {} = '{}'\n".format(param,attrs[param]))

            if (type(attrs[param])==int):
               f.write("   integer(ip), parameter :: {} = {}_ip\n".format(param,attrs[param]))

            if (type(attrs[param])==bool):
               if attrs[param]==True:
                  f.write("   logical, parameter :: {} = .true.\n".format(param))
               if attrs[param]==False:
                  f.write("   logical, parameter :: {} = .false.\n".format(param))

            if (type(attrs[param])==float):
               f.write("   real(rp), parameter :: {} = {}_rp\n".format(param,attrs[param]))

            if (type(attrs[param])==complex):
               f.write("   complex(rp), parameter :: {} = ({}_rp,{}_rp)\n".format(param,attrs[param].real,attrs[param].imag))

            if ((type(attrs[param])==list)
            and  type(attrs[param][0])==int
            ):
               l= list(set(attrs[param])) ## to get rid of duplicates in list 
               n= len(l)
               lstr= "[{}_ip".format(l[0])
               for item in l[1:]:
                  lstr+= ",{}_ip".format(item)
               lstr+= "]"
               f.write("   integer(ip), dimension({}), parameter :: {} = {}\n".format(n,param,lstr))

         f.write('end module mod_params\n')
      subprocess.call('make clean_obj',shell=True)
      subprocess.call('make '+self.bin_name,shell=True)
