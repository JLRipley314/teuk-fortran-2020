#!/usr/bin/env python
#=============================================================================
## parameter file for evolution
## usage:
## ./setup.py [run_type] [debug]
## example:
## ./setup.py basic_run
#=============================================================================
import sys, time
from sim_class import Sim
#=============================================================================
args= sys.argv
sim= Sim(args)
#=============================================================================
sim.black_hole_mass= float(0.5)	
sim.black_hole_spin= round(0.0*sim.black_hole_mass,6)
sim.compactification_length= float(1)
#=============================================================================
sim.evolve_time= float(10) ## units of black hole mass
sim.num_saved_times= int(100)
#=============================================================================
sim.nx= 64 ## num radial pts 
sim.nl= 20  ## num angular values
#=============================================================================
## evolution and write: take boolean values 

sim.metric_recon=     True
sim.scd_order=        True
sim.constrained_evo = True

sim.write_indep_res=           True
sim.write_metric_recon_fields= False
sim.write_scd_order_source=    True
sim.write_coefs=               False
sim.write_sphere_coefs=        True 
#=============================================================================
sim.computer= 'home'#'della'#
sim.della_out_stem= '/tigress/jripley/tf-out/'

## for della cluster/slurm script

sim.walltime= '144:00:00' ## (hh:mm:ss)
sim.memory=  '2048' ## MB 
sim.email=  'lloydripley@gmail.com' ## for slurm notification
#=============================================================================
## we can only do metric reconstruction starting from psi4 for now.
## For pure first order Teukolsky evolution we can consider other
## spin weighted fields though.
## psi4 is spin -2, boost -2
## psi3 is spin -1, boost -1
## psi2 is spin  0, boost  0 

sim.psi_spin=  int(-2)
sim.psi_boost= int(-2)
#=============================================================================
## start multiple for second order metric evolution 

sim.start_multiple= float(1.0)
#=============================================================================
## Initial data
#=============================================================================
## p(m)n:                  +(-) m angular number
## l_ang:                  initial data is a particular swal function
## initial_data_direction: which way pulse is approximately "heading"
## amp_re(im)_pm:          initial amplitude of real/imaginary parts of psi4
## rl(ru)_pm_0:            lower(upper) bounds of initial data as a multiple
##                         of the black hole horizon
#=============================================================================
## initial data for mode m1
#=============================================================================
sim.pm1_ang =  int(2) 
#-----------------------------------------------------------------------------
sim.l_ang_pm1= int(2)

sim.initial_data_direction_pm1= "ingoing"#"outgoing"#"time_symmetric"#

sim.amp_re_pm1= float(0.1)
sim.amp_im_pm1= float(0.0)

sim.rl_pm1_0= float( 1.1)
sim.ru_pm1_0= float( 2.5)
#-----------------------------------------------------------------------------
sim.l_ang_nm1= int(2)

sim.initial_data_direction_nm1= "ingoing"#"time_symmetric"#"outgoing"#

sim.amp_re_nm1= float( 0.0)
sim.amp_im_nm1= float( 0.0)

sim.rl_nm1_0= float( 1.1)
sim.ru_nm1_0= float( 2.5)
#=============================================================================
## initial data for mode m2
#=============================================================================
sim.pm2_ang =  int(3) ## m_ang is preserved by time evolution
#-----------------------------------------------------------------------------
sim.l_ang_pm2= int(3)

sim.initial_data_direction_pm2= "ingoing"#"time_symmetric"#"outgoing"#

sim.amp_re_pm2= float(0.0)
sim.amp_im_pm2= float(0.0)

sim.rl_pm2_0= float(-1.5) 
sim.ru_pm2_0= float( 1.5) 
#-----------------------------------------------------------------------------
sim.l_ang_nm2= int(3)
sim.l_ang_nm2= int(3)

sim.initial_data_direction_nm2= "ingoing"#"time_symmetric"#"outgoing"#

sim.amp_re_nm2= float(0.0)
sim.amp_im_nm2= float(0.0)

sim.rl_nm2_0= float(-1.5) 
sim.ru_nm2_0= float( 1.5)
#=============================================================================
## which m angular values to evolve

sim.lin_m= [
   -sim.pm1_ang,
    sim.pm1_ang
]

sim.scd_m= [
   -2*sim.pm1_ang,
    2*sim.pm1_ang,
   0
]
#=============================================================================
## which m angular values to write to file

sim.lin_write_m= [
   -sim.pm1_ang,
    sim.pm1_ang
]

sim.scd_write_m= [
   -2*sim.pm1_ang,
    2*sim.pm1_ang,
   0
]
#=============================================================================
if (sim.run_type == "basic_run"):
   sim.launch_run()
#=============================================================================
elif (sim.run_type == "multiple_runs"):
#-----------------------------------------------------------------------------
   default_sm  = 1

   default_nx_07 = 176
   default_nl_07 = 32
#   default_nx_07 = 144
#   default_nl_07 = 24

   default_nx_0998 = 214
   default_nl_0998 = 44
#   default_nx_0998 = 160
#   default_nl_0998 = 28
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
   sim.black_hole_spin= round(0.7*sim.black_hole_mass,6)
   sim.start_multiple = default_sm
   sim.nx = default_nx_07
   sim.nl = default_nl_07

   nxs = [160, 176, 192]
   nls = [ 28,  32,  36]
#   nxs = [128, 144, 160]
#   nls = [ 20,  24,  28]
   for i in range(len(nxs)):
      sim.nx = nxs[i]
      sim.nl = nls[i]
      sim.launch_run() 
#-----------------------------------------------------------------------------
   sim.black_hole_spin= round(0.7*sim.black_hole_mass,6)
   sim.start_multiple = default_sm
   sim.nx = default_nx_07
   sim.nl = default_nl_07

   sms = [2,3]
   for sm in sms:
      sim.start_multiple= sm
      sim.launch_run()
#-----------------------------------------------------------------------------
   sys.exit()
#-----------------------------------------------------------------------------
   sim.black_hole_spin= round(0.99998*sim.black_hole_mass,6)
   sim.start_multiple = default_sm
   sim.nx = default_nx_0998
   sim.nl = default_nl_0998

   nxs = [208, 214, 230]
   nls = [ 40,  44,  48]
#   nxs = [144, 160, 176]
#   nls = [ 24,  28,  32]
   for i in range(len(nxs)):
      sim.nx = nxs[i]
      sim.nl = nls[i]
      sim.launch_run() 
#-----------------------------------------------------------------------------
   sim.black_hole_spin= round(0.99998*sim.black_hole_mass,6)
   sim.start_multiple = default_sm
   sim.nx = default_nx_0998
   sim.nl = default_nl_0998

   sms = [2,3]
   for sm in sms:
      sim.start_multiple= sm
      sim.launch_run()

#=============================================================================
elif (sim.run_type == "spin_ramp"):
   for bhs in [0,0.01,0.02,0.04,0.08,0.12,0.16,0.2,0.24,0.28,0.32]:
      sim.black_hole_spin= bhs
      sim.launch_run()
#=============================================================================
else:
   raise ValueError("run_type = "+str(sim.run_type)) 
