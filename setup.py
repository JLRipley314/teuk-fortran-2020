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
sim.black_hole_spin= round(0.7*sim.black_hole_mass,4)
sim.compactification_length= float(1)
#=============================================================================
sim.evolve_time= float(125) ## units of black hole mass
sim.num_saved_times= int(1800)
#=============================================================================
sim.nx= 96 ## num radial pts 
sim.nl= 24  ## num angular values
#=============================================================================
## evolution and write: take boolean values 
#=============================================================================
sim.metric_recon= True
sim.scd_order=    True

sim.write_indep_res=           True
sim.write_metric_recon_fields= False
sim.write_scd_order_source=    True
sim.write_coefs=               False

sim.constrained_evo = False
#=============================================================================
sim.computer= 'feynman'#'home'#
sim.feyn_out_stem= '/mnt/grtheory/tf-out/'
#=============================================================================
## for feynman cluster/slurm script
sim.walltime= '24:00:00' ## (hh:mm:ss)
sim.memory=  '512' ## MB 
sim.num_nodes= '1'
sim.num_tasks_per_node= '1'
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
## initial data for mode m1
#=============================================================================
sim.pm1_ang =  int(2) ## m_ang is preserved by time evolution
#-----------------------------------------------------------------------------
sim.l_ang_pm1= int(2) ## l_ang: support of initial swal 

sim.initial_data_direction_pm1= "ingoing"#"outgoing"#"time_symmetric"#

sim.amp_pm1= float( 0.1)  ## amplitude of the initial perturbation

sim.rl_pm1_0= float( 1.1)  ## lower r value as a multiple of horizon
sim.ru_pm1_0= float( 3.0)  ## upper r value as a multiple of horizon 
#-----------------------------------------------------------------------------
sim.l_ang_nm1= int(2) ## support over single spin weighted spherical harmonic

sim.initial_data_direction_nm1= "ingoing"#"time_symmetric"#"outgoing"#

sim.amp_nm1= float( 0.1)  ## amplitude of the initial perturbation

sim.rl_nm1_0= float( 1.1)  ## lower r value as multiple of horizon
sim.ru_nm1_0= float( 3.0)  ## upper r value as multiple of horizon
#=============================================================================
## initial data for mode m2
#=============================================================================
sim.pm2_ang =  int(2) ## m_ang is preserved by time evolution
#-----------------------------------------------------------------------------
sim.l_ang_pm2= int(3) ## support over single spin weighted spherical harmonic

sim.initial_data_direction_pm2= "ingoing"#"time_symmetric"#"outgoing"#

sim.amp_pm2= float(0.0)  ## amplitude of the initial perturbation

sim.rl_pm2_0= float(-1.5)  ## lower r value as multiple of horizon 
sim.ru_pm2_0= float( 1.5)  ## upper r value as multiple of horizon 

sim.l_ang_pm2= int(3)     ## support over single spin weighted spherical harmonic
#-----------------------------------------------------------------------------
sim.l_ang_nm2= int(2) ## support over single spin weighted spherical harmonic

sim.initial_data_direction_nm2= "ingoing"#"time_symmetric"#"outgoing"#

sim.amp_nm2= float(0.0)  ## amplitude of the initial perturbation

sim.rl_nm2_0= float(-1.5)  ## lower r value as multiple of horizon 
sim.ru_nm2_0= float( 1.5)  ## upper r value as multiple of horizon
#=============================================================================
## which m angular values to evolve
#=============================================================================
sim.lin_m = [  -sim.pm1_ang,   sim.pm1_ang]
sim.scd_m = [-2*sim.pm1_ang, 2*sim.pm1_ang]
#=============================================================================
## which m angular values to write to file
#=============================================================================
sim.lin_write_m = [  -sim.pm1_ang,   sim.pm1_ang]
sim.scd_write_m = [-2*sim.pm1_ang, 2*sim.pm1_ang]
#=============================================================================
if (sim.run_type == "basic_run"):
   sim.launch_run()
#=============================================================================
elif (sim.run_type == "multiple_runs"):
#-----------------------------------------------------------------------------
   default_amp = 0.1
   default_sm  = 1

   default_nx_07 = 96
   default_nl_07 = 24

   default_nx_099 = 112
   default_nl_099 = 24

   default_nx_0998 = 112
   default_nl_0998 = 28
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
   sim.black_hole_spin= round(0.7*sim.black_hole_mass,4)
   sim.start_multiple = default_sm
   sim.amp_pm1 = default_amp
   sim.amp_nm1 = default_amp
   sim.nx = default_nx_07
   sim.nl = default_nl_07

   nxs = [80,  96, 112]
   nls = [20,  24,  28]
   for i in range(len(nxs)):
      sim.nx = nxs[i]
      sim.nl = nls[i]
      sim.launch_run() 
      time.sleep(60)
#-----------------------------------------------------------------------------
   sim.black_hole_spin= round(0.7*sim.black_hole_mass,4)
   sim.start_multiple = default_sm
   sim.amp_pm1 = default_amp
   sim.amp_nm1 = default_amp
   sim.nx = default_nx_07
   sim.nl = default_nl_07

   amps = [0.3, 0.5]
   for amp in amps:
      sim.amp_pm1= amp
      sim.amp_nm1= amp
      sim.launch_run() 
      time.sleep(60)
#-----------------------------------------------------------------------------
   sim.black_hole_spin= round(0.7*sim.black_hole_mass,4)
   sim.start_multiple = default_sm
   sim.amp_pm1 = default_amp
   sim.amp_nm1 = default_amp
   sim.nx = default_nx_07
   sim.nl = default_nl_07

   sms = [2,3]
   for sm in sms:
      sim.start_multiple= sm
      sim.launch_run()
      time.sleep(60)
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
   sim.black_hole_spin= round(0.998*sim.black_hole_mass,4)
   sim.start_multiple = default_sm
   sim.amp_pm1 = default_amp
   sim.amp_nm1 = default_amp
   sim.nx = default_nx_0998
   sim.nl = default_nl_0998

   nxs = [96, 112, 128]
   nls = [24,  28,  32]
   for i in range(len(nxs)):
      sim.nx = nxs[i]
      sim.nl = nls[i]
      sim.launch_run() 
      time.sleep(60)
#-----------------------------------------------------------------------------
   sim.black_hole_spin= round(0.998*sim.black_hole_mass,4)
   sim.start_multiple = default_sm
   sim.amp_pm1 = default_amp
   sim.amp_nm1 = default_amp
   sim.nx = default_nx_0998
   sim.nl = default_nl_0998

   amps = [0.3, 0.5]
   for amp in amps:
      sim.amp_pm1= amp
      sim.amp_nm1= amp
      sim.launch_run() 
      time.sleep(60)
#-----------------------------------------------------------------------------
   sim.black_hole_spin= round(0.998*sim.black_hole_mass,4)
   sim.start_multiple = default_sm
   sim.amp_pm1 = default_amp
   sim.amp_nm1 = default_amp
   sim.nx = default_nx_0998
   sim.nl = default_nl_0998

   sms = [2,3]
   for sm in sms:
      sim.start_multiple= sm
      sim.launch_run()
      time.sleep(60)

   sys.exit()
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
   sim.black_hole_spin= round(0.99*sim.black_hole_mass,4)
   sim.start_multiple = default_sm
   sim.amp_pm1 = default_amp
   sim.amp_nm1 = default_amp
   sim.nx = default_nx_099
   sim.nl = default_nl_099

   nxs = [96, 112, 128]
   nls = [20,  24,  28]
   for i in range(len(nxs)):
      sim.nx = nxs[i]
      sim.nl = nls[i]
      sim.launch_run() 
      time.sleep(60)
#-----------------------------------------------------------------------------
   sim.black_hole_spin= round(0.99*sim.black_hole_mass,4)
   sim.start_multiple = default_sm
   sim.amp_pm1 = default_amp
   sim.amp_nm1 = default_amp
   sim.nx = default_nx_099
   sim.nl = default_nl_099

   amps = [0.3, 0.5]
   for amp in amps:
      sim.amp_pm1= amp
      sim.amp_nm1= amp
      sim.launch_run() 
      time.sleep(60)
#-----------------------------------------------------------------------------
   sim.black_hole_spin= round(0.99*sim.black_hole_mass,4)
   sim.start_multiple = default_sm
   sim.amp_pm1 = default_amp
   sim.amp_nm1 = default_amp
   sim.nx = default_nx_099
   sim.nl = default_nl_099

   sms = [2,3]
   for sm in sms:
      sim.start_multiple= sm
      sim.launch_run()
      time.sleep(60)
#=============================================================================
elif (sim.run_type == "start_times"):
   sms = [1,2,3]

   for sm in sms:
      sim.start_multiple= sm
      sim.launch_run()
      time.sleep(60)
#=============================================================================
elif (sim.run_type == "spin_ramp"):
   for bhs in [0,0.01,0.02,0.04,0.08,0.12,0.16,0.2,0.24,0.28,0.32]:
      sim.black_hole_spin= bhs
      sim.launch_run()

      time.sleep(60)
#=============================================================================
else:
   raise ValueError("run_type = "+str(sim.run_type)) 
