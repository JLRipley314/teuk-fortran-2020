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
sim.evolve_time= float(100) ## units of black hole mass
sim.num_saved_times= int(200)
#=============================================================================
sim.nx= int(pow(2,4)*pow(3,1)*pow(5,0)*pow(7,0)) ## num radial pts 
sim.nl= int(pow(2,2)*pow(3,1)*pow(5,0)*pow(7,0)) ## num swaL angular pts 
#=============================================================================
## evolution and write: take boolean values 
#=============================================================================
sim.metric_recon= True
sim.scd_order=    True

sim.write_indep_res=           True
sim.write_metric_recon_fields= True
sim.write_scd_order_source=    True
sim.write_coefs=               False
#=============================================================================
## change start time
sim.start_multiple= float(1.0)
#=============================================================================
sim.computer= 'feynman'#'home'#
#=============================================================================
## for feynman cluster/slurm script
sim.walltime= '4:00:00' ## (hh:mm:ss)
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
## initial data for mode m1
#=============================================================================
sim.pm1_ang =  int(0) ## m_ang is preserved by time evolution
#-----------------------------------------------------------------------------
sim.l_ang_pm1= int(2) ## l_ang: support of initial swal 

sim.initial_data_direction_pm1= "ingoing"#"outgoing"#"time_symmetric"#

sim.amp_pm1= float(10.0)  ## amplitude of the initial perturbation
sim.rl_pm1=  float(-1.5)  ## compact support: lower r value
sim.ru_pm1=  float( 1.5)  ## compact support: upper r value 
#-----------------------------------------------------------------------------
sim.l_ang_nm1= int(2) ## support over single spin weighted spherical harmonic

sim.initial_data_direction_nm1= "ingoing"#"time_symmetric"#"outgoing"#

sim.amp_nm1= float( 0.0)  ## amplitude of the initial perturbation
sim.rl_nm1=  float(-1.5)  ## compact support: lower r value
sim.ru_nm1=  float( 1.5)  ## compact support: upper r value 
#=============================================================================
## initial data for mode m2
#=============================================================================
sim.pm2_ang =  int(2) ## m_ang is preserved by time evolution
#-----------------------------------------------------------------------------
sim.l_ang_pm2= int(3) ## support over single spin weighted spherical harmonic

sim.initial_data_direction_pm2= "time_symmetric"#"outgoing"#"ingoing"#

sim.amp_pm2= float( 0.0)  ## amplitude of the initial perturbation
sim.rl_pm2=  float(-1.5)  ## compact support: lower r value
sim.ru_pm2=  float( 1.5)  ## compact support: upper r value 

sim.l_ang_pm2= int(3)     ## support over single spin weighted spherical harmonic
#-----------------------------------------------------------------------------
sim.l_ang_nm2= int(2) ## support over single spin weighted spherical harmonic

sim.initial_data_direction_nm2= "time_symmetric"#"outgoing"#"ingoing"#

sim.amp_nm2= float( 0.0)  ## amplitude of the initial perturbation
sim.rl_nm2=  float(-1.5)  ## compact support: lower r value
sim.ru_nm2=  float( 1.5)  ## compact support: upper r value 
#=============================================================================
## which m angular values to evolve
#=============================================================================
sim.lin_m = [  -sim.pm1_ang,   sim.pm1_ang]
sim.scd_m = [-2*sim.pm1_ang, 2*sim.pm1_ang, 0]
#=============================================================================
## which m angular values to write to file
#=============================================================================
sim.lin_write_m = [  -sim.pm1_ang,   sim.pm1_ang]
sim.scd_write_m = [-2*sim.pm1_ang, 2*sim.pm1_ang, 0]
#=============================================================================
if (sim.run_type == "basic_run"):
   sim.launch_run()
#=============================================================================
elif (sim.run_type == "multiple_runs"):
   nxs = [48, 64, 80]
   nls = [12, 16, 20]

   bhss = [round(0.7*sim.black_hole_mass,4), round(0.998*sim.black_hole_mass,4)]
   sms = [1,2,3]

   for bhs in bhss:
      for i in range(len(nxs)):
         sim.black_hole_spin = bhs 
         sim.nx = nxs[i]
         sim.nl = nls[i]
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
