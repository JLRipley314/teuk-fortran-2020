#!/usr/bin/env python

import sys, time
from sim_class import Sim
#=============================================================================
args= sys.argv
#=============================================================================
## input paramters: set these by hand 
#=============================================================================
sim= Sim(args)
#-----------------------------------------------------------------------------
sim.computer= 'home'#'feynman'#
#-----------------------------------------------------------------------------
sim.black_hole_mass= float(0.5)	
sim.black_hole_spin= float(0.35)
sim.compactification_length= float(1)
#-----------------------------------------------------------------------------
## evolve time: in units of black hole mass
#-----------------------------------------------------------------------------
sim.evolve_time= float(30.0)
sim.num_saved_times= int(100)
#-----------------------------------------------------------------------------
## initial data
#-----------------------------------------------------------------------------
sim.initial_data_direction= str("ingoing")#str("time_symmetric")#
sim.initial_data_type= str("compact_bump")
sim.amp_pm= float(10.0)  ## amplitude of the initial perturbation
sim.rl_pm=  float(-1.5)  ## compact support: lower r value
sim.ru_pm=  float( 1.5)  ## compact support: upper r value 
sim.l_ang_pm= int(2)     ## support over single spin weighted spherical harmonic

sim.amp_nm= float( 0.0)  ## amplitude of the initial perturbation
sim.rl_nm=  float(-1.5)  ## compact support: lower r value
sim.ru_nm=  float( 1.5)  ## compact support: upper r value 
sim.l_ang_nm= int(2)     ## support over single spin weighted spherical harmonic
#-----------------------------------------------------------------------------
##  Teukolsky equation preserves m 
sim.pm_ang= int(2)
assert(sim.pm_ang>=0)
#-----------------------------------------------------------------------------
## psi_4 is spin -2, psi_0 is spin +2 (code only reconstructs for psi_4) 
sim.spin= int(-2)
#-----------------------------------------------------------------------------
sim.nx= int(pow(2,4)*pow(3,1)*pow(5,0)*pow(7,0)) ## num radial pts 
sim.nl= int(pow(2,2)*pow(3,0)*pow(5,1)*pow(7,0)) ## num swaL angular pts 
#-----------------------------------------------------------------------------
## further diagnostics
sim.save_indep_res= "true"#"false"#
sim.save_coefs= "true"#"false"#
sim.save_metric= "true"#"false"#
sim.save_source= "true"#"false"#
#-----------------------------------------------------------------------------
## change start time
sim.start_multiple= float(1.0)
#=============================================================================
## for feynman cluster/slurm script
sim.walltime= '12:00:00' ### (hh:mm:ss)
sim.memory=  '512' ### MB 
sim.num_nodes= '1'
sim.num_tasks_per_node= '1'
#=============================================================================
if (sim.run_type == "basic_run"):
   sim.launch_run()
#=============================================================================
elif (sim.run_type == "multiple_runs"):
   m_vals = [0,1,2]

   sim.nx = 64
   sim.ny = 12
   for m in m_vals:
      sim.pm_ang = m
      sim.launch_run() 
      time.sleep(120)

#   sim.nx = 80
#   sim.ny = 16
#   for m in m_vals:
#      sim.pm_ang = m
#      sim.launch_run() 
#      time.sleep(180)
#=============================================================================
elif (sim.run_type == "convergence_test"):
   nx_vals= [64, 80, 96]
   nl_vals= [16, 20, 24]

   for i in range(len(nx_vals)):
      sim.nx= nx_vals[i]
      sim.nl= nl_vals[i]

      sim.launch_run()

      time.sleep(1)
#=============================================================================
else:
   raise ValueError("run_type = "+str(sim.run_type)) 
