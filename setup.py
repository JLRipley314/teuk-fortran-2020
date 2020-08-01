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
sim.evolve_time= float(60)
sim.num_saved_times= int(200)
#-----------------------------------------------------------------------------
## initial data
#-----------------------------------------------------------------------------
sim.initial_data_direction= str("time_symmetric")#str("outgoing")#str("ingoing")#
sim.initial_data_type= str("compact_bump")
sim.amp_pm= float( 0.1)  ## amplitude of the initial perturbation
sim.rl_pm=  float(-1.5)  ## compact support: lower r value
sim.ru_pm=  float( 1.5)  ## compact support: upper r value 
sim.l_ang_pm= int(2)     ## support over single spin weighted spherical harmonic

sim.amp_nm= float(10.0)  ## amplitude of the initial perturbation
sim.rl_nm=  float(-1.5)  ## compact support: lower r value
sim.ru_nm=  float( 1.5)  ## compact support: upper r value 
sim.l_ang_nm= int(2)     ## support over single spin weighted spherical harmonic
#-----------------------------------------------------------------------------
##  Teukolsky equation preserves m 
sim.pm_ang= int(2)
assert(sim.pm_ang>=0)
#-----------------------------------------------------------------------------
## we can only do metric reconstruction starting from psi4 for now.
## For pure first order Teukolsky evolution we can consider other
## spin weighted fields though.
## psi4 is spin -2, boost -2
## psi3 is spin -1, boost -1
## psi2 is spin  0, boost  0 
sim.psi_spin=  int(-2)
sim.psi_boost= int(-2)
#-----------------------------------------------------------------------------
sim.nx= int(pow(2,4)*pow(3,1)*pow(5,0)*pow(7,0)) ## num radial pts 
sim.nl= int(pow(2,4)*pow(3,0)*pow(5,0)*pow(7,0)) ## num swaL angular pts 
#-----------------------------------------------------------------------------
sim.metric_recon= True#False# 
sim.scd_order=    True#False# 

sim.write_indep_res=           True#False#
sim.write_metric_recon_fields= False#True#
sim.write_scd_order_source=    True#False#
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
elif (sim.run_type == "spin_ramp"):
   for bhs in [0,0.01,0.02,0.04,0.08,0.12,0.16,0.2,0.24,0.28,0.32]:
      sim.black_hole_spin= bhs
      sim.launch_run()

      time.sleep(60)
#=============================================================================
else:
   raise ValueError("run_type = "+str(sim.run_type)) 
