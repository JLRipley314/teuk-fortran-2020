!
! automatically generated from sim_class.py
!
module mod_params
use mod_prec
implicit none
   character(*), parameter :: home_dir = '/home/jripley/teuk-fortran'
   character(*), parameter :: run_type = 'basic_run'
   character(*), parameter :: computer = 'home'
   character(*), parameter :: output_dir = 'output/Wed_Jul_15_13_59_09_2020_a0.7_nx64_ny28_nl16_s-2_lpm2_lnm2_pm2'
   real(rp), parameter :: black_hole_mass = 0.5_rp
   real(rp), parameter :: black_hole_spin = 0.35_rp
   real(rp), parameter :: compactification_length = 1.0_rp
   real(rp), parameter :: evolve_time = 20.0_rp
   integer(ip), parameter :: num_saved_times = 50_ip
   character(*), parameter :: initial_data_direction = 'ingoing'
   character(*), parameter :: initial_data_type = 'compact_bump'
   real(rp), parameter :: amp_pm = 10.0_rp
   real(rp), parameter :: rl_pm = -1.5_rp
   real(rp), parameter :: ru_pm = 1.5_rp
   integer(ip), parameter :: l_ang_pm = 2_ip
   real(rp), parameter :: amp_nm = 0.0_rp
   real(rp), parameter :: rl_nm = -1.5_rp
   real(rp), parameter :: ru_nm = 1.5_rp
   integer(ip), parameter :: l_ang_nm = 2_ip
   integer(ip), parameter :: pm_ang = 2_ip
   integer(ip), parameter :: spin = -2_ip
   integer(ip), parameter :: nx = 64_ip
   integer(ip), parameter :: nl = 16_ip
   character(*), parameter :: save_indep_res = 'true'
   character(*), parameter :: save_coefs = 'true'
   character(*), parameter :: save_metric = 'true'
   character(*), parameter :: save_source = 'true'
   real(rp), parameter :: start_multiple = 1.0_rp
   real(rp), parameter :: horizon = 0.8570714214271424_rp
   real(rp), parameter :: R_max = 1.1667639067172042_rp
   integer(ip), parameter :: lmin = 2_ip
   integer(ip), parameter :: ny = 28_ip
   real(rp), parameter :: dt = 0.001953125_rp
   integer(ip), parameter :: nt = 5120_ip
   integer(ip), parameter :: t_step_save = 102_ip
   integer(ip), parameter :: max_m = 2_ip
   integer(ip), parameter :: min_m = -2_ip
   integer(ip), parameter :: max_s = 2_ip
   integer(ip), parameter :: min_s = -3_ip
   character(*), parameter :: walltime = '168:00:00'
   character(*), parameter :: memory = '512'
   character(*), parameter :: output_stem = 'Wed_Jul_15_13_59_09_2020_a0.7_nx64_ny28_nl16_s-2_lpm2_lnm2_pm2'
   character(*), parameter :: bin = 'Wed_Jul_15_13_59_09_2020_a0.7_nx64_ny28_nl16_s-2_lpm2_lnm2_pm2.run'
   character(*), parameter :: tables_dir = 'output/Wed_Jul_15_13_59_09_2020_a0.7_nx64_ny28_nl16_s-2_lpm2_lnm2_pm2/tables'
end module mod_params