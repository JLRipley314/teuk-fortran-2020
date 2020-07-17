!
! automatically generated from sim_class.py
!
module mod_params
use mod_prec
implicit none
   character(*), parameter :: home_dir = '/home/jripley/teuk-fortran'
   character(*), parameter :: run_type = 'multiple_runs'
   character(*), parameter :: computer = 'home'
   character(*), parameter :: output_dir = 'output/Fri_Jul_17_13_19_06_2020_a0.7_nx80_ny28_nl16_s-2_lpm2_lnm2_pm2'
   real(rp), parameter :: black_hole_mass = 0.5_rp
   real(rp), parameter :: black_hole_spin = 0.35_rp
   real(rp), parameter :: compactification_length = 1.0_rp
   real(rp), parameter :: evolve_time = 100.0_rp
   integer(ip), parameter :: num_saved_times = 500_ip
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
   integer(ip), parameter :: nx = 80_ip
   integer(ip), parameter :: nl = 16_ip
   character(*), parameter :: save_indep_res = 'true'
   character(*), parameter :: save_coefs = 'true'
   character(*), parameter :: save_metric = 'true'
   character(*), parameter :: save_source = 'true'
   real(rp), parameter :: start_multiple = 1.0_rp
   character(*), parameter :: walltime = '168:00:00'
   character(*), parameter :: memory = '512'
   integer(ip), parameter :: ny = 28_ip
   integer(ip), parameter :: lmax = 15_ip
   real(rp), parameter :: horizon = 0.8570714214271424_rp
   real(rp), parameter :: R_max = 1.1667639067172042_rp
   integer(ip), parameter :: lmin = 2_ip
   real(rp), parameter :: dt = 0.00125_rp
   integer(ip), parameter :: nt = 40000_ip
   integer(ip), parameter :: t_step_save = 80_ip
   integer(ip), parameter :: max_m = 2_ip
   integer(ip), parameter :: min_m = -2_ip
   integer(ip), parameter :: max_s = 2_ip
   integer(ip), parameter :: min_s = -3_ip
   character(*), parameter :: output_stem = 'Fri_Jul_17_13_19_06_2020_a0.7_nx80_ny28_nl16_s-2_lpm2_lnm2_pm2'
   character(*), parameter :: bin = 'Fri_Jul_17_13_19_06_2020_a0.7_nx80_ny28_nl16_s-2_lpm2_lnm2_pm2.run'
   character(*), parameter :: tables_dir = 'output/Fri_Jul_17_13_19_06_2020_a0.7_nx80_ny28_nl16_s-2_lpm2_lnm2_pm2/tables'
end module mod_params
