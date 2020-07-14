!
! automatically generated from sim_class.py
!
module mod_params
use mod_prec
implicit none
   character(*), parameter :: home_dir = '/home/jripley/teuk-fortran'
   character(*), parameter :: run_type = 'basic_run'
   character(*), parameter :: computer = 'home'
   character(*), parameter :: output_dir = 'output/Tue_Jul_14_15_16_32_2020_a0.7_nx64_ny28_nl16_s-2_lpm2_lnm2_pm2'
   integer(ip), parameter :: num_saved_times = 50_ip
   character(*), parameter :: initial_data_direction = 'ingoing'
   character(*), parameter :: initial_data_type = 'compact_bump'
   integer(ip), parameter :: l_ang_pm = 2_ip
   integer(ip), parameter :: l_ang_nm = 2_ip
   integer(ip), parameter :: pm_ang = 2_ip
   integer(ip), parameter :: spin = -2_ip
   integer(ip), parameter :: nx = 64_ip
   integer(ip), parameter :: nl = 16_ip
   character(*), parameter :: save_indep_res = 'true'
   character(*), parameter :: save_coefs = 'true'
   character(*), parameter :: save_metric = 'true'
   character(*), parameter :: save_source = 'true'
   integer(ip), parameter :: lmin = 2_ip
   integer(ip), parameter :: ny = 28_ip
   integer(ip), parameter :: nt = 5120_ip
   integer(ip), parameter :: t_step_save = 102_ip
   character(*), parameter :: walltime = '168:00:00'
   character(*), parameter :: memory = '512'
   character(*), parameter :: tables_dir = '/home/jripley/teuk-fortran/output/swaL_tables'
   character(*), parameter :: output_stem = 'Tue_Jul_14_15_16_32_2020_a0.7_nx64_ny28_nl16_s-2_lpm2_lnm2_pm2'
   character(*), parameter :: bin = 'Tue_Jul_14_15_16_32_2020_a0.7_nx64_ny28_nl16_s-2_lpm2_lnm2_pm2.run'
end module mod_params
