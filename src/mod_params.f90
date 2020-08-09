!
! automatically generated from sim_class.py
!
module mod_params
use mod_prec
implicit none
   character(*), parameter :: home_dir = '/home/jripley/teuk-fortran'
   character(*), parameter :: run_type = 'basic_run'
   logical, parameter :: debug = .false.
   real(rp), parameter :: black_hole_mass = 0.5_rp
   real(rp), parameter :: black_hole_spin = 0.35_rp
   real(rp), parameter :: compactification_length = 1.0_rp
   real(rp), parameter :: evolve_time = 30.0_rp
   integer(ip), parameter :: num_saved_times = 150_ip
   integer(ip), parameter :: nx = 64_ip
   integer(ip), parameter :: nl = 12_ip
   logical, parameter :: metric_recon = .true.
   logical, parameter :: scd_order = .false.
   logical, parameter :: write_indep_res = .true.
   logical, parameter :: write_metric_recon_fields = .true.
   logical, parameter :: write_scd_order_source = .false.
   real(rp), parameter :: start_multiple = 1.0_rp
   character(*), parameter :: computer = 'home'
   character(*), parameter :: walltime = '12:00:00'
   character(*), parameter :: memory = '512'
   character(*), parameter :: num_nodes = '1'
   character(*), parameter :: num_tasks_per_node = '1'
   integer(ip), parameter :: psi_spin = -2_ip
   integer(ip), parameter :: psi_boost = -2_ip
   integer(ip), parameter :: pm1_ang = 2_ip
   integer(ip), parameter :: l_ang_pm1 = 2_ip
   character(*), parameter :: initial_data_direction_pm1 = 'ingoing'
   real(rp), parameter :: amp_pm1 = 10.0_rp
   real(rp), parameter :: rl_pm1 = -1.5_rp
   real(rp), parameter :: ru_pm1 = 1.5_rp
   integer(ip), parameter :: l_ang_nm1 = 2_ip
   character(*), parameter :: initial_data_direction_nm1 = 'time_symmetric'
   real(rp), parameter :: amp_nm1 = 0.0_rp
   real(rp), parameter :: rl_nm1 = -1.5_rp
   real(rp), parameter :: ru_nm1 = 1.5_rp
   integer(ip), parameter :: pm2_ang = 1_ip
   integer(ip), parameter :: l_ang_pm2 = 3_ip
   character(*), parameter :: initial_data_direction_pm2 = 'time_symmetric'
   real(rp), parameter :: amp_pm2 = 0.0_rp
   real(rp), parameter :: rl_pm2 = -1.5_rp
   real(rp), parameter :: ru_pm2 = 1.5_rp
   integer(ip), parameter :: l_ang_nm2 = 2_ip
   character(*), parameter :: initial_data_direction_nm2 = 'time_symmetric'
   real(rp), parameter :: amp_nm2 = 0.0_rp
   real(rp), parameter :: rl_nm2 = -1.5_rp
   real(rp), parameter :: ru_nm2 = 1.5_rp
   integer(ip), parameter :: max_l = 11_ip
   real(rp), parameter :: horizon = 0.8570714214271424_rp
   real(rp), parameter :: R_max = 1.1667639067172042_rp
   real(rp), parameter :: scd_order_start_time = 4.810510846598895_rp
   integer(ip), parameter :: ny = 24_ip
   real(rp), parameter :: dt = 0.002197265625_rp
   integer(ip), parameter :: nt = 6826_ip
   integer(ip), parameter :: t_step_save = 45_ip
   integer(ip), parameter :: max_m = 4_ip
   integer(ip), parameter :: min_m = -4_ip
   integer(ip), parameter :: max_s = 3_ip
   integer(ip), parameter :: min_s = -3_ip
   character(*), parameter :: output_stem = 'Sun_12_28_bhm0.5_bhs0.35_nx64_ny24_nl12_s-2_pm12_pm21'
   character(*), parameter :: output_dir = 'output/Sun_12_28_bhm0.5_bhs0.35_nx64_ny24_nl12_s-2_pm12_pm21'
   character(*), parameter :: bin_name = 'Sun_12_28_bhm0.5_bhs0.35_nx64_ny24_nl12_s-2_pm12_pm21.run'
   character(*), parameter :: tables_dir = 'output/Sun_12_28_bhm0.5_bhs0.35_nx64_ny24_nl12_s-2_pm12_pm21/tables'
end module mod_params
