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
   real(rp), parameter :: evolve_time = 20.0_rp
   integer(ip), parameter :: num_saved_times = 100_ip
   integer(ip), parameter :: nx = 64_ip
   integer(ip), parameter :: nl = 16_ip
   logical, parameter :: metric_recon = .true.
   logical, parameter :: scd_order = .true.
   logical, parameter :: constrained_evo = .true.
   logical, parameter :: write_indep_res = .true.
   logical, parameter :: write_metric_recon_fields = .false.
   logical, parameter :: write_scd_order_source = .true.
   logical, parameter :: write_coefs = .false.
   character(*), parameter :: computer = 'home'
   character(*), parameter :: feyn_out_stem = '/mnt/grtheory/tf-out/'
   character(*), parameter :: walltime = '168:00:00'
   character(*), parameter :: memory = '512'
   character(*), parameter :: email = 'lloydripley@gmail.com'
   integer(ip), parameter :: psi_spin = -2_ip
   integer(ip), parameter :: psi_boost = -2_ip
   real(rp), parameter :: start_multiple = 1.0_rp
   integer(ip), parameter :: pm1_ang = 2_ip
   integer(ip), parameter :: l_ang_pm1 = 2_ip
   character(*), parameter :: initial_data_direction_pm1 = 'ingoing'
   real(rp), parameter :: amp_pm1 = 0.1_rp
   real(rp), parameter :: rl_pm1_0 = -2.0_rp
   real(rp), parameter :: ru_pm1_0 = 2.0_rp
   integer(ip), parameter :: l_ang_nm1 = 2_ip
   character(*), parameter :: initial_data_direction_nm1 = 'ingoing'
   real(rp), parameter :: amp_nm1 = 0.0_rp
   real(rp), parameter :: rl_nm1_0 = 1.1_rp
   real(rp), parameter :: ru_nm1_0 = 3.0_rp
   integer(ip), parameter :: pm2_ang = 3_ip
   integer(ip), parameter :: l_ang_pm2 = 3_ip
   character(*), parameter :: initial_data_direction_pm2 = 'ingoing'
   real(rp), parameter :: amp_pm2 = 0.0_rp
   real(rp), parameter :: rl_pm2_0 = -1.5_rp
   real(rp), parameter :: ru_pm2_0 = 1.5_rp
   integer(ip), parameter :: l_ang_nm2 = 3_ip
   character(*), parameter :: initial_data_direction_nm2 = 'ingoing'
   real(rp), parameter :: amp_nm2 = 0.0_rp
   real(rp), parameter :: rl_nm2_0 = -1.5_rp
   real(rp), parameter :: ru_nm2_0 = 1.5_rp
   integer(ip), dimension(2), parameter :: lin_m = [2_ip,-2_ip]
   integer(ip), dimension(3), parameter :: scd_m = [0_ip,-4_ip,4_ip]
   integer(ip), dimension(2), parameter :: lin_write_m = [2_ip,-2_ip]
   integer(ip), dimension(3), parameter :: scd_write_m = [0_ip,-4_ip,4_ip]
   integer(ip), parameter :: max_l = 15_ip
   real(rp), parameter :: horizon = 0.8567143500057153_rp
   real(rp), parameter :: R_max = 1.167250204217227_rp
   real(rp), parameter :: rl_pm1 = -1.7134287000114305_rp
   real(rp), parameter :: ru_pm1 = 1.7134287000114305_rp
   real(rp), parameter :: rl_nm1 = 0.9423857850062869_rp
   real(rp), parameter :: ru_nm1 = 2.570143050017146_rp
   real(rp), parameter :: rl_pm2 = -1.285071525008573_rp
   real(rp), parameter :: ru_pm2 = 1.285071525008573_rp
   real(rp), parameter :: rl_nm2 = -1.285071525008573_rp
   real(rp), parameter :: ru_nm2 = 1.285071525008573_rp
   real(rp), parameter :: constraint_damping = 36.51483710615301_rp
   real(rp), parameter :: scd_order_start_time = 11.248163954718162_rp
   integer(ip), dimension(1), parameter :: lin_pos_m = [2_ip]
   integer(ip), parameter :: len_lin_pos_m = 1_ip
   integer(ip), parameter :: len_lin_m = 2_ip
   integer(ip), parameter :: len_scd_m = 3_ip
   integer(ip), parameter :: len_lin_write_m = 2_ip
   integer(ip), parameter :: len_scd_write_m = 3_ip
   integer(ip), parameter :: num_threads = 1_ip
   integer(ip), parameter :: ny = 30_ip
   real(rp), parameter :: dt = 0.002197265625_rp
   integer(ip), parameter :: nt = 4551_ip
   integer(ip), parameter :: t_step_save = 45_ip
   integer(ip), parameter :: max_m = 6_ip
   integer(ip), parameter :: min_m = -6_ip
   integer(ip), parameter :: max_s = 3_ip
   integer(ip), parameter :: min_s = -3_ip
   character(*), parameter :: output_stem = 'Tue_15_02_bhm0.5_bhs0.35_nx64_nl16_m_2_3'
   character(*), parameter :: output_dir = 'output/Tue_15_02_bhm0.5_bhs0.35_nx64_nl16_m_2_3'
   character(*), parameter :: bin_name = 'Tue_15_02_bhm0.5_bhs0.35_nx64_nl16_m_2_3.run'
   character(*), parameter :: tables_dir = 'output/Tue_15_02_bhm0.5_bhs0.35_nx64_nl16_m_2_3/tables'
end module mod_params
