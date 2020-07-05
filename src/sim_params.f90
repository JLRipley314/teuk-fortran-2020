module mod_sim_params
!-----------------------------------------------------------------------------
    use, intrinsic :: iso_fortran_env, only: &
        ip => int64, rp => real64

    implicit none
!-----------------------------------------------------------------------------
    integer(ip), parameter :: nt = 10_ip
    integer(ip), parameter :: nx = 4_ip 
    integer(ip), parameter :: nl = 4_ip 
    integer(ip), parameter :: ny = 4_ip
    integer(ip), parameter :: num_t_save = 10_ip

    integer(ip), parameter :: spin = -2_ip

    real(rp), parameter :: start_multiple = 1.0_rp

    logical, parameter :: save_coefs = .false.
    logical, parameter :: save_indep_res = .false.
    logical, parameter :: save_metric = .false.
    logical, parameter :: save_source = .false.

    real(rp), parameter :: bh_mass = 0.5_rp
    real(rp), parameter :: bh_spin = 0.0_rp
    real(rp), parameter :: com_len = 1_rp
    !-----------------
    ! for initial data
    !-----------------
    complex(rp), parameter :: amp_pm = 0_rp
    complex(rp), parameter :: amp_nm = 0_rp
    complex(rp), parameter ::  rl_pm = 0_rp
    complex(rp), parameter ::  rl_nm = 0_rp
    complex(rp), parameter ::  ru_pm = 0_rp
    complex(rp), parameter ::  ru_nm = 0_rp

    integer(ip), parameter :: pm_ang = 2_rp
    integer(ip), parameter :: l_ang_pm = 2_rp
    integer(ip), parameter :: l_ang_nm = 2_rp

    character(*), parameter :: initial_data_direction = 'ingoing'
    character(*), parameter :: initial_data_type = 'bump'
    !-------------------
    ! derived parameters
    !-------------------
    integer(ip), parameter :: t_step_save = nt / num_t_save

    real(rp), parameter :: dt = 6_rp / ((nx ** 2) * (ny ** 2))

    real(rp), parameter :: r_max = 1_rp 
    real(rp), parameter :: dx_over_dr = 1_rp

    real(rp), parameter :: r_pts(nx) = [ 0, 0, 0, 0] 
    real(rp), parameter :: y_pts(ny) = [ 0, 0, 0, 0] 
    real(rp), parameter :: cs(ny) = cos(y_pts)
    real(rp), parameter :: sn(ny) = sin(y_pts)
!-----------------------------------------------------------------------------
end module mod_sim_params
