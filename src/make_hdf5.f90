! This code generates .h5 file of heating rates, compositions, and so on.
! (C) Sho Fujibayashi
program main

  use module_trajectory
  use module_tdata

  implicit  none
  
  type(traj) :: p
  integer :: it,ntraj

  ! dir_traj = "/scratch/sfujibayashi/Particle_trace_data/data_3D_SFHo_120_150-150mstg_Fugaku/data_sk1_th24_r1.0e+09_test"
  ! dir_traj = "/scratch/sfujibayashi/DD2Tim326_135_135_0028_12.5mstg_B15.5_HLLD_CT_GS/Analysis_ptr/data"

  block
    use inputparser

    character(256) :: dir_traj
    character(256) :: fn_para

    call getarg(1,fn_para)
    write(6,*) "parameter file = ", trim(fn_para)

    call get_string_parameter(fn_para,"dir_traj",dir_traj)

    call traj_set_dir(dir_traj)
  end block


  call traj_read_info
  !call traj_allocate
  !call traj_show_info
  
  call traj_make_hdf5

  stop

  call traj_read_info_hdf
  ntraj = traj_get_ntraj()
  write(6,*) ntraj
  stop

  call traj_read_hdf(1000,p)

  write(6,*) get_ye_5gk(p)
  write(6,*) get_crossing_time_radius(p,1d8)
  write(6,*) get_mass(p)
  ! do it=1,p%ntime
  !    write(6,'(i6,99es12.4)') it,p%time(it),p%tem(it),p%ye(it), sqrt(p%x(it)**2 + p%y(it)**2 + p%z(it)**2)
  ! enddo

end program main
