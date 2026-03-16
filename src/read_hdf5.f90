program read

  use module_set_tdata

  implicit none

  block
    use inputparser

    character(256) :: fn_hdf
    character(256) :: fn_para

    call getarg(1,fn_para)
    write(6,*) "parameter file = ", trim(fn_para)

    call get_string_parameter(fn_para,"fn_hdf",fn_hdf)
    
    call read_set_tdata_from_hdf5(fn_hdf)
  end block

  block
    integer :: ip
    real(8) :: t, x,y,z

    t = 0.5d0
    
    do ip=1,ntraj
       call get_coordinates_from_time(particles(ip), t, x, y, z)
       write(6,'(i5,99es12.4)') ip,x,y,z,sqrt(x*x+y*y+z*z)
    end do
  end block
  ! call traj_read_info_hdf
  ! ntraj = traj_get_ntraj()
  ! write(6,*) ntraj
  ! stop

  ! call traj_read_hdf(1000,p)

end program read
