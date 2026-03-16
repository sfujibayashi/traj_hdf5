module module_tdata
  implicit none
  
  type traj
     ! data saved in the trajectory file
     integer :: ip
     integer :: ntime
     real(8) :: mass, ut1_final, hut_final
     real(8),allocatable :: time(:)
     real(8),allocatable :: x(:),y(:),z(:),&
          vlx (:),&
          vly (:),&
          vlz (:),&
          qrho(:),&
          tem (:),&
          ye  (:),&
          sen (:),&
          rne (:),&
          rae (:),&
          deptn(:),&
          depta(:),&
          ut  (:),&
          hhh (:),&
          eta_nue  (:),&
          eta_nub  (:),&
          b2  (:),&
          pres(:)
     
  end type traj

contains

  subroutine tdata_allocate(p,ntime_in)
    type(traj),intent(inout) :: p
    integer,intent(in) :: ntime_in

    p%ntime = ntime_in

    call tdata_deallocate(p)

    allocate(p%time(p%ntime),&
         p%x(p%ntime),p%y(p%ntime),p%z(p%ntime),&
         p%vlx (p%ntime),&
         p%vly (p%ntime),&
         p%vlz (p%ntime),&
         p%qrho(p%ntime),&
         p%tem (p%ntime),&
         p%ye  (p%ntime),&
         p%sen (p%ntime),&
         p%rne (p%ntime),&
         p%rae (p%ntime),&
         p%deptn(p%ntime),&
         p%depta(p%ntime),&
         p%ut  (p%ntime),&
         p%hhh (p%ntime),&
         p%eta_nue  (p%ntime),&
         p%eta_nub  (p%ntime),&
         p%b2  (p%ntime),&
         p%pres(p%ntime))
    
  end subroutine tdata_allocate

  subroutine tdata_deallocate(p)
    type(traj),intent(inout) :: p
    
    if(allocated(p%time))then
       
       deallocate(p%time,&
            p%x,p%y,p%z,&
            p%vlx ,&
            p%vly ,&
            p%vlz ,&
            p%qrho,&
            p%tem ,&
            p%ye  ,&
            p%sen ,&
            p%rne ,&
            p%rae ,&
            p%deptn,&
            p%depta,&
            p%ut  ,&
            p%hhh ,&
            p%eta_nue  ,&
            p%eta_nub  ,&
            p%b2  ,&
            p%pres)
    endif
           
  end subroutine tdata_deallocate

  real(8) function get_ye_5GK(p) result(ye_5gk)
    type(traj),intent(in) :: p

    integer :: it

    real(8) :: temp
    integer :: it_5gk
    
    real(8) :: s1,s0

    it_5gk = 0

    do it = 1,p%ntime-1
       temp = p%tem(it)

       ! time at which the particle crosses T = 5 GK
       if( 5.d9 <= p%tem(it) .and. 5.d9 > p%tem(it+1) ) then
          if(it_5gk==0)it_5gk=it
       endif
    
    enddo

    it = it_5gk
    if(it>0)then
       s1 = (5.d9-p%tem(it))/(p%tem(it+1)-p%tem(it))
       s0 = 1.d0-s1
       ye_5gk= s1* p%ye(it+1) + s0* p%ye(it)
    else
       ye_5gk   = 0.d0
    endif
    
  end function get_ye_5GK
  
  real(8) function get_crossing_time_radius(p,r_ext) result(t)
    type(traj),intent(in) :: p
    real(8),intent(in) :: r_ext

    integer :: it

    real(8) :: r,rp
    integer :: it_cross
    
    real(8) :: s1,s0

    it_cross = 0

    do it = 1,p%ntime-1
       r  = sqrt(p%x(it)**2 + p%y(it)**2 + p%z(it)**2)
       rp = sqrt(p%x(it+1)**2 + p%y(it+1)**2 + p%z(it+1)**2)
       
       if( r <= r_ext .and. r_ext < rp ) then
          it_cross=it
       endif
    enddo

    it=it_cross
    if(it>0)then
       r  = sqrt(p%x(it)**2 + p%y(it)**2 + p%z(it)**2)
       rp = sqrt(p%x(it+1)**2 + p%y(it+1)**2 + p%z(it+1)**2)

       s1 = (r_ext-r)/(rp-r)
       s0 = 1.d0-s1
       t = s1* p%time(it+1) + s0*p%time(it)
    else
       t = 0d0
    endif
    
  end function get_crossing_time_radius

  real(8) function get_mass(p) result(m)
    type(traj),intent(in) :: p
    m=p%mass
  end function get_mass

  subroutine get_coordinates_from_time(p, t, x, y, z)
    type(traj),intent(in) :: p
    real(8),intent(in) :: t
    real(8),intent(out) :: x,y,z

    if (p%ntime <= 0) then
       x = 0.d0
       y = 0.d0
       z = 0.d0
       return
    elseif (p%ntime == 1) then
       x = p%x(1)
       y = p%y(1)
       z = p%z(1)
       return
    endif
    
    if    (t <= p%time(1))then

       x = p%x(1)
       y = p%y(1)
       z = p%z(1)
       return

    elseif(t >= p%time(p%ntime))then

       x = p%x(p%ntime)
       y = p%y(p%ntime)
       z = p%z(p%ntime)
       return

    else

       ! block
       !   integer :: it
       !   real(8) :: w0, w1
       !   do it = 2, p%ntime
       !      if (p%time(it) >= t) then
       !         w1 = (t - p%time(it-1)) / (p%time(it) - p%time(it-1))
       !         w0 = 1.d0 - w1
               
       !         x = w0*p%x(it-1) + w1*p%x(it)
       !         y = w0*p%y(it-1) + w1*p%y(it)
       !         z = w0*p%z(it-1) + w1*p%z(it)
       !         return
       !      endif
       !   enddo
       ! end block

       block
         integer :: lo, hi, mid
         real(8) :: w0, w1, t0, t1

         lo = 1
         hi = p%ntime
         
         do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (p%time(mid) <= t) then
               lo = mid
            else
               hi = mid
            endif
            ! write(6,*) lo,hi,p%time(lo),p%time(hi)
         enddo
         
         t0 = p%time(lo)
         t1 = p%time(hi)
         w1 = (t - t0) / (t1 - t0)
         w0 = 1.d0 - w1
         
         x = w0*p%x(lo) + w1*p%x(hi)
         y = w0*p%y(lo) + w1*p%y(hi)
         z = w0*p%z(lo) + w1*p%z(hi)
         
         return
         
       end block

    endif
    
  end subroutine get_coordinates_from_time
  
end module module_tdata


module module_set_tdata
  use module_tdata
  use hdf5
  implicit none

  integer :: ntraj = 0
  type(traj), allocatable :: particles(:)

contains

  subroutine read_set_tdata_from_hdf5(fn)
    character(*), intent(in) :: fn

    integer :: hdf_err, ip
    integer(HID_T) :: file_id

    call h5open_f(hdf_err)
    if (hdf_err /= 0) then
       write(6,*) "ERROR: h5open_f failed."
       stop
    endif

    call h5fopen_f(trim(fn), H5F_ACC_RDONLY_F, file_id, hdf_err)
    if (hdf_err /= 0) then
       write(6,*) "ERROR: h5fopen_f failed. file = ", trim(fn)
       stop
    endif

    call read_scalar_integer_hdf5(file_id, "/ntraj", ntraj)
    
    if (allocated(particles)) call deallocate_particles()
    allocate(particles(ntraj))

    do ip = 1, ntraj
       call read_one_tdata_from_hdf5(file_id, ip, particles(ip))
    enddo
    
    call h5fclose_f(file_id, hdf_err)
    call h5close_f(hdf_err)
  end subroutine read_set_tdata_from_hdf5


  subroutine read_one_tdata_from_hdf5(file_id, ip, p)
    integer(HID_T), intent(in) :: file_id
    integer, intent(in) :: ip
    type(traj), intent(inout) :: p

    integer :: hdf_err, ntime_p
    integer(HID_T) :: group_id
    character(256) :: gname

    write(gname,'(i8.8)') ip

    call h5gopen_f(file_id, trim(gname), group_id, hdf_err)
    if (hdf_err /= 0) then
       write(6,*) "ERROR: h5gopen_f failed. group = ", trim(gname)
       stop
    endif

    call get_1d_dataset_size_hdf5(group_id, "time", ntime_p)
    call tdata_allocate(p, ntime_p)

    p%ip        = ip
    p%mass      = 0.d0
    p%ut1_final = 0.d0
    p%hut_final = 0.d0

    p%time    = 0.d0
    p%x       = 0.d0
    p%y       = 0.d0
    p%z       = 0.d0
    p%vlx     = 0.d0
    p%vly     = 0.d0
    p%vlz     = 0.d0
    p%qrho    = 0.d0
    p%tem     = 0.d0
    p%ye      = 0.d0
    p%sen     = 0.d0
    p%rne     = 0.d0
    p%rae     = 0.d0
    p%deptn   = 0.d0
    p%depta   = 0.d0
    p%ut      = 0.d0
    p%hhh     = 0.d0
    p%eta_nue = 0.d0
    p%eta_nub = 0.d0
    p%b2      = 0.d0
    p%pres    = 0.d0

    ! scalar datasets
    call read_scalar_real4_to_real8_hdf5(group_id, "mass", p%mass)
    call read_scalar_real4_to_real8_hdf5(group_id, "ut1_final", p%ut1_final)
    call read_scalar_real4_to_real8_hdf5(group_id, "hut_final", p%hut_final)

    ! required 1D datasets
    call read_1d_real4_to_real8_hdf5(group_id, "time", p%time)
    call read_1d_real4_to_real8_hdf5(group_id, "x",    p%x)
    call read_1d_real4_to_real8_hdf5(group_id, "y",    p%y)
    call read_1d_real4_to_real8_hdf5(group_id, "z",    p%z)
    call read_1d_real4_to_real8_hdf5(group_id, "rho",  p%qrho)
    call read_1d_real4_to_real8_hdf5(group_id, "temp", p%tem)
    call read_1d_real4_to_real8_hdf5(group_id, "Ye",   p%ye)
    call read_1d_real4_to_real8_hdf5(group_id, "entr", p%sen)
    call read_1d_real4_to_real8_hdf5(group_id, "v^x",     p%vlx)
    call read_1d_real4_to_real8_hdf5(group_id, "v^y",     p%vly)
    call read_1d_real4_to_real8_hdf5(group_id, "v^z",     p%vlz)
    call read_1d_real4_to_real8_hdf5(group_id, "Enue",     p%rne)
    call read_1d_real4_to_real8_hdf5(group_id, "Enub",     p%rae)
    call read_1d_real4_to_real8_hdf5(group_id, "tau_nue",   p%deptn)
    call read_1d_real4_to_real8_hdf5(group_id, "tau_nub",   p%depta)
    call read_1d_real4_to_real8_hdf5(group_id, "u_t",      p%ut)
    call read_1d_real4_to_real8_hdf5(group_id, "h",     p%hhh)
    call read_1d_real4_to_real8_hdf5(group_id, "eta_nue", p%eta_nue)
    call read_1d_real4_to_real8_hdf5(group_id, "eta_nub", p%eta_nub)
    call read_1d_real4_to_real8_hdf5(group_id, "b^2",      p%b2)
    call read_1d_real4_to_real8_hdf5(group_id, "pres",    p%pres)

    call h5gclose_f(group_id, hdf_err)
  end subroutine read_one_tdata_from_hdf5


  subroutine read_scalar_integer_hdf5(loc_id, name, val)
    integer(HID_T), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer, intent(out) :: val

    integer :: hdf_err
    integer(HID_T) :: dset_id
    integer(HSIZE_T) :: dims1(1)

    dims1(1) = 1

    call h5dopen_f(loc_id, trim(name), dset_id, hdf_err)
    if (hdf_err /= 0) then
       write(6,*) "ERROR: h5dopen_f failed for ", trim(name)
       stop
    endif

    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, val, dims1, hdf_err)
    if (hdf_err /= 0) then
       write(6,*) "ERROR: h5dread_f failed for ", trim(name)
       stop
    endif

    call h5dclose_f(dset_id, hdf_err)
  end subroutine read_scalar_integer_hdf5


  subroutine read_scalar_real4_to_real8_hdf5(loc_id, name, val)
    integer(HID_T), intent(in) :: loc_id
    character(*), intent(in) :: name
    real(8), intent(out) :: val

    integer :: hdf_err
    integer(HID_T) :: dset_id
    integer(HSIZE_T) :: dims1(1)
    real(4) :: tmp

    dims1(1) = 1

    call h5dopen_f(loc_id, trim(name), dset_id, hdf_err)
    if (hdf_err /= 0) then
       write(6,*) "ERROR: h5dopen_f failed for ", trim(name)
       stop
    endif

    call h5dread_f(dset_id, H5T_NATIVE_REAL, tmp, dims1, hdf_err)
    if (hdf_err /= 0) then
       write(6,*) "ERROR: h5dread_f failed for ", trim(name)
       stop
    endif

    val = real(tmp,8)

    call h5dclose_f(dset_id, hdf_err)
  end subroutine read_scalar_real4_to_real8_hdf5


  subroutine try_read_scalar_real4_to_real8_hdf5(loc_id, name, val)
    integer(HID_T), intent(in) :: loc_id
    character(*), intent(in) :: name
    real(8), intent(inout) :: val

    integer :: hdf_err
    integer(HID_T) :: dset_id
    integer(HSIZE_T) :: dims1(1)
    real(4) :: tmp

    dims1(1) = 1

    call h5dopen_f(loc_id, trim(name), dset_id, hdf_err)
    if (hdf_err /= 0) return

    call h5dread_f(dset_id, H5T_NATIVE_REAL, tmp, dims1, hdf_err)
    if (hdf_err == 0) val = real(tmp,8)

    call h5dclose_f(dset_id, hdf_err)
  end subroutine try_read_scalar_real4_to_real8_hdf5


  subroutine read_1d_real4_to_real8_hdf5(loc_id, name, arr)
    integer(HID_T), intent(in) :: loc_id
    character(*), intent(in) :: name
    real(8), intent(out) :: arr(:)

    integer :: hdf_err, n
    integer(HID_T) :: dset_id
    integer(HSIZE_T) :: dims1(1)
    real(4), allocatable :: tmp(:)

    n = size(arr,1)
    dims1(1) = n
    allocate(tmp(n))

    call h5dopen_f(loc_id, trim(name), dset_id, hdf_err)
    if (hdf_err /= 0) then
       write(6,*) "ERROR: h5dopen_f failed for ", trim(name)
       stop
    endif

    call h5dread_f(dset_id, H5T_NATIVE_REAL, tmp, dims1, hdf_err)
    if (hdf_err /= 0) then
       write(6,*) "ERROR: h5dread_f failed for ", trim(name)
       stop
    endif

    arr = real(tmp,8)

    deallocate(tmp)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine read_1d_real4_to_real8_hdf5


  subroutine try_read_1d_real4_to_real8_hdf5(loc_id, name, arr)
    integer(HID_T), intent(in) :: loc_id
    character(*), intent(in) :: name
    real(8), intent(inout) :: arr(:)

    integer :: hdf_err, n
    integer(HID_T) :: dset_id
    integer(HSIZE_T) :: dims1(1)
    real(4), allocatable :: tmp(:)

    n = size(arr,1)
    dims1(1) = n

    call h5dopen_f(loc_id, trim(name), dset_id, hdf_err)
    if (hdf_err /= 0) return

    allocate(tmp(n))
    call h5dread_f(dset_id, H5T_NATIVE_REAL, tmp, dims1, hdf_err)
    if (hdf_err == 0) arr = real(tmp,8)
    deallocate(tmp)

    call h5dclose_f(dset_id, hdf_err)
  end subroutine try_read_1d_real4_to_real8_hdf5


  subroutine get_1d_dataset_size_hdf5(loc_id, name, n)
    integer(HID_T), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer, intent(out) :: n

    integer :: hdf_err, rank
    integer(HID_T) :: dset_id, dspace_id
    integer(HSIZE_T) :: dims(1), maxdims(1)

    call h5dopen_f(loc_id, trim(name), dset_id, hdf_err)
    if (hdf_err /= 0) then
       write(6,*) "ERROR: h5dopen_f failed for ", trim(name)
       stop
    endif

    call h5dget_space_f(dset_id, dspace_id, hdf_err)
    call h5sget_simple_extent_ndims_f(dspace_id, rank, hdf_err)
    if (rank /= 1) then
       write(6,*) "ERROR: dataset rank is not 1 for ", trim(name), " rank=", rank
       stop
    endif

    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, hdf_err)
    n = int(dims(1))

    call h5sclose_f(dspace_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine get_1d_dataset_size_hdf5


  integer function get_ntraj_set() result(n)
    n = ntraj
  end function get_ntraj_set


  subroutine deallocate_particles()
    integer :: ip

    if (allocated(particles)) then
       do ip = 1, size(particles)
          call tdata_deallocate(particles(ip))
       enddo
       deallocate(particles)
    endif

    ntraj = 0
  end subroutine deallocate_particles


  subroutine copy_particle(ip, p)
    integer, intent(in) :: ip
    type(traj), intent(inout) :: p

    if (.not. allocated(particles)) then
      write(6,*) "ERROR: particles is not allocated."
      stop
    endif

    if (ip < 1 .or. ip > ntraj) then
       write(6,*) "ERROR: ip out of range. ip=", ip
       stop
    endif

    call tdata_allocate(p, particles(ip)%ntime)

    p%ip        = particles(ip)%ip
    p%mass      = particles(ip)%mass
    p%ut1_final = particles(ip)%ut1_final
    p%hut_final = particles(ip)%hut_final

    p%time    = particles(ip)%time
    p%x       = particles(ip)%x
    p%y       = particles(ip)%y
    p%z       = particles(ip)%z
    p%vlx     = particles(ip)%vlx
    p%vly     = particles(ip)%vly
    p%vlz     = particles(ip)%vlz
    p%qrho    = particles(ip)%qrho
    p%tem     = particles(ip)%tem
    p%ye      = particles(ip)%ye
    p%sen     = particles(ip)%sen
    p%rne     = particles(ip)%rne
    p%rae     = particles(ip)%rae
    p%deptn   = particles(ip)%deptn
    p%depta   = particles(ip)%depta
    p%ut      = particles(ip)%ut
    p%hhh     = particles(ip)%hhh
    p%eta_nue = particles(ip)%eta_nue
    p%eta_nub = particles(ip)%eta_nub
    p%b2      = particles(ip)%b2
    p%pres    = particles(ip)%pres
  end subroutine copy_particle

end module module_set_tdata
