module module_trajectory
  implicit none

  private

  public :: traj_set_dir, traj_read_info, traj_show_info, traj_make_hdf5

  ! from trajectries
  integer :: ntraj
  integer :: ntime
  character(256) :: dir_traj
  
  real(8),parameter :: pi = 4d0*atan(1d0)
  real(8),parameter :: day = 24d0*60d0*60d0
  real(8),parameter :: clight = 2.99792458d10
  
contains
  
  
  subroutine traj_show_info
    write(6,'(a)') trim(dir_traj)
    write(6,*) ntraj
    write(6,*) ntime
    
  end subroutine traj_show_info
  
  subroutine traj_set_dir(dir_traj_in)
    character(*),intent(in) :: dir_traj_in
    dir_traj = trim(dir_traj_in)

  end subroutine traj_set_dir

  subroutine traj_read_info
    
    !$ use omp_lib
    
    integer :: iostat
    
    ! openmp
    integer :: my_thr, max_thr
    !

    integer :: access

    character(256) :: fn
    integer :: ip, it

    block
      integer,parameter :: ntraj_limit = 1000000
      integer :: i
      logical :: exists

      ip=0
      find_n:do i=1,ntraj_limit
         write(fn,'(a,"/traj_",i8.8,".dat")') trim(dir_traj), i
         inquire(file=fn, exist=exists)
         if(exists)then
            ip=i
         else
            exit find_n
         endif
      enddo find_n
      ntraj=ip
    end block
!     fn = trim(dir_traj)//"/ana_traj.dat"
!     open(10,file=fn,iostat=iostat,status="old",action="read")
!     read(10,*);read(10,*)
!     ip = 0
!     do
!        read(10,*,end=99)
!        ip=ip+1
!     enddo
! 99  close(10)
!     ntraj = ip

    write(6,*) "ntraj = ",ntraj

    ip=1
    write(fn,'(a,"/traj_",i8.8,".dat")') trim(dir_traj), ip
    ! write(6,'(a)') fn

    open(10,file=fn,iostat=iostat,status="old",action="read")
    read(10,*);read(10,*);read(10,*);read(10,*)
    it = 0
    do
       read(10,*,end=999)
       it=it+1
    enddo
999  close(10)
    ntime = it
 
    write(6,*) "ntime = ",ntime
    
    
  end subroutine traj_read_info



  subroutine traj_make_hdf5
    !$ use omp_lib
    use hdf5

    real(4), allocatable :: time_p(:)
    real(4), allocatable :: x_p(:), y_p(:), z_p(:)
    real(4), allocatable :: vlx_p(:), vly_p(:), vlz_p(:)
    real(4), allocatable :: qrho_p(:), tem_p(:), ye_p(:), sen_p(:)
    real(4), allocatable :: rne_p(:), rae_p(:)
    real(4), allocatable :: deptn_p(:), depta_p(:)
    real(4), allocatable :: ut_p(:), hhh_p(:)
    real(4), allocatable :: eta_nue_p(:), eta_nub_p(:), b2_p(:), pres_p(:)

    real(4) :: mass_p
    real(4) :: ut1_final_p, hut_final_p

    integer :: my_thr, max_thr
    integer :: ip, it, ntime_p
    integer :: nunit, iostat
    integer :: hdf_err

    character(256) :: fn, fn_hdf, str1

    integer(HID_T) :: file_id, group_id

    fn_hdf = trim(dir_traj) // "/all_traj.h5"

    call h5open_f(hdf_err)
    call h5fcreate_f(fn_hdf, H5F_ACC_TRUNC_F, file_id, hdf_err)
    if (hdf_err /= 0) then
       write(6,'(a)') "ERROR: HDF5 file creation failed."
       stop
    endif

    call write_scalar_int_hdf5_root(file_id, "ntraj", ntraj)

    my_thr  = 0
    max_thr = 1

    !$omp parallel default(none) &
    !$omp shared(ntraj,ntime,dir_traj,file_id) &
    !$omp private(ip,it,nunit,iostat,fn,str1,group_id,hdf_err,my_thr,max_thr,ntime_p, &
    !$omp         mass_p,ut1_final_p,hut_final_p, &
    !$omp         time_p,x_p,y_p,z_p,vlx_p,vly_p,vlz_p, &
    !$omp         qrho_p,tem_p,ye_p,sen_p,rne_p,rae_p, &
    !$omp         deptn_p,depta_p,ut_p,hhh_p,eta_nue_p,eta_nub_p,b2_p,pres_p)

    !$ my_thr  = omp_get_thread_num()
    !$ max_thr = omp_get_num_threads()

    !$omp do schedule(dynamic)
    do ip = 1, ntraj

       allocate(time_p(ntime), x_p(ntime), y_p(ntime), z_p(ntime), &
            vlx_p(ntime), vly_p(ntime), vlz_p(ntime), &
            qrho_p(ntime), tem_p(ntime), ye_p(ntime), sen_p(ntime), &
            rne_p(ntime), rae_p(ntime), &
            deptn_p(ntime), depta_p(ntime), &
            ut_p(ntime), hhh_p(ntime), &
            eta_nue_p(ntime), eta_nub_p(ntime), b2_p(ntime), pres_p(ntime))


       mass_p      = 0.0
       ut1_final_p = 0.0
       hut_final_p = 0.0
       ntime_p     = 0

       time_p    = 0.0
       x_p       = 0.0
       y_p       = 0.0
       z_p       = 0.0
       vlx_p     = 0.0
       vly_p     = 0.0
       vlz_p     = 0.0
       qrho_p    = 0.0
       tem_p     = 0.0
       ye_p      = 0.0
       sen_p     = 0.0
       rne_p     = 0.0
       rae_p     = 0.0
       deptn_p   = 0.0
       depta_p   = 0.0
       ut_p      = 0.0
       hhh_p     = 0.0
       eta_nue_p = 0.0
       eta_nub_p = 0.0
       b2_p      = 0.0
       pres_p    = 0.0

       write(fn,'(a,"/traj_",i8.8,".dat")') trim(dir_traj), ip
       open(newunit=nunit, file=fn, iostat=iostat, status="old", action="read")
       if (iostat /= 0) then
          write(6,*) "ERROR: failed to open ", trim(fn), " iostat=", iostat

          deallocate(time_p, x_p, y_p, z_p, &
               vlx_p, vly_p, vlz_p, &
               qrho_p, tem_p, ye_p, sen_p, &
               rne_p, rae_p, deptn_p, depta_p, &
               ut_p, hhh_p, eta_nue_p, eta_nub_p, b2_p, pres_p)
          
          cycle
       endif

       read(nunit,*)
       read(nunit,*)
       read(nunit,'(16x,es13.5,20x,2es13.5)') mass_p, ut1_final_p, hut_final_p
       read(nunit,*)

       do it = 1, ntime
          read(nunit,*,end=99) &
               time_p(it),    &
               x_p(it),       &
               y_p(it),       &
               z_p(it),       &
               vlx_p(it),     &
               vly_p(it),     &
               vlz_p(it),     &
               qrho_p(it),    &
               tem_p(it),     &
               ye_p(it),      &
               sen_p(it),     &
               rne_p(it),     &
               rae_p(it),     &
               deptn_p(it),   &
               depta_p(it),   &
               ut_p(it),      &
               hhh_p(it),     &
               eta_nue_p(it), &
               eta_nub_p(it), &
               b2_p(it),      &
               pres_p(it)
          ntime_p = it
       enddo

99     close(nunit)

       if (ntime_p <= 0) then
          write(6,*) "WARNING: no trajectory data for ip=", ip
          
          deallocate(time_p, x_p, y_p, z_p, &
               vlx_p, vly_p, vlz_p, &
               qrho_p, tem_p, ye_p, sen_p, &
               rne_p, rae_p, deptn_p, depta_p, &
               ut_p, hhh_p, eta_nue_p, eta_nub_p, b2_p, pres_p)

          cycle
       endif

       write(str1,'(i8.8)') ip

       !$omp critical(hdf5_io)

       call h5gcreate_f(file_id, trim(str1), group_id, hdf_err)
       if (hdf_err /= 0) then
          write(6,*) "ERROR: h5gcreate_f failed for group ", trim(str1)
       else
          call write_scalar_int_hdf5(group_id,  "ntime",      ntime_p)
          call write_scalar_real_hdf5(group_id, "mass",       mass_p)
          call write_scalar_real_hdf5(group_id, "ut1_final",  ut1_final_p)
          call write_scalar_real_hdf5(group_id, "hut_final",  hut_final_p)

          call write_1d_real_hdf5(group_id, "time",   time_p(1:ntime_p),   ntime_p)
          call write_1d_real_hdf5(group_id, "x",      x_p(1:ntime_p),      ntime_p)
          call write_1d_real_hdf5(group_id, "y",      y_p(1:ntime_p),      ntime_p)
          call write_1d_real_hdf5(group_id, "z",      z_p(1:ntime_p),      ntime_p)
          call write_1d_real_hdf5(group_id, "v^x",    vlx_p(1:ntime_p),    ntime_p)
          call write_1d_real_hdf5(group_id, "v^y",    vly_p(1:ntime_p),    ntime_p)
          call write_1d_real_hdf5(group_id, "v^z",    vlz_p(1:ntime_p),    ntime_p)
          call write_1d_real_hdf5(group_id, "rho",    qrho_p(1:ntime_p),   ntime_p)
          call write_1d_real_hdf5(group_id, "temp",   tem_p(1:ntime_p),    ntime_p)
          call write_1d_real_hdf5(group_id, "Ye",     ye_p(1:ntime_p),     ntime_p)
          call write_1d_real_hdf5(group_id, "entr",   sen_p(1:ntime_p),    ntime_p)
          call write_1d_real_hdf5(group_id, "Enue",    rne_p(1:ntime_p),    ntime_p)
          call write_1d_real_hdf5(group_id, "Enub",    rae_p(1:ntime_p),    ntime_p)
          call write_1d_real_hdf5(group_id, "tau_nue",  deptn_p(1:ntime_p),  ntime_p)
          call write_1d_real_hdf5(group_id, "tau_nub",  depta_p(1:ntime_p),  ntime_p)
          call write_1d_real_hdf5(group_id, "u_t",     ut_p(1:ntime_p),     ntime_p)
          call write_1d_real_hdf5(group_id, "h",    hhh_p(1:ntime_p),    ntime_p)

          ! 現在のtype(traj)には無いが、HDF5には保存しておく
          call write_1d_real_hdf5(group_id, "eta_nue", eta_nue_p(1:ntime_p), ntime_p)
          call write_1d_real_hdf5(group_id, "eta_nub", eta_nub_p(1:ntime_p), ntime_p)
          call write_1d_real_hdf5(group_id, "b^2",      b2_p(1:ntime_p),      ntime_p)
          call write_1d_real_hdf5(group_id, "pres",    pres_p(1:ntime_p),    ntime_p)

          call h5gclose_f(group_id, hdf_err)
       endif

       !$omp end critical(hdf5_io)

       if (my_thr == 0) write(6,*) ip, ntraj/max_thr

       deallocate(time_p, x_p, y_p, z_p, &
            vlx_p, vly_p, vlz_p, &
            qrho_p, tem_p, ye_p, sen_p, &
            rne_p, rae_p, deptn_p, depta_p, &
            ut_p, hhh_p, eta_nue_p, eta_nub_p, b2_p, pres_p)
       

    enddo
    !$omp end do
    !$omp end parallel


    call h5fclose_f(file_id, hdf_err)
    call h5close_f(hdf_err)

  contains

    subroutine write_scalar_real_hdf5(loc_id, name, val)
      integer(HID_T), intent(in) :: loc_id
      character(*), intent(in)   :: name
      real(4), intent(in)        :: val

      integer :: ierr, rank0
      integer(HID_T) :: dsid, spid
      integer(HSIZE_T) :: dd(1)

      rank0 = 1
      dd(1) = 1
      call h5screate_simple_f(rank0, dd, spid, ierr)
      call h5dcreate_f(loc_id, trim(name), H5T_NATIVE_REAL, spid, dsid, ierr)
      call h5dwrite_f(dsid, H5T_NATIVE_REAL, val, dd, ierr)
      call h5dclose_f(dsid, ierr)
      call h5sclose_f(spid, ierr)
    end subroutine write_scalar_real_hdf5

    subroutine write_scalar_int_hdf5(loc_id, name, val)
      integer(HID_T), intent(in) :: loc_id
      character(*), intent(in)   :: name
      integer, intent(in)        :: val

      integer :: ierr, rank0
      integer(HID_T) :: dsid, spid
      integer(HSIZE_T) :: dd(1)

      rank0 = 1
      dd(1) = 1
      call h5screate_simple_f(rank0, dd, spid, ierr)
      call h5dcreate_f(loc_id, trim(name), H5T_NATIVE_INTEGER, spid, dsid, ierr)
      call h5dwrite_f(dsid, H5T_NATIVE_INTEGER, val, dd, ierr)
      call h5dclose_f(dsid, ierr)
      call h5sclose_f(spid, ierr)
    end subroutine write_scalar_int_hdf5

    subroutine write_scalar_int_hdf5_root(loc_id, name, val)
      integer(HID_T), intent(in) :: loc_id
      character(*), intent(in)   :: name
      integer, intent(in)        :: val

      call write_scalar_int_hdf5(loc_id, name, val)
    end subroutine write_scalar_int_hdf5_root

    subroutine write_1d_real_hdf5(loc_id, name, arr, n)
      integer(HID_T), intent(in) :: loc_id
      character(*), intent(in)   :: name
      real(4), intent(in)        :: arr(:)
      integer, intent(in)        :: n

      integer :: ierr, rank0
      integer(HID_T) :: dsid, spid
      integer(HSIZE_T) :: dd(1)

      rank0 = 1
      dd(1) = n
      call h5screate_simple_f(rank0, dd, spid, ierr)
      call h5dcreate_f(loc_id, trim(name), H5T_NATIVE_REAL, spid, dsid, ierr)
      call h5dwrite_f(dsid, H5T_NATIVE_REAL, arr, dd, ierr)
      call h5dclose_f(dsid, ierr)
      call h5sclose_f(spid, ierr)
    end subroutine write_1d_real_hdf5

  end subroutine traj_make_hdf5

end module module_trajectory
