module module_trajectory
  implicit none

  private

  public :: traj_set_dir, traj_read_info, traj_show_info, traj_make_hdf5, traj_read_hdf, traj_read_info_hdf,traj_get_ntraj

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

    fn = trim(dir_traj)//"/ana_traj.dat"
    open(10,file=fn,iostat=iostat,status="old",action="read")
    read(10,*);read(10,*)
    ip = 0
    do
       read(10,*,end=99)
       ip=ip+1
    enddo
99  close(10)
    ntraj = ip
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
    
    real(4),allocatable :: time_p(:)
    real(4),allocatable :: x_p(:),y_p(:),z_p(:),&
         qrho_p(:),&
         ye_p  (:),&
         tem_p (:),&
         ut_p  (:),&
         qb_p  (:),&
         sen_p (:),&
         vlx_p (:),&
         vly_p (:),&
         vlz_p (:),&
         hhh_p (:),&
         dt_p  (:),&
         rne_p (:),&
         rae_p (:),&
         deptn_p(:),&
         depta_p(:)

    real(4) :: mass_p, ut1, hut, ut

    ! openmp
    integer :: my_thr, max_thr
    !
    character(256) :: fn

    integer :: it, ip, ntime_p
    integer :: nunit, iostat


    
    character(256) :: fn_hdf
    integer :: hdf_err,rank
    integer(HID_T) :: file_id, dset_id, dspace_id, group_id
    integer(HSIZE_T) :: dims1(1), dims2(2), dims3(3)


    character(256) :: str1

    !
    fn_hdf = trim(dir_traj) // "/all_traj.h5"
    ! fn_hdf = "all_traj.h5"

    call h5open_f(hdf_err)
    call h5fcreate_f(fn_hdf, H5F_ACC_TRUNC_F, file_id, hdf_err)
    if(hdf_err/=0)then
       write(6,'(a)') "File creation failed."
       stop
    endif
    rank=1
    dims1(1) = 1
    call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
    call h5dcreate_f(file_id, "/ntraj", H5T_NATIVE_INTEGER, &
         dspace_id,dset_id, hdf_err)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, ntraj, dims1, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    
    !
    allocate(time_p(ntime),&
         x_p(ntime),y_p(ntime),z_p(ntime),&
         qrho_p(ntime),&
         ye_p  (ntime),&
         tem_p (ntime),&
         ut_p  (ntime),&
         qb_p  (ntime),&
         sen_p (ntime),&
         vlx_p (ntime),&
         vly_p (ntime),&
         vlz_p (ntime),&
         hhh_p (ntime),&
         dt_p  (ntime),&
         rne_p (ntime),&
         rae_p (ntime),&
         deptn_p(ntime),&
         depta_p(ntime) )

    my_thr = 0
    max_thr= 1

    !$omp parallel default(none) &
    !$omp shared(ntraj,ntime,dir_traj, file_id, H5T_NATIVE_REAL) &
    !$omp private(nunit,fn,iostat, my_thr, max_thr, group_id, dspace_id, dset_id, hdf_err, rank, dims1, &
    !$omp   mass_p, ut1, hut, ut, ntime_p, str1, &
    !$omp   time_p,x_p,y_p,z_p,vlx_p,vly_p,vlz_p,qrho_p,tem_p,ye_p,sen_p,rne_p,rae_p)

    !$ my_thr  = omp_get_thread_num()
    !$ max_thr = omp_get_num_threads()

    !$omp do schedule(dynamic)
    do ip=1,ntraj,1
       
       write(fn,'(a,"/traj_",i8.8,".dat")') trim(dir_traj), ip
       open(newunit=nunit,file=fn,iostat=iostat,status="old",action="read")
       read(nunit,*)
       read(nunit,*)
       read(nunit,'(16x,es13.5,20x,2es13.5)') mass_p, ut1, hut
       read(nunit,*)
       
       ntime_p = 0
       do it=1,ntime
          
          read(nunit,*,end=99) &
               time_p(it), &
               x_p(it), &
               y_p(it), &
               z_p(it), &
               vlx_p(it), &
               vly_p(it), &
               vlz_p(it), &
               qrho_p(it), &
               tem_p(it), &
               ye_p(it), &
               sen_p(it), &
               rne_p(it), &
               rae_p(it)

          ntime_p = it
          
       enddo
       
99     close(nunit)

       
       ! do it=1,ntime_p
       !    write(6,'(99es12.4)') time_p(it), x_p(it)
       ! enddo

       ut = ut1-1d0

       write(str1,'(i8.8)') ip
       
       !$omp critical(hdf5_io)

       call h5gcreate_f(file_id, trim(str1) , group_id, hdf_err)

       rank=1
       dims1(1) = 1
       call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
       call h5dcreate_f(group_id, "mass", H5T_NATIVE_REAL, &
            dspace_id,dset_id, hdf_err)
       call h5dwrite_f(dset_id, H5T_NATIVE_REAL, mass_p, dims1, hdf_err)
       call h5dclose_f(dset_id, hdf_err)
       call h5sclose_f(dspace_id, hdf_err)

       rank=1
       dims1(1) = 1
       call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
       call h5dcreate_f(group_id, "u_t", H5T_NATIVE_REAL, &
            dspace_id,dset_id, hdf_err)
       call h5dwrite_f(dset_id, H5T_NATIVE_REAL, ut, dims1, hdf_err)
       call h5dclose_f(dset_id, hdf_err)
       call h5sclose_f(dspace_id, hdf_err)

       rank=1
       dims1(1) = ntime_p
       call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
       call h5dcreate_f(group_id, "time", H5T_NATIVE_REAL, &
            dspace_id,dset_id, hdf_err)
       call h5dwrite_f(dset_id, H5T_NATIVE_REAL, time_p(1:ntime_p), dims1, hdf_err)
       call h5dclose_f(dset_id, hdf_err)
       call h5sclose_f(dspace_id, hdf_err)

       call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
       call h5dcreate_f(group_id, "rho", H5T_NATIVE_REAL, &
            dspace_id,dset_id, hdf_err)
       call h5dwrite_f(dset_id, H5T_NATIVE_REAL, qrho_p(1:ntime_p), dims1, hdf_err)
       call h5dclose_f(dset_id, hdf_err)
       call h5sclose_f(dspace_id, hdf_err)

       call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
       call h5dcreate_f(group_id, "temp", H5T_NATIVE_REAL, &
            dspace_id,dset_id, hdf_err)
       call h5dwrite_f(dset_id, H5T_NATIVE_REAL, tem_p(1:ntime_p), dims1, hdf_err)
       call h5dclose_f(dset_id, hdf_err)
       call h5sclose_f(dspace_id, hdf_err)

       call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
       call h5dcreate_f(group_id, "Ye", H5T_NATIVE_REAL, &
            dspace_id,dset_id, hdf_err)
       call h5dwrite_f(dset_id, H5T_NATIVE_REAL, ye_p(1:ntime_p), dims1, hdf_err)
       call h5dclose_f(dset_id, hdf_err)
       call h5sclose_f(dspace_id, hdf_err)

       call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
       call h5dcreate_f(group_id, "entr", H5T_NATIVE_REAL, &
            dspace_id,dset_id, hdf_err)
       call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sen_p(1:ntime_p), dims1, hdf_err)
       call h5dclose_f(dset_id, hdf_err)
       call h5sclose_f(dspace_id, hdf_err)

       call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
       call h5dcreate_f(group_id, "x", H5T_NATIVE_REAL, &
            dspace_id,dset_id, hdf_err)
       call h5dwrite_f(dset_id, H5T_NATIVE_REAL, x_p(1:ntime_p), dims1, hdf_err)
       call h5dclose_f(dset_id, hdf_err)
       call h5sclose_f(dspace_id, hdf_err)

       call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
       call h5dcreate_f(group_id, "y", H5T_NATIVE_REAL, &
            dspace_id,dset_id, hdf_err)
       call h5dwrite_f(dset_id, H5T_NATIVE_REAL, y_p(1:ntime_p), dims1, hdf_err)
       call h5dclose_f(dset_id, hdf_err)
       call h5sclose_f(dspace_id, hdf_err)

       call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
       call h5dcreate_f(group_id, "z", H5T_NATIVE_REAL, &
            dspace_id,dset_id, hdf_err)
       call h5dwrite_f(dset_id, H5T_NATIVE_REAL, z_p(1:ntime_p), dims1, hdf_err)
       call h5dclose_f(dset_id, hdf_err)
       call h5sclose_f(dspace_id, hdf_err)


       call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
       call h5dcreate_f(group_id, "v^x", H5T_NATIVE_REAL, &
            dspace_id,dset_id, hdf_err)
       call h5dwrite_f(dset_id, H5T_NATIVE_REAL, vlx_p(1:ntime_p), dims1, hdf_err)
       call h5dclose_f(dset_id, hdf_err)
       call h5sclose_f(dspace_id, hdf_err)

       call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
       call h5dcreate_f(group_id, "v^y", H5T_NATIVE_REAL, &
            dspace_id,dset_id, hdf_err)
       call h5dwrite_f(dset_id, H5T_NATIVE_REAL, vly_p(1:ntime_p), dims1, hdf_err)
       call h5dclose_f(dset_id, hdf_err)
       call h5sclose_f(dspace_id, hdf_err)

       call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
       call h5dcreate_f(group_id, "v^z", H5T_NATIVE_REAL, &
            dspace_id,dset_id, hdf_err)
       call h5dwrite_f(dset_id, H5T_NATIVE_REAL, vlz_p(1:ntime_p), dims1, hdf_err)
       call h5dclose_f(dset_id, hdf_err)
       call h5sclose_f(dspace_id, hdf_err)

       
    
       call h5gclose_f(group_id, hdf_err)
       
       !$omp end critical(hdf5_io)
       
       ! exit
       if(my_thr==0) write(6,*) ip,ntraj/max_thr
       ! stop

    enddo
    !$omp end do
    !$omp end parallel

    call h5fclose_f(file_id, hdf_err)    
    
  end subroutine traj_make_hdf5

  subroutine traj_read_info_hdf
    use hdf5

    character(256) :: fn_hdf
    integer :: hdf_err

    integer(HID_T) :: file_id, dset_id, dspace_id, group_id
    integer(HSIZE_T) :: dims1(1), dims2(2), dims3(3)

    fn_hdf = trim(dir_traj) // "/all_traj.h5"
    
    call h5open_f(hdf_err)

    call h5fopen_f(fn_hdf, H5F_ACC_RDONLY_F, file_id, hdf_err)
    if(hdf_err/=0)then
       write(6,*) hdf_err
       write(6,'(a)') fn_hdf
       write(6,'(a)') "File open failed."
       stop
    endif

    dims1(1) = 1
    call h5dopen_f(file_id, "ntraj", dset_id, hdf_err)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ntraj, dims1, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    
    write(6,'(a,i10)') "# number of trajectory read:", ntraj

  end subroutine traj_read_info_hdf

  real(8) function traj_get_ntraj() result(n)
    n=ntraj
  end function traj_get_ntraj
  
  subroutine traj_read_hdf(ip,p)
    use hdf5
    use module_tdata
    
    ! particle id
    integer,intent(in) :: ip
    type(traj),intent(inout) :: p

    character(256) :: fn_hdf
    integer :: hdf_err

    integer(HID_T) :: file_id, dset_id, dspace_id, group_id
    integer(HSIZE_T) :: dims1(1), dims2(2), dims3(3)
    real(4),allocatable :: buf_real4(:)

    character(256) :: group
    integer(HSIZE_T) :: npoints

    fn_hdf = trim(dir_traj) // "/all_traj.h5"
    
    call h5open_f(hdf_err)

    call h5fopen_f(fn_hdf, H5F_ACC_RDONLY_F, file_id, hdf_err)
    if(hdf_err/=0)then
       write(6,*) hdf_err
       write(6,'(a)') fn_hdf
       write(6,'(a)') "File open failed."
       stop
    endif

    dims1(1) = 1
    call h5dopen_f(file_id, "ntraj", dset_id, hdf_err)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ntraj, dims1, hdf_err)
    call h5dclose_f(dset_id, hdf_err)

    if(ip<=0.or.ntraj<ip)then
       write(6,'(a)') "particle id out of range."
       stop
    endif
    
    write(group,'(i8.8)') ip
    
    call h5dopen_f(file_id, "/"//trim(group)//"/time", dset_id, hdf_err)
        
    call h5dget_space_f(dset_id, dspace_id, hdf_err) 
    call h5sget_simple_extent_npoints_f(dspace_id, npoints, hdf_err) 
    ntime = npoints
    call h5sclose_f(dspace_id,hdf_err)
    call h5dclose_f(dset_id,hdf_err)
    
    call tdata_allocate(p,ntime)


    dims1(1) = 1
    allocate(buf_real4(1))
    call h5dopen_f(file_id, "/"//trim(group)//"/mass", dset_id, hdf_err)
    call h5dread_f(dset_id, H5T_NATIVE_REAL, buf_real4, dims1, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    p%mass = dble(buf_real4(1))
    deallocate(buf_real4)


    allocate(buf_real4(ntime))

    dims1(1) = ntime

    call h5dopen_f(file_id, "/"//trim(group)//"/time", dset_id, hdf_err)
    call h5dread_f(dset_id, H5T_NATIVE_REAL, buf_real4, dims1, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    p%time(:) = dble(buf_real4(:))

    call h5dopen_f(file_id, "/"//trim(group)//"/x", dset_id, hdf_err)
    call h5dread_f(dset_id, H5T_NATIVE_REAL, buf_real4, dims1, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    p%x(:) = dble(buf_real4(:))

    call h5dopen_f(file_id, "/"//trim(group)//"/y", dset_id, hdf_err)
    call h5dread_f(dset_id, H5T_NATIVE_REAL, buf_real4, dims1, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    p%y(:) = dble(buf_real4(:))

    call h5dopen_f(file_id, "/"//trim(group)//"/z", dset_id, hdf_err)
    call h5dread_f(dset_id, H5T_NATIVE_REAL, buf_real4, dims1, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    p%z(:) = dble(buf_real4(:))

    call h5dopen_f(file_id, "/"//trim(group)//"/Ye", dset_id, hdf_err)
    call h5dread_f(dset_id, H5T_NATIVE_REAL, buf_real4, dims1, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    p%ye(:) = dble(buf_real4(:))
    
    call h5dopen_f(file_id, "/"//trim(group)//"/temp", dset_id, hdf_err)
    call h5dread_f(dset_id, H5T_NATIVE_REAL, buf_real4, dims1, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    p%tem(:) = dble(buf_real4(:))

    deallocate(buf_real4)

  end subroutine traj_read_hdf

  

!   subroutine qtab_read_hdf(fn)
    
!     use hdf5
!     use h5lt
    
!     character(*),intent(in) :: fn
    
!     ! hdf5 output
!     integer :: hdf_err,rank
!     integer(HID_T) :: file_id, dset_id, dspace_id
!     integer(HSIZE_T) :: dims1(1), dims2(2), dims3(3)

!     character(256) :: str1

!     call h5fopen_f(fn, H5F_ACC_RDONLY_F, file_id, hdf_err)

!     call H5LTread_dataset_string_f(file_id, "/original_model_name", original_model_name, hdf_err)

!     call H5LTread_dataset_string_f(file_id, "/directory_trajectory", directory_trajectory, hdf_err)

!     call H5LTread_dataset_string_f(file_id, "/directory_heating_rate", directory_heating_rate, hdf_err)

!     dims1(1) = 1
!     call h5dopen_f(file_id, "ntraj", dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ntraj, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5dopen_f(file_id, "ntraj_inside", dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ntraj_inside, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
    
!     call h5dopen_f(file_id, "ntnc", dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ntnc, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5dopen_f(file_id, "nz", dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nz, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call qtab_allocate

!     ! allocate(data_exists(ntraj), dm(ntraj), ye(ntraj), entr(ntraj), texp(ntraj), t_ncl0(ntraj))

!     dims1(1) = ntraj
!     call h5dopen_f(file_id, "data_exists",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data_exists, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
    
!     call h5dopen_f(file_id, "dm",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dm, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
    
!     call h5dopen_f(file_id, "ye",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, ye, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
    
!     call h5dopen_f(file_id, "entr",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, entr, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
    
!     call h5dopen_f(file_id, "texp",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, texp, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
 
!     call h5dopen_f(file_id, "t_ncl0",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, t_ncl0, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5dopen_f(file_id, "ye_ext",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, ye_ext, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
    
!     call h5dopen_f(file_id, "vr_ext",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, vr_ext, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
        
!     !allocate(t_ncl(ntnc))
!     dims1(1) = ntnc
!     call h5dopen_f(file_id, "t_ncl",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, t_ncl, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
    
!     ! allocate( qtot(ntnc,ntraj),qgam(ntnc,ntraj),qele(ntnc,ntraj),qnu(ntnc,ntraj),qfis(ntnc,ntraj),qalp(ntnc,ntraj),mexc(ntnc,ntraj), &
!     !      xmass_1day(0:nz,ntraj), xmass_10day(0:nz,ntraj) )
    
!     dims2(1) = ntnc; dims2(2) = ntraj
!     call h5dopen_f(file_id, "qtot",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, qtot, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
    
!     call h5dopen_f(file_id, "qgam",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, qgam, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5dopen_f(file_id, "qele",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, qele, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5dopen_f(file_id, "qnu",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, qnu, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5dopen_f(file_id, "qfis",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, qfis, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5dopen_f(file_id, "qalp",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, qalp, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5dopen_f(file_id, "mexc",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, mexc, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     dims2(1) = nz+1; dims2(2) = ntraj
!     call h5dopen_f(file_id, "xmass_1day",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, xmass_1day, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
    
!     call h5dopen_f(file_id, "xmass_10day",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, xmass_10day, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

    
!     ! extraction
!     if(flux_based==1)then
       
!        dims1(1) = 1
!        call h5dopen_f(file_id, "n_time_ext", dset_id, hdf_err)
!        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, n_time_ext, dims1, hdf_err)
!        call h5dclose_f(dset_id, hdf_err)
       
!        call h5dopen_f(file_id, "n_theta_ext", dset_id, hdf_err)
!        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, n_theta_ext, dims1, hdf_err)
!        call h5dclose_f(dset_id, hdf_err)
       
!        call h5dopen_f(file_id, "r_ext", dset_id, hdf_err)
!        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, r_ext, dims1, hdf_err)
!        call h5dclose_f(dset_id, hdf_err)
       
!        call h5dopen_f(file_id, "dtheta_ext", dset_id, hdf_err)
!        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dtheta_ext, dims1, hdf_err)
!        call h5dclose_f(dset_id, hdf_err)
       
!        allocate(n_phi_ext(n_theta_ext))
!        allocate(time_ext(n_time_ext))
              
!        dims1(1) = n_theta_ext
!        call h5dopen_f(file_id, "n_phi_ext", dset_id, hdf_err)
!        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, n_phi_ext, dims1, hdf_err)
!        call h5dclose_f(dset_id, hdf_err)
       
!        dims1(1) = n_time_ext
!        call h5dopen_f(file_id, "time_ext", dset_id, hdf_err)
!        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, time_ext, dims1, hdf_err)
!        call h5dclose_f(dset_id, hdf_err)

!        dims1(1) = ntraj
!        ! call h5dopen_f(file_id, "time_ext", dset_id, hdf_err)
!        ! call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, time_ext, dims1, hdf_err)
!        ! call h5dclose_f(dset_id, hdf_err)

!        call h5dopen_f(file_id, "theta_ext", dset_id, hdf_err)
!        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, theta_ext, dims1, hdf_err)
!        call h5dclose_f(dset_id, hdf_err)

!        call h5dopen_f(file_id, "phi_ext", dset_id, hdf_err)
!        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, phi_ext, dims1, hdf_err)
!        call h5dclose_f(dset_id, hdf_err)

!        call h5dopen_f(file_id, "itime_ext", dset_id, hdf_err)
!        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, itime_ext, dims1, hdf_err)
!        call h5dclose_f(dset_id, hdf_err)

!        ! call h5dopen_f(file_id, "itheta_ext", dset_id, hdf_err)
!        ! call h5dread_f(dset_id, H5T_NATIVE_INTEGER, itheta_ext, dims1, hdf_err)
!        ! call h5dclose_f(dset_id, hdf_err)

!        ! call h5dopen_f(file_id, "iphi_ext", dset_id, hdf_err)
!        ! call h5dread_f(dset_id, H5T_NATIVE_INTEGER, iphi_ext, dims1, hdf_err)
!        ! call h5dclose_f(dset_id, hdf_err)
       
!     endif
    
!     call h5fclose_f(file_id, hdf_err)
    
    

!     write(6,'("fraction of existing data:",i7,"/",i7)') sum(data_exists(:)),  ntraj


!   end subroutine qtab_read_hdf


!   subroutine qtab_read_ref(fn)
    
!     use hdf5
    
!     character(*),intent(in) :: fn
    
!     ! hdf5 output
!     integer :: hdf_err,rank
!     integer(HID_T) :: file_id, dset_id, dspace_id
!     integer(HSIZE_T) :: dims1(1), dims2(2), dims3(3)

!     write(6,'("read reference qtable : ",a)') trim(fn)

!     file_used_to_fill_data = trim(fn)
    
!     call h5fopen_f(fn, H5F_ACC_RDONLY_F, file_id, hdf_err)
    
!     dims1(1) = 1
!     call h5dopen_f(file_id, "ntraj", dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ntraj_ref, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5dopen_f(file_id, "ntraj_inside", dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ntraj_inside_ref, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

    
!     allocate( data_exists_ref(ntraj_ref) )
!     allocate( ye_ref(ntraj_ref), entr_ref(ntraj_ref), texp_ref(ntraj_ref) )
!     allocate( qtot_ref(ntnc,ntraj_ref),qgam_ref(ntnc,ntraj_ref),qele_ref(ntnc,ntraj_ref),qnu_ref(ntnc,ntraj_ref),qfis_ref(ntnc,ntraj_ref),qalp_ref(ntnc,ntraj_ref), mexc_ref(ntnc,ntraj_ref) )
!     allocate( xmass_1day_ref(0:nz,ntraj_ref), xmass_10day_ref(0:nz,ntraj_ref) )


!     dims1(1) = ntraj
!     call h5dopen_f(file_id, "data_exists",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data_exists_ref, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
    
!     call h5dopen_f(file_id, "ye",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, ye_ref, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
    
!     call h5dopen_f(file_id, "entr",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, entr_ref, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
    
!     call h5dopen_f(file_id, "texp",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, texp_ref, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
 
    
!     dims2(1) = ntnc; dims2(2) = ntraj
!     call h5dopen_f(file_id, "qtot",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, qtot_ref, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
    
!     call h5dopen_f(file_id, "qgam",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, qgam_ref, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5dopen_f(file_id, "qele",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, qele_ref, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5dopen_f(file_id, "qnu",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, qnu_ref, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5dopen_f(file_id, "qfis",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, qfis_ref, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5dopen_f(file_id, "qalp",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, qalp_ref, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5dopen_f(file_id, "mexc",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, mexc_ref, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     dims2(1) = nz+1; dims2(2) = ntraj
!     call h5dopen_f(file_id, "xmass_1day",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, xmass_1day_ref, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5dopen_f(file_id, "xmass_10day",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, xmass_10day_ref, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5fclose_f(file_id, hdf_err)
    

!   end subroutine qtab_read_ref




!   subroutine qtab_write_hdf(fn_hdf)

!     use hdf5
!     use h5lt
    
!     character(*),intent(in) :: fn_hdf

!     ! hdf5 output
!     integer :: hdf_err,rank
!     integer(HID_T) :: file_id, dset_id, dspace_id
!     integer(HSIZE_T) :: dims1(1), dims2(2), dims3(3)
    
!     call h5open_f(hdf_err)
!     call h5fcreate_f(fn_hdf, H5F_ACC_TRUNC_F, file_id, hdf_err)

!     if(hdf_err/=0) stop "error"
    
!     ! original model name
!     call h5ltmake_dataset_string_f(file_id, "/original_model_name", original_model_name, hdf_err)

!     ! path to the directory of tracers
!     call h5ltmake_dataset_string_f(file_id, "/directory_trajectory", directory_trajectory, hdf_err)
    
!     ! path to the directory of heating rate
!     call h5ltmake_dataset_string_f(file_id, "/directory_heating_rate", directory_heating_rate, hdf_err)

!     ! name of the file used to fill data
!     call h5ltmake_dataset_string_f(file_id, "/file_used_to_fill_data", file_used_to_fill_data, hdf_err)

!     ! name of the file for dynamical ejecta (single-data)
!     call h5ltmake_dataset_string_f(file_id, "/file_dynamical", file_dynamical, hdf_err)
        
!     ! 1-dim
!     rank=1
    
!     dims1(1) = 1
!     ! ntraj
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/ntraj", H5T_NATIVE_INTEGER, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, ntraj, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/ntraj_inside", H5T_NATIVE_INTEGER, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, ntraj_inside, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! ntnc
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/ntnc", H5T_NATIVE_INTEGER, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, ntnc, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! nz
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/nz", H5T_NATIVE_INTEGER, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, nz, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     dims1(1) = ntnc
!     ! t - t_ncu0 array
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/t_ncl", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, t_ncl, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)


!     dims1(1) = ntraj
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/dm", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dm, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/ye", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, ye, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/entr", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, entr, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/texp", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, texp, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/ye_ext", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, ye_ext, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)


!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/t_ncl0", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, t_ncl0, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! extra. 
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/vr_ext", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vr_ext, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! if data exists
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/data_exists", H5T_NATIVE_INTEGER, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data_exists, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! if dynamical
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/is_dynamical", H5T_NATIVE_INTEGER, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, is_dynamical, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)


!     dims1(1) = 1
!     ! n_wing for dmdt
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/n_wing", H5T_NATIVE_INTEGER, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, n_wing, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! 2-dim

!     rank=2

!     dims2(1) = ntnc
!     dims2(2) = ntraj
!     ! qtot
!     call h5screate_simple_f(rank, dims2, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/qtot", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, qtot, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! qgam
!     call h5screate_simple_f(rank, dims2, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/qgam", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, qgam, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! qele
!     call h5screate_simple_f(rank, dims2, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/qele", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, qele, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! qnu
!     call h5screate_simple_f(rank, dims2, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/qnu", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, qnu, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! qfis
!     call h5screate_simple_f(rank, dims2, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/qfis", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, qfis, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! qalp
!     call h5screate_simple_f(rank, dims2, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/qalp", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, qalp, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! mexc
!     call h5screate_simple_f(rank, dims2, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/mexc", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mexc, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)
    
!     ! qtot_dmdt
!     call h5screate_simple_f(rank, dims2, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/qtot_dmdt", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, qtot_dmdt, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)


!     dims2(1) = nz+1
!     dims2(2) = ntraj
!     ! xmass_1day
!     call h5screate_simple_f(rank, dims2, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/xmass_1day", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, xmass_1day, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! xmass_10day
!     call h5screate_simple_f(rank, dims2, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/xmass_10day", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, xmass_10day, dims2, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)


!     ! additional information
!     if(flux_based==1)then
!        ! extraction
!        rank = 1
!        dims1(1) = 1
!        ! n_time_ext
!        call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!        call h5dcreate_f(file_id, "/n_time_ext", H5T_NATIVE_INTEGER, &
!             dspace_id,dset_id, hdf_err)
!        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, n_time_ext, dims1, hdf_err)
!        call h5dclose_f(dset_id, hdf_err)
!        call h5sclose_f(dspace_id, hdf_err)
!        ! n_theta_ext
!        call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!        call h5dcreate_f(file_id, "/n_theta_ext", H5T_NATIVE_INTEGER, &
!             dspace_id,dset_id, hdf_err)
!        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, n_theta_ext, dims1, hdf_err)
!        call h5dclose_f(dset_id, hdf_err)
!        call h5sclose_f(dspace_id, hdf_err)

!        ! r_ext
!        call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!        call h5dcreate_f(file_id, "/r_ext", H5T_NATIVE_DOUBLE, &
!             dspace_id,dset_id, hdf_err)
!        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, r_ext, dims1, hdf_err)
!        call h5dclose_f(dset_id, hdf_err)
!        call h5sclose_f(dspace_id, hdf_err)

!        ! dtheta
!        call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!        call h5dcreate_f(file_id, "/dtheta_ext", H5T_NATIVE_DOUBLE, &
!             dspace_id,dset_id, hdf_err)
!        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dtheta_ext, dims1, hdf_err)
!        call h5dclose_f(dset_id, hdf_err)
!        call h5sclose_f(dspace_id, hdf_err)

!        ! dims1(1) = n_theta_ext
!        ! ! n_phi_ext (it is a function of ith)
!        ! call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!        ! call h5dcreate_f(file_id, "/n_phi_ext", H5T_NATIVE_INTEGER, &
!        !      dspace_id,dset_id, hdf_err)
!        ! call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, n_phi_ext, dims1, hdf_err)
!        ! call h5dclose_f(dset_id, hdf_err)
!        ! call h5sclose_f(dspace_id, hdf_err)

!        dims1(1) = n_time_ext
!        ! extraction time
!        call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!        call h5dcreate_f(file_id, "/time_ext", H5T_NATIVE_DOUBLE, &
!             dspace_id,dset_id, hdf_err)
!        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, time_ext, dims1, hdf_err)
!        call h5dclose_f(dset_id, hdf_err)
!        call h5sclose_f(dspace_id, hdf_err)
       
!        dims1(1) = ntraj
!        ! ! time_ext
!        ! call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!        ! call h5dcreate_f(file_id, "/time_ext", H5T_NATIVE_DOUBLE, &
!        !      dspace_id,dset_id, hdf_err)
!        ! call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, time_ext, dims1, hdf_err)
!        ! call h5dclose_f(dset_id, hdf_err)
!        ! call h5sclose_f(dspace_id, hdf_err)
!        ! theta
!        call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!        call h5dcreate_f(file_id, "/theta_ext", H5T_NATIVE_DOUBLE, &
!             dspace_id,dset_id, hdf_err)
!        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, theta_ext, dims1, hdf_err)
!        call h5dclose_f(dset_id, hdf_err)
!        call h5sclose_f(dspace_id, hdf_err)
!        ! phi
!        call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!        call h5dcreate_f(file_id, "/phi_ext", H5T_NATIVE_DOUBLE, &
!             dspace_id,dset_id, hdf_err)
!        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, phi_ext, dims1, hdf_err)
!        call h5dclose_f(dset_id, hdf_err)
!        call h5sclose_f(dspace_id, hdf_err)

!        ! itime_ext
!        call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!        call h5dcreate_f(file_id, "/itime_ext", H5T_NATIVE_INTEGER, &
!             dspace_id,dset_id, hdf_err)
!        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, itime_ext, dims1, hdf_err)
!        call h5dclose_f(dset_id, hdf_err)
!        call h5sclose_f(dspace_id, hdf_err)
!        ! ! itheta
!        ! call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!        ! call h5dcreate_f(file_id, "/itheta_ext", H5T_NATIVE_INTEGER, &
!        !      dspace_id,dset_id, hdf_err)
!        ! call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, itheta_ext, dims1, hdf_err)
!        ! call h5dclose_f(dset_id, hdf_err)
!        ! call h5sclose_f(dspace_id, hdf_err)
!        ! ! iphi
!        ! call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!        ! call h5dcreate_f(file_id, "/iphi_ext", H5T_NATIVE_INTEGER, &
!        !      dspace_id,dset_id, hdf_err)
!        ! call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, iphi_ext, dims1, hdf_err)
!        ! call h5dclose_f(dset_id, hdf_err)
!        ! call h5sclose_f(dspace_id, hdf_err)

!        ! mass flux
!        call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!        call h5dcreate_f(file_id, "/mflux_ext", H5T_NATIVE_DOUBLE, &
!             dspace_id,dset_id, hdf_err)
!        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mflux_ext, dims1, hdf_err)
!        call h5dclose_f(dset_id, hdf_err)
!        call h5sclose_f(dspace_id, hdf_err)

!     endif
    
!     call h5fclose_f(file_id, hdf_err)
!     call h5close_f(hdf_err)
    
!     write(6,'("saved with name ",a)') trim(fn_hdf)
!   end subroutine qtab_write_hdf

!   subroutine qtab_mass_average_nucltime

!     integer :: ip
!     real(8) :: dm_tmp, mass_tot
    
!     logical,allocatable :: flag_avail(:)
    
!     call qtab_allocate_single
    
!     allocate(flag_avail(ntraj))

!     flag_avail(:) = .false.

!     if(flux_based==1)then
!        do ip=ntraj_inside+1,ntraj
!           if(data_exists(ip)==1.and.is_dynamical(ip)==0)then
!              flag_avail(ip) = .true.
!           endif
!        enddo
!     else
!        do ip=ntraj_inside+1,ntraj
!           if(data_exists(ip)==1.and.is_dynamical(ip)==1)then
!              flag_avail(ip) = .true.
!           endif
!        enddo
!     endif
    
!     mass_tot = sum(dm(:)*dble(data_exists(:)))
!     t_ncl0_single = minval(t_ncl0(:))

!     mass_avail = 0d0
!     mass_post  = 0d0
!     mass_dyn   = 0d0
!     mass_total = 0d0
!     t_ncl0_single = 1d99
    
!     do ip=ntraj_inside+1,ntraj
!        ! write(6,*) ip, data_exists(ip), is_dynamical(ip)
!        mass_total = mass_total + dm(ip)
          
!        if(is_dynamical(ip)==0)then
!           mass_post = mass_post + dm(ip)
!        else
!           mass_dyn = mass_dyn + dm(ip)
!        endif
       
!        if(flag_avail(ip))then
!           mass_avail = mass_avail + dm(ip)
!        endif
       
!     end do

!     ! mass-average (only the trajectory with data)

!     qtot_single(:)=0d0
!     qgam_single(:)=0d0
!     qele_single(:)=0d0
!     qnu_single (:)=0d0
!     qfis_single(:)=0d0
!     qalp_single(:)=0d0
!     mexc_single(:)=0d0
!     xmass_1day_single (:)=0d0
!     xmass_10day_single(:)=0d0
    
!     qtot_dmdt_single(:) = 0d0
!     ! t_ncl0_single = sum(dm(:)*t_ncl0(:))
!     do ip=1,ntraj
!        if(flag_avail(ip))then
!           qtot_single(:) = qtot_single(:) + dm(ip)*qtot(:,ip)
!           qgam_single(:) = qgam_single(:) + dm(ip)*qgam(:,ip)
!           qele_single(:) = qele_single(:) + dm(ip)*qele(:,ip)
!           qnu_single (:) = qnu_single (:) + dm(ip)*qnu (:,ip)
!           qfis_single(:) = qfis_single(:) + dm(ip)*qfis(:,ip)
!           qalp_single(:) = qalp_single(:) + dm(ip)*qalp(:,ip)
!           mexc_single(:) = mexc_single(:) + dm(ip)*mexc(:,ip)
!           xmass_1day_single (:) = xmass_1day_single (:) + dm(ip)*xmass_1day (:,ip)
!           xmass_10day_single(:) = xmass_10day_single(:) + dm(ip)*xmass_10day(:,ip)
!           ye_single = ye_single + dm(ip)*ye(ip)
!           entr_single = entr_single + dm(ip)*entr(ip)

!           qtot_dmdt_single(:) = qtot_dmdt_single(:) + dm(ip)*qtot_dmdt(:,ip)
!        endif
!     enddo
!     ! t_ncl0_single = t_ncl0_single / mass_tot
    
!     qtot_single(:) = qtot_single(:) / mass_avail
!     qgam_single(:) = qgam_single(:) / mass_avail
!     qele_single(:) = qele_single(:) / mass_avail
!     qnu_single (:) = qnu_single (:) / mass_avail
!     qfis_single(:) = qfis_single(:) / mass_avail
!     qalp_single(:) = qalp_single(:) / mass_avail
!     mexc_single(:) = mexc_single(:) / mass_avail
!     xmass_1day_single (:) = xmass_1day_single (:) / mass_avail
!     xmass_10day_single(:) = xmass_10day_single(:) / mass_avail

!     ye_single = ye_single / mass_avail
!     entr_single = entr_single / mass_avail

!     qtot_dmdt_single(:) = qtot_dmdt_single(:) / mass_avail

!     ! remove nz=110 and renormalize
!     xmass_1day_single (nz) = 0d0
!     xmass_10day_single(nz) = 0d0
!     xmass_1day_single (:) = xmass_1day_single (:) / sum(xmass_1day_single (:))
!     xmass_10day_single(:) = xmass_10day_single(:) / sum(xmass_10day_single(:))

!     deallocate(flag_avail)

!   end subroutine qtab_mass_average_nucltime

  
!   subroutine qtab_mass_average_hydrotime
    
!     integer :: ip, itnc, itnc_prime
!     real(8) :: mass_tot
!     real(8) :: tt0, tt1, t_ip

!     logical,allocatable :: flag_avail(:)

!     call qtab_allocate_single

!     allocate(flag_avail(ntraj))

!     flag_avail(:) = .false.

!     if(flux_based==1)then
!        do ip=ntraj_inside+1,ntraj
!           if(data_exists(ip)==1.and.is_dynamical(ip)==0)then
!              flag_avail(ip) = .true.
!           endif
!        enddo
!     else
!        do ip=ntraj_inside+1,ntraj
!           if(data_exists(ip)==1.and.is_dynamical(ip)==1)then
!              flag_avail(ip) = .true.
!           endif
!        enddo
!     endif
       
    
!     ! mass_tot = sum(dm(:)*dble(data_exists(:))*dble(1-is_dynamical(:)))
    
!     mass_avail = 0d0
!     mass_post  = 0d0
!     mass_dyn   = 0d0
!     mass_total = 0d0
!     t_ncl0_single = 1d99
    
!     do ip=ntraj_inside+1,ntraj
!        ! write(6,*) ip, data_exists(ip), is_dynamical(ip)
!        mass_total = mass_total + dm(ip)
          
!        if(is_dynamical(ip)==0)then
!           mass_post = mass_post + dm(ip)
!        else
!           mass_dyn = mass_dyn + dm(ip)
!        endif

!        if(flag_avail(ip))then
!           mass_avail = mass_avail + dm(ip)
!           t_ncl0_single = min(t_ncl0_single, t_ncl0(ip))
!        endif
!     end do

!     qtot_single(:)=0d0
!     qgam_single(:)=0d0
!     qele_single(:)=0d0
!     qnu_single (:)=0d0
!     qfis_single(:)=0d0
!     qalp_single(:)=0d0
!     mexc_single(:)=0d0
!     xmass_1day_single (:)=0d0
!     xmass_10day_single(:)=0d0
!     ye_single = 0d0
!     entr_single = 0d0
    
!     qtot_dmdt_single(:)=0d0

!     ! time origin is the minimum t_ncl0 available.
!     ! t_sim = t_ncl(itnc) + t_ncl0_min is the simulation time.
!     ! t_sim - t_ncl0(ip) is then the nucleosynthesis time for individual trajectory
!     do ip=ntraj_inside+1,ntraj
!        if(flag_avail(ip))then
!           ! write(6,*) ip,t_ncl0(ip)
!           do itnc=1,ntnc
             
!              ! (time - t_ncl0_single) + t_ncl0
!              t_ip = t_ncl(itnc) + t_ncl0_single - t_ncl0(ip)

!              ! this is buggy expression !
!              ! t_ip = t_ncl(itnc) - t_ncl0_single + t_ncl0(ip)
             
!              if    ( t_ip < t_ncl(1) )then
!                 ! contribution of the trajectory for which the nucleosynthesis is not yet started is zero.
!                 itnc_prime = 1
!                 tt1 = 0d0
!                 tt0 = 0d0
!              elseif( t_ip > t_ncl(ntnc) )then
!                 itnc_prime = ntnc-1
!                 tt1 = 1d0
!                 tt0 = 0d0
!              else
!                 itnc_prime = 1
!                 search_itnc_prime:do
!                    if( t_ncl(itnc_prime) <= t_ip .and. t_ip <= t_ncl(itnc_prime+1) )then
!                       exit search_itnc_prime
!                    endif
!                    itnc_prime = itnc_prime + 1
!                 enddo search_itnc_prime
                
!                 ! log interpolation
!                 tt1 = ( log10(t_ip) - log10(t_ncl(itnc_prime)) ) / ( log10(t_ncl(itnc_prime+1)) - log10(t_ncl(itnc_prime)) )
!                 tt0 = 1d0 - tt1
!              endif
             
!              qtot_single(itnc) = qtot_single(itnc) + (tt0*qtot(itnc_prime,ip)+tt1*qtot(itnc_prime+1,ip))*dm(ip)
!              qgam_single(itnc) = qgam_single(itnc) + (tt0*qgam(itnc_prime,ip)+tt1*qgam(itnc_prime+1,ip))*dm(ip)
!              qele_single(itnc) = qele_single(itnc) + (tt0*qele(itnc_prime,ip)+tt1*qele(itnc_prime+1,ip))*dm(ip)
!              qnu_single (itnc) = qnu_single (itnc) + (tt0*qnu (itnc_prime,ip)+tt1*qnu (itnc_prime+1,ip))*dm(ip)
!              qfis_single(itnc) = qfis_single(itnc) + (tt0*qfis(itnc_prime,ip)+tt1*qfis(itnc_prime+1,ip))*dm(ip)
!              qalp_single(itnc) = qalp_single(itnc) + (tt0*qalp(itnc_prime,ip)+tt1*qalp(itnc_prime+1,ip))*dm(ip)
!              mexc_single(itnc) = mexc_single(itnc) + (tt0*mexc(itnc_prime,ip)+tt1*mexc(itnc_prime+1,ip))*dm(ip)
             
!              qtot_dmdt_single(itnc) = qtot_dmdt_single(itnc) + (tt0*qtot_dmdt(itnc_prime,ip)+tt1*qtot_dmdt(itnc_prime+1,ip))*dm(ip)

!              ! write(6,'(i5,99es15.7)') itnc, t_ip + t_ncl0(ip), &
!              !      (tt0*qtot(itnc_prime,ip)+tt1*qtot(itnc_prime+1,ip)), &
!              !      (tt0*qtot_dmdt(itnc_prime,ip)+tt1*qtot_dmdt(itnc_prime+1,ip))
!           enddo
          
!           ye_single = ye_single + ye(ip)*dm(ip)
!           entr_single = entr_single + entr(ip)*dm(ip)
!           xmass_1day_single  (:) = xmass_1day_single  (:) + xmass_1day  (:,ip)*dm(ip)
!           xmass_10day_single (:) = xmass_10day_single (:) + xmass_10day (:,ip)*dm(ip)

!        endif
!     enddo

!     qtot_single(:) = qtot_single(:) / mass_avail
!     qgam_single(:) = qgam_single(:) / mass_avail
!     qele_single(:) = qele_single(:) / mass_avail
!     qnu_single (:) = qnu_single (:) / mass_avail
!     qfis_single(:) = qfis_single(:) / mass_avail
!     qalp_single(:) = qalp_single(:) / mass_avail
!     mexc_single(:) = mexc_single(:) / mass_avail
!     xmass_1day_single (:) = xmass_1day_single (:) / mass_avail
!     xmass_10day_single(:) = xmass_10day_single(:) / mass_avail

!     ye_single = ye_single / mass_avail
!     entr_single = entr_single / mass_avail

!     qtot_dmdt_single(:) = qtot_dmdt_single(:) / mass_avail

!     ! remove nz=110 and renormalize
!     xmass_1day_single (nz) = 0d0
!     xmass_10day_single(nz) = 0d0
!     xmass_1day_single (:) = xmass_1day_single (:) / sum(xmass_1day_single (:))
!     xmass_10day_single(:) = xmass_10day_single(:) / sum(xmass_10day_single(:))

!     deallocate(flag_avail)
    
!   end subroutine qtab_mass_average_hydrotime


!   subroutine qtab_write_hdf_single(fn_hdf)

!     use hdf5
!     use h5lt
    
!     character(*),intent(in) :: fn_hdf

!     ! hdf5 output
!     integer :: hdf_err,rank
!     integer(HID_T) :: file_id, dset_id, dspace_id
!     integer(HSIZE_T) :: dims1(1), dims2(2), dims3(3)
    
!     call h5open_f(hdf_err)
!     call h5fcreate_f(fn_hdf, H5F_ACC_TRUNC_F, file_id, hdf_err)
    
!     ! original model name
!     call h5ltmake_dataset_string_f(file_id, "/original_model_name", original_model_name, hdf_err)

!     ! path to the directory of tracers
!     call h5ltmake_dataset_string_f(file_id, "/directory_trajectory", directory_trajectory, hdf_err)

!     ! path to the directory of heating rate
!     call h5ltmake_dataset_string_f(file_id, "/directory_heating_rate", directory_heating_rate, hdf_err)

!     ! 1-dim
!     rank=1

!     dims1(1) = 1
!     ! ntnc
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/ntnc", H5T_NATIVE_INTEGER, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, ntnc, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! nz
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/nz", H5T_NATIVE_INTEGER, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, nz, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! t_nucl0
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/t_ncl0", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, t_ncl0_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! mass
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/mass_avail", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mass_avail, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/mass_post", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mass_post, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/mass_dyn", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mass_dyn, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/mass_total", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mass_total, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! Ye
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/ye", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, ye_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)
!     ! S
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/entr", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, entr_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! n_wing for dmdt
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/n_wing", H5T_NATIVE_INTEGER, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, n_wing, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)


!     dims1(1) = ntnc
!     ! t - t_ncu0 array
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/t_ncl", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, t_ncl, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! heating rate
!     ! qtot
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/qtot", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, qtot_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! qgam
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/qgam", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, qgam_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! qele
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/qele", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, qele_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! qnu
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/qnu", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, qnu_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! qfis
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/qfis", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, qfis_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! qalp
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/qalp", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, qalp_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! mexc
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/mexc", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, mexc_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! qtot_dmdt
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/qtot_dmdt", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, qtot_dmdt_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)


!     dims1(1) = nz+1
!     ! xmass_1day
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/xmass_1day", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, xmass_1day_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

!     ! xmass_10day
!     call h5screate_simple_f(rank, dims1, dspace_id, hdf_err)
!     call h5dcreate_f(file_id, "/xmass_10day", H5T_NATIVE_DOUBLE, &
!          dspace_id,dset_id, hdf_err)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, xmass_10day_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5sclose_f(dspace_id, hdf_err)

    
!     call h5fclose_f(file_id, hdf_err)
!     call h5close_f(hdf_err)
    
!     write(6,'("saved with name ",a)') trim(fn_hdf)
    
!   end subroutine qtab_write_hdf_single
  
!   subroutine qtab_read_single(fn)
    
!     use hdf5
    
!     character(*),intent(in) :: fn
    
!     ! hdf5 output
!     integer :: hdf_err,rank
!     integer(HID_T) :: file_id, dset_id, dspace_id
!     integer(HSIZE_T) :: dims1(1), dims2(2), dims3(3)

!     write(6,'("read qtable with single data : ",a)') trim(fn)

!     file_dynamical = trim(fn)
    
!     call qtab_allocate_single
    

!     call h5fopen_f(fn, H5F_ACC_RDONLY_F, file_id, hdf_err)
    
    
!     dims1(1) = 1
!     call h5dopen_f(file_id, "t_ncl0",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, t_ncl0_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

    
!     dims1(1) = ntnc
!     call h5dopen_f(file_id, "qtot",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, qtot_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
    
!     call h5dopen_f(file_id, "qgam",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, qgam_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5dopen_f(file_id, "qele",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, qele_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5dopen_f(file_id, "qnu",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, qnu_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5dopen_f(file_id, "qfis",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, qfis_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5dopen_f(file_id, "qalp",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, qalp_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5dopen_f(file_id, "mexc",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, mexc_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     dims1(1) = nz+1
!     call h5dopen_f(file_id, "xmass_1day",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, xmass_1day_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)
!     call h5dopen_f(file_id, "xmass_10day",dset_id, hdf_err)
!     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, xmass_10day_single, dims1, hdf_err)
!     call h5dclose_f(dset_id, hdf_err)

!     call h5fclose_f(file_id, hdf_err)
    

!   end subroutine qtab_read_single



  
!   subroutine qtab_check_dynamical
!     integer :: ip
    
!     do ip=1,ntraj
!        if(ye_ext(ip) < 0.08d0 )then
!           is_dynamical(ip) = 1
!        else
!           is_dynamical(ip) = 0
!        endif
!        ! if(data_exists(ip)==1) write(6,*) ip, is_dynamical(ip), data_exists(ip), ye_ext(ip), ye(ip)
!     enddo
    
!   end subroutine qtab_check_dynamical



!   subroutine qtab_fill_data_dynamical
    
!     use utils

!     integer :: ip

!     ! block
!     !   integer :: itnc
!     !   do itnc=1,ntnc
!     !      write(6,'(i5,99es12.4)') itnc, t_ncl(itnc), qtot_single(itnc)
!     !   enddo
      
!     ! end block
    
!     do ip=1,ntraj
!        ! write(6,*) ip, is_dynamical(ip)
!        if(is_dynamical(ip) == 1)then
          
!           qtot(:,ip) = qtot_single(:)
!           qgam(:,ip) = qgam_single(:)
!           qele(:,ip) = qele_single(:)
!           qnu (:,ip) = qnu_single (:)
!           qfis(:,ip) = qfis_single(:)
!           qalp(:,ip) = qalp_single(:)
!           mexc(:,ip) = mexc_single(:)
!           xmass_1day(:,ip) = xmass_1day_single(:)
!           xmass_10day(:,ip) = xmass_10day_single(:)
          
!        endif
!     enddo

!   end subroutine qtab_fill_data_dynamical




!   subroutine qtab_fill_data_self
    
!     use utils

!     integer :: ip, ip_closest
!     real(8) :: dist_min
!     integer,allocatable :: ip_closest_traj(:)

!     file_used_to_fill_data="self"

!     allocate( ip_closest_traj(ntraj) )

!     do ip=1,ntraj
!        if(is_dynamical(ip) == 0 .and. data_exists(ip) == 0)then
          
!           call find_closest_p(ye(ip), entr(ip), texp(ip), &
!                ye, entr, texp, data_exists, ntraj, &
!                ip_closest, dist_min)
          
!           ip_closest_traj(ip) = ip_closest
!           write(99,'(2i6,99es12.4)') ip, ip_closest, dist_min, ye(ip), entr(ip), texp(ip), &
!                ye(ip_closest), entr(ip_closest), texp(ip_closest), dble(data_exists(ip_closest))
!           if(dist_min > 1d0) write(99,'(2i6,99es12.4)') ip, ip_closest, dist_min, ye(ip), entr(ip), texp(ip), &
!                ye(ip_closest), entr(ip_closest), texp(ip_closest), dble(data_exists(ip_closest))

!           if(ip == ip_closest) write(6,*) "itself!", ip
!        else
!           ip_closest_traj(ip) = 0
!        endif
!     enddo

!     do ip=1,ntraj
!        if(ip_closest_traj(ip) > 0)then ! i.e., post-merger & data does not exists

!           ip_closest = ip_closest_traj(ip)
          
!           qtot(:,ip) = qtot(:,ip_closest)
!           qgam(:,ip) = qgam(:,ip_closest)
!           qele(:,ip) = qele(:,ip_closest)
!           qnu (:,ip) = qnu (:,ip_closest)
!           qfis(:,ip) = qfis(:,ip_closest)
!           qalp(:,ip) = qalp(:,ip_closest)
!           mexc(:,ip) = mexc(:,ip_closest)
!           xmass_1day(:,ip) = xmass_1day(:,ip_closest)
!           xmass_10day(:,ip) = xmass_10day(:,ip_closest)
          
!        endif
!     enddo
    
!     deallocate( ip_closest_traj )
    
!   end subroutine qtab_fill_data_self




!   subroutine qtab_fill_data_ref
    
!     use utils

!     integer :: ip, ip_closest
!     real(8) :: dist_min
!     integer,allocatable :: ip_closest_traj(:)
    
!     allocate( ip_closest_traj(ntraj) )

!     do ip=1,ntraj
!        if(is_dynamical(ip) == 0 .and. data_exists(ip) == 0)then
          
!           call find_closest_p(ye(ip), entr(ip), texp(ip), &
!                ye_ref, entr_ref, texp_ref, data_exists_ref, ntraj_ref, &
!                ip_closest, dist_min)
          
!           ip_closest_traj(ip) = ip_closest
!           write(99,'(2i6,99es12.4)') ip, ip_closest, dist_min, ye(ip), entr(ip), texp(ip), &
!                ye_ref(ip_closest), entr_ref(ip_closest), texp_ref(ip_closest), dble(data_exists_ref(ip_closest))
!           if(dist_min > 1d0) write(99,'(2i6,99es12.4)') ip, ip_closest, dist_min, ye(ip), entr(ip), texp(ip), &
!                ye_ref(ip_closest), entr_ref(ip_closest), texp_ref(ip_closest), dble(data_exists_ref(ip_closest))
          
!        else
!           ip_closest_traj(ip) = 0
!        endif
!     enddo

!     do ip=1,ntraj
!        if(ip_closest_traj(ip) > 0)then ! i.e., post-merger & data does not exists

!           ip_closest = ip_closest_traj(ip)
          
!           qtot(:,ip) = qtot_ref(:,ip_closest)
!           qgam(:,ip) = qgam_ref(:,ip_closest)
!           qele(:,ip) = qele_ref(:,ip_closest)
!           qnu (:,ip) = qnu_ref (:,ip_closest)
!           qfis(:,ip) = qfis_ref(:,ip_closest)
!           qalp(:,ip) = qalp_ref(:,ip_closest)
!           mexc(:,ip) = mexc_ref(:,ip_closest)
!           xmass_1day(:,ip) = xmass_1day_ref(:,ip_closest)
!           xmass_10day(:,ip) = xmass_10day_ref(:,ip_closest)
          
!        endif
!     enddo
    
!     deallocate( ip_closest_traj )
    
!   end subroutine qtab_fill_data_ref

  
!   subroutine qtab_test_output
!     integer :: ip
    
!     real(8),allocatable :: time_edge(:), time(:)
!     real(8) :: t1, dt1, ratio
    
!     integer :: itime, ntime

!     integer :: itheta_ip, iphi_ip
!     real(8) :: theta, phi

!     ! t1 = 2.011231000000000d-002
!     ! dt1 = 5.188599999999998d-003
!     ! ratio = 1.04540815936347d0

!     ! ntime = 55
!     ! allocate(time_edge(ntime+1), time(ntime))
    

!     ! do itime=1,ntime
!     !    time(itime) = t1 + dt1*(ratio**(itime-1)-1d0)/(ratio-1d0)
!     ! enddo
    
!     ! time_edge(1) = 0d0
!     ! do itime=2,ntime
!     !    time_edge(itime) = 0.5d0*(time(itime) + time(itime-1))
!     ! enddo
    
!     ! do itime=1,ntime
!     !    write(6,'(99es12.4)') time(itime), time_edge(itime), time_edge(itime+1)
!     ! enddo
!     ! stop

!     open(10,file="ye_qtable_orig.dat", status="replace", action="write")
    
!     do ip=ntraj_inside+1,ntraj
!        if(ip>ntraj_inside+1)then
!           if(itime_ext(ip) /= itime_ext(ip-1))then
!              write(10,*)
!              write(10,*)
!           endif
!        endif
!        ! itheta_ip = itheta_ext(ip)
!        ! iphi_ip   = iphi_ext(ip)

!        ! theta = dtheta_ext*(dble(itheta_ip-1)+0.5d0)
!        ! phi   = 2d0*pi*dble(iphi_ip-1)/dble(n_phi_ext(itheta_ip))
       
!        theta = theta_ext(ip)
!        phi   = phi_ext(ip)

!        write(10,'(99es12.4)') time_ext(itime_ext(ip)), theta, phi, ye(ip)
!     enddo
!     close(10)

!   end subroutine qtab_test_output


!   subroutine qtab_interp_distance(x_in, y_in, z_in, time_in, n_rank, &
!        ip_used, x_dyn_out, t_ncl0_out, qtot_out, qgam_out, qele_out, qnu_out , qfis_out, qalp_out, mexc_out, xmass_1day_out, xmass_10day_out, ye_out, flag_test)
    
!     use utils

!     real(8),intent(in) :: x_in, y_in, z_in, time_in
!     integer,intent(in) :: n_rank
!     integer,intent(out) :: ip_used(n_rank)
!     real(8),intent(out) ::  x_dyn_out, t_ncl0_out, &
!          qtot_out(ntnc), qgam_out(ntnc), qele_out(ntnc), qnu_out (ntnc), qfis_out(ntnc), qalp_out(ntnc), mexc_out(ntnc), &
!          xmass_1day_out(0:nz), xmass_10day_out(0:nz)
!     real(8),intent(out) :: ye_out
    
!     logical,intent(in),optional :: flag_test

!     real(8),allocatable :: dist(:)

!     integer :: ip, itheta, itime_prim, itime_ip!, itheta_ip, iphi_ip
!     real(8) :: vr_ip, dt_ip, theta_ip, phi_ip, dr_ip, x_ip, y_ip, z_ip, rho_ip, flux_ip
  
!     ! real(8),allocatable :: theta_ext(:), phi_ext(:)

!     integer :: irank
!     integer,allocatable :: ind_rank(:)
!     real(8),allocatable :: weight(:)
  
!     real(8) :: t_ncl0_min
!     integer :: itnc, itnc_prime
!     real(8) :: t_ip, tt1, tt0
!     ! allocate( theta_ext(n_theta_ext) )

!     ! do itheta=1,n_theta_ext
!     !    theta_ext(itheta) = dtheta_ext*(dble(itheta_ip-1)+0.5d0)
!     ! enddo

!     itime_prim = 1
!     do while(time_ext(itime_prim+1) < time_in .and. itime_prim+1<n_time_ext)
!        itime_prim = itime_prim + 1
!     enddo
       
!     allocate(dist(ntraj))
!     dist(:) = 1d99
!     do ip=ntraj_inside+1,ntraj

!        itime_ip  = itime_ext(ip)
!        !itheta_ip = itheta_ext(ip)
!        !iphi_ip   = iphi_ext(ip)
!        vr_ip     = vr_ext(ip)

!        theta_ip = theta_ext(ip)
!        phi_ip   = phi_ext(ip)
       
!        if( abs(itime_ip-itime_prim)<=1 )then
          
!           dt_ip    = time_in - time_ext(itime_ip)
          
!           !theta_ip = dtheta_ext*(dble(itheta_ip-1)+0.5d0)
!           !phi_ip   = 2d0*pi*dble(iphi_ip-1)/dble(n_phi_ext(itheta_ip))
          
!           dr_ip    = vr_ip*dt_ip
          
!           x_ip = (r_ext-dr_ip)*sin(theta_ip)*cos(phi_ip)
!           y_ip = (r_ext-dr_ip)*sin(theta_ip)*sin(phi_ip)
!           z_ip = (r_ext-dr_ip)*cos(theta_ip)
          
!           dist(ip) = sqrt( (x_in-x_ip)**2 + (y_in-y_ip)**2 + (z_in-z_ip)**2 )
!           ! write(97,'(4i5,99es12.4)') ip, itime_ip, itheta_ip, iphi_ip, dist(ip), dr_ip, vr_ip, x_ip, y_ip, z_ip, theta_ip, phi_ip
          
!           ! if( time_ext(itime_ip) < 0.1d0 .and. is_dynamical(ip)==0)then
!           !    dist(ip) = 1d99
!           ! endif
          
!        endif
       
       
!     enddo

!     allocate(ind_rank(n_rank), weight(n_rank))

!     call index_rank(n_rank, ntraj, dist, ind_rank, -1d0)
    
!     do irank=1,n_rank
!        weight(irank) = 1d0/(dist(ind_rank(irank))+1d-99)
!     enddo

!     ! do irank=1,n_rank
!     !    ip = ind_rank(irank)

!     !    itime_ip  = itime_ext(ip)
!     !    !itheta_ip = itheta_ext(ip)
!     !    !iphi_ip   = iphi_ext(ip)
!     !    vr_ip     = vr_ext(ip)

!     !    theta_ip = theta_ext(ip)
!     !    phi_ip   = phi_ext(ip)
       
       
!     !    dt_ip    = time_in - time_ext(itime_ip)
       
!     !    dr_ip    = vr_ip*dt_ip
       
!     !    flux_ip = r_ext**2/(r_ext+dr_ip)**2 * mflux_ext(ip)
!     !    rho_ip = flux_ip/vr_ip
       
!     !    ! weight(irank) = weight(irank) * rho_ip
!     !    ! weight(irank) = weight(irank) * flux_ip
!     ! enddo

!     weight(:) = weight(:)/sum(weight(:))
    
!     do irank=1,n_rank
!        ip_used(irank) = ind_rank(irank)
!     enddo
!     if(present(flag_test))then
!        if(flag_test)then
!           write(6,'(99i6)') ip_used(:)
!           do irank=1,n_rank
!              write(6,'(i6,99es14.6)') ip_used(irank), weight(irank), t_ncl0(ip_used(irank))
!           enddo
!        endif
!     endif
    

!     x_dyn_out   = 0d0
!     t_ncl0_out  = 0d0
!     qtot_out(:) = 0d0
!     qgam_out(:) = 0d0
!     qele_out(:) = 0d0
!     qnu_out (:) = 0d0
!     qfis_out(:) = 0d0
!     qalp_out(:) = 0d0
!     mexc_out(:) = 0d0
!     xmass_1day_out (:) = 0d0
!     xmass_10day_out(:) = 0d0
!     ye_out = 0d0

!     t_ncl0_min = 1d99
!     do irank=1,n_rank
!        ip = ind_rank(irank)
!        t_ncl0_min = min(t_ncl0_min, t_ncl0(ip))
!        x_dyn_out = x_dyn_out + dble(is_dynamical(ip))*weight(irank)
!     enddo

!     if(present(flag_test))then
!        if(flag_test)then
!           write(6,'("t_ncl0_min=",99es12.4)') t_ncl0_min
!        endif
!     endif
    
!     t_ncl0_out = t_ncl0_min
   
!     do irank=1,n_rank
!        ip = ind_rank(irank)

!        do itnc=1,ntnc
             
!           ! (time - t_ncl0_single) + t_ncl0
!           t_ip = t_ncl(itnc) + t_ncl0_min - t_ncl0(ip)

          

!           if    ( t_ip < t_ncl(1) )then
!              ! contribution of the trajectory for which the nucleosynthesis is not yet started is zero.
!              itnc_prime = 1
!              tt1 = 0d0
!              tt0 = 0d0
!           elseif( t_ip > t_ncl(ntnc) )then
!              itnc_prime = ntnc-1
!              tt1 = 1d0
!              tt0 = 0d0
!           else
!              itnc_prime = 1
!              search_itnc_prime:do
!                 if( t_ncl(itnc_prime) <= t_ip .and. t_ip <= t_ncl(itnc_prime+1) )then
!                    exit search_itnc_prime
!                 endif
!                 itnc_prime = itnc_prime + 1
!              enddo search_itnc_prime
             
!              ! log interpolation
!              ! tt1 = ( log10(t_ip) - log10(t_ncl(itnc_prime)) ) / ( log10(t_ncl(itnc_prime+1)) - log10(t_ncl(itnc_prime)) )
!              ! linear
!              tt1 = ( (t_ip) - (t_ncl(itnc_prime)) ) / ( (t_ncl(itnc_prime+1)) - (t_ncl(itnc_prime)) )
             
!              tt0 = 1d0 - tt1
!           endif
          
!           qtot_out(itnc) = qtot_out(itnc) + (tt0*qtot(itnc_prime,ip)+tt1*qtot(itnc_prime+1,ip))*weight(irank)
!           qgam_out(itnc) = qgam_out(itnc) + (tt0*qgam(itnc_prime,ip)+tt1*qgam(itnc_prime+1,ip))*weight(irank)
!           qele_out(itnc) = qele_out(itnc) + (tt0*qele(itnc_prime,ip)+tt1*qele(itnc_prime+1,ip))*weight(irank)
!           qnu_out (itnc) = qnu_out (itnc) + (tt0*qnu (itnc_prime,ip)+tt1*qnu (itnc_prime+1,ip))*weight(irank)
!           qfis_out(itnc) = qfis_out(itnc) + (tt0*qfis(itnc_prime,ip)+tt1*qfis(itnc_prime+1,ip))*weight(irank)
!           qalp_out(itnc) = qalp_out(itnc) + (tt0*qalp(itnc_prime,ip)+tt1*qalp(itnc_prime+1,ip))*weight(irank)
!           mexc_out(itnc) = mexc_out(itnc) + (tt0*mexc(itnc_prime,ip)+tt1*mexc(itnc_prime+1,ip))*weight(irank)
          
!           if(present(flag_test))then
!              if(flag_test.and.1d0 < t_ncl(itnc) + t_ncl0_min .and. t_ncl(itnc) + t_ncl0_min < 3d0)then
!                 block
!                   integer :: itnc_d, itnc_u
!                   integer :: n_w
!                   real(8) :: dmdt
!                   n_w=1
!                   itnc_d = max(itnc_prime-n_w, 1)
!                   itnc_u = min(itnc_prime+n_w, ntnc)
!                   ! dmdt = (mexc (itnc_u,ip) - mexc(itnc_d,ip)) / (t_ncl(itnc_u) - t_ncl(itnc_d))
!                   dmdt = ( (tt0*mexc(itnc_u,ip)+tt1*mexc(itnc_u+1,ip)) - (tt0*mexc(itnc_d,ip)+tt1*mexc(itnc_d+1,ip))) / (t_ncl(itnc_u) - t_ncl(itnc_d))
                  
!                   write(6,'(3i6,99es14.6)') ip, itnc, itnc_prime, tt1, t_ncl(itnc) + t_ncl0_min, t_ncl(itnc), &
!                        tt0*qtot(itnc_prime,ip)+tt1*qtot(itnc_prime+1,ip), -dmdt * clight**2, tt0*mexc(itnc_prime,ip)+tt1*mexc(itnc_prime+1,ip)
!                 end block
!              endif
!           endif

!        enddo

!        ! qtot_out(:) = qtot_out(:) + qtot(:,ip)*weight(irank)
!        ! qgam_out(:) = qgam_out(:) + qgam(:,ip)*weight(irank)
!        ! qele_out(:) = qele_out(:) + qele(:,ip)*weight(irank)
!        ! qnu_out (:) = qnu_out (:) + qnu (:,ip)*weight(irank)
!        ! qfis_out(:) = qfis_out(:) + qfis(:,ip)*weight(irank)
!        ! qalp_out(:) = qalp_out(:) + qalp(:,ip)*weight(irank)
!        ! mexc_out(:) = mexc_out(:) + mexc(:,ip)*weight(irank)
       
!        xmass_1day_out(:) = xmass_1day_out(:) + xmass_1day(:,ip)*weight(irank)
!        xmass_10day_out(:) = xmass_10day_out(:) + xmass_10day(:,ip)*weight(irank)

!        ye_out = ye_out + ye(ip)*weight(irank)

!     if(present(flag_test))then
!        if(flag_test)then

!           do itnc=1,ntnc
!              if(1d0 < t_ncl(itnc) + t_ncl0_min .and. t_ncl(itnc) + t_ncl0_min < 3d0)then
!                 block
!                   integer :: itnc_d, itnc_u
!                   integer :: n_w
!                   real(8) :: dmdt
!                   n_w=10
!                   itnc_d = max(itnc-n_w, 1)
!                   itnc_u = min(itnc+n_w, ntnc)
!                   dmdt = (mexc_out (itnc_u) - mexc_out(itnc_d)) / (t_ncl(itnc_u) - t_ncl(itnc_d))

!                   write(6,'(99es12.4)') t_ncl(itnc) + t_ncl0_min, t_ncl(itnc), qtot_out(itnc), -dmdt * clight**2, mexc_out (itnc)
!                 end block
!              endif
!           enddo
!        endif
!     endif

!  enddo

!     if(present(flag_test))then
!        if(flag_test)then
!           write(6,*)
!           do itnc=1,ntnc
!              if(1d0 < t_ncl(itnc) + t_ncl0_min .and. t_ncl(itnc) + t_ncl0_min < 3d0)then
!                 block
!                   integer :: itnc_d, itnc_u
!                   integer :: n_w
!                   real(8) :: dmdt
!                   n_w=10
!                   itnc_d = max(itnc-n_w, 1)
!                   itnc_u = min(itnc+n_w, ntnc)
!                   dmdt = (mexc_out (itnc_u) - mexc_out(itnc_d)) / (t_ncl(itnc_u) - t_ncl(itnc_d))

!                   write(6,'(99es12.4)') t_ncl(itnc) + t_ncl0_min, t_ncl(itnc), qtot_out(itnc), -dmdt * clight**2, mexc_out (itnc)
!                 end block
!              endif
!           enddo
!        endif
!     endif

!     xmass_1day_out(:) = xmass_1day_out(:) / sum(xmass_1day_out(:))
!     xmass_10day_out(:) = xmass_10day_out(:) / sum(xmass_10day_out(:))
    
!     ! write(6,'(4i5,5f12.4,4es12.4)') ip_used(:), ye(ip_used(:)), ye_out, mflux_ext(ip_used(:))
!     ! do irank=1,n_rank
!     !    write(6,*) ind_rank(irank), is_dynamical(ind_rank(irank))
!     ! enddo
    
!     deallocate(ind_rank)
!     deallocate(dist)
    
!   end subroutine qtab_interp_distance

!   subroutine qtab_interp_timeangle(theta_in, phi_in, time_in, &
!        ip_used, x_dyn_out, t_ncl0_out, qtot_out, qgam_out, qele_out, qnu_out , qfis_out, qalp_out, mexc_out, xmass_1day_out, xmass_10day_out, ye_out)
!     ! interpolate the q and so on with respect to theta, phi and time of the injection.
!     real(8),intent(in) :: theta_in, phi_in, time_in
!     integer,intent(out) ::  ip_used(8)
!     real(8),intent(out) ::  x_dyn_out, t_ncl0_out, &
!          qtot_out(ntnc), qgam_out(ntnc), qele_out(ntnc), qnu_out (ntnc), qfis_out(ntnc), qalp_out(ntnc), mexc_out(ntnc), &
!          xmass_1day_out(0:nz), xmass_10day_out(0:nz)
!     real(8),intent(out) :: ye_out
    
!     integer :: ip

!     integer :: it_0, it_1
!     real(8) :: tt_0, tt_1
    
!     ip_used(:) = 0
!     x_dyn_out = 1d-99; t_ncl0_out = 1d99
!     qtot_out(:)=0d0; qgam_out(:)=0d0; qele_out(:)=0d0; qnu_out(:)=0d0; qfis_out(:)=0d0; qalp_out(:)=0d0; mexc_out(:)=0d0
!     xmass_1day_out(0:nz)=0d0; xmass_10day_out(0:nz)=0d0
!     ye_out=0d0

!     ! search theta and phi indices
!     ! assume HEALPix discretization
!     ! 1) search time
!     block
!       integer :: it_ext
!       if    (time_in < time_ext(1))then
!          it_0 = 1
!          it_1 = 2
!       elseif(time_ext(n_time_ext) <= time_in)then
!          it_0 = n_time_ext-1
!          it_1 = n_time_ext
!       else
!          do it_ext = 1,n_time_ext-1
!             if(time_ext(it_ext) <= time_in .and. time_in < time_ext(it_ext+1))then
!                it_0 = it_ext
!                it_1 = it_ext + 1
!                exit
!             endif
!          enddo
!       endif
!       tt_1 = min(1d0, max(0d0, (time_in-time_ext(it_0))/(time_ext(it_1)-time_ext(it_0))))
!       tt_0 =  1d0 - tt_1
!     end block

!     ! 2) search angle
!     block
!       integer :: ip_t_1, ip_t_2
!       integer :: ip_t_theta_1, ip_t_theta_2
!       integer :: ip_000,ip_001,ip_010,ip_011,ip_100,ip_101,ip_110,ip_111
      
!       real(8) :: theta_prime

!       ! for t=t0 slice
!       ip_t_1 = 0
!       do ip=ntraj_inside+1,ntraj
!          if(itime_ext(ip)==it_0)then
!             ip_t_1 = ip
!             exit
!          endif
!       enddo
!       ip_t_2 = 0
!       do ip=ntraj,ntraj_inside+1,-1
!          if(itime_ext(ip)==it_0)then
!             ip_t_2 = ip
!             exit
!          endif
!       enddo

!       write(6,*) "theta=", theta_in
!       write(6,*) theta_ext(ip_t_1), theta_ext(ip_t_2)
!       write(6,*)

!       ! get theta index
!       theta_prime = max(theta_in, theta_ext(ip_t_1)) * (1d0-1d-5)
!       do ip=ip_t_1,ip_t_2-1
!          !write(6,*) theta_prime, theta_ext(ip)
!          if(theta_ext(ip) <= theta_prime .and. theta_prime < theta_ext(ip+1))then
!             ip_t_theta_1 = ip
!             exit
!          endif
!       enddo

!       theta_prime = min(theta_in, theta_ext(ip_t_2)) * (1d0+1d-5)
!       do ip=ip_t_2,ip_t_1,-1
!          !write(6,*) theta_prime, theta_ext(ip)
!          if(theta_ext(ip) <= theta_prime .and. theta_prime < theta_ext(ip+1))then
!             ip_t_theta_2 = ip
!             exit
!          endif
!       enddo

!       ! in the same theta, find phi
!       if(phi_in < phi_ext(ip_t_theta_1) .or. phi_ext(ip_t_theta_2) <= phi_in )then
!          ip_000 = ip_t_theta_1
!          ip_001 = ip_t_theta_2
!       else
!          do ip=ip_t_theta_1,ip_t_theta_2-1
!             if(phi_ext(ip) <= phi_in .and. phi_in < phi_ext(ip+1) )then
!                ip_000 = ip
!                ip_001 = ip+1
!             endif
!          enddo
!       endif

!       stop
      
!       theta_prime = min(theta_in, theta_ext(ip_t_2)) * (1d0+1d-5)
!       do ip=ip_t_2,ip_t_1,-1
!          write(6,*) theta_prime, theta_ext(ip)
!          if( theta_prime > theta_ext(ip))then
!             ip_t_theta_2 = ip
!             exit
!          endif
!       enddo

!       do ip=ip_t_theta_1, ip_t_theta_2
!          write(6,'(i5,i5,99es12.4)') ip, itime_ext(ip), theta_ext(ip), phi_ext(ip)
!       enddo
      
!       stop
      
!     end block
    
    
!   end subroutine qtab_interp_timeangle
  
!   subroutine qtab_index_passer(ntnc_out, nz_out, n_theta_out, n_phi_out, n_time_out)
  
!     integer,intent(out) :: ntnc_out, nz_out, n_theta_out, n_phi_out, n_time_out

!     ntnc_out = ntnc
!     nz_out = nz
!     n_theta_out = n_theta_ext
!     n_phi_out  = 72!n_phi_ext(n_theta_ext)
!     n_time_out = n_time_ext
    
!   end subroutine qtab_index_passer
  
!   subroutine qtab_array_passer(r_ext_out, dtheta_out, t_ncl_out, time_ext_out)
!     real(8),intent(out) :: r_ext_out, dtheta_out, t_ncl_out(ntnc), time_ext_out(n_time_ext)
    
!     r_ext_out = r_ext
!     dtheta_out = dtheta_ext
!     t_ncl_out(:) = t_ncl(:)
!     time_ext_out(:) = time_ext(:)
    
!   end subroutine qtab_array_passer
!   ! subroutine qtab_add_info(dir_traj)

!   !   use utils

!   !   character(*),intent(in) :: dir_traj


!   !   character(256) :: fn

!   !   fn = trim(dir_traj) // "/parameters.dat"
!   !   call param_parser(fn, &
!   !        n_theta = n_theta_ext)
    
!   !   write(6,*) n_theta_ext
!   !   allocate(n_phi_ext(n_theta_ext))
    
!   !   call pset_angle(n_theta_ext, pi/2d0, n_phi_ext)
!   !   dtheta_ext = pi/2d0/dble(n_theta_ext)
    
!   ! end subroutine qtab_add_info
  
! !   subroutine read_qtable(fn,if_uni2d)
! !     character(*),intent(in) :: fn
! !     logical,intent(in) :: if_uni2d
    
! ! !!! read data file
! !     write(6,'(a)') fn
! !     open(10,file=fn,status='old',action="read",form='binary')
! !     read(10) ntraj,ntnc,nz
    
! !     write(6,*) ntraj,ntnc,nz

! !     allocate(t_ncl(ntnc))
! !     allocate(qtot(ntnc,ntraj),qgam(ntnc,ntraj),qele(ntnc,ntraj),qnu(ntnc,ntraj),qfis(ntnc,ntraj),qalp(ntnc,ntraj),t_ncl0(ntraj),ye(ntraj),entr(ntraj), &
! !          mexc(ntnc,ntraj),xmass_1day(0:nz,ntraj),xmass_10day(0:nz,ntraj))
    
! !     read(10) t_ncl
! !     read(10) qtot,qgam,qele,qnu,qfis,qalp,t_ncl0,ye,entr
! !     read(10) xmass_1day,xmass_10day,mexc

! !     if(if_uni2d)then
! !        read(10) qtable_n_r,qtable_n_theta
! !        allocate(qtable_r(qtable_n_r),qtable_theta(qtable_n_theta),qtable_ip_rt(qtable_n_r,qtable_n_theta))
! !        read(10) qtable_r, qtable_theta
! !        read(10) qtable_ip_rt
! !     endif

! !     close(10)
    
! !   end subroutine read_qtable

!   ! subroutine output_qtable
!   !   integer :: i_r, i_theta
!   !   integer :: ip

!   !   do i_r=1,qtable_n_r
!   !      do i_theta=1,qtable_n_theta
!   !         ip = qtable_ip_rt(i_r,i_theta)
          
!   !         write(6,'(99es12.4)') qtable_r(i_r),qtable_theta(i_theta), ye(ip), entr(ip), sum(xmass_1day(57:71,ip)), sum(xmass_1day(:,ip))
          
!   !      enddo
!   !   enddo

!   ! end subroutine output_qtable

!   subroutine qtab_qtot_from_dmdt(n_wing_in)
!     integer,intent(in) :: n_wing_in
    
!     integer :: itnc, ip, itnc_d, itnc_u
!     real(8) :: dlogm, dlogt, dmdt

!     n_wing = n_wing_in

!     do ip=1,ntraj
       
!        do itnc=1,ntnc
          
!           itnc_d = max(itnc-n_wing, 1)                
!           itnc_u = min(itnc+n_wing, ntnc)
          
          
!           dmdt = (mexc (itnc_u,ip) - mexc(itnc_d,ip)) / (t_ncl(itnc_u) - t_ncl(itnc_d))
          
!           qtot_dmdt(itnc,ip) = -dmdt * clight**2
!           qtot(itnc,ip) = max(qtot(itnc,ip), qtot_dmdt(itnc,ip))
          
!        enddo
       
!     enddo


 
!   end subroutine qtab_qtot_from_dmdt


end module module_trajectory
