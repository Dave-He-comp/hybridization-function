!
! routine for computing vector cross-product
!
subroutine vcross(a,b,c)
  implicit real*8(a-h,o-z)
  dimension a(3), b(3), c(3)
  
  a(1) = b(2) * c(3) - b(3) * c(2)
  a(2) = b(3) * c(1) - b(1) * c(3)
  a(3) = b(1) * c(2) - b(2) * c(1)
  return
end subroutine vcross

subroutine self_overlap(ispin,iwk,sys_nband,&
           sys_nplane,env_nband,env_nplane)
  implicit none
  complex*8,  allocatable     :: sys_coeff(:,:), env_coeff(:,:)
  real*8,     allocatable     :: freqx(:), freqy(:), freqz(:)
  integer,    allocatable     :: sys_igall(:,:,:), env_igall(:,:,:)
  integer                     :: i,j,k, ndimx, ndimy, ndimz, nrecl
  integer                     :: g_mesh(3), sys_file
  integer                     :: nthread
  real*8                      :: dtmp1, dtmp2, dtmp3, dtmp4
  real*8                      :: sys_occ, env_occ
  integer                     :: ntmp1, ntmp2, m, n
  integer                     :: sys_nband, sys_nplane, env_nband
  integer                     :: env_nplane
  integer                     :: istat, gall_tmp(3), sys_igall_tmp(3)
  integer                     :: ispin, iwk, sys_iband, sys_iplane
  integer                     :: env_iband, env_iplane
  complex*16                  :: coupling_tmp1, coupling_tmp2, cdtmp1
  complex*16                  :: cdtmp2
  complex*8                   :: sys_eig, env_eig
  !
  integer :: isp, nk, nb, ni, nl, itemp
  integer :: ix, iy, iz, nx, ny, nz, tid
  integer :: nspin, nband, nwk, nionp, lm, env_file
  real*8  :: sys_fermi, env_fermi, delta_scaling
  real*8, allocatable :: vkpt(:,:), wtkpt(:)
  real*8, allocatable :: sys_eig_tot(:,:), sys_occ_tot(:,:)
  real*8, allocatable :: env_eig_tot(:,:), env_occ_tot(:,:)

  character(len=2)  :: ch_spin
  character(len=5)  :: ch_wk
  character(len=5)  :: ch_sys_band, ch_env_band
  character(len=75) :: sys_title, env_title, title
  character(len=256) :: ch_tmp, buffer
  character(len=80) ::  dump

  logical :: l_chosen_sys_bands, l_sys_file_exist, lfile_exist
  logical :: lfile_open, lz_skip, ltmp, file_exists
  integer :: sys_nband_save, env_nband_save, nz_start, nz_end
  integer :: num_chosen_sys_bands, num_chosen_env_bands
  integer, allocatable :: chosen_sys_bands(:), chosen_env_bands(:)
  !
  real*8  :: a1(3), a2(3), a3(3)
  real*8  :: b1(3), b2(3), b3(3)
  real*8  :: vtmp(3), sumkg(3), Vcell, z_up_limit, z_low_limit
  real*8  :: delx, dely, delz, wkx, wky, wkz, x, y, z
  real*8, allocatable :: space_mesh(:,:,:,:)
  complex*16, allocatable :: fft_wavefunc(:,:,:), coupling(:,:)
  !
  namelist / wavefunc / nx, ny, nz, nz_start, nz_end, &
                        lz_skip, z_up_limit, z_low_limit, sys_iband
  !dimension a1(3), a2(3), a3(3)
  !dimension b1(3), b2(3), b3(3)
  !dimension vtmp(3), sumkg(3), wk(3)
       
  ! constant 'c' below is 2m/hbar**2 in units of 1/eV Ang^2 
  ! (value is adjusted in final decimal places to agree with VASP value;
  ! program checks for discrepancy of any results between this and VASP
  ! values)
  
  !data c/0.262465831d0/ 
  ! data c/0.26246582250210965422d0/ 
  !pi = 4.0d0 * atan(1.0d0)
        
  ! read input of system
       
  nrecl = 24
  nspin = 2
  delta_scaling = 1.d0
  lm = 16 ! number of s+p+d orbtials
  open(unit=64,file='G_MAX')
  rewind(64)
  read(64,*) (g_mesh(i), i=1,3)
  close(64)
  do i=1,3
  !  if(mod(g_mesh(i),2) .eq. 0) then
  !    g_mesh(i) = 2 * g_mesh(i)
  !  else
      g_mesh(i) = 2 * g_mesh(i) - 1
  !  end if
  end do
  if (sys_nplane .ne. env_nplane) then
    write(6,*) 'fatal error!', &
               ' wrong size of sys_nplane and env_nplane'
    write(6,*)  sys_nplane, env_nplane
    write(6,*) 'please cheak the cell of POSCAR'
    stop
  end if
  !
  inquire(file='fermi_sys', exist=file_exists)
  if (file_exists) then
     title = 'fermi_sys'
     open(unit=63, file=title)
     rewind(63)
     read(63,'(A)') buffer
     read(63,'(A)') buffer
     read(buffer(1:23),'(A10,F12.7)') ch_tmp, sys_fermi
     close(63)
     write(6,*) 'Fermi level of system: ', sys_fermi
     write(6,*)
  end if
  !
  inquire(file='fermi_sub', exist=file_exists)
  if (file_exists) then
      title = 'fermi_sub'
      open(unit=64, file=title)
      rewind(64)
      read(64,'(A)') buffer
      read(64,'(A)') buffer
      read(buffer(1:23),'(A10,F12.7)') ch_tmp, env_fermi
      close(64)
      write(6,*) 'Fermi level of sub: ', env_fermi
      write(6,*)
  end if
  !
  inquire(file='fermi_tip', exist=file_exists)
  if (file_exists) then
      title = 'fermi_tip'
      open(unit=64, file=title)
      rewind(64)
      read(64,'(A)') buffer
      read(64,'(A)') buffer
      read(buffer(1:23),'(A10,F12.7)') ch_tmp, env_fermi
      close(64)
      write(6,*) 'Fermi level of tip: ', env_fermi
      write(6,*)
  end if
  !
  if (ispin .eq. 1) then
     inquire(file='to_be_read_sys_bands_sp1', exist=lfile_exist)
  else if (ispin .eq. 2) then
     inquire(file='to_be_read_sys_bands_sp2', exist=lfile_exist)
  end if
  !
  if (lfile_exist) then
     num_chosen_sys_bands = -1
     istat = 0
     write(6,*) 'read the chosen system band file'
     if (ispin .eq. 1) then
        open(unit=29, file='to_be_read_sys_bands_sp1')
     else if (ispin .eq. 2) then
        open(unit=29, file='to_be_read_sys_bands_sp2')
     end if
     rewind(30)
     rewind(29)
     do while(.not. is_iostat_end(istat))
       read(29,'(a)', iostat=istat) ch_tmp
       num_chosen_sys_bands = num_chosen_sys_bands + 1
     end do
     write(6,*) 'total number of chosen system bands',&
                 num_chosen_sys_bands
     !
     sys_nband_save = sys_nband
     if (&!num_chosen_sys_bands .le. sys_nband .and. &
         num_chosen_sys_bands .gt. 0) then
         sys_nband = num_chosen_sys_bands
     end if
     allocate(chosen_sys_bands(sys_nband))
     allocate(sys_eig_tot(sys_nband,2), sys_occ_tot(sys_nband,2))
     rewind(29)
     do i=1, num_chosen_sys_bands
       read(29,620) chosen_sys_bands(i),&
                    sys_occ_tot(i,1),sys_eig_tot(i,1),&
                    sys_occ_tot(i,2),sys_eig_tot(i,2)
       write(6,*) 'read number of system band',&
                   chosen_sys_bands(i)
     end do
     close(29)
     write(6,*)
  end if
  !
  if (ispin .eq. 1) then
     inquire(file='to_be_read_env_bands_sp1', exist=lfile_exist)
  else if (ispin .eq. 2) then
     inquire(file='to_be_read_env_bands_sp2', exist=lfile_exist)
  end if
  !
  if (lfile_exist) then
     num_chosen_env_bands = -1
     istat = 0
     write(6,*) 'read the chosen environment band file'
     if (ispin .eq. 1) then
        open(unit=30, file='to_be_read_env_bands_sp1')
     else if (ispin .eq. 2) then
        open(unit=30, file='to_be_read_env_bands_sp2')
     end if
     rewind(30)
     do while(.not. is_iostat_end(istat))
       read(30,'(a)', iostat=istat) ch_tmp
       num_chosen_env_bands = num_chosen_env_bands + 1
     end do
     write(6,*) 'total number of chosen environment bands',&
                 num_chosen_env_bands
     if (&!num_chosen_env_bands .le. env_nband .and. &
         num_chosen_env_bands .gt. 0) then
        env_nband = num_chosen_env_bands
     end if
     allocate(chosen_env_bands(env_nband))
     allocate(env_eig_tot(env_nband,2), env_occ_tot(env_nband,2))
     rewind(30)
     do i=1, num_chosen_env_bands
       read(30,620) chosen_env_bands(i),&
                    env_occ_tot(i,1),env_eig_tot(i,1),&
                    env_occ_tot(i,2),env_eig_tot(i,2)
       write(6,*) 'read number of environment bands',&
                   chosen_env_bands(i)
     end do
     close(30)
     write(6,*)
  end if
  !
  620 format('No. of band',1x,I5,2x,&
             ',occ. spin-1',2x,F9.6,2x,&
             ',eigenergy',2x,F11.6,2x,&
             ',occ. spin-2',2x,F9.6,2x,&
             ',eigenergy',2x,F11.6)
  !
  coupling_tmp1 = 0.d0
  coupling_tmp2 = 0.d0
  ch_spin = ''
  ch_wk = ''
  ch_sys_band = ''
  sys_title = ''
  ch_env_band = ''
  env_title = ''
  sys_file = 101
  env_file = 100
  cdtmp1 = 0.d0
  cdtmp2 = 0.d0
  write(ch_spin,'(I1)') ispin
  write(ch_wk,'(I5)') iwk
  write(6,*) 'to be calculate number of system bands', sys_nband
  write(6,*) 'to be calculate number of enviornment bands', env_nband
  allocate(sys_coeff(sys_nplane,sys_nband))
  allocate(sys_igall(3,sys_nplane,sys_nband))
  allocate(env_coeff(env_nplane,env_nband))
  allocate(env_igall(3,env_nplane,env_nband))
  allocate(coupling(sys_nband,env_nband))
  sys_coeff = 0
  env_coeff = 0
  sys_igall = 0
  env_igall = 0
  !
  write(6,*) 'reading environment wavecar'
  do m=1, env_nband
    env_iband = chosen_env_bands(m)
    ! reading environment wavecar
    write(ch_env_band,'(I5)') env_iband
    env_title = ('spin-'//trim(adjustl(ch_spin))//'-wk-'//&
                 trim(adjustl(ch_wk))//'-band-'//trim(adjustl(ch_env_band)))
    call system('cp ./ENV/'//trim(adjustl(env_title))//&
                '   ENV-'//trim(adjustl(env_title)))
    !
    write(6,40) env_iband
    open(unit=env_file,file=('ENV-'//trim(adjustl(env_title))),&
         iostat=istat)
    if (istat .ne. 0) then
       write(6,*) 'Error when reading ', env_title
       write(6,*)
       continue
    end if
    rewind(env_file)
    ! skip 
    do i=1, 31
       select case (i)
          case (28) ! read eig-energy
             read(env_file,*) dtmp1, dtmp2
             env_eig = dcmplx(dtmp1,dtmp2)
          case (30) ! read occupancy
             read(env_file,*) env_occ
          case default
             read(env_file,'(a)') ch_tmp
       end select
    end do
    write(6,*) 'distance from fermi level &
                of environment band', env_iband
    write(6,*) dreal(env_eig - env_fermi)
    do i=1, env_nplane
       read(env_file,*) env_igall(1,i,m),&
                        env_igall(2,i,m),&
                        env_igall(3,i,m),&
                        dtmp1,dtmp2
       env_coeff(i,m) = dcmplx(dtmp1,dtmp2)
    end do
    close(env_file)
    write(6,41) env_iband
    write(6,*)
    !
    call system('rm ENV-'//trim(adjustl(env_title)))
    env_file = env_file + 2
  end do
  !
  write(6,*) 'reading system wavecar'
  do n = 1, sys_nband
     ! reading system wavecar
     sys_iband = chosen_sys_bands(n)
     write(ch_sys_band,'(I5)') sys_iband
     sys_title = ('spin-'//trim(adjustl(ch_spin))//'-wk-'//&
                  trim(adjustl(ch_wk))//'-band-'//trim(adjustl(ch_sys_band)))
     call system('cp ./SYS/'//trim(adjustl(sys_title))//&
                 '   SYS-'//trim(adjustl(sys_title)))
     write(6,42) sys_iband
     !end if
     open(unit=sys_file,file=('SYS-'//trim(adjustl(sys_title))),&
          iostat=istat)
     if (istat .ne. 0) then
        write(6,*) 'Error when reading ', sys_title
        write(6,*)
        continue
     end if
     rewind(sys_file)
     ! skip 
     do i=1, 31
        select case (i)
           case (28) ! read eig-energy
              read(sys_file,*) dtmp1, dtmp2
              sys_eig = dcmplx(dtmp1,dtmp2)
           case (30) ! read occupancy
              read(sys_file,*) sys_occ
           case (13) ! read a1
              read(sys_file,*) (a1(j),j=1,3)
           case (15) ! read a2
              read(sys_file,*) (a2(j),j=1,3)
           case (17) ! read a3
              read(sys_file,*) (a3(j),j=1,3)
           case (20) ! read b1
              read(sys_file,*) (b1(j),j=1,3)
           case (22) ! read b2  
              read(sys_file,*) (b2(j),j=1,3)
           case (24) ! read b3  
              read(sys_file,*) (b3(j),j=1,3)
           case default
              read(sys_file,'(a)') ch_tmp
        end select
     end do
     write(6,*) 'distance from fermi level of system band', sys_iband
     write(6,*) dreal(sys_eig - sys_fermi)
     do i=1, sys_nplane
        read(sys_file,*) sys_igall(1,i,n),&
                         sys_igall(2,i,n),&
                         sys_igall(3,i,n),&
                         dtmp1,dtmp2
        sys_coeff(i,n) = dcmplx(dtmp1,dtmp2)
     end do
     close(sys_file)
     write(6,43) sys_iband
     write(6,*)
     sys_file = sys_file + 2
     call system('rm SYS-'//trim(adjustl(sys_title)))
  end do
  !
  40 format('read ENV bands  ',I6)
  41 format('close ENV bands  ',I6)
  42 format('read SYS bands  ',I6)
  43 format('close SYS bands  ',I6)
  !
  write(6,*) 'real space lattice vectors:'
  write(6,*) 'a1 =',((a1(j)),j=1,3)
  write(6,*) 'a2 =',((a2(j)),j=1,3)
  write(6,*) 'a3 =',((a3(j)),j=1,3)
  write(6,*) ' '
  call vcross(vtmp,a2,a3)
  Vcell=a1(1)*vtmp(1)+a1(2)*vtmp(2)+a1(3)*vtmp(3)
  write(6,*) 'volume unit cell =',(Vcell)
  write(6,*) ' '
  write(6,*) 'reciprocal lattice vectors:'
  write(6,*) 'b1 =',((b1(j)),j=1,3)
  write(6,*) 'b2 =',((b2(j)),j=1,3)
  write(6,*) 'b3 =',((b3(j)),j=1,3)
  write(6,*) ' '
  !
  write(6,*) 'finish reading'
  write(6,*)
  !
  call flush(6)
  x = 0
  y = 0
  z = 0
  nx = g_mesh(1)
  ny = g_mesh(2)
  nz = g_mesh(3)
  nz_start = 1
  nz_end = nz
  lz_skip = .False.
  sys_iband = 1
  ni = 1
  z_up_limit = (a1(3) + a2(3) + a3(3))/dble(g_mesh(3)) &
             * (0.01d0+dble(g_mesh(3)))
  z_low_limit = -0.01d0*(a1(3) + a2(3) + a3(3))/dble(g_mesh(3))
  !
  open(unit=5,file='input')
  rewind(5)
  read(5, wavefunc, end=103)
  103 continue
  close(5)
  delx = (a1(1) + a2(1) + a3(1)) / (nx - 1)
  dely = (a1(2) + a2(2) + a3(2)) / (ny - 1)
  delz = (a1(3) + a2(3) + a3(3)) / (nz - 1)
  !
  !
  !ltmp1 = .False.
  !ltmp2 = .False.
  !if (nz_end .gt. nz) then
  !  ltmp1 = .True.
  !  write(6,*) 'does not define or define wrong nz_end'
  !  write(6,*) 'program will automatically check nz_end'
  !end if
  if (nz_start .lt. 1) then
  !  ltmp2 = .True.
    nz_start = 1
    write(6,*) 'does not define or define wrong nz_start'
    write(6,*) 'program will automatically check nz_start'
  end if
  !
!  ltmp1 = .False.
!  ltmp2 = .False.
  if (z_up_limit .gt. (a1(3)+a2(3)+a3(3)+delz)) then
!    ltmp1 = .True.
    z_up_limit = (a1(3)+a2(3)+a3(3)+0.01d0*delz)
    write(6,*) 'does not define or define wrong z_up_limit'
    write(6,*) 'program will automatically check z_up_limit'
  end if
  if (z_low_limit .lt. -delz) then
!    ltmp2 = .True.
    z_low_limit = -0.01d0 * delz
    write(6,*) 'does not define or define wrong z_low_limit'
    write(6,*) 'program will automatically check z_low_limit'
  end if
  !
  if (sys_nband .eq. 1) then
    sys_iband = chosen_sys_bands(1)
    ni = 1
  else
    do ni=1, sys_nband
      itemp = chosen_sys_bands(ni)
      if (itemp .eq. sys_iband) then
        exit
      end if
    end do
  end if
  write(6,'("x varying from 0.d0 to ",E15.6," by step ", I5, " with length ", E15.6)')&
  (a1(1) + a2(1) + a3(1)), nx, delx
  write(6,'("y varying from 0.d0 to ",E15.6," by step ", I5, " with length ", E15.6)')&
  (a1(2) + a2(2) + a3(2)), ny, dely
  write(6,'("z varying from 0.d0 to ",E15.6," by step ", I5, " with length ", E15.6)')&
  (a1(3) + a2(3) + a3(3)), nz, delz
  write(6,'("upper limit in the z dierction ",F12.6)') z_up_limit
  write(6,'("lower limit in the z dierction ",F12.6)') z_low_limit
  !write(6,*)
  !
  !if (ltmp) then
  !  if (nz_end .gt. nz) then
  !    nz_end = nz
  !  end if
  !  if (nz_start .lt. 1) then
  !   nz_start = 1
  !  end if
  !  z_up_limit = (nz_end-1) * delz
  !  z_low_limit = (nz_start-1) * delz
  !  write(6,'("upper limit in the z dierction ",F12.6)') z_up_limit
  !  write(6,'("lower limit in the z dierction ",F12.6)') z_low_limit
  !else
  !  z_up_limit = (nz_end-1) * delz
  !  z_low_limit = (nz_start-1) * delz
  !  write(6,'("upper limit in the z dierction ",F12.6)') z_up_limit
  !  write(6,'("lower limit in the z dierction ",F12.6)') z_low_limit
  !end if
  !write(6,*)
  !
  open(unit=23,file='freqx')
  rewind(23)
  allocate(freqx(nx))
  read(23,*) freqx
  close(23)
  write(6,*) 'read freqx done'
  !
  open(unit=24,file='freqy')
  rewind(24)
  allocate(freqy(ny))
  read(24,*) freqy
  close(24)
  write(6,*) 'read freqy done'
  !
  open(unit=25,file='freqz')
  rewind(25)
  allocate(freqz(nz))
  read(25,*) freqz
  close(25)
  write(6,*) 'read freqz done'
  !
  allocate(fft_wavefunc(nx,ny,nz))
  fft_wavefunc = 0
  !
  do n=1, sys_nband
    sys_iband = chosen_sys_bands(n)
    sys_file = sys_file + 2
    write(ch_sys_band,'(I5)') sys_iband
    sys_title = ('wavefunc-fft-spin-'//trim(adjustl(ch_spin))//'-wk-'//&
                 trim(adjustl(ch_wk))//'-band-'//trim(adjustl(ch_sys_band)))
    open(unit=sys_file,file=sys_title)
    rewind(sys_file)
    fft_wavefunc = 0.d0
    do iz=1, nz
      do iy=1, ny
        do ix=1, nx
          read(sys_file,*) dtmp1,dtmp2
          fft_wavefunc(ix,iy,iz) = dcmplx(dtmp1,dtmp2)
        end do
      end do
    end do
    close(sys_file)
    do m=1, env_nband
      env_iband = chosen_env_bands(m)
      cdtmp1 = 0
      !
      do i=1, sys_nplane
        wkx = sys_igall(1,i,n)
        wky = sys_igall(2,i,n)
        wkz = sys_igall(3,i,n)
        !
        ix = -1
        do nk = 1, nx
          if (int(freqx(nk)) .eq. wkx) then
            ix = nk
            exit
          end if
        end do
        if (ix < 0) then
          write(6,*) 'wkx', wkx, ' do not find'
          stop
        end if
        !
        iy = -1
        do nk = 1, ny
          if (int(freqy(nk)) .eq. wky) then
            iy = nk
            exit
          end if
        end do
        if (iy < 0) then
          write(6,*) 'wky', wky, ' do not find'
          stop
        end if
        !
        iz = -1
        do nk = 1, nz
          if (int(freqz(nk)) .eq. wkz) then
            iz = nk
            exit
          end if
        end do
        if (iz < 0) then
          write(6,*) 'wkz', wkz, ' do not find'
          stop
        end if
        !
        cdtmp2 =  conjg(env_coeff(i,m))/dsqrt(Vcell) * fft_wavefunc(ix,iy,iz) 
        cdtmp1 = cdtmp1 + cdtmp2 * (sys_eig_tot(n,ispin))
        !write(6,*)'env_coeff',conjg(env_coeff(i,m))
        !write(6,*)'Vcell',dsqrt(Vcell)
        !write(6,'(3(I6,1x),2(E14.5))')ix,iy,iz,real(fft_wavefunc(ix,iy,iz)),imag(fft_wavefunc(ix,iy,iz))
        !write(6,*) (sys_eig_tot(n,ispin))
        !stop
      end do
      coupling(n,m) = cdtmp1
      write(6,'("No. sys_band ",I4,1x,"No. env_band ",I4,1x,"overlap ",E14.5,2x)') &
                        sys_nband, env_iband, &
                        cdabs(cdtmp1)!,cdabs(cdtmp1)**2
    end do
  end do
  !
  570 format('computation information', 3i10)
  580 format(3F12.7, F20.10, ' s')
  !
  env_file = env_file + 2
  open(unit=env_file,file=('spin-'//trim(adjustl(ch_spin))//'-coupling-details'))
  write(env_file,*) 'total hopping coefficient is ',&
                     dreal(sum(coupling)), dimag(sum(coupling))
  write(env_file,*) 'its absolute value is ', cdabs(sum(coupling))
  write(env_file,*) 'sum the abs values of hopping coefficients ',&
                     sum(cdabs(coupling))
  write(env_file,*)
  write(env_file,*) 'details of hopping between system and &
                     evnironment bands'
  write(env_file,*) '  No. sys       No. env ',&
                      '       eig       real',&
                      '       imag      abs'
  dtmp1 = 0.d0
  dtmp2 = 0.d0
  dtmp3 = 0.d0
  do m=1, env_nband
    env_iband = chosen_env_bands(m)
    do n=1, sys_nband
      sys_iband = chosen_sys_bands(n)
      dtmp1 = dtmp1 + dreal(coupling(n,m))
      dtmp2 = dtmp2 + dimag(coupling(n,m))
      dtmp3 = dtmp3 + cdabs(coupling(n,m))
      write(env_file,590) sys_iband, env_iband, &
                          (env_eig_tot(m,ispin)-env_fermi),&
                          dreal(coupling(n,m)), &
                          dimag(coupling(n,m)), &
                          cdabs(coupling(n,m))
    end do
  end do
  write(env_file,591) dtmp1, dtmp2, dtmp3 
  write(env_file,*) 
  close(env_file)
  !
  590 format(I10, 2X, I10, 4E18.10)
  591 format(5X,"summation",6X,2X,3E18.10)
end subroutine self_overlap

program main
  implicit none
  integer  :: ncpus, ispin, nspin, iwk
  integer  :: sys_band, sys_plane
  integer  :: sub_band, sub_plane
  integer  :: tip_band, tip_plane
  logical  :: lfile_exist
  
  nspin = 2
  iwk = 1
  ispin = 1
  ncpus = 32
  !open(unit=10,file='number_cpus')
  !rewind(10)
  !read(10,*) ncpus
  !close(10)
  !
  inquire(file='sys_para', exist=lfile_exist)
  if (lfile_exist) then
     open(unit=7,file='sys_para')
     read(7,*) sys_band 
     read(7,*) sys_plane
     close(7)
  end if
  !
  inquire(file='tip_para', exist=lfile_exist)
  if (lfile_exist) then
     open(unit=8,file='tip_para')
     read(8,*) tip_band
     read(8,*) tip_plane
     close(8)
  end if
  !
  inquire(file='sub_para', exist=lfile_exist)
  if (lfile_exist) then
     open(unit=9,file='sub_para')
     read(9,*) sub_band
     read(9,*) sub_plane
     close(9)
  end if
  !
  do ispin=1, nspin
     write(6,*)
     write(6,*) '************************'
     write(6,*) 'calculating for spin', ispin
     inquire(file='OUTCAR_TIP', exist=lfile_exist)
     if (lfile_exist) then
       call system('mv ./TIP ./ENV')
       call self_overlap(ispin,iwk,sys_band,&
            sys_plane,tip_band,tip_plane)
       call system('mv ./ENV ./TIP')
     end if
     inquire(file='OUTCAR_SUB', exist=lfile_exist)
     if (lfile_exist) then
       call system('mv ./SUB ./ENV')
       call self_overlap(ispin,iwk,sys_band,&
            sys_plane,sub_band,sub_plane)
       call system('mv ./ENV ./SUB')
     end if
     write(6,*) '************************'
     write(6,*) 
  end do
  !
end program main
