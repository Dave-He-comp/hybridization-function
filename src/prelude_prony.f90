program main
  implicit none
  integer  :: i, n, m, ntmp1, ntmp2, istat
  integer  :: sys_nband, env_nband, num_ww
  real*8   :: dtmp1, dtmp2
  real*8   :: sys_fermi, sub_fermi
  real*8   :: ww_i, ww_f, ww_step, width, heom_width
  integer,allocatable :: chosen_env_bands(:)
  real*8, allocatable :: ww(:), hyb(:)
  real*8, parameter   :: pi = 3.1415926d0
  !real*8, allocatable :: sys_eig(:,:), sys_occ(:,:)
  real*8, allocatable :: env_eig(:,:), env_occ(:,:)
  real*8, allocatable :: hopping(:,:), sum_hopping(:)
  character(len=128)  :: buffer, ch_tmp, title
  logical :: lfile_exist1, lfile_exist2
  namelist / hybridization / ww_i, ww_f, num_ww, width
  !
  integer :: j, npade, num_tt
  real*8  :: T, sigma(2), crit
  real*8  :: tt_f, tt_step, tt
  real*8, parameter :: hbar = 0.658211928d0
  complex*16 :: ctmp1, ctmp2
  complex*16, parameter :: cunit = (1.0d0, 0.0d0), &
                           eye = (0.0d0, 1.0d0), &
                           czero = (0.0d0, 0.0d0)
  complex*16, allocatable :: band_center(:), coupling(:)
  complex*16, allocatable :: zpeigv(:), zpcoef(:), zppole(:,:)
  complex*16, allocatable :: eta_tmp(:,:,:), gama_tmp(:,:,:)
  complex*16, allocatable :: ct_tot(:,:), ct(:,:,:)
  data sigma /1, -1/
  namelist / pade / npade, T, tt_f, num_tt, crit

  ww_i = -10
  ww_f =  10  
  num_ww = 10001
  allocate(ww(num_ww), hyb(num_ww))
  ww = 0.d0
  hyb = 0.d0
  width = 0.03d0
  open(unit=5,file='input')
  rewind(5)
  read(5, hybridization, end=83)
  83 continue
  close(5)
  heom_width = 0.3d0
  !
  npade = 1000
  T = 0.001d0
  tt_f = 200.d0
  num_tt = 10001
  crit = 1.d-8
  open(unit=5,file='input')
  rewind(5)
  read(5, pade, end=93)
  93 continue
  close(5)
  if (mod(num_tt,2) .eq. 0) then
    write(6,*) 'warning: num_tt should be odd'
    num_tt = num_tt + 1
  end if
  tt_step = tt_f / (num_tt - 1)
  allocate(zpeigv(npade), zpcoef(npade), zppole(npade,2))
  zpeigv = czero
  zpcoef = czero
  zppole = czero
  call zpout_psd(npade, 2, zpeigv, zpcoef, zppole)
  !
  write(6,*) '*********************************************************'
  write(6,*) '          generate hybdridization funtion     '
  write(6,*) '*********************************************************'
  write(6,'("Width of Lorentzian functions to approximate delta functions",E13.5)') width
  write(6,'("Width of Lorentzian function in HEOM calculations", E13.5)') heom_width
  write(6,'("Number of pade poles", E13.5)') npade
  write(6,'("Temperature of slab", E13.5)') T
  write(6,'("Final time of C(t)", E13.5)') tt_f
  write(6,'("Number of time stpe ", I7)') num_tt
  write(6,'("Time stpe of C(t)", E13.5)') tt_step
  write(6,'("Crition for minimal hopping strength", E13.5)') crit
  write(6,*)
  !
  sys_fermi = 0.d0
  sub_fermi = 0.d0
  inquire(file='fermi_sys', exist=lfile_exist2)
  if (lfile_exist2) then
     title = 'fermi_sys'
     open(unit=63, file=title)
     rewind(63)
     read(63,'(A)') buffer
     read(63,'(A)') buffer
     read(buffer(1:23),'(A10,F12.7)') ch_tmp, sys_fermi
     close(63)
     write(6,*) 'Fermi level of system: ', sys_fermi
  end if
  !
  inquire(file='fermi_sub', exist=lfile_exist1)
  if (lfile_exist1) then
      title = 'fermi_sub'
      open(unit=64, file=title)
      rewind(64)
      read(64,'(A)') buffer
      read(64,'(A)') buffer
      read(buffer(1:23),'(A10,F12.7)') ch_tmp, sub_fermi
      close(64)
      write(6,*) 'Fermi level of slab: ', sub_fermi
  end if
  !
  inquire(file='to_be_read_env_bands_sp1', exist=lfile_exist1)
  inquire(file='to_be_read_sys_bands_sp1', exist=lfile_exist2)
  if (lfile_exist1 .and. lfile_exist2) then
    open(unit=10,file='to_be_read_sys_bands_sp1')
    rewind(10)
    sys_nband = -1
    istat = 0
    write(6,*) 'read to_be_read_sys_bands_sp1 file'
    do while(.not. is_iostat_end(istat))
      read(10,'(a)', iostat=istat) buffer
      sys_nband = sys_nband + 1
    end do
    close(10)
    !
    write(6,*) 'number of system bands', sys_nband
    !
    open(unit=11,file='to_be_read_env_bands_sp1')
    rewind(11)
    env_nband = -1
    istat = 0
    write(6,*) 'read to_be_read_env_bands_sp1 file'
    do while(.not. is_iostat_end(istat))
      read(11,'(a)', iostat=istat) buffer
      env_nband = env_nband + 1
    end do
    !
    write(6,*) 'number of environment bands', env_nband
    !
    allocate(chosen_env_bands(env_nband))
    allocate(env_eig(env_nband,2), env_occ(env_nband,2))
    rewind(11)
    do i=1, env_nband
      read(11,100) chosen_env_bands(i),&
                   env_occ(i,1),env_eig(i,1),&
                   env_occ(i,2),env_eig(i,2)
      !write(6,*) 'read number of environment bands',&
      !            chosen_env_bands(i)
    end do
    close(11)
    !
    ww_i = minval(env_eig(:,1)) - 15.d0
    ww_f = maxval(env_eig(:,1)) + 15.d0
    open(unit=5,file='input')
    rewind(5)
    read(5, hybridization, end=113)
    113 continue
    close(5)
    ww_step = (ww_f - ww_i) / num_ww
    do m=1, num_ww
       ww(m) = (m-1) * ww_step + ww_i
    end do
    !
    allocate(hopping(sys_nband,env_nband))
    allocate(sum_hopping(env_nband))
    hopping = 0.d0
    sum_hopping = 0.d0
    open(unit=12, file='spin-1-coupling-details')
    rewind(12)
    do m=1, 6
       read(12,'(a)') buffer
    end do
    !
    do m=1, env_nband 
      do n=1, sys_nband
        read(12,101) ntmp1, ntmp2, &
                     dtmp1, dtmp2, &
                     hopping(n,m)
        sum_hopping(m) = sum_hopping(m) + &
                         hopping(n,m)**2
      end do  
    end do
    close(12)
    !
    open(unit=13,file='hyb-sp1')
    rewind(13)
    do m=1, num_ww
      do n=1, env_nband
        hyb(m) = hyb(m) + sum_hopping(n) * &
        width / ((ww(m)-env_eig(n,1))**2+(width*width))
      end do
      write(13,102) ww(m) - sub_fermi, hyb(m)
    end do
    close(13)
    write(6,103) (1.0d0/heom_width)*sum(sum_hopping)!/&
                !(atan(ww_f/heom_width)-atan(ww_i/heom_width))
    !
    ! Pade decomposition
    allocate(band_center(env_nband), coupling(env_nband))
    allocate(eta_tmp(npade+1,env_nband,2), gama_tmp(npade+1,env_nband,2))
    allocate(ct_tot(num_tt,2), ct(num_tt,env_nband,2))
    band_center = env_eig(:,1) - sub_fermi
    coupling = sum_hopping
    eta_tmp = czero
    gama_tmp = czero
    ct_tot = czero
    ct = czero
    !
    do i=1, 2
      do m=1, env_nband
        if(cdabs(coupling(m)) .lt. crit) then
          cycle
        end if
        !if(cdabs(coupling(m) - maxval(dreal(coupling))) .le. crit) then
        !  ntmp1 = m
        !end if
        ! generate spectral results
        gama_tmp(1,m,i) = (width - eye * sigma(i) * band_center(m)) / hbar
        ctmp1 = eye * gama_tmp(1,m,i) / T * hbar
        ctmp2 = dcmplx(5.d-1, 0.d0)
        do n=1, npade
          ctmp2 = ctmp2 + zpcoef(n) / (ctmp1 + zppole(n,i)) + &
                     zpcoef(n) / (ctmp1 - zppole(n,i))
        end do
        eta_tmp(1,m,i) = 0.5d0 * coupling(m) * width * ctmp2 / hbar**2
        do n=2, npade + 1
          gama_tmp(n,m,i) = ( -1.d0 * eye * sigma(i) * T * zppole(n-1,i) ) / hbar
          eta_tmp(n,m,i) = 2.d0 * eye * zpcoef(n-1) * coupling(m) * 0.5d0 * width**2 * T &
                        / ( (zppole(n-1,i) * T - band_center(m))**2 + width**2 ) / hbar**2
        end do
        ! calculate C(t)
        do j=1, num_tt
          tt = (j-1) * tt_step
          do n=1, npade+1
            ct(j,m,i) = ct(j,m,i) + &
             eta_tmp(n,m,i) * cdexp(-1.d0 * gama_tmp(n,m,i) * tt)
          end do
          ct_tot(j,i) = ct_tot(j,i) + ct(j,m,i)
        end do
        write(6,'("calculating for sigma = ",I2, " and No. ", I5," band with hopping ",E12.5)')&
        i, chosen_env_bands(m), dreal(coupling(m))
      end do
    end do
    !
    open(unit=21, file='ct-sp1.data')
    rewind(21)
    !write(21,'(4E18.10)') dble(coupling(ntmp1)), dble(band_center(ntmp1)), width, T
    !do j=1, num_tt
    !  tt = (j-1) * tt_step
    !  write(21,'(5E18.10)') tt, & 
    !                        dreal(ct(j,ntmp1,1)), dimag(ct(j,ntmp1,1)),&
    !                        dreal(ct(j,ntmp1,2)), dimag(ct(j,ntmp1,2))
    !end do
    do j=1, num_tt
      tt = (j-1) * tt_step
      write(21,'(5E18.10)') tt, & 
                            dreal(ct_tot(j,1)), dimag(ct_tot(j,1)),&
                            dreal(ct_tot(j,2)), dimag(ct_tot(j,2))
    end do
    close(21)
    !
    deallocate(sum_hopping,hopping,chosen_env_bands)
    deallocate(env_occ,env_eig)
    deallocate(band_center,coupling,eta_tmp,gama_tmp,ct_tot,ct)
  else
    write(6,*) 'Do not find to_be_read_sys_bands_sp1 '&
               'or to_be_read_env_bands_sp1'
  end if
  write(6,*)
  write(6,*)
  !
  ww = 0.d0
  hyb = 0.d0
  !
  zpeigv = czero
  zpcoef = czero
  zppole = czero
  call zpout_psd(npade, 2, zpeigv, zpcoef, zppole)
  !
  inquire(file='to_be_read_env_bands_sp2', exist=lfile_exist1)
  inquire(file='to_be_read_sys_bands_sp2', exist=lfile_exist2)
  if (lfile_exist1 .and. lfile_exist2) then
    open(unit=14,file='to_be_read_sys_bands_sp2')
    rewind(14)
    sys_nband = -1
    istat = 0
    write(6,*) 'read to_be_read_sys_bands_sp2 file'
    do while(.not. is_iostat_end(istat))
      read(14,'(a)', iostat=istat) buffer
      sys_nband = sys_nband + 1
    end do
    close(14)
    !
    write(6,*) 'number of system bands', sys_nband
    !
    open(unit=15,file='to_be_read_env_bands_sp2')
    rewind(15)
    env_nband = -1
    istat = 0
    write(6,*) 'read to_be_read_env_bands_sp2 file'
    do while(.not. is_iostat_end(istat))
      read(15,'(a)', iostat=istat) buffer
      env_nband = env_nband + 1
    end do
    !
    write(6,*) 'number of environment bands', env_nband
    !
    allocate(chosen_env_bands(env_nband))
    allocate(env_eig(env_nband,2), env_occ(env_nband,2))
    rewind(15)
    do i=1, env_nband
      read(15,100) chosen_env_bands(i),&
                   env_occ(i,1),env_eig(i,1),&
                   env_occ(i,2),env_eig(i,2)
      !write(6,*) 'read number of environment bands',&
      !            chosen_env_bands(i)
    end do
    close(15)
    !
    ww_i = minval(env_eig(:,2)) - 15.d0
    ww_f = maxval(env_eig(:,2)) + 15.d0
    open(unit=5,file='input')
    rewind(5)
    read(5, hybridization, end=123)
    123 continue
    close(5)
    ww_step = (ww_f - ww_i) / num_ww
    do m=1, num_ww
       ww(m) = (m-1) * ww_step + ww_i
    end do
    !
    allocate(hopping(sys_nband,env_nband))
    allocate(sum_hopping(env_nband))
    hopping = 0.d0
    sum_hopping = 0.d0
    open(unit=16, file='spin-2-coupling-details')
    rewind(16)
    do m=1, 6
       read(16,'(a)') buffer
    end do
    do m=1, env_nband 
      do n=1, sys_nband
        read(16,101) ntmp1, ntmp2, &
                     dtmp1, dtmp2, &
                     hopping(n,m)
        sum_hopping(m) = sum_hopping(m) + &
                         hopping(n,m)**2
      end do  
    end do
    close(16)
    !
    open(unit=17,file='hyb-sp2')
    rewind(17)
    do m=1, num_ww
      do n=1, env_nband
        hyb(m) = hyb(m) + sum_hopping(n) * &
        width / ((ww(m)-env_eig(n,2))**2+(width*width))
      end do
      write(17,102) ww(m) - sub_fermi, hyb(m)
    end do
    close(17)
    write(6,104) (1.0d0/heom_width)*sum(sum_hopping)!/&
                !(atan(ww_f/heom_width)-atan(ww_i/heom_width))
    !
    ! Pade decomposition
    allocate(band_center(env_nband), coupling(env_nband))
    allocate(eta_tmp(npade+1,env_nband,2), gama_tmp(npade+1,env_nband,2))
    allocate(ct_tot(num_tt,2), ct(num_tt,env_nband,2))
    band_center = env_eig(:,1) - sub_fermi
    coupling = sum_hopping
    eta_tmp = czero
    gama_tmp = czero
    ct_tot = czero
    ct = czero
    !
    do i=1, 2
      do m=1, env_nband
        if(cdabs(coupling(m)) .lt. crit) then
          cycle
        end if
        !if(cdabs(coupling(m) - maxval(dreal(coupling))) .le. crit) then
        !  ntmp1 = m
        !end if
        ! generate spectral results
        gama_tmp(1,m,i) = (width - eye * sigma(i) * band_center(m)) / hbar
        ctmp1 = eye * gama_tmp(1,m,i) / T * hbar
        ctmp2 = dcmplx(5.d-1, 0.d0)
        do n=1, npade
          ctmp2 = ctmp2 + zpcoef(n) / (ctmp1 + zppole(n,i)) + &
                     zpcoef(n) / (ctmp1 - zppole(n,i))
        end do
        eta_tmp(1,m,i) = 0.5d0 * coupling(m) * width * ctmp2 / hbar**2
        do n=2, npade + 1
          gama_tmp(n,m,i) = ( -1.d0 * eye * sigma(i) * T * zppole(n-1,i) ) / hbar
          eta_tmp(n,m,i) = 2.d0 * eye * zpcoef(n-1) * coupling(m) * 0.5d0 * width**2 * T &
                        / ( (zppole(n-1,i) * T - band_center(m))**2 + width**2 ) / hbar**2
        end do
        ! calculate C(t)
        do j=1, num_tt
          tt = (j-1) * tt_step
          do n=1, npade+1
            ct(j,m,i) = ct(j,m,i) + &
             eta_tmp(n,m,i) * cdexp(-1.d0 * gama_tmp(n,m,i) * tt)
          end do
          ct_tot(j,i) = ct_tot(j,i) + ct(j,m,i)
        end do
        write(6,'("calculating for sigma = ",I2, " and No. ", I5," band with hopping ",E12.5)')&
        i, chosen_env_bands(m), dreal(coupling(m))
      end do
    end do
    open(unit=22, file='ct-sp2.data')
    rewind(22)
!    write(22,'(4E18.10)') dble(coupling(ntmp1)), dble(band_center(ntmp1)), width, T
!    do j=1, num_tt
!      tt = (j-1) * tt_step
!      write(22,'(5E18.10)') tt, & 
!                            dreal(ct(j,ntmp1,1)), dimag(ct(j,ntmp1,1)),&
!                            dreal(ct(j,ntmp1,2)), dimag(ct(j,ntmp1,2)) 
!    end do
    do j=1, num_tt
      tt = (j-1) * tt_step
      write(22,'(5E18.10)') tt, & 
                            dreal(ct_tot(j,1)), dimag(ct_tot(j,1)),&
                            dreal(ct_tot(j,2)), dimag(ct_tot(j,2))
    end do
    close(22)
    !
    deallocate(sum_hopping,hopping,chosen_env_bands)
    deallocate(env_occ,env_eig)
  else
    write(6,*) 'Do not find to_be_read_sys_bands_sp2 '&
               'or to_be_read_env_bands_sp2'
  end if
  !
  100 format('No. of band',1x,I5,2x,&
             ',occ. spin-1',2x,F9.6,2x,&
             ',eigenergy',2x,F11.6,2x,&
             ',occ. spin-2',2x,F9.6,2x,&
             ',eigenergy',2x,F11.6)
  101 format(I10, 2X, I10, 3E18.10)
  102 format(2(E18.10,2x))
  103 format("For Lorentzian type spectral function,",&
             " effective spin-1 coupling is", 2x, E18.10)
  104 format("For Lorentzian type spectral function,",&
             " effective spin-2 coupling is", 2x, E18.10)
  !
end program main

subroutine zpout_psd(mpole, msgn, zeig, zcoef, zpole)
!
! purpose : find poles and coefficients for PSD scheme
! 
! f(z) = 0.5 + sum_i zcoef * 2 * z / (z^2 - zpole^2)
! with z = beta * w
!
implicit none
!
integer, intent (in)     :: mpole, msgn
complex*16, intent (out)  :: zeig(*), zcoef(*), zpole(mpole,*)
integer              :: istat, info, ni, nj, nk
integer              :: dim0, dim1, n0, n1
integer              :: lwork
real*8               :: dtmp0
real*8, allocatable     :: delta(:,:), dmat0(:,:), dmat1(:,:)
real*8, allocatable     :: dval0(:), dval1(:), dwork(:)
real*8, allocatable     :: eval0(:), eval1(:), mattmp1(:,:)
real*8, allocatable     :: tmpa(:), dcoefeta(:), mattmp0(:,:)
complex*16            :: cunity
!
n0   = mpole
n1   = mpole - 1
dim0 = 2 * mpole
dim1 = dim0 - 1
cunity = (1.d0, 0.d0)
!
lwork = 5 * dim0
allocate(delta(dim0,dim0), dmat0(dim0,dim0), dmat1(dim1,dim1), mattmp0(dim0,dim0), mattmp1(dim1,dim1), STAT=istat)
allocate(dval0(dim0), dval1(dim1), dwork(lwork), STAT=istat)
!
! Kronecker delta function
!
delta = 0.d0
do ni=1,dim0
   delta(ni,ni) = 1.d0
end do
!
! construct real symmetric matrices
!
dmat0 = 0.d0
do nj=1,dim0
   do ni=1,dim0
      if (ni .eq. nj + 1 .or. ni .eq. nj - 1) then
         dtmp0 = (2.d0 * dble(ni) - 1.d0) * (2.d0 * dble(nj) - 1.d0)
         dmat0(ni,nj) = dmat0(ni,nj) + 1.d0 / dsqrt(dtmp0)
      end if
   end do
end do
!
dmat1 = 0.d0
do nj=1,dim1
   do ni=1,dim1
      if (ni .eq. nj + 1 .or. ni .eq. nj - 1) then
         dtmp0 = (2.d0 * dble(ni) + 1.d0) * (2.d0 * dble(nj) + 1.d0)
         dmat1(ni,nj) = dmat1(ni,nj) + 1.d0 / dsqrt(dtmp0)
      end if
   end do
end do
! eigenvalues 
!
call dsyev('N', 'U', dim0, dmat0, dim0, dval0, dwork, lwork, info)
if (info .ne. 0) then
   write(6,*)'zpout_psd: error! diagonalization failed for dmat0', info
   stop
end if
call dsyev('N', 'U', dim1, dmat1, dim1, dval1, dwork, lwork, info)
if (info .ne. 0) then
   write(6,*)'zpout_psd: error! diagonalization failed for dmat1', info
   stop
end if
!
allocate(eval0(n0), eval1(n1), STAT=istat)
allocate(tmpa(n0), dcoefeta(n0), STAT=istat)
!
eval0(1:n0) = dval0(1:n0)
eval1(1:n1) = dval1(1:n1)
do ni=1,n0
   eval0(ni) = 4.d0 / eval0(ni)**2  ! eigXi
end do
do ni=1,n1
   eval1(ni) = 4.d0 / eval1(ni)**2  ! eigZeta
end do
!
zeig(1:n0) = dcmplx(eval0(1:n0), 0.d0)
!
tmpa = 1.d0
do nj=1,n0
   do nk=1,nj-1
      tmpa(nj) = tmpa(nj) * (eval1(nk) - eval0(nj)) / (eval0(nk) - eval0(nj)) 
   end do
   do nk=nj+1,n0
      tmpa(nj) = tmpa(nj) * (eval1(nk - 1) - eval0(nj)) / (eval0(nk) - eval0(nj))
   end do
end do

!
!write(6,*)'zpout_psd: show dcoefeta'
do ni=1,n0
   dcoefeta(ni) = 5.d-1 * dble(n0 * (2 * n0 + 1)) * tmpa(ni)
   !write(6,100) ni, dcoefeta(ni)
   if (isnan(dcoefeta(ni))) then
      write(6,*)'zpout_psd: error! NAN found for dcoefeta(', ni, ')'
      stop
   end if
end do
100 format('zpout_psd: dcoefeta(', I4, ') = ', e14.6e3)
!
zcoef(1:n0) = -cunity * dcoefeta(1:n0)
!
do ni=1,n0
   dtmp0 = dsqrt(eval0(ni))
   !zpole(ni) = dcmplx(0.d0, dabs(dtmp0))
   zpole(ni,1) = dcmplx(0.d0, dabs(dtmp0))
   zpole(ni,2) = dconjg(zpole(ni,1))
end do

!
deallocate(delta, dmat0, dmat1, STAT=istat)
deallocate(dval0, dval1, dwork, STAT=istat)
deallocate(eval0, eval1, STAT=istat)
deallocate(tmpa, dcoefeta, STAT=istat)
!
end subroutine


