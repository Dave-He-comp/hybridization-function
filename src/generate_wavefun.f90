!$**************************** WaveTrans *******************************
!!$*
!!$*   input the WAVECAR file in binary format from VASP, and output 
!!$*   in text format a file GCOEFF of G values and corresponding plane
!!$*   wave coefficients, together with energy eigenvalues and
!occupations
!!$*
!!$*   Compile with gfortran or ifort. Flag "-assume byterecl" is
!required
!!$*   for ifort.
!!$*
!!$*   version 2.0 - July 5, 2012 - R. M. Feenstra and M. Widom
!!$*   version 2.1 - Sept 30, 2012 - changed estimator for max. no. of
!!$*                                 plane waves
!!$*   version 1.3 - Sept 17, 2014 - updated 'c' value
!!$*
!!$*   format of GCOEFF file:
!!$*     no. wavevector values (integer)
!!$*     no. bands (integer)
!!$*     real space lattice vector (a1) x,y,z coefficients (real)
!!$*     real space lattice vector (a2) x,y,z coefficients (real)
!!$*     real space lattice vector (a3) x,y,z coefficients (real)
!!$*     recip. lattice vector 1 (b1) x,y,z coefficients (real)
!!$*     recip. lattice vector 2 (b2) x,y,z coefficients (real)
!!$*     recip. lattice vector 3 (b3) x,y,z coefficients (real)
!!$*     loop over spins
!!$*        loop through wavevector values:
!!$*           wavevector x,y,z components (real)
!!$*           energy (complex), occupation (real)
!!$*           loop through bands:
!!$*              index of band (integer), no. plane waves (integer)
!!$*              loop through plane waves:
!!$*                 ig1,ig2,ig3 values (integer) and coefficient
!(complex)
!!$*              end loop through plane waves
!!$*           end loop through bands
!!$*         end loop through wavevector values
!!$*     end loop over spins
!!$*
!!$*   the x,y,z components of each G value are given in terms of the
!!$*   ig1,ig2,ig3 values and the components of the recip. lattice
!vectors
!!$*   according to:
!!$*   ig1*b1_x + ig2*b2_x + ig3*b3_x,
!!$*   ig1*b1_y + ig2*b2_y + ig3*b3_y, and
!!$*   ig1*b1_z + ig2*b2_z + ig3*b3_z, respectively
!!$*
!!$*   note that the energy eigenvalues are complex, as provided in the
!!$*   WAVECAR file, but the imaginary part is zero (at least for all
!cases
!!$*   investigated thus far)
!!$*     

subroutine read_wavecar(file_name,energy_windows_i,energy_windows_f,lread_para)
                                                                        
implicit real*8             (a-h, o-z)
complex*8,  allocatable     :: coeff(:,:)
complex*16, allocatable     :: cener(:,:)
real*8,     allocatable     :: occ(:,:)
real*8                      :: energy_windows_i, energy_windows_f
real*8                      :: fermi_level
integer,    allocatable     :: igall(:,:)
integer,    allocatable     :: chosen_sys_bands(:), chosen_env_bands(:)
integer,    allocatable     :: read_band_tmp(:,:)
integer,    allocatable     :: cnt_sys_band_read(:)
integer,    allocatable     :: cnt_env_band_read(:)
integer                     :: num_file, lread_para
integer                     :: num_chosen_sys_bands
integer                     :: num_chosen_env_bands
integer                     :: gall_tmp(3)
logical                     :: file_exists
logical                     :: l_chosen_sys_bands, l_sys_bands_chosen
logical                     :: l_chosen_env_bands, l_env_bands_chosen
logical                     :: l_sys_file_exist, l_env_file_exist
!
character(len=2)  :: ch_spin
character(len=5)  :: ch_wk
character(len=5)  :: ch_band
character(len=25) :: title
character(len=*) :: file_name
character(len=128) :: buffer, ch_tmp
!
dimension a1(3), a2(3), a3(3)
dimension b1(3), b2(3), b3(3)
dimension vtmp(3), sumkg(3), wk(3)
     
! constant 'c' below is 2m/hbar**2 in units of 1/eV Ang^2 
! (value is adjusted in final decimal places to agree with VASP value;
! program checks for discrepancy of any results between this and VASP
! values)

data c/0.262465831d0/ 
! data c/0.26246582250210965422d0/ 
pi = 4.0d0 * atan(1.0d0)
      
! input
     
nrecl = 24
open(unit = 10, file = file_name, access = 'direct', &
!open(unit = 10, file = 'WAVECAR', access = 'direct', &
     recl = nrecl, iostat = iost, status = 'old')
!
if (iost .ne. 0) then
   write(6,*) 'open error - iostat =', iost 
   write(6,*) 'file name is ', file_name
   stop
end if
!
rewind(10)
read(unit = 10, rec = 1) xnrecl, xnspin, xnprec
close(unit = 10)
!
nrecl = nint(xnrecl)
nspin = nint(xnspin)
nprec = nint(xnprec)
!
if(nprec .eq. 45210) then
   write(6,*) '*** error - WAVECAR_double requires complex*16'
   stop
endif
!
write(6,*) 
write(6,*) 'record length  =', nrecl, &
           ' spins =', nspin, ' prec flag ', nprec
open(unit = 10, file = file_name, access = 'direct', &
!open(unit = 10, file = 'WAVECAR', access = 'direct',&
     recl = nrecl, iostat = iost, status = 'old')
!
if (iost .ne. 0) then
   write(6,*) 'open error - iostat =', iost
   stop
end if
!           
!open(unit = 11, file = 'GCOEFF.txt')
read(unit = 10, rec = 2) xnwk, xnband, ecut, &
    (a1(j), j=1, 3), (a2(j), j=1, 3), (a3(j), j=1, 3)
nwk = nint(xnwk)
nband = nint(xnband)
!
allocate(occ(nband,nspin))
allocate(cener(nband,nspin))
!
write(6,*) 'no. k points =', nwk
write(6,*) 'no. bands =', nband
write(6,*) 'max. energy =', ecut
write(6,*) 'real space lattice vectors:'
write(6,*) 'a1 =', ((a1(j)), j=1, 3)
write(6,*) 'a2 =', ((a2(j)), j=1, 3)
write(6,*) 'a3 =', ((a3(j)), j=1, 3)
write(6,*) ' '
!write(11,*) nspin
!write(11,*) nwk
!write(11,*) nband

! compute reciprocal properties

write(6,*) ' '
call vcross(vtmp,a2,a3)
Vcell = a1(1) * vtmp(1) + a1(2) * vtmp(2) + a1(3) * vtmp(3)
write(6,*) 'volume unit cell =', (Vcell)
!
call vcross(b1,a2,a3)
call vcross(b2,a3,a1)
call vcross(b3,a1,a2)
!
do j=1, 3
   b1(j) = 2.0d0 * pi * b1(j) / Vcell
   b2(j) = 2.0d0 * pi * b2(j) / Vcell
   b3(j) = 2.0d0 * pi * b3(j) / Vcell
end do
!
b1mag = dsqrt(b1(1)**2 + b1(2)**2 + b1(3)**2)
b2mag = dsqrt(b2(1)**2 + b2(2)**2 + b2(3)**2)
b3mag = dsqrt(b3(1)**2 + b3(2)**2 + b3(3)**2)
!
write(6,*) 'reciprocal lattice vectors:'
write(6,*) 'b1 =', ((b1(j)), j=1,3)
write(6,*) 'b2 =', ((b2(j)), j=1,3)
write(6,*) 'b3 =', ((b3(j)), j=1,3)
write(6,*) 'reciprocal lattice vector magnitudes:'
write(6,*) (b1mag), (b2mag), (b3mag)
write(6,*) ' '
!

phi12 = dacos((b1(1) * b2(1) + &
               b1(2) * b2(2) + &
               b1(3) * b2(3)) / (b1mag * b2mag))
call vcross(vtmp, b1, b2)
vmag = dsqrt(vtmp(1)**2 + vtmp(2)**2 + vtmp(3)**2)
sinphi123 = (b3(1) * vtmp(1) + &
             b3(2) * vtmp(2) + &
             b3(3) * vtmp(3)) / (vmag * b3mag)
nb1maxA = (dsqrt(ecut * c) / (b1mag * dabs(dsin(phi12)))) + 1
nb2maxA = (dsqrt(ecut * c) / (b2mag * dabs(dsin(phi12)))) + 1
nb3maxA = (dsqrt(ecut * c) / (b3mag * dabs(sinphi123))) + 1
npmaxA = nint(4.0d0 * pi * nb1maxA * nb2maxA * nb3maxA / 3.0d0)
      
phi13 = dacos((b1(1) * b3(1) + &
               b1(2) * b3(2) + &
               b1(3) * b3(3)) / (b1mag * b3mag))
call vcross(vtmp, b1, b3)
vmag = dsqrt(vtmp(1)**2 + vtmp(2)**2 + vtmp(3)**2)
sinphi123 =( b2(1) * vtmp(1) + &
             b2(2) * vtmp(2) + &
             b2(3) * vtmp(3)) / (vmag * b2mag)
phi123 = dabs(dasin(sinphi123))
nb1maxB = (dsqrt(ecut * c) / (b1mag * dabs(dsin(phi13)))) + 1
nb2maxB = (dsqrt(ecut * c) / (b2mag * dabs(sinphi123))) + 1
nb3maxB = (dsqrt(ecut * c) / (b3mag * dabs(dsin(phi13)))) + 1
npmaxB = nint(4.0d0 * pi * nb1maxB * nb2maxB * nb3maxB / 3.0d0)
      
phi23 = dacos((b2(1) * b3(1) + &
               b2(2) * b3(2) + &
               b2(3) * b3(3)) / (b2mag * b3mag))
call vcross(vtmp, b2, b3)
vmag = dsqrt(vtmp(1)**2 + vtmp(2)**2 + vtmp(3)**2)
sinphi123 = (b1(1) * vtmp(1) + &
             b1(2) * vtmp(2) + &
             b1(3) * vtmp(3)) / (vmag * b1mag)
phi123 = dabs(dasin(sinphi123))
nb1maxC = (dsqrt(ecut * c) / (b1mag * dabs(sinphi123))) + 1
nb2maxC = (dsqrt(ecut * c) / (b2mag * dabs(dsin(phi23)))) + 1
nb3maxC = (dsqrt(ecut * c) / (b3mag * dabs(dsin(phi23)))) + 1 
npmaxC = nint(4.0d0 * pi * nb1maxC * nb2maxC * nb3maxC / 3.0d0)

nb1max = max0(nb1maxA, nb1maxB, nb1maxC)
nb2max = max0(nb2maxA, nb2maxB, nb2maxC)
nb3max = max0(nb3maxA, nb3maxB, nb3maxC)
npmax = min0(npmaxA, npmaxB, npmaxC)

inquire(file='G_MAX', exist=file_exists)
if (file_exists) then
   write(6,*) 'max. no. G values; 1,2,3 =', nb1max, nb2max, nb3max
   write(6,*) ' '
   gall_tmp = 0
   open(unit=30,file='G_MAX')
   read(30,'(3I8)') (gall_tmp(i), i=1,3)
   if ((gall_tmp(1) .ne. nb1max) .and. &
       (gall_tmp(2) .ne. nb2max) .and. &
       (gall_tmp(3) .ne. nb3max)) then
       rewind(30)
       write(30,'(3I8)') nb1max, nb2max, nb3max
   end if
   close(30)
else
   write(6,*) 'max. no. G values; 1,2,3 =', nb1max, nb2max, nb3max
   write(6,*) ' '
   open(unit=30,file='G_MAX')
   write(30,'(3I8)') nb1max, nb2max, nb3max
   close(30)
end if
!
allocate (igall(3, npmax))
igall = 0
write(6,*) 'estimated max. no. plane waves =', npmax
allocate (coeff(npmax, nband))
coeff = 0
!write(11,*) npmax

! Begin loops over spin, k-points and bands
! output the files that specifies spin, k-vector and band
irec = 2
num_file = 100
ch_spin = ''
ch_wk = ''
ch_band = ''
title = ''
nplane = 0
max_nplane = 0
allocate(read_band_tmp(nband,nspin))
allocate(cnt_sys_band_read(nspin))
allocate(cnt_env_band_read(nspin))
l_chosen_sys_bands = .False.
l_chosen_env_bands = .False.
num_chosen_sys_bands = -1
num_chosen_env_bands = -1
cnt_sys_band_read = 0
cnt_env_band_read = 0
read_band_tmp = 0
!
! search the chosen system bands
inquire(file='chosen_sys_bands', exist=l_sys_file_exist)
if (l_sys_file_exist .and. file_name(1:11) .eq. "WAVECAR_SYS") then
   l_chosen_sys_bands = .True.
   num_chosen_sys_bands = -1
   write(6,*) 'find the chosen band file'
   open(unit=29, file='chosen_sys_bands')
   rewind(29)
   do while(.not. is_iostat_end(istat))
     read(29,'(a)', iostat=istat) ch_tmp
     num_chosen_sys_bands = num_chosen_sys_bands + 1
   end do
   write(6,*) 'total number of chosen system bands',&
               num_chosen_sys_bands
   allocate(chosen_sys_bands(num_chosen_sys_bands))
   rewind(29)
   do i=1, num_chosen_sys_bands
     read(29,*) chosen_sys_bands(i)
     write(6,*) 'read number of system band', chosen_sys_bands(i)
   end do
   close(29)
   write(6,*)
end if
!
inquire(file='chosen_env_bands', exist=l_env_file_exist)
if (l_env_file_exist .and. &
    (file_name(1:11) .eq. "WAVECAR_SUB" .or. &
     file_name(1:11) .eq. "WAVECAR_TIP" .or. &
     file_name(1:11) .eq. "WAVECAR_ENV")) then
   l_chosen_env_bands = .True.
   num_chosen_env_bands = -1
   write(6,*) 'find the chosen band file'
   open(unit=30, file='chosen_env_bands')
   rewind(30)
   do while(.not. is_iostat_end(istat))
     read(30,'(a)', iostat=istat) ch_tmp
     num_chosen_env_bands = num_chosen_env_bands + 1
   end do
   write(6,*) 'total number of chosen environment bands',&
               num_chosen_env_bands
   allocate(chosen_sys_bands(num_chosen_env_bands))
   rewind(30)
   do i=1, num_chosen_env_bands
     read(30,*) chosen_env_bands(i)
     write(6,*) 'read number of environment band',&
                 chosen_env_bands(i)
   end do
   close(30)
   write(6,*)
end if
!
write(6,*)
if (l_chosen_sys_bands) then
  write(6,*) 'number of the output files supposed to be genrated:', &
              num_chosen_sys_bands
else
  write(6,*) 'number of the output files supposed to be genrated:', &
             (nband * nspin * nwk)
end if
!
if (l_chosen_env_bands) then
  write(6,*) 'number of the output files supposed to be genrated:', &
              num_chosen_env_bands
else
  write(6,*) 'number of the output files supposed to be genrated:', &
             (nband * nspin * nwk)
end if
!
fermi_level = 0.d0
if (file_name(1:11) .eq. "WAVECAR_SYS") then
  !call system('grep Fermi OUTCAR_SYS > outmp')
  open(unit=63, file='fermi_sys')
  rewind(63)
  read(63,'(A)') buffer
  read(63,'(A)') buffer
  !read(63,'(A14, F16.10)') buffer, fermi_level
  read(buffer(1:23),'(A10,F12.7)') ch_tmp, fermi_level
  close(63)
  write(6,*) 'Fermi level of system: ', fermi_level
  write(6,*)
end if
!
if (file_name(1:11) .eq. "WAVECAR_SUB" .or. &
    file_name(1:11) .eq. "WAVECAR_TIP" .or. &
    file_name(1:11) .eq. "WAVECAR_ENV") then  
  !
  !call system('grep Fermi OUTCAR_ENV > outmp')
  inquire(file='fermi_sub', exist=file_exists)
  if (file_name(1:11) .eq. "WAVECAR_SUB" .and. &
      file_exists) then
     title = 'fermi_sub'
  end if
  inquire(file='fermi_tip', exist=file_exists)
  if (file_name(1:11) .eq. "WAVECAR_TIP" .and. &
      file_exists) then
     title = 'fermi_tip'
  end if
  open(unit=64, file=title)
  rewind(64)
  read(64,'(A)') buffer
  read(64,'(A)') buffer
  read(buffer(1:23),'(A10,F12.7)') ch_tmp, fermi_level
  !read(64,'(A14, F16.10)') buffer, fermi_level
  close(64)
  write(6,*) 'Fermi level of environment: ', fermi_level
  write(6,*)
end if
title=''
!
do isp = 1, nspin
   write(ch_spin,'(I1)') isp
   write(6,*) ' '
   write(6,*) '******'
   write(6,*) 'reading spin ', isp
   do iwk=1, nwk
      write(ch_wk,'(I4)') iwk
      irec = irec + 1
      read(unit=10, rec=irec) xnplane, (wk(i), i=1,3), &
      (cener(iband,isp), occ(iband,isp), iband=1, nband)
      nplane = nint(xnplane)
      if (nplane .gt. max_nplane) then
        max_nplane = nplane 
      end if 
      write(6,*) 'k point #', iwk, '  input no. of plane waves =',&
                 nplane
      write(6,*) 'k value =', ((wk(j)), j=1,3)
      !write(11,*) ((wk(j)), j=1,3) 
      ! Calculate plane waves
      ncnt = 0
      do ig3 = 0, 2 * nb3max
         ig3p = ig3
         if (ig3 .gt. nb3max) then
            ig3p = ig3 - 2 * nb3max - 1
         end if
         do ig2 = 0, 2 * nb2max
            ig2p = ig2
            if (ig2 .gt. nb2max) then
               ig2p = ig2 - 2 * nb2max - 1
            end if
            do ig1 = 0, 2*nb1max
               ig1p = ig1
               if (ig1 .gt. nb1max) then
                  ig1p = ig1 - 2 * nb1max - 1
               end if
               do j = 1, 3
                  sumkg(j) = (wk(1) + ig1p) * b1(j) + &
                             (wk(2) + ig2p) * b2(j) + &
                             (wk(3) + ig3p) * b3(j)
               end do
               gtot = dsqrt(sumkg(1)**2 + sumkg(2)**2 + sumkg(3)**2)
               etot = gtot**2 / c
               if (etot .lt. ecut) then
                  ncnt = ncnt + 1
                  igall(1, ncnt) = ig1p
                  igall(2, ncnt) = ig2p
                  igall(3, ncnt) = ig3p
               end if
            end do
         end do
      end do
      if (ncnt .ne. nplane) then
         write(6,*) '*** error - computed no. != input no.'
         stop
      end if
      if (ncnt .gt. npmax) then
         write(6,*) '*** error - plane wave count exceeds estimate'
         stop
      end if
      do iband = 1, nband
         irec = irec + 1
         read(unit=10, rec=irec) (coeff(iplane, iband), &
                                  iplane=1,nplane)
      end do
      ! output G values and coefficients
      do iband = 1, nband
         write(ch_band,'(I4)') iband
         !write(11,*) iband, nplane
         !write(11,560) cener(iband,isp), occ(iband,isp)
         ! output files that specify spin, wk and band
         !
         title = ('spin-'//trim(adjustl(ch_spin))//&
                  '-wk-'//trim(adjustl(ch_wk))//&
                  '-band-'//trim(adjustl(ch_band)))
         if (lread_para .eq. 1) then
            num_file = 99
            open(unit = num_file, name = 'wavecar_output',&
                 iostat = iost)
            rewind(num_file)
         else
            ! check system band
            if (file_name(1:11) .eq. "WAVECAR_SYS") then
               if (l_sys_file_exist) then
                  l_sys_bands_chosen = .True.
                  do i=1, num_chosen_sys_bands
                    if (iband .eq. chosen_sys_bands(i)) then
                      l_sys_bands_chosen = .False.
                      exit
                    end if
                  end do
                  !
                  if (l_sys_bands_chosen) then
                    cycle
                  end if
                  cnt_sys_band_read(isp) = cnt_sys_band_read(isp) + 1
                  read_band_tmp(iband,isp) = iband
               else 
                  if (&!(occ(iband,isp) .eq. 0.d0) .or. & ! null state
                      !(occ(iband,isp) .eq. 1.d0) .or. & ! full state
                      !(dreal(cener(iband,isp)) .gt. fermi_level) .or. & !
                      !higher than fermi level
                      (dreal(cener(iband,isp)-fermi_level) .lt.&
                      energy_windows_i) .or. &
                      (dreal(cener(iband,isp)-fermi_level) .gt.&
                      energy_windows_f)) then  ! near fermi level
                      cycle
                  end if
                  cnt_sys_band_read(isp) = cnt_sys_band_read(isp) + 1
                  read_band_tmp(iband,isp) = iband
               end if
            ! check environment band
            else if  (file_name(1:11) .eq. "WAVECAR_SUB" .or. &
                      file_name(1:11) .eq. "WAVECAR_TIP" .or. &
                      file_name(1:11) .eq. "WAVECAR_ENV") then
               if (l_env_file_exist) then
                  l_env_bands_chosen = .True.
                  do i=1, num_chosen_env_bands
                     if (iband .eq. chosen_env_bands(i)) then
                        l_env_bands_chosen = .False.
                        exit
                     end if
                  end do
                  !
                  if (l_env_bands_chosen) then
                    cycle
                  end if
                  cnt_env_band_read(isp) = cnt_env_band_read(isp) + 1
                  read_band_tmp(iband,isp) = iband
               else
                  if (&!(occ(iband,isp) .eq. 0.d0) .or. & ! null state
                     !(occ(iband,isp) .eq. 1.d0) .or. & ! full state
                     !(dreal(cener(iband,isp)) .gt. fermi_level) .or. & !
                     !higher than fermi level
                     (dreal(cener(iband,isp)-fermi_level) .lt.&
                      energy_windows_i) .or. &
                     (dreal(cener(iband,isp)-fermi_level) .gt.&
                      energy_windows_f)) then  ! near fermi level
                     cycle
                  end if
                  cnt_env_band_read(isp) = cnt_env_band_read(isp) + 1
                  read_band_tmp(iband,isp) = iband
               end if
            end if
            !
            num_file = num_file + 1
            open(unit = num_file, name = title, iostat = iost)
         end if
         if (iost .ne. 0) then
            write(6,*) 'fatal error ecountered in opening a new file',&
                        title
            stop
         end if
         write(num_file,*) '# number of spin:'
         write(num_file,*) nspin
         write(num_file,*) '# number of k points:'
         write(num_file,*) nwk
         write(num_file,*) '# number of bands:'
         write(num_file,*) nband
         write(num_file,*) '# number of plane waves:'
         write(num_file,*) nplane
         write(num_file,*) '# cut energy:'
         write(num_file,*) ecut
         write(num_file,*) '# real space lattice vectors:'
         write(num_file,*) '# a1'
         write(num_file,*) ((a1(j)), j=1, 3)
         write(num_file,*) '# a2'
         write(num_file,*) ((a2(j)), j=1, 3)
         write(num_file,*) '# a3'
         write(num_file,*) ((a3(j)), j=1, 3)  
         write(num_file,*) '# recip. lattice vectors:'
         write(num_file,*) '# b1'
         write(num_file,*) ((b1(j)), j=1, 3)
         write(num_file,*) '# b2'
         write(num_file,*) ((b2(j)), j=1, 3)
         write(num_file,*) '# b3'
         write(num_file,*) ((b3(j)), j=1, 3)
         write(num_file,*) '# k value' 
         write(num_file,*) ((wk(j)), j=1, 3)
         write(num_file,*) '# eigen-energy of this band state'
         write(num_file,*) dble(cener(iband,isp)), &
                           dimag(cener(iband,isp))
         write(num_file,*) '# occupancy number of this band state'
         write(num_file,*) occ(iband,isp)
         write(num_file,*)
         !
         if (lread_para .eq. 1) then
             close(num_file)
             close(10)
             return
         end if
         do iplane = 1, nplane
            !write(11,570) (igall(j, iplane), j=1,3), coeff(iplane,
            !iband)
            write(num_file,590) (igall(j, iplane), j=1,3),& 
                                coeff(iplane, iband)
         end do
         close(num_file)
      end do
   end do
end do
!
if (lread_para .eq. 0 .and. file_name(1:11) .eq. "WAVECAR_SYS") then
   if (l_chosen_sys_bands) then
      if (cnt_sys_band_read(1) .ne. num_chosen_sys_bands .and.&
          cnt_sys_band_read(2) .ne. num_chosen_sys_bands) then
         write(6,*) 'Error'
         write(6,600) cnt_sys_band_read(1),& 
                      cnt_sys_band_read(2),&
                      num_chosen_sys_bands
         stop
      end if
   end if
   open(unit=31,file='to_be_read_sys_bands_sp1')
   rewind(31)
   do i=1, nband
      if (read_band_tmp(i,1) .ne. 0) then
         write(31,610) read_band_tmp(i,1),&
                       occ(i,1), dreal(cener(i,1)),&
                       occ(i,2), dreal(cener(i,2))
      end if
   end do
   close(31)
   open(unit=32,file='to_be_read_sys_bands_sp2')
   !
   rewind(32)
   do i=1, nband
      if (read_band_tmp(i,2) .ne. 0) then
         write(32,610) read_band_tmp(i,2),&
                       occ(i,1), dreal(cener(i,1)),&
                       occ(i,2), dreal(cener(i,2))
      end if
   end do
   close(32)
else if (file_name(1:11) .eq. "WAVECAR_SUB" .or. &
         file_name(1:11) .eq. "WAVECAR_TIP" .or. &
         file_name(1:11) .eq. "WAVECAR_ENV") then
   if (l_chosen_env_bands) then
      if (cnt_env_band_read(1) .ne. num_chosen_env_bands .or.&
          cnt_env_band_read(2) .ne. num_chosen_env_bands) then
         write(6,*) 'Error'
         write(6,600) cnt_env_band_read, num_chosen_env_bands
         stop
      end if
   end if
   open(unit=33,file='to_be_read_env_bands_sp1')
   rewind(33)
   do i=1, nband
      if (read_band_tmp(i,1) .ne. 0) then
         write(33,610) read_band_tmp(i,1), &
                       occ(i,1), dreal(cener(i,1)),&
                       occ(i,2), dreal(cener(i,2))
      end if
   end do
   close(33)  
   open(unit=34,file='to_be_read_env_bands_sp2')
   rewind(34)
   do i=1, nband
      if (read_band_tmp(i,2) .ne. 0) then
         write(34,610) read_band_tmp(i,2), &
                       occ(i,1), dreal(cener(i,1)),&
                       occ(i,2), dreal(cener(i,2))
      end if
   end do
   close(34)
end if

!call flush(11)
!close(11)

write(6,*)
write(6,*) 'number of the generated output files:', (num_file - 100)

560 format('( ', g14.6,' , ', g14.6,' ) ', g14.6)
570 format(3i6, '  ( ', g14.6, ' , ', g14.6, ' )')
580 format(g20.12, ' ', g20.12, ' ', g20.12, ' ', g20.12)
590 format(3i6, ' ', g21.13, ' ', g21.13)
600 format('number of reading bands:',2x,2(I6,4x),&
           'number of given bands:',2x,I6)
610 format('No. of band',1x,I5,2x,&
           ',occ. spin-1',2x,F9.6,2x,&
           ',eigenergy',2x,F11.6,2x,&
           ',occ. spin-2',2x,F9.6,2x,&
           ',eigenergy',2x,F11.6)
return

end subroutine read_wavecar

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

subroutine generate_wavefunc(ispin,iwk,sys_nband,&
           sys_nplane,nthread)
  use omp_lib
  implicit none
  complex*8,  allocatable     :: sys_coeff(:,:)
  real*8,     allocatable     :: freqx(:), freqy(:), freqz(:)
  integer,    allocatable     :: sys_igall(:,:,:)
  integer                     :: i,j,k, ndimx, ndimy, ndimz, nrecl
  integer                     :: g_mesh(3), sys_file
  integer                     :: nthread
  real*8                      :: dtmp1, dtmp2, dtmp3, dtmp4
  real*8                      :: sys_occ
  integer                     :: ntmp1, ntmp2, m, n
  integer                     :: sys_nband, sys_nplane
  integer                     :: istat, gall_tmp(3), sys_igall_tmp(3)
  integer                     :: ispin, iwk, sys_iband, sys_iplane
  complex*16                  :: coupling_tmp1, coupling_tmp2, cdtmp1
  complex*16                  :: cdtmp2
  complex*8                   :: sys_eig
  !
  integer :: isp, nk, nb, ni, nl, itemp
  integer :: ix, iy, iz, nx, ny, nz, tid
  integer :: nspin, nband, nwk, nionp, lm
  real*8  :: sys_fermi, delta_scaling
  real*8, allocatable :: vkpt(:,:), wtkpt(:)
  real*8, allocatable :: sys_eig_tot(:,:), sys_occ_tot(:,:)

  character(len=2)  :: ch_spin
  character(len=5)  :: ch_wk
  character(len=5)  :: ch_sys_band
  character(len=50) :: sys_title
  character(len=256) :: ch_tmp
  character(len=80) ::  dump

  logical :: l_chosen_sys_bands, l_sys_file_exist, lfile_exist
  logical :: lfile_open, lz_skip, ltmp, ltmp2
  integer :: sys_nband_save, nz_start, nz_end
  integer :: num_chosen_sys_bands
  integer, allocatable :: chosen_sys_bands(:)
  !
  real*8  :: a1(3), a2(3), a3(3)
  real*8  :: b1(3), b2(3), b3(3)
  real*8  :: vtmp(3), sumkg(3), Vcell, z_up_limit, z_low_limit
  real*8  :: delx, dely, delz, wkx, wky, wkz, x, y, z
  real*8, allocatable :: space_mesh(:,:,:,:)
  complex*16, allocatable :: sys_func(:,:,:)
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
  sys_file = 101
  cdtmp1 = 0.d0
  cdtmp2 = 0.d0
  write(ch_spin,'(I1)') ispin
  write(ch_wk,'(I5)') iwk
  write(6,*) 'to be calculate number of system bands', sys_nband
  allocate(sys_coeff(sys_nplane,sys_nband))
  allocate(sys_igall(3,sys_nplane,sys_nband))
  sys_coeff = 0
  sys_igall = 0
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
  !write(6,*) 'begin calculating coupling ...'
  !write(6,*)
  !
  !call system_clock(count1, count_rate, count_max)
  !write(6,*) 'count running time setting'
  !write(6,*) 'count rate =', count_rate
  !write(6,*) 'count max  =', count_max
  !write(6,*)
  !
  write(6,*) 'output wavefunction in 3D space'
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
  z_up_limit = (a1(3) + a2(3) + a3(3))/dble(g_mesh(3)) * (0.01d0+dble(g_mesh(3))) 
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
  write(6,'("band is ", I5)') chosen_sys_bands(ni)
  !
  write(6,'("x varying from 0.d0 to ",E15.6," by step ", I5, " with length ", E15.6)')&
  (a1(1) + a2(1) + a3(1)), nx, delx
  write(6,'("y varying from 0.d0 to ",E15.6," by step ", I5, " with length ", E15.6)')&
  (a1(2) + a2(2) + a3(2)), ny, dely
  write(6,'("z varying from 0.d0 to ",E15.6," by step ", I5, " with length ", E15.6)')&
  (a1(3) + a2(3) + a3(3)), nz, delz
  write(6,'("upper limit in the z dierction ",F12.6)') z_up_limit
  write(6,'("lower limit in the z dierction ",F12.6)') z_low_limit
  !
  !if (ltmp1 .and. ltmp2 .and. ltmp3 .and. ltmp4) then
  !  if (nz_end .gt. nz) then
  !    nz_end = nz
  !  end if
  !  if (nz_start .lt. 1) then
  !   nz_start = 1
  !  end if
  !  z_up_limit = (nz_end-1) * delz + 0.01d0 * delz
  !  z_low_limit = (nz_start-1) * delz - 0.01d0 * delz
  !  write(6,'("upper limit in the z dierction ",F12.6)') z_up_limit
  !  write(6,'("lower limit in the z dierction ",F12.6)') z_low_limit
  !else if (.not. ltmp1 .and. .not. ltmp2 .and. ltmp3 .and. ltmp4) then
  !  z_up_limit = (nz_end-1) * delz + 0.01d0 * delz
  !  z_low_limit = (nz_start-1) * delz - 0.01d0 * delz
  !  write(6,'("upper limit in the z dierction ",F12.6)') z_up_limit
  !  write(6,'("lower limit in the z dierction ",F12.6)') z_low_limit
  !else if (ltmp1 .and. ltmp2 .and. .not. ltmp3 .and. .not. ltmp4) then
  !  if (nz_end .gt. nz) then
  !    nz_end = nz
  !  end if
  !  if (nz_start .lt. 1) then
  !   nz_start = 1
  !  end if
  !  write(6,'("upper limit in the z dierction ",F12.6)') z_up_limit
  !  write(6,'("lower limit in the z dierction ",F12.6)') z_low_limit 
  !else if (ltmp1 .and. .not. ltmp2 .and. .not. ltmp3 .and. .not. ltmp4) then
  !end if
  write(6,*)
  allocate(sys_func(nx,ny,nz))!,space_mesh(nx,ny,nz,3))
  !
  sys_func = 0
  !$omp parallel default(shared)
  !nthread = omp_get_num_threads()
  if (nthread >= 48) then
     nthread = 48
  end if
  !$omp end parallel
  call omp_set_num_threads(nthread)
  write(6,*) 'the number of parallel threads,', nthread
  write(6,*)
  !
  !$omp parallel do private(tid,ix,iy,iz,x,y,z,lz_skip,wkx,wky,wkz,sumkg,cdtmp1,nb,j) default(shared)
  do iz=nz_start, nz_end
    z = (iz - 1) * delz
    !$omp critical
    lz_skip = .False.
    tid = omp_get_thread_num()
    if (z .gt. z_up_limit .or. z .lt. z_low_limit) then
      write(6,'("set zero-value for wavefunction at z=", F12.6," on tid ",I4)') z,tid
      lz_skip = .True.
    else
      write(6,'("calculate value for wavefunction at z = ",F12.6," on tid ", I4)') z,tid
    end if
    !$omp end critical
    if (.not. lz_skip) then
    do iy=1, ny
      y = (iy - 1) * dely
      !$omp critical
      if (mod(iy,100) .eq. 0 ) then
        write(6,'("(y, z) = (",F12.6," ,",F12.6,") on tid ",I4)') y, z, tid
      end if
      !$omp end critical
      call flush(6)
      do ix=1, nx
        x = (ix - 1) * delx
        !space_mesh(ix,iy,iz,1) = x
        !space_mesh(ix,iy,iz,2) = y
        !space_mesh(ix,iy,iz,3) = z
        !
        cdtmp1 = dcmplx(0.0d0, 0.0d0)
        do nb=1, sys_nplane
          wkx = sys_igall(1,nb,ni)
          wky = sys_igall(2,nb,ni)
          wkz = sys_igall(3,nb,ni)
          do j=1,3
            sumkg(j) = (wkx) * b1(j) + (wky) * b2(j) + (wkz) * b3(j)
           end do
          cdtmp1 = cdtmp1 + sys_coeff(nb,ni) * &
                     cdexp(dcmplx(0.0d0, 1.0d0) *&
                   (sumkg(1) * x + sumkg(2) * y + sumkg(3) * z))
        end do
        !write(6,'("(ix, iy, iz) = ("I4," ,",I4," ,",I4,") on tid ",I4)') ix, iy, iz, tid
        sys_func(ix,iy,iz) = cdtmp1 / dsqrt(Vcell)
      end do
    end do
    end if
  end do
  !$omp end parallel do
  !
  ch_sys_band = ''
  !sys_file = sys_file + 2
  write(ch_sys_band,'(I5)') sys_iband
  goto 114
  sys_title = ('wavefunc-spin-'//trim(adjustl(ch_spin))//'-wk-'//&
          trim(adjustl(ch_wk))//'-band-'//trim(adjustl(ch_sys_band)))
  open(unit=sys_file,file=sys_title)
  rewind(sys_file)
  do iz=1, nz
    z = (iz - 1) * delz
    do iy=1, ny
      y = (iy - 1) * dely
      do ix=1, nx
        x = (ix - 1) * delx
        write(sys_file,590) x,y,z,&
                            dreal(sys_func(ix,iy,iz)),&
                            dimag(sys_func(ix,iy,iz))
      end do
      call flush(sys_file)
    end do
  end do 
  close(sys_file)
  write(6,*) 'output wavefunction done'
  !
  114 continue
  sys_file = sys_file + 1
  sys_title = ('wavefunc-real-spin-'//trim(adjustl(ch_spin))//'-wk-'//&
            trim(adjustl(ch_wk))//'-band-'//trim(adjustl(ch_sys_band))//'.vasp')
  open(unit=sys_file,file=sys_title)
  rewind(sys_file)
  write(sys_file,*) nx, ny, nz
  write(sys_file,'(5(E14.5))') (((dreal(sys_func(ix,iy,iz)),&
                        ix=1,nx),iy=1,ny),iz=1,nz)
  call flush(sys_file)
  close(sys_file)
  !
  goto 115
  sys_file = sys_file + 1
  sys_title = ('wavefunc-real-spin-'//trim(adjustl(ch_spin))//'-wk-'//&
            trim(adjustl(ch_wk))//'-band-'//trim(adjustl(ch_sys_band))//'-check.vasp')
  open(unit=sys_file,file=sys_title)
  rewind(sys_file)
  write(sys_file,*) 'nx'
  write(sys_file,'(5(E14.5))') (((dreal(sys_func(ix,iy,iz)),&
                        ix=1,nx),iy=1,1),iz=1,1)
  write(sys_file,*) 'ny'
  write(sys_file,'(5(E14.5))') (((dreal(sys_func(ix,iy,iz)),&
                        ix=1,1),iy=1,ny),iz=1,1)
  write(sys_file,*) 'nz'
  write(sys_file,'(5(E14.5))') (((dreal(sys_func(ix,iy,iz)),&
                        ix=1,1),iy=1,1),iz=1,nz)
  call flush(sys_file)
  close(sys_file)
  115 continue
  write(6,*) 'output wavefunction real part done'
  !
  sys_file = sys_file + 1
  sys_title = ('wavefunc-imag-spin-'//trim(adjustl(ch_spin))//'-wk-'//&
            trim(adjustl(ch_wk))//'-band-'//trim(adjustl(ch_sys_band))//'.vasp')
  open(unit=sys_file,file=sys_title)
  rewind(sys_file)
  write(sys_file,*) nx, ny, nz
  write(sys_file,'(5(E14.5))') (((dimag(sys_func(ix,iy,iz)),&
                    ix=1,nx),iy=1,ny),iz=1,nz) 
  call flush(sys_file)
  close(sys_file)
  !
  goto 116
  sys_file = sys_file + 1
  sys_title = ('wavefunc-imag-spin-'//trim(adjustl(ch_spin))//'-wk-'//&
            trim(adjustl(ch_wk))//'-band-'//trim(adjustl(ch_sys_band))//'-check.vasp')
  open(unit=sys_file,file=sys_title)
  rewind(sys_file)
  write(sys_file,*) 'nx'
  write(sys_file,'(5(E14.5))') (((dimag(sys_func(ix,iy,iz)),&
                        ix=1,nx),iy=1,1),iz=1,1)
  write(sys_file,*) 'ny'
  write(sys_file,'(5(E14.5))') (((dimag(sys_func(ix,iy,iz)),&
                        ix=1,1),iy=1,ny),iz=1,1)
  write(sys_file,*) 'nz'
  write(sys_file,'(5(E14.5))') (((dimag(sys_func(ix,iy,iz)),&
                        ix=1,1),iy=1,1),iz=1,nz)
  call flush(sys_file)
  close(sys_file)
  !
  116 continue
  write(6,*) 'output wavefunction imaginary part done'
  !
  write(6,*)
  ! output wavefunction
  !
  570 format('computation information', 3i10)
  580 format(3F12.7, F20.10, ' s')
  590 format(3(F16.8,2x), 2(E20.11,2x))
end subroutine generate_wavefunc


program main
  implicit none
  integer  :: ncpus, ispin, nspin, iwk
  integer  :: sys_band, sys_plane
  logical  :: lfile_exist
  
  nspin = 2
  iwk = 1
  ispin = 1
  ncpus = 16
  !
  write(6,*) 'Version 1.0'
  open(unit=10,file='number_cpus')
  rewind(10)
  read(10,*) ncpus
  close(10)
  !
  inquire(file='sys_para', exist=lfile_exist)
  if (lfile_exist) then
     open(unit=7,file='sys_para')
     read(7,*) sys_band 
     read(7,*) sys_plane
     close(7)
  end if
  !
  do ispin=1, nspin
     write(6,*)
     write(6,*) '************************'
     write(6,*) 'calculating for spin', ispin
     inquire(file='OUTCAR_SYS', exist=lfile_exist)
     if (lfile_exist) then
       call generate_wavefunc(ispin,iwk,sys_band,&
            sys_plane,ncpus)
     end if
     write(6,*) '************************'
     write(6,*) 
  end do
  !
end program main
