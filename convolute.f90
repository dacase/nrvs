program convolute

! take an nrvs.stk file from the output file of orca_vib, here called
!   <basename>.ovout
! convolute with an 8 cm**-1 line width, create a corresonding stick file

   implicit none
   integer i,k,nmodes
   real freq(1000), int(1000), gau, f, s, int_comp(1000), int_vdos(1000)
   character(len=7) label
   character(len=256) line,basename

   call get_command_argument(1, basename)
   open(unit=5,file=trim(basename)//'.ovout')
   open(unit=6,file=trim(basename)//'.nrvs.dat')
   open(unit=7,file=trim(basename)//'.nrvs.stk')

   do i=1,99999
      read(5,'(a)', end=99) line
         if( line(1:38) .eq. 'Fe(57) NORMAL MODE COMPOSITION FACTORS' ) exit
   end do
   do i=1,3
      read(5,*)
   end do
   do i=1,248
      read(5,'(9x,f11.6)') int_comp(i)
   end do

   do i=1,99999
      read(5,'(a)', end=99) line
         if( line(1:38) .eq. 'WEIGHTED VIBRATIONAL DENSITY OF STATES' ) exit
   end do
   do i=1,4
      read(5,*)
   end do
   do i=1,248
      read(5,'(30x,f11.6)') int_vdos(i)
   end do

   do i=1,99999
      read(5,'(a)', end=99) line
         if( line(1:38) .eq. 'NUCLEAR RESONANCE VIBRATIONAL SPECTRUM' ) exit
   end do
   do i=1,4
      read(5,*)
   end do
   do i=1,1000
      read(5,'(a7,5x,f10.2,18x,f11.6)') label,freq(i),int(i)
      if( label .eq. '       ' ) exit
      write(7,'(f7.2,3f10.5)') freq(i), int(i), int_comp(i), int_vdos(i)
   end do
   nmodes = i-1
         
!  now, convolute this with an 8 cm**-1 Gaussian: naive implementation for now

   do k=1,650
      f = k
      s = 0.0
      do i=1,nmodes
         gau = exp( -(f - freq(i))**2/128.)
         s = s + int_vdos(i)*gau
      end do
      write(6,'(f6.1,f10.4)' ) f, s
   end do
   stop

99 write(0,*) 'end of file found on input'
   stop

end program convolute
