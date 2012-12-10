program nrvs

   implicit none
   integer nat,iat
   parameter(nat=36)
   real freq(3*nat), vec(3*nat,3*nat), mass(3*nat), num, den, int(3*nat), &
         gau, f, s, rmass(3*nat), x(3*nat), massr(nat)
   integer iset,i,j,jat,l,k,nmodes,il,imass(3*nat),nat3,mode
   character(len=18) label(3*nat)
   character(len=256) line
   character(len=2) atsymb
   logical adf, gaussian, jaguar, orca

adf = .false.
gaussian = .false.
jaguar = .false.
orca = .true.

if (adf) then
   i = 0
   do iset=1,nat-2

      read(5,*) freq(i+1),freq(i+2),freq(i+3)
      read(5,*)
      j = 0
      do jat=1,nat
         read(5,*) label(jat),vec(j+1,i+1),vec(j+2,i+1),vec(j+3,i+1),  &
                              vec(j+1,i+2),vec(j+2,i+2),vec(j+3,i+2),  &
                              vec(j+1,i+3),vec(j+2,i+3),vec(j+3,i+3)
         j = j + 3
      end do
      read(5,*)
      read(5,*)
      i = i + 3
   end do
   write(6,'(6f10.3)') (vec(j,1), j=1,12)
   write(6,'(6f10.3)') (vec(j,6), j=1,12)
   write(6,'(6f10.3)') (vec(j,13), j=1,12)

   j = 0
   do jat=1,nat
      l = len_trim( label(jat) )
      if( label(jat)(l-1:l) == "Fe" ) then
         mass(j+1:j+3) = 55.85
      else if( label(jat)(l:l) == "S" ) then
         mass(j+1:j+3) = 32.06
      else if( label(jat)(l:l) == "C" ) then
         mass(j+1:j+3) = 12.01
      else if( label(jat)(l:l) == "H" ) then
         mass(j+1:j+3) =  1.008
      else
         write(6,*) 'unknown element: ',label(jat)
         stop
      end if
      write(6,*) jat, mass(j+1)
      j = j + 3
   end do

    do i=1,3*nat-6
       s = 0.0
       do j=1,3*nat
          s = s + mass(j)*vec(j,i)*vec(j,i)
       end do
       f = 1.d0/sqrt(s)
       vec(1:3*nat,i) = vec(1:3*nat,i)*f
   end do
   nmodes = 3*nat-6

else if( gaussian ) then
   i = 0
   do iset=1,(3*nat-6)/5

      read(5,*)
      read(5,*)
      read(5,'(23x,5f10.4)') freq(i+1),freq(i+2),freq(i+3),freq(i+4),freq(i+5)
      read(5,'(23x,5f10.4)') rmass(i+1),rmass(i+2),rmass(i+3),rmass(i+4),rmass(i+5)
      read(5,*)
      read(5,*)
      read(5,*)
      do j=1,3*nat
         read(5,*) il,jat,imass(j), vec(j,i+1),vec(j,i+2),vec(j,i+3),  &
                           vec(j,i+4),vec(j,i+5)
      end do
      i = i + 5
   end do
   nmodes = i

   do j=1,3*nat
      if( imass(j) == 26 ) then
         mass(j) = 55.85
      else if( imass(j) == 16 ) then
         mass(j) = 32.06
      else if( imass(j) ==  6 ) then
         mass(j) = 12.01
      else if( imass(j) ==  1 ) then
         mass(j) = 1.008
      else
         write(6,*) 'unknown element: ',j,imass(j)
         stop
      end if

      do i=1,nmodes
         vec(j,i) = vec(j,i)/sqrt(rmass(i))
      end do
   end do

else if( orca ) then

   do i=1,99999
      read(5,'(a)', end=99) line
         if( line(1:23) .eq. 'MULLIKEN ATOMIC CHARGES' ) exit
   end do
   read(5,*)
   j = 0
   do i=1,nat
      read(5,'(5x,a2)') atsymb
      if( atsymb == 'H ') then
         mass(j+1:j+3) = 1.008
      else if( atsymb == 'C ') then
         mass(j+1:j+3) = 12.010
      else if( atsymb == 'N ') then
         mass(j+1:j+3) = 14.001
      else if( atsymb == 'O ') then
         mass(j+1:j+3) = 16.00
      else if( atsymb == 'Rh') then
         mass(j+1:j+3) = 102.91
      else
         write(0,*) 'unknown element: ',atsymb
         stop
      end if
!     write(0,*) i, mass(i)
      j = j + 3
   end do
   do i=1,99999
      read(5,'(a)', end=99) line
         if( line(1:23) .eq. 'VIBRATIONAL FREQUENCIES' ) exit
   end do
   read(5,*)
   read(5,*)
   do i=1,3*nat
      read(5,'(8x,f10.2)') freq(i)
      if( i > 6 ) freq(i-6) = freq(i)
   end do
   do i=1,10
      read(5,*)
   end do
   i = 0
   do iset=1,(3*nat)/6
      read(5,*)
      do j=1,3*nat
         read(5,*) jat, vec(j,i+1),vec(j,i+2),vec(j,i+3),  &
                        vec(j,i+4),vec(j,i+5),vec(j,i+6)
      end do
      i = i + 6
      if( iset == 1 ) i=0
   end do
   nmodes = i
   do i=1,nmodes
      num = 0.0
      do k = 1,3*nat
         num = num + mass(k)*vec(k,i)*vec(k,i)
      end do
      num = 1./sqrt(num)
      do k = 1,3*nat
         vec(k,i) = num*vec(k,i)
      end do
   end do

else if (jaguar) then

   j=0
   do iat=1,nat
      read(5,*) label(iat),x(j+1),x(j+2),x(j+3)
      j = j + 3
   end do

   i = 0
   do iset=1,nat/2-1

      read(5,'(18x,6f9.2)') &
          freq(i+1),freq(i+2),freq(i+3),freq(i+4),freq(i+5),freq(i+6)
      read(5,'(18x,6f9.2)') &
          rmass(i+1),rmass(i+2),rmass(i+3),rmass(i+4),rmass(i+5),rmass(i+6)
      read(5,*)
      do j=1,3*nat
         read(5,'(a18,6f9.5)') label(j), &
                               vec(j,i+1),vec(j,i+2),vec(j,i+3),  &
                               vec(j,i+4),vec(j,i+5),vec(j,i+6)
      end do
      read(5,*)
      i = i + 6
   end do

   do j=1,3*nat
      if( label(j)(5:6) == "Fe" ) then
         mass(j) = 55.85
      else if( label(j)(5:6) == "Mo" ) then
         mass(j) = 95.94
      else if( label(j)(5:5) == "S" ) then
         mass(j) = 32.06
      else if( label(j)(5:5) == "C" ) then
         mass(j) = 12.01
      else if( label(j)(5:5) == "N" ) then
         mass(j) = 14.01
      else if( label(j)(5:5) == "O" ) then
         mass(j) = 16.00
      else if( label(j)(5:5) == "H" ) then
         mass(j) =  1.008
      else
         write(6,*) 'unknown element: ',label(jat)
         stop
      end if
   end do
   nmodes = 6*(nat/2-1)

else   ! amber modes

   read(5,'(5e16.8)') (massr(k), k=1,nat)
   i = 0
   do k=1,nat
      mass(i+1) = massr(k)
      mass(i+2) = massr(k)
      mass(i+3) = massr(k)
      i = i + 3 
   end do

   read(5,*)  ! title
   read(5,*) nat3
   read(5,'(7f11.5)') (x(i), i=1,nat3)

   nmodes = nat3-6
   do i=1,6
      read(5,*) !asterisks
      read(5,*) mode, freq(i)
      read(5,'(7f11.5)') (vec(j,i),j=1,nat3)
   end do

   do i=1,nmodes
      read(5,*) !asterisks
      read(5,*) mode, freq(i)
      read(5,'(7f11.5)') (vec(j,i),j=1,nat3)
   end do

end if

!  check orthonormalization:

   do i=1,nmodes
      do j=i,nmodes
         num = 0.0
         do k = 1,3*nat
            num = num + mass(k)*vec(k,i)*vec(k,j)
         end do
         if ( j==i ) then
            if (abs(num - 1.0) > 0.001) write(6,'(i5,f10.6)' ) j,num
         else
            if (abs(num) > 0.001) write(6,'(10x, 2i5,f10.6)' ) i,j,num
         end if
      end do
   end do
         
!  get the intensity factor for the Fe atoms

   do i=1,nmodes
      num = 0.0
      do j=1,3*nat
         if( abs(mass(j)-102.91)<0.1 ) num = num + mass(j)*vec(j,i)*vec(j,i)
      end do
      int(i) = num
      write(6,'(i4,2f10.3)') i, freq(i), num
   end do

!  now, convolute this with an 8 cm**-1 Gaussian: naive implementation for now

   do k=1,800
      f = k
      s = 0.0
      do i=1,nmodes
         gau = exp( -(f - freq(i))**2/128.)
         s = s + int(i)*gau
      end do
      write(6,'(f6.1,f10.4)' ) f, s
   end do
   stop

99 write(0,*) 'end of file found on input'
   stop

end program nrvs
