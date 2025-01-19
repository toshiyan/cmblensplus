!////////////////////////////////////////////////////!
! * Module to generate random number
! * This code is part of forutils/RandUtils.f90
!////////////////////////////////////////////////////!

module random
  implicit none
  integer :: Feedback = 1, rand_inst = 0 
  logical, parameter :: use_ziggurat = .false.

contains


subroutine poisson(mu,f)
! generate random number along with the Poisson distribution
  implicit none
  double precision, intent(in) :: mu
  integer, intent(out) :: f
  double precision :: p, rnd

  p = mu
  f = -1
  do while (p.ge.0d0)
    call random_number(rnd)
    if (rnd.le.0d0)  call random_number(rnd)
    p = p + dlog(rnd)
    f = f + 1 
  end do

end subroutine poisson


subroutine initRandom(i,i2)
  implicit none
  integer, optional, intent(IN) :: i
  integer, optional, intent(IN) :: i2
  integer seed_in,kl,ij
  character(len=10) :: fred
  real :: klr
  
  if (present(i)) then
    seed_in = i
  else
    seed_in = -1
  end if

  if (seed_in /=-1) then
    if (present(i2)) then
      kl = i2
    else
      kl = 9373
    end if
    ij = i
  else
    call system_clock(count=ij)
    ij = mod(ij + rand_inst*100, 31328)
    call date_and_time(time=fred)
    read (fred,'(e10.3)') klr
    kl = mod(int(klr*1000), 30081)       
  end if

  call rmarin(ij,kl)

end subroutine initRandom


subroutine RandIndices(indices, nmax, n)
  integer, intent(in) :: nmax, n
  integer indices(n),i, ix
  integer tmp(nmax)
 
  do i=1, nmax
    tmp(i)=i
  end do
  do i=1, n
    ix = int(ranmar()*(nmax +1 -i)) + 1
    indices(i) = tmp(ix)
    tmp(ix) = tmp(nmax+1-i)
  end do

end subroutine RandIndices


subroutine RandRotation(R, N)
  integer, intent(in) :: N
  real R(N,N), vec(N), norm
  integer i,j
    
  do j = 1, N
    do
      do i = 1, N
        vec(i) = Gaussian1()
      end do
      do i = 1, j-1
        vec = vec - sum(vec*R(i,:))*R(i,:)
      end do
      norm = sum(vec**2)
      if (norm > 1e-3) exit
    end do
    R(j,:) = vec / sqrt(norm)
  end do
    
end subroutine RandRotation


double precision function GAUSSIAN1()
  implicit none
  double precision R, V1, V2, FAC
  integer, save :: iset = 0
  double precision, save :: gset

  !Box muller
  if (ISET==0) then
    R=2
    do while (R >= 1.d0)
      V1=2.d0*ranmar()-1.d0
      V2=2.d0*ranmar()-1.d0
      R=V1**2+V2**2
    end do
    FAC=sqrt(-2.d0*log(R)/R)
    GSET=V1*FAC
    GAUSSIAN1=V2*FAC
    ISET=1
  else
    GAUSSIAN1=GSET
    ISET=0
  endif

end function GAUSSIAN1


real FUNCTION RANDEXP1()
! Random-number generator for the exponential distribution
! Algorithm EA from J. H. Ahrens and U. Dieter,
! Communications of the ACM, 31 (1988) 1330--1337.
! Coded by K. G. Hamilton, December 1996, with corrections.
  real u, up, g, y
  real, parameter ::   alog2= 0.6931471805599453
  real, parameter ::      a = 5.7133631526454228
  real, parameter ::      b = 3.4142135623730950
  real, parameter ::     c = -1.6734053240284925
  real, parameter ::      p = 0.9802581434685472
  real, parameter ::     aa = 5.6005707569738080
  real, parameter ::     bb = 3.3468106480569850
  real, parameter ::     hh = 0.0026106723602095
  real, parameter ::     dd = 0.0857864376269050

  u = ranmar()
  do while (u.le.0)                 ! Comment out this block 
    u = ranmar()                    ! if your RNG can never
  end do                            ! return exact zero
  g = c
  u = u+u
  do while (u.lt.1.0)
    g = g + alog2
    u = u+u
  end do
  u = u-1.0
  if (u.le.p) then
    randexp1 = g + aa/(bb-u)
    return
  end if
  do
    u = ranmar()
    y = a/(b-u)
    up = ranmar()
    if ((up*hh+dd)*(b-u)**2 .le. exp(-(y+c))) then
      randexp1 = g+y
      return
    end if
  end do

end function randexp1


subroutine RMARIN(IJ,KL)
! This is the initialization routine for the random number generator RANMAR()
! NOTE: The seed variables can have values between:    0 <= IJ <= 31328
!                                                      0 <= KL <= 30081
!The random number sequences created by these two seeds are of sufficient 
! length to complete an entire calculation with. For example, if sveral 
! different groups are working on different parts of the same calculation,
! each group could be assigned its own IJ seed. This would leave each group
! with 30000 choices for the second seed. That is to say, this random 
! number generator can create 900 million different subsequences -- with 
! each subsequence having a length of approximately 10^30.
!
! Use IJ = 1802 & KL = 9373 to test the random number generator. The
! subroutine RANMAR should be used to generate 20000 random numbers.
! Then display the next six random numbers generated multiplied by 4096*4096
! If the random number generator is working properly, the random numbers
!    should be:
!           6533892.0  14220222.0  7275067.0
!           6172232.0  8354498.0   10633180.0
      double precision U(97), C, CD, CM, S, T
      integer I97, J97,i,j,k,l,m
      integer ij,kl
      integer ii,jj
           
    
!      INTEGER IRM(103)
      
      common /RASET1/ U, C, CD, CM, I97, J97
      if( IJ .lt. 0  .or.  IJ .gt. 31328  .or. &
         KL .lt. 0  .or.  KL .gt. 30081 ) then
          print '(A)', ' The first random number seed must have a value  between 0 and 31328'
          print '(A)',' The second seed must have a value between 0 and   30081'
            stop
      endif
      I = mod(IJ/177, 177) + 2
      J = mod(IJ    , 177) + 2
      K = mod(KL/169, 178) + 1
      L = mod(KL,     169) 
      do 2 II = 1, 97
         S = 0.0
         T = 0.5
         do 3 JJ = 1, 24
            M = mod(mod(I*J, 179)*K, 179)
            I = J
            J = K
            K = M
            L = mod(53*L+1, 169)
            if (mod(L*M, 64) .ge. 32) then
               S = S + T
            endif
            T = 0.5 * T
3        continue
         U(II) = S
2     continue
      C = 362436.0 / 16777216.0
      CD = 7654321.0 / 16777216.0
      CM = 16777213.0 /16777216.0
      I97 = 97
      J97 = 33
end subroutine RMARIN


double precision function RANMAR()
! This is the random number generator proposed by George Marsaglia in 
! Florida State University Report: FSU-SCRI-87-50
! It was slightly modified by F. James to produce an array of pseudorandom numbers.
  double precision U(97), C, CD, CM
  integer I97, J97
  double precision uni
   
  common /RASET1/ U, C, CD, CM, I97, J97
!      INTEGER IVEC
  UNI = U(I97) - U(J97)
  if( UNI .lt. 0.0 ) UNI = UNI + 1.0
  U(I97) = UNI
  I97 = I97 - 1
  if(I97 .eq. 0) I97 = 97
  J97 = J97 - 1
  if(J97 .eq. 0) J97 = 97
  C = C - CD
  if( C .lt. 0.d0 ) C = C + CM
  UNI = UNI - C
  if( UNI .lt. 0.d0 ) UNI = UNI + 1.0 ! bug?
  RANMAR = UNI
      
end function RANMAR


end module random
  
