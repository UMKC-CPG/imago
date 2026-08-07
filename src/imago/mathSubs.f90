!! SPDX-License-Identifier: ECL-2.0
!! Copyright (c) 2026 Paul Rulis

module O_MathSubs

   ! Import necessary modules.
   use O_Kinds
   use O_Constants

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   ! Define module variables.
   integer(8), allocatable, dimension(:) :: preCompFactorial

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

! This is a calculation of an integral over the error function.
function error_fn(x)

   ! Use precision parameters
   use O_Kinds

   ! Make certain that no implicit variables are accidently declared.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   real (kind=double) :: x
   real (kind=double) :: error_fn

   ! Define the local variables used in this subroutine.
   real (kind=double) :: absX
   real (kind=double) :: absXSqrd
   real (kind=double) :: invAbsXSqrd
   real (kind=double) :: invSqrtPi
   real (kind=double) :: a
   real (kind=double), dimension(6) :: p, q

!-----------------------------------------------
!
!     SANDIA MATHEMATICAL PROGRAM LIBRARY
!     APPLIED MATHEMATICS DIVISION 2613
!     SANDIA LABORATORIES
!     ALBUQUERQUE, NEW MEXICO  87185
!     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!                    ISSUED BY SANDIA LABORATORIES                     *
!  *                   A PRIME CONTRACTOR TO THE                       *
!  *                UNITED STATES DEPARTMENT OF ENERGY                 *
!  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * * *
!  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE   *
!  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE      *
!  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,    *
!  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES  *
!  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL      *
!  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR     *
!  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS      *
!  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE          *
!  * OWNED RIGHTS.                                                     *
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS     *
!  * PART IS SAND77-1441.                                              *
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!     WRITTEN BY J.E. VOGEL FROM APPROXIMATIONS DERIVED BY W.J. CODY .
!
!     ABSTRACT
!
!          ERF(X) COMPUTES 2.0/SQRT(PI) TIMES THE INTEGRAL FROM 0 TO X
!          OF EXP(-X**2). THIS IS DONE USING RATIONAL APPROXIMATIONS.
!          ELEVEN CORRECT SIGNIFICANT FIGURES ARE PROVIDED.
!
!     DESCRIPTION OF PARAMETERS
!
!          X MAY BE ANY REAL VALUE
!
!     ERF IS DOCUMENTED COMPLETELY IN SC-M-70-275
!
! Note that this procedure has been modified to use if-then-else-endif
!   statements instead of GO TO statements.  It has also been changed so that
!   the appropriate constants will be defined only when explicitly needed
!   for each case.  For the Imago purposes, that great majority of the time
!   is spent in the > 6.0 case where there is no need for any of the constants
!   defined previously above.  So, those assignments have been moved into the
!   body of the function.  Finally, to improve the efficiency even more, the
!   number of arrays that need to be allocated for this routine has been
!   reduced from 6 to 2.  The p1, p2, and p3 are now just p, and the same with
!   the q1, q2, and q3.  This is because only one of each are used for any
!   given case.


   absX = abs(x)
   if (absX <= 6.0_double) then
      absXSqrd = absX * absX
      if (absX <= 4.0_double) then
         if (absX <= 0.46875_double) then
            ! Assign parameters where x is between 0.46875 and 4.
            p(1) = 242.6679552305318_double
            p(2) = 21.97926161829415_double
            p(3) = 6.996383488619136_double
            p(4) = -0.03560984370181539_double
            q(1) = 215.0588758698612_double
            q(2) = 91.16490540451490_double
            q(3) = 15.08279763040779_double
            q(4) = 1.000000000000000_double

            a = absX*(p(1) + absXSqrd * (p(2) + &
              & absXSqrd * (p(3) + absXSqrd * p(4))))
            a = a  / (q(1) + absXSqrd * (q(2) + &
              & absXSqrd * (q(3) + absXSqrd * q(4))))
            if (x < 0.0_double) a = -a
            error_fn = a
         else
            ! Assign parameters where x is between 0.46875 and 4.
            p(1) = 22.898992851659_double
            p(2) = 26.094746956075_double
            p(3) = 14.571898596926_double
            p(4) = 4.2677201070898_double
            p(5) = 0.56437160686381_double
            p(6) = -0.0000060858151959688_double
            q(1) = 22.898985749891_double
            q(2) = 51.933570687552_double
            q(3) = 50.273202863803_double
            q(4) = 26.288795758761_double
            q(5) = 7.5688482293618_double
            q(6) = 1.0000000000000_double

            a = exp((-absXSqrd)) * (p(1) + absX * (p(2) + &
              & absX * (p(3) + absX * (p(4) + absX * (p(5) + &
              & absX * p(6))))))
            a = a/(q(1) + absX * (q(2) + absX * (q(3) + absX * &
              & (q(4) + absX * (q(5) + absX * q(6))))))
            error_fn = sign(1.0_double - a,x)
         endif
      else
         ! Assign parameters where x is between 4 and 6.
         invSqrtPi = 0.564189583547756_double
         p(1) = -0.0121308276389978_double
         p(2) = -0.1199039552681460_double
         p(3) = -0.243911029488626_double
         p(4) = -0.0324319519277746_double
         q(1) = 0.0430026643452770_double
         q(2) = 0.489552441961437_double
         q(3) = 1.43771227937118_double
         q(4) = 1.00000000000000_double

         invAbsXSqrd = 1.0_double / absXSqrd
         a = invAbsXSqrd * (p(1) + invAbsXSqrd * (p(2) + invAbsXSqrd * &
           & (p(3) + invAbsXSqrd * p(4)))) / (q(1) + invAbsXSqrd * (q(2) + &
           & invAbsXSqrd * (q(3) + invAbsXSqrd * q(4))))
         a = exp(-absXSqrd) * (invSqrtPi + a) / absX
         error_fn = sign(1.0_double - a,x)
      endif
   else
      error_fn = x/absX
   endif

end function error_fn


function stepFunction (x,stepFnRange)
   use O_Kinds
   implicit none
   real (kind=double), intent(in) :: x
   real (kind=double), intent(in) :: stepFnRange
   real (kind=double) :: stepFunction

   ! Note regarding the value of the step function range. The default input
   !   value is currently 11.5. This will cause the fermi function to cut
   !   off (and not evaluate the exponential) when the fermi function would
   !   evaluate to something any closer to 1.0 than 0.99999 or any closer to
   !   0.0 than 0.00001. This can be seen by solving the stepFunction for x
   !   x=ln(1/stepFn - 1) and putting 0.99999 or 0.00001 in for stepFn.

   if     (x >  stepFnRange) then
      stepFunction = 0.0_double
   elseif (x < -stepFnRange) then
      stepFunction = 1.0_double
   else
      stepFunction = 1.0_double / (1.0_double + exp(x))
   endif
end function stepFunction


subroutine crossProduct (answer, vector1, vector2)

   ! Import the necessary kinds definitions.
   use O_Kinds

   ! Make sure that no variables are accidentally defined.
   implicit none

   ! Define the passed parameters.
   real (kind=double), dimension(3) :: answer
   real (kind=double), dimension(3) :: vector1  ! ax + ay + az
   real (kind=double), dimension(3) :: vector2  ! bx + by + bz

   ! Compute the cross product.  Note the correctly reversed sign when
   !   computing answer(2).
   answer(1) = vector1(2)*vector2(3) - vector1(3)*vector2(2) ! ay*bz - az*by
   answer(2) = vector1(3)*vector2(1) - vector1(1)*vector2(3) ! az*bx - ax*bz
   answer(3) = vector1(1)*vector2(2) - vector1(2)*vector2(1) ! ax*by - ay*bx

end subroutine crossProduct


subroutine computeFactorials(maxFact)

   implicit none

   ! Define the passed parameters.
   integer :: maxFact ! Highest factorial parameter.

   ! Define local variables.
   integer(8) :: i

   allocate (preCompFactorial(0:maxFact))

   preCompFactorial(0) = 1 ! This is 0!
   preCompFactorial(1) = 1 ! This is 1!
   do i = 2, maxFact
      preCompFactorial(i) = preCompFactorial(i-1) * i
   enddo

end subroutine computeFactorials


function wignerD(twoj, twom, twomp, eulerCoords)

   ! Use necessary modules.
   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! Define function definition variable.
   complex (kind=double) :: wignerD

   ! Define passed parameters.
   integer :: twoj, twom, twomp
   real(kind=double), dimension(3) :: eulerCoords ! alpha, beta, gamma
         ! Pulled from theta and phi of the ttpCoords.

   ! Define local variables.
   real (kind=double) :: term1Exp, term3Exp
   complex (kind=double) :: term1, term3
   real (kind=double) :: term2


   ! Condition 1
!   if ((j < 0) .or. (.not. (modulo(j,1.0) == 0)) .or. &
!         & (modulo(j,1.0) == 0.5 .and. ((modulo(m,1.0) /= 0) &
!         & .or. (modulo(mp,1.0) /= 0)))) then
!      write (20,*) "Invalid input parameters:"
!      write (20,fmt="(a,3e10.3)") "j, m, mp: ", j, m, mp
!      write (20,*) "Parameter j must be non-negative integer or half-integer."
!      write (20,*) "Parameters m and mp must be between -j and j."
!      stop
!   endif

!   ! Condition 2
!   if ((eulerCoords(1) < 0) .or. (eulerCoords(1) > 2.0 * pi) &
!         & .or. (eulerCoords(2) < 0) .or. (eulerCoords(2) > pi) &
!         & .or. (eulerCoords(3) < 0) .or. (eulerCoords(3) > 2.0 * pi)) then
!      write (20,*) "Invalid input parameters:"
!      write (20,fmt="(a,3e12.3)") "phi theta phi: ", eulerCoords(1), &
!            & eulerCoords(2), eulerCoords(3)
!      write (20,fmt="(a,3e12.3)") "phi theta phi: ",2.0*pi,2.0*pi,2.0*pi
!   endif
!   if ((theta_0 < 0) .or. (theta_0 > 2.0 * pi) .or. (theta < 0) &
!         & .or. (theta > 2.0 * pi) .or. (phi < 0) .or. (phi > 2.0 * pi)) then
!      write (20,*) "Invalid input parameters:"
!      write (20,fmt="(a,3e10.3)") "j, m, mp: ", j, m, mp
!      write (20,*) "Parameter j must be non-negative integer or half-integer."
!   endif
!write (6,fmt="(a,3i3,3e14.5)") "twoj twom twomp EulerCoords",twoj,&
!   & twom,twomp,eulerCoords(:)

   ! Compute term1.
   term1Exp = real(twom) / 2.0_double * eulerCoords(1)
   term1 = cmplx(cos(term1Exp), -sin(term1Exp), double)

   ! Compute term2.
   term2 = smalld(twoj, twom, twomp, eulerCoords(2))

   ! Compute term3.
   term3Exp = real(twomp) / 2.0_double * eulerCoords(3)
   term3 = cmplx(cos(term3Exp), -sin(term3Exp), double)

!write (6,fmt="(a,2e15.4)") "WignerD term1", term1
!write (6,fmt="(a,2e15.4)") "WignerD term2", term2
!write (6,fmt="(a,2e15.4)") "WignerD term3", term3

   wignerD = term1 * term2 * term3

end function wignerD


function smalld (twoj, twom, twomp, theta)

   ! Use necessary modules
   use O_Kinds

   implicit none

   ! Define function return variable.
   real(kind=double) :: smalld

   ! Define passed parameters.
   integer :: twoj, twom, twomp
   real(kind=double) :: theta

   ! Define local variables.
   integer :: twok, twokStart, twokEnd, halfk
   real(kind=double) :: dsum
!   real(kind=double) :: coeff
   real(kind=double) :: cos_theta_2
   real(kind=double) :: sin_theta_2
   real(kind=double) :: numerator
   real(kind=double) :: denominator


   ! Following equation (4) in section 4.3.1 of "Quantum Theory of Angular
   !   Momentum: Irreducible Tensors, Spherical Harmonics, Vector Coupling
   !   Coefficients, 3nj Symbols" by Varshalovich DA, Moskalev AN, and
   !   Khersonski VK.; Singapore; Teaneck, NJ, USA: World Scientific Pub;
   !   1988. 514 p.; Equation found at the bottom of page 76.
   

   cos_theta_2 = cos(0.5_double * theta)
   sin_theta_2 = sin(0.5_double * theta)
!write (6,*) "cos sin",cos_theta_2,sin_theta_2

   twokStart = max(0, twom - twomp)
   twokEnd = min(twoj + twom, twoj - twomp)

!write (6,*) "twoj, twom, twomp = ", twoj, twom, twomp
!write (6,*) "twokStart twokEnd", twokStart, twokEnd
   smalld = sqrt( &
         &   real(preCompFactorial((twoj+twom)/2),double) &
         & * real(preCompFactorial((twoj-twom)/2),double) &
         & * real(preCompFactorial((twoj+twomp)/2),double) &
         & * real(preCompFactorial((twoj-twomp)/2),double))
!write (6,*) "smalld", smalld

   dsum = 0.0_double
   do twok = twokStart, twokEnd, 2
      halfk = real(twok) / 2.0_double
      numerator = (-1.0d0)**halfk &
            & * cos_theta_2**(twoj - twok + (twom - twomp)/2) &
            & * sin_theta_2**(twok + (twomp - twom)/2)
!write (6,*) "(twom - twomp)/2",(twom - twomp)/2
!write (6,*) "(twomp - twom)/2",(twomp - twom)/2
!write (6,fmt="(a,i3,3e14.5)") "twok 1 2 3", twok, (-1.0d0)**halfk, &
!      & cos_theta_2**(twoj - twok + (twom - twomp)/2), &
!      & sin_theta_2**(twok + (twomp - twom)/2)

      denominator = real(preCompFactorial(halfk),double) &
            & * real(preCompFactorial((twoj + twom - twok)/2),double) &
            & * real(preCompFactorial((twoj - twomp - twok)/2),double) &
            & * real(preCompFactorial((twomp - twom + twok)/2),double)

!write (6,*) "numerator denom ", numerator, denominator
      dsum = dsum + numerator / denominator
   enddo

   smalld = smalld * dsum

end function smalld


function hypersphericalHarmonic4D(twoj,twom,twomp,ttpCoords)

   ! Use necessary modules
   use O_Kinds

   implicit none

   ! Define function return variable.
   complex(kind=double) :: hypersphericalHarmonic4D

   ! Define passed parameters.
   integer :: twoj, twom, twomp
   real(kind=double), dimension(3) :: ttpCoords ! theta_0, theta, phi

   ! Define local variables.
   integer :: i
   real(kind=double) :: term2Param
   real(kind=double), dimension(3) :: ttpTempCoords
   complex(kind=double) :: term1, term2, term3

   ! Initialize the accumulation variable.
   hypersphericalHarmonic4D = cmplx(0.0, 0.0, double)
   do i = -twoj, twoj, 2
      ttpTempCoords(1) = ttpCoords(3)  !  Phi
      ttpTempCoords(2) = ttpCoords(2)  !  Theta
      ttpTempCoords(3) = -ttpCoords(3) ! -Phi

      ! Compute term 1.
      term1 = wignerD(twoj, twom, i, ttpTempCoords)

      ! Compute term 2.
      term2Param = real(i)/2.0_double * ttpCoords(1)
      term2 = cmplx(cos(term2Param), -sin(term2Param), double)

      ! Compute term 3.
      ttpTempCoords(2) = -ttpCoords(2) ! -Theta
      term3 = wignerD(twoj, i, twomp, ttpTempCoords)

      ! HSH = product(term1,term2,term3).
      hypersphericalHarmonic4D = hypersphericalHarmonic4D &
            & + term1 * term2 * term3

!write (6,*) "i HSH = ", i, hypersphericalHarmonic4D
!write (6,*) "term1 = ", term1
!write (6,*) "term2 = ", term2
!write (6,*) "term3 = ", term3
   enddo
   
end function hypersphericalHarmonic4D


function clebschGordan(twoj1, twoj2, twoj, twom1, twom2, twom)

   ! Use necessary modules
   use O_Kinds

   implicit none

   ! Define function return variable.
   real (kind=double) :: clebschGordan

   ! Define passed paramters.
   integer :: twoj1, twoj2, twoj, twom1, twom2, twom

   ! Define local variables.
   integer :: twoz ! Index in the summation for computing CGC values.
   integer :: twocgcMin ! Min value of twoz.
   integer :: twocgcMax ! Max value of twoz.
   real(kind=double) :: preFactor
   real(kind=double) :: coefficient
   real(kind=double) :: numerator
   real(kind=double) :: denominator
!write (6,fmt="(a,i)") "(twoj1 + twoj2 - twoj)/2", (twoj1 + twoj2 - twoj)/2
!write (6,fmt="(a,i)") "(twoj1 - twoj2 + twoj)/2", (twoj1 - twoj2 + twoj)/2
!write (6,fmt="(a,i)") "(-twoj1 + twoj2 + twoj)/2", (-twoj1 + twoj2 + twoj)/2
!write (6,fmt="(a,i)") "(twoj + twoj1 + twoj2 + 2)/2", (twoj + twoj1 + twoj2 + 2)/2

   ! Following QTAM chapter 8, we see in sub-section 8.1.1 the constraints on
   !   parameters (j, j1, j2, m, m1, m2) to the CGCs, and in section 8.2 we
   !   see the explicit definition of the CGC values. Specifically, we use
   !   equations (1) from section 8.2 of QTAM (page 237) to define the
   !   preFactor (Delta(abc) with a = j, b = j1, c = j2).
   preFactor = sqrt( &
         & real(preCompFactorial((twoj1 + twoj2 - twoj)/2),double) &
         & * real(preCompFactorial((twoj1 - twoj2 + twoj)/2),double) &
         & * real(preCompFactorial((-twoj1 + twoj2 + twoj)/2),double) &
         & / real(preCompFactorial((twoj + twoj1 + twoj2 + 2)/2),double))
!write (6,fmt="(a,e15.4)") "preFactor", preFactor
!
!write (6,fmt="(a,i)") "(twoj + twom)/2", (twoj + twom)/2
!write (6,fmt="(a,i)") "(twoj - twom)/2", (twoj - twom)/2
!write (6,fmt="(a,i)") "(twoj + 1)", twoj + 1
!write (6,fmt="(a,i)") "(twoj1 + twom1)/2", (twoj1 + twom1)/2
!write (6,fmt="(a,i)") "(twoj1 - twom1)/2", (twoj1 - twom1)/2
!write (6,fmt="(a,i)") "(twoj2 + twom2)/2", (twoj2 + twom2)/2
!write (6,fmt="(a,i)") "(twoj2 - twom2)/2", (twoj2 - twom2)/2

   ! Then, we use equation (5) in section 8.2 of QTAM (page 238) to compute
   !   explicit values for the CGCs.
   coefficient = sqrt( &
      & real(preCompFactorial((twoj + twom)/2),double) &
      & * real(preCompFactorial((twoj - twom)/2) * (twoj + 1),double) &
      & / (real(preCompFactorial((twoj1 + twom1)/2),double) &
      & * real(preCompFactorial((twoj1 - twom1)/2),double) &
      & * real(preCompFactorial((twoj2 + twom2)/2),double) &
      & * real(preCompFactorial((twoj2 - twom2)/2),double)))
!write (6,fmt="(a,e15.4)") "coefficient", coefficient

   clebschGordan = 0.0_double
   twocgcMin = max(0, twom1 - twoj1, twoj2 - twoj1 + twom)
   twocgcMax = min(twoj2 + twoj + twom1, twoj - twoj1 + twoj2,&
         & twoj + twom)
!write (6,fmt="(a,2i)") "cgc Min,Max", twocgcMin, twocgcMax
!   cgcMin = (max(0, twom1 - twoj1, twoj2 - twoj1 + twom))/2
!   cgcMax = (min(twoj2 + twoj + twom1, twoj - twoj1 + twoj2,&
!         & twoj + twom))/2
!!   clebschGordan = 0.0_double
!!   cgcMin = max(0, int(m1 - j1), int(j2 - j1 + m))
!!   cgcMax = min(int(j2 + j + m1), int(j - j1 + j2), int(j + m))

   do twoz = twocgcMin, twocgcMax, 2
!write (6,*) "twoz = ", twoz
!write (6,fmt="(a,i)") "(twoj2 + twom2 + i)/2", (twoj2 + twom2 + twoz)

      numerator = real(((-1)**((twoj2 + twom2 + twoz)/2)),double) &
         & * real(preCompFactorial((twoj + twoj2 + twom1 - twoz)/2),double) &
         & * real(preCompFactorial((twoj1 - twom1 + twoz)/2),double)
!      numerator = real(((-1)**(j2 + m2 + i)) &
!            & * preCompFactorial(int(j2 + j + m1 - i)) &
!            & * preCompFactorial(int(j1 - m1 + i)),double)

      denominator = real(preCompFactorial(twoz/2),double) &
         & * real(preCompFactorial((twoj - twoj1 + twoj2 - twoz)/2),double) &
         & * real(preCompFactorial((twoj + twom - twoz)/2),double) &
         & * real(preCompFactorial((twoj1 - twoj2 - twom + twoz)/2),double)
!      denominator = real(preCompFactorial(i) &
!            & * preCompFactorial(int(j - j1 + j2 - i)) &
!            & * preCompFactorial(int(j + m - i)) &
!            & * preCompFactorial(int(j1 - j2 - m + i)),double)

      clebschGordan = clebschGordan + numerator / denominator
   enddo

   clebschGordan = preFactor * coefficient * clebschGordan

end function clebschGordan



! ------------------------------------------------------------
!  Bloechl linear-analytic-tetrahedron corner weights
!  (DESIGN 1.3 / 1.5 / 1.6; PSEUDOCODE 2, 3, 3a).
!
!  These are pure functions of one tetrahedron's four sorted
!    corner eigenvalues and an energy, so they live here rather
!    than in any one consumer.  Three modules need them and the
!    dependency graph forbids sharing them from O_DOS: the DOS
!    path integrates them over energy, the bond/effective-charge
!    path evaluates them once at the Fermi level, and the SCF
!    occupation path (O_Populate) evaluates them at each trial
!    energy of its Fermi search -- and O_DOS already uses
!    O_Populate, so the reverse import would be a module cycle.
!    Duplicating them instead would leave two copies of the
!    Bloechl formulas to drift apart.
! ------------------------------------------------------------
! Compute the four Bloechl corner integration weights for an arbitrary energy E
!   and sorted eigenvalues e1 <= e2 <= e3 <= e4. These weights determine how
!   much of each corner's partial-DOS projection to include at energy E. They
!   are the LAT replacement for Gaussian broadening in the energy-resolved PDOS
!   (DESIGN 1.4, PSEUDOCODE 3a).
!
!   The formulas follow from the vertex-averaging property of linear functions
!   over tetrahedra: the integral of barycentric coordinate lambda_i over any
!   sub-tetrahedron equals V_sub/4 times the sum of lambda_i at the sub-tet's
!   four vertices.
!
!   Four cases arise depending on where E falls relative to the sorted corner
!   eigenvalues:
!     Case 0: E < e1 or E >= e4 (trivial bounds) Case 1: e1 <= E < e2 (small
!     sub-tet near e1) Case 2: e2 <= E < e3 (pentahedral middle) Case 3: e3 <= E
!     < e4 (complement of Case 1)
!
!   Same physics as populate.F90's corner weight computation, but parameterized
!   on an arbitrary energy E rather than E_Fermi = 0.
subroutine bloechlCornerWeights(energy, sortedEps, &
      & cornerWt)

   use O_Kinds

   implicit none

   ! Passed parameters.
   real (kind=double), intent(in) :: energy
   real (kind=double), dimension(4), intent(in) :: &
         & sortedEps
   real (kind=double), dimension(4), intent(out) :: &
         & cornerWt

   ! Local variables.
   real (kind=double) :: e1, e2, e3, e4
   real (kind=double) :: e31, e41, e32, e42
   real (kind=double) :: t2, t3, t4, f
   real (kind=double) :: a, b, c, d
   real (kind=double) :: v_I, v_II, v_III
   real (kind=double) :: s1, s2, s3, f_un
   real (kind=double) :: denom
   real (kind=double), parameter :: tol = 1.0d-12

   e1 = sortedEps(1)
   e2 = sortedEps(2)
   e3 = sortedEps(3)
   e4 = sortedEps(4)

   ! Case 0a: energy below all corners. No spectral weight from this tetrahedron
   !   at this energy.
   if (energy < e1) then
      cornerWt(:) = 0.0_double
      return
   endif

   ! Case 0b: energy above all corners. By vertex averaging over the full
   !   tetrahedron, each corner gets an equal share of 1/4.
   if (energy >= e4) then
      cornerWt(:) = 0.25_double
      return
   endif

   ! Case 1: e1 <= E < e2. The iso-energy surface cuts the three edges from
   !   corner 1, forming a small sub-tetrahedron with apex at corner 1.
   if (energy < e2) then
      denom = (e2-e1) * (e3-e1) * (e4-e1)
      if (abs(denom) < tol) then
         cornerWt(:) = 0.0_double
         return
      endif
      t2 = (energy - e1) / (e2 - e1)
      t3 = (energy - e1) / (e3 - e1)
      t4 = (energy - e1) / (e4 - e1)
      f = t2 * t3 * t4
      cornerWt(2) = f * t2 / 4.0_double
      cornerWt(3) = f * t3 / 4.0_double
      cornerWt(4) = f * t4 / 4.0_double
      cornerWt(1) = f - cornerWt(2) &
            & - cornerWt(3) - cornerWt(4)
      return
   endif

   ! Case 2: e2 <= E < e3 (middle range). Corners 1 and 2 lie below; corners 3
   !   and 4 lie above. The occupied region is a pentahedron decomposed into
   !   three sub-tetrahedra (T_I, T_II, T_III).
   if (energy < e3) then
      e31 = e3 - e1
      e41 = e4 - e1
      e32 = e3 - e2
      e42 = e4 - e2
      if (e31 * e41 < tol .or. &
            & e32 * e42 < tol) then
         cornerWt(:) = 0.0_double
         return
      endif

      ! Intersection parameters: fractional positions where the iso-energy
      !   surface cuts each edge.
      a = (energy - e1) / e31
      b = (energy - e1) / e41
      c = (energy - e2) / e32
      d = (energy - e2) / e42

      ! Sub-tetrahedra volume ratios (as fractions of the full tetrahedron
      !   volume).
      v_I   = a * b
      v_II  = a * d * (1.0_double - b)
      v_III = (1.0_double - a) * c * d

      ! Corner weights from vertex averaging over the three sub-tetrahedra.
      cornerWt(1) = ( &
            & v_I * (3.0_double - a - b) &
            & + v_II * (2.0_double - a - b) &
            & + v_III * (1.0_double - a) &
            & ) / 4.0_double
      cornerWt(2) = ( &
            & v_I &
            & + v_II * (2.0_double - d) &
            & + v_III * (3.0_double - c - d) &
            & ) / 4.0_double
      cornerWt(3) = ( &
            & v_I * a + v_II * a &
            & + v_III * (a + c) &
            & ) / 4.0_double
      cornerWt(4) = ( &
            & v_I * b &
            & + v_II * (b + d) &
            & + v_III * d &
            & ) / 4.0_double
      return
   endif

   ! Case 3: e3 <= E < e4. Only corner 4 lies above the energy. The unoccupied
   !   region is a small sub-tet near corner 4 (complement of Case 1).
   denom = (e4-e1) * (e4-e2) * (e4-e3)
   if (abs(denom) < tol) then
      cornerWt(:) = 0.25_double
      return
   endif
   s1 = (e4 - energy) / (e4 - e1)
   s2 = (e4 - energy) / (e4 - e2)
   s3 = (e4 - energy) / (e4 - e3)
   f_un = s1 * s2 * s3
   cornerWt(1) = 0.25_double &
         & - f_un * s1 / 4.0_double
   cornerWt(2) = 0.25_double &
         & - f_un * s2 / 4.0_double
   cornerWt(3) = 0.25_double &
         & - f_un * s3 / 4.0_double
   cornerWt(4) = (1.0_double - f_un) &
         & - cornerWt(1) - cornerWt(2) &
         & - cornerWt(3)

end subroutine bloechlCornerWeights


! Compute the four per-corner DOS density weights (cornerDOSWt_LAT) for one
!   tetrahedron at a given energy. These are the energy derivatives of the
!   cumulative corner integration weights returned by bloechlCornerWeights:
!
!     cornerDOSWt_LAT(c) = d/dE [cornerIntgWt(c)]
!
!   Units: 1/energy (same as eigenvalue units). Their sum equals the total
!   per-tetrahedron DOS (the dosContrib from the original TDOS code).
!
!   Used by both the TDOS (sum only) and the PDOS (per-corner, to weight
!   Mulliken projections). See DESIGN 1.3 and PSEUDOCODE 2a.
subroutine bloechlCornerDOSWt(energy, sortedEps, &
      & cornerDOSWt)

   use O_Kinds

   implicit none

   ! Passed parameters.
   real (kind=double), intent(in) :: energy
   real (kind=double), dimension(4), intent(in) :: &
         & sortedEps
   real (kind=double), dimension(4), intent(out) :: &
         & cornerDOSWt

   ! Local variables.
   real (kind=double) :: e1, e2, e3, e4
   real (kind=double) :: e31, e41, e32, e42
   real (kind=double) :: t2, t3, t4, f, gTotal
   real (kind=double) :: a, b, cv, dv
   real (kind=double) :: da, db, dc, dd
   real (kind=double) :: v_I, v_II, v_III
   real (kind=double) :: dv_I, dv_II, dv_III
   real (kind=double) :: s1, s2, s3, f_un
   real (kind=double) :: denom
   real (kind=double), parameter :: tol = 1.0d-12

   e1 = sortedEps(1)
   e2 = sortedEps(2)
   e3 = sortedEps(3)
   e4 = sortedEps(4)

   ! Case 0: energy outside eigenvalue range. No spectral density from this
   !   tetrahedron.
   if (energy < e1 .or. energy >= e4) then
      cornerDOSWt(:) = 0.0_double
      return
   endif

   ! ------------------------------------------------- Case 1: e1 <= E < e2
   ! ------------------------------------------------- The cumulative weights
   ! are w(j) = f*t_j/4
   !   for j=2,3,4, and w(1) = f - w(2) - w(3) - w(4), where f = t2*t3*t4 and
   !   t_j = (E-e1) / (e_j - e1). Apply the product rule:
   !     d(f*t_j)/dE = gTotal*t_j + f/(e_j - e1)
   !   where gTotal = df/dE = 3*(E-e1)^2 / denom.
   if (energy < e2) then
      denom = (e2-e1) * (e3-e1) * (e4-e1)
      if (abs(denom) < tol) then
         cornerDOSWt(:) = 0.0_double
         return
      endif
      t2 = (energy - e1) / (e2 - e1)
      t3 = (energy - e1) / (e3 - e1)
      t4 = (energy - e1) / (e4 - e1)
      f = t2 * t3 * t4
      gTotal = 3.0_double &
            & * (energy - e1)**2 / denom

      cornerDOSWt(2) = (gTotal * t2 &
            & + f / (e2 - e1)) / 4.0_double
      cornerDOSWt(3) = (gTotal * t3 &
            & + f / (e3 - e1)) / 4.0_double
      cornerDOSWt(4) = (gTotal * t4 &
            & + f / (e4 - e1)) / 4.0_double
      cornerDOSWt(1) = gTotal &
            & - cornerDOSWt(2) &
            & - cornerDOSWt(3) &
            & - cornerDOSWt(4)
      return
   endif

   ! ------------------------------------------------- Case 2: e2 <= E < e3
   ! (middle range) ------------------------------------------------- The
   ! cumulative weights use three sub-tetrahedra
   !   volumes v_I, v_II, v_III with intersection parameters a, b, c, d. We
   !   compute volume derivatives dv_I, dv_II, dv_III and parameter derivatives
   !   da, db, dc, dd, then apply the product rule to each corner weight
   !   expression.
   if (energy < e3) then
      e31 = e3 - e1
      e41 = e4 - e1
      e32 = e3 - e2
      e42 = e4 - e2
      if (e31 * e41 < tol .or. &
            & e32 * e42 < tol) then
         cornerDOSWt(:) = 0.0_double
         return
      endif

      ! Intersection parameters.
      a  = (energy - e1) / e31
      b  = (energy - e1) / e41
      cv = (energy - e2) / e32
      dv = (energy - e2) / e42

      ! Sub-tetrahedra volume ratios.
      v_I   = a * b
      v_II  = a * dv * (1.0_double - b)
      v_III = (1.0_double - a) * cv * dv

      ! Parameter derivatives (d/dE).
      da = 1.0_double / e31
      db = 1.0_double / e41
      dc = 1.0_double / e32
      dd = 1.0_double / e42

      ! Volume derivatives (d/dE).
      dv_I   = b / e31 + a / e41
      dv_II  = dv * (1.0_double - b) / e31 &
            & + a * (1.0_double - b) / e42 &
            & - a * dv / e41
      dv_III = -cv * dv / e31 &
            & + (1.0_double - a) * dv / e32 &
            & + (1.0_double - a) * cv / e42

      ! Corner 1: w(1) = [v_I*(3-a-b) + v_II*(2-a-b) + v_III*(1-a)] / 4
      cornerDOSWt(1) = ( &
            & dv_I * (3.0_double - a - b) &
            & + v_I * (-da - db) &
            & + dv_II * (2.0_double - a - b) &
            & + v_II * (-da - db) &
            & + dv_III * (1.0_double - a) &
            & + v_III * (-da) &
            & ) / 4.0_double

      ! Corner 2: w(2) = [v_I + v_II*(2-d) + v_III*(3-c-d)] / 4
      cornerDOSWt(2) = ( &
            & dv_I &
            & + dv_II * (2.0_double - dv) &
            & + v_II * (-dd) &
            & + dv_III &
            & * (3.0_double - cv - dv) &
            & + v_III * (-dc - dd) &
            & ) / 4.0_double

      ! Corner 3: w(3) = [v_I*a + v_II*a + v_III*(a+c)] / 4
      cornerDOSWt(3) = ( &
            & dv_I * a + v_I * da &
            & + dv_II * a + v_II * da &
            & + dv_III * (a + cv) &
            & + v_III * (da + dc) &
            & ) / 4.0_double

      ! Corner 4: w(4) = [v_I*b + v_II*(b+d) + v_III*d] / 4
      cornerDOSWt(4) = ( &
            & dv_I * b + v_I * db &
            & + dv_II * (b + dv) &
            & + v_II * (db + dd) &
            & + dv_III * dv &
            & + v_III * dd &
            & ) / 4.0_double
      return
   endif

   ! ------------------------------------------------- Case 3: e3 <= E < e4
   ! ------------------------------------------------- The unoccupied region is
   ! a small sub-tet near
   !   corner 4 with fraction f_un = s1*s2*s3 where s_j = (e4-E)/(e4-e_j).
   !   Derivatives:
   !     ds_j/dE = -1/(e4-e_j) df_un/dE = -gTotal
   !   For j=1,2,3: dw(j)/dE = (gTotal*s_j
   !       + f_un/(e4-e_j)) / 4
   denom = (e4-e1) * (e4-e2) * (e4-e3)
   if (abs(denom) < tol) then
      cornerDOSWt(:) = 0.0_double
      return
   endif
   s1 = (e4 - energy) / (e4 - e1)
   s2 = (e4 - energy) / (e4 - e2)
   s3 = (e4 - energy) / (e4 - e3)
   f_un = s1 * s2 * s3
   gTotal = 3.0_double &
         & * (e4 - energy)**2 / denom

   cornerDOSWt(1) = (gTotal * s1 &
         & + f_un / (e4 - e1)) / 4.0_double
   cornerDOSWt(2) = (gTotal * s2 &
         & + f_un / (e4 - e2)) / 4.0_double
   cornerDOSWt(3) = (gTotal * s3 &
         & + f_un / (e4 - e3)) / 4.0_double
   cornerDOSWt(4) = gTotal &
         & - cornerDOSWt(1) &
         & - cornerDOSWt(2) &
         & - cornerDOSWt(3)

end subroutine bloechlCornerDOSWt


! Compute the four Bloechl curvature correction terms for one tetrahedron at a
!   given energy (DESIGN 1.3.1; PSEUDOCODE 3a).
!
!   Linear interpolation of the band energy inside a tetrahedron misplaces the
!   iso-energy surface, and the resulting error falls off only slowly as the
!   k-point mesh is made denser. Bloechl's correction (PRB 49, 16223 (1994),
!   equation 22) compensates for it by shifting integration weight between the
!   four corners:
!
!!     dw_i = (1/40) * D_T(E_F) * sum_{j=1..4} (eps_j - eps_i)
!!          = (1/10) * D_T(E_F) * (epsBar - eps_i)
!
!   where D_T(E_F) is THIS tetrahedron's density of states at the energy in
!   question and epsBar is the mean of its four corner eigenvalues. The two
!   forms are the same quantity, since sum_j (eps_j - eps_i) = 4*(epsBar -
!   eps_i); the second is implemented because it makes the central property
!   visible rather than leaving it to be derived.
!
!   That property is that the four terms SUM TO ZERO, identically, for any
!   energy and any corner energies. The correction therefore moves weight
!   between corners and never changes a tetrahedron's total, which is what
!   makes it safe to add to a converged SCF path: the electron count is
!   unchanged, so the Fermi search is unchanged, so the calibration of the
!   total occupation against numElectrons cannot be broken by adding this.
!   What does change is WHICH k-points hold the occupation, and therefore the
!   charge density and the band energy.
!
!   The correction also vanishes wherever no iso-energy surface passes through
!   the tetrahedron, since D_T is then zero. Insulators are untouched and only
!   straddling tetrahedra move.
!
!   The caller adds these terms to the weights from bloechlCornerWeights rather
!   than receiving them already folded in. That keeps that routine a literal
!   transcription of the paper's uncorrected expressions, so a reader can check
!   it against the reference without mentally subtracting a correction, and it
!   gives this routine a self-contained test: sum(cornerCorrWt) must be zero to
!   rounding, which needs no reference values and no other quantity.
!
!   Equation 22 is the whole correction. The reference list in DESIGN once
!   cited "equations 22-24" together, which reads as though two further terms
!   were owed here. They are not: 23 and 24 compare the true Fermi surface
!   against the linearly interpolated polyhedral one, an assessment the paper
!   makes rather than a formula anything computes. This routine is complete as
!   written, and a reader should not go hunting for the rest of it.
subroutine bloechlCornerCorrection(energy, sortedEps, &
      & cornerCorrWt)

   use O_Kinds

   implicit none

   ! Passed parameters.
   real (kind=double), intent(in) :: energy
   real (kind=double), dimension(4), intent(in) :: &
         & sortedEps
   real (kind=double), dimension(4), intent(out) :: &
         & cornerCorrWt

   ! Local variables.
   real (kind=double), dimension(4) :: cornerDOSWt
   real (kind=double) :: tetraDOSAtEnergy, meanEps
   integer :: corner

   ! No iso-energy surface passes through a tetrahedron lying wholly below or
   !   wholly above the energy, so its DOS there is zero and so is the
   !   correction. Both bounds are tested explicitly rather than left to the
   !   arithmetic below: they are the common case (every fully occupied and
   !   every empty tetrahedron in the mesh), and this is the test that keeps
   !   gapped systems untouched.
   if ((energy < sortedEps(1)) .or. &
         & (energy >= sortedEps(4))) then
      cornerCorrWt(:) = 0.0_double
      return
   endif

   ! This tetrahedron's density of states at the energy is the sum of its four
   !   corner DOS weights, taken before any tetrahedron-volume factor. Reusing
   !   that routine rather than writing the per-case DOS expressions again is
   !   what keeps the case logic and the degenerate-corner guards in one place.
   call bloechlCornerDOSWt (energy, sortedEps, cornerDOSWt)
   tetraDOSAtEnergy = sum(cornerDOSWt(:))

   meanEps = sum(sortedEps(:)) / 4.0_double

   do corner = 1, 4
      cornerCorrWt(corner) = tetraDOSAtEnergy &
            & * (meanEps - sortedEps(corner)) / 10.0_double
   enddo

end subroutine bloechlCornerCorrection


! Average an atom-resolved spectrum over the crystal point group
!   (DESIGN 1.7; PSEUDOCODE 20).
!
! Why this is needed at all. The tetrahedron decomposition of DESIGN 1.2 is
!   carried onto itself by the point group for lattices whose operations
!   permute the mesh axes up to sign -- cubic, tetragonal, orthorhombic --
!   and not for hexagonal or rhombohedral ones, where a six-fold rotation
!   sends one mesh axis onto the SUM of two and so does not map a grid box
!   onto a grid box at all. Wherever the decomposition is not invariant,
!   k-points related by symmetry receive unrelated integration weights, and a
!   quantity resolved onto individual atoms inherits the difference: atoms
!   that must be equivalent come out unequal. Averaging the finished result
!   over the group removes exactly that part of the error.
!
! Why it is not a cosmetic patch. Averaging an atom-resolved result over the
!   orbit its atom belongs to is EXACTLY equal to replacing each integration
!   weight by its average over the star of its k-point. It is therefore the
!   projection of the quadrature onto the symmetric subspace rather than an
!   adjustment applied to a finished number, and two useful things follow.
!   Totals do not move at all, because summing the averaged weights over
!   k-points returns the original sum. And it is exact on every lattice,
!   including the ones the decomposition alone cannot fix, because it
!   averages over the star that actually exists.
!
! No orbit is enumerated here, and none needs to be. Summing over EVERY
!   operation of the group is the orbit average, because the operations
!   carrying a channel onto a given orbit member form a coset and therefore
!   contribute that member equally often. For the same reason the direction
!   of the permutation does not matter: a table built from the forward atom
!   map and one built from its inverse give the same average, since the
!   operation index runs over the whole group either way.
!
! The channel index is whatever the caller resolves by -- one channel per
!   atom, or per atom and radial function. The permutation table says which
!   channel each one becomes under each operation, and the caller owns it.
! The routine also reports how much it changed. An imposed equality that
!   leaves no trace cannot be told apart from an earned one, so the caller is
!   handed the largest disagreement found within any group of channels the
!   averaging merged, along with the largest value in the spectrum to measure
!   it against. Measuring it here rather than in a separate pass keeps the
!   reported number describing the averaging that actually happened.
subroutine symmetrizeChannels(values, permTable, numOps, &
      & numChannels, numEnergyPoints, largestSpread, &
      & largestValue)

   use O_Kinds

   implicit none

   ! Define passed parameters.
   real (kind=double), dimension(:,:), intent(inout) :: values
         !   (numChannels, numEnergyPoints). Modified in place.
   integer, dimension(:,:), intent(in) :: permTable
         !   (numOps, numChannels). permTable(R, alpha) is the channel that
         !   alpha becomes under operation R.
   integer, intent(in) :: numOps
   integer, intent(in) :: numChannels
   integer, intent(in) :: numEnergyPoints
   real (kind=double), intent(out) :: largestSpread
         !   The largest gap between channels that had to be made equal.
   real (kind=double), intent(out) :: largestValue
         !   The largest magnitude anywhere in the spectrum, so the caller can
         !   report the spread as a fraction of something meaningful.

   ! Define local variables.
   integer :: energyPoint  ! Energy grid loop index.
   integer :: opIndex      ! Point group operation loop index.
   integer :: channel      ! Channel loop index.
   real (kind=double) :: thisValue ! The value being folded in.
   real (kind=double), allocatable, dimension(:) :: averaged
         !   One energy point's averaged channels, built separately because
         !   the average must be formed from the UNSYMMETRIZED values: writing
         !   into `values` as we go would feed already-averaged numbers back
         !   into later operations. One energy point at a time keeps this to a
         !   single channel vector rather than a copy of the whole spectrum.
   real (kind=double), allocatable, dimension(:) :: smallestSeen
   real (kind=double), allocatable, dimension(:) :: largestSeen
         !   The extremes within each channel's group, tracked alongside the
         !   sum so that the disagreement being averaged away is measured in
         !   the same pass that averages it.

   allocate (averaged(numChannels))
   allocate (smallestSeen(numChannels))
   allocate (largestSeen(numChannels))

   largestSpread = 0.0_double
   largestValue = 0.0_double

   do energyPoint = 1, numEnergyPoints
      averaged(:) = 0.0_double
      smallestSeen(:) = huge(1.0_double)
      largestSeen(:) = -huge(1.0_double)

      do opIndex = 1, numOps
         do channel = 1, numChannels
            thisValue = values(permTable(opIndex, channel), &
                  & energyPoint)
            averaged(channel) = averaged(channel) + thisValue
            if (thisValue < smallestSeen(channel)) then
               smallestSeen(channel) = thisValue
            endif
            if (thisValue > largestSeen(channel)) then
               largestSeen(channel) = thisValue
            endif
         enddo
      enddo

      largestSpread = max(largestSpread, &
            & maxval(largestSeen(:) - smallestSeen(:)))
      largestValue = max(largestValue, &
            & maxval(abs(values(:, energyPoint))))

      values(:, energyPoint) = averaged(:) &
            & / real(numOps, double)
   enddo

   deallocate (averaged)
   deallocate (smallestSeen)
   deallocate (largestSeen)

end subroutine symmetrizeChannels


! Average a PAIR-resolved spectrum over the crystal point group
!   (DESIGN 1.7; PSEUDOCODE 20).
!
! The partial optical properties resolve each transition by the partial the
!   initial state belongs to AND the partial the final state belongs to, so
!   the quantity carries two channel indices rather than one. Both are
!   permuted by the same operation, which is what preserves the meaning of a
!   pair: the partials at one end of a transition and at the other are
!   carried to their images together, exactly as the atoms they sit on are.
!
! Everything said about symmetrizeChannels above applies unchanged --
!   why the averaging is needed, why it is equivalent to symmetrizing the
!   integration weights, and why no orbit has to be enumerated.
subroutine symmetrizePairs(values, permTable, numOps, &
      & numPartials, numComponents, numEnergyPoints, &
      & largestSpread, largestValue)

   use O_Kinds

   implicit none

   ! Define passed parameters.
   real (kind=double), dimension(:,:,:,:), intent(inout) :: values
         !   (numPartials, numPartials, numComponents, numEnergyPoints).
         !   Modified in place.
   integer, dimension(:,:), intent(in) :: permTable
         !   (numOps, numPartials).
   integer, intent(in) :: numOps
   integer, intent(in) :: numPartials
   integer, intent(in) :: numComponents
   integer, intent(in) :: numEnergyPoints
   real (kind=double), intent(out) :: largestSpread
         !   The largest gap between pairs that had to be made equal.
   real (kind=double), intent(out) :: largestValue
         !   The largest magnitude anywhere in the spectra.

   ! Define local variables.
   integer :: energyPoint  ! Energy grid loop index.
   integer :: opIndex      ! Point group operation loop index.
   integer :: initPartial  ! Initial-state partial loop index.
   integer :: finPartial   ! Final-state partial loop index.
   integer :: initRotated  ! Image of initPartial under this operation.
   integer :: finRotated   ! Image of finPartial under this operation.
   real (kind=double), allocatable, dimension(:,:,:) :: averaged
         !   One energy point's averaged pair matrix, built separately for the
         !   same reason as in symmetrizeChannels. One energy point at a time
         !   matters more here: a copy of the whole array would be the pair
         !   count squared times the component count times the energy grid.
   real (kind=double), allocatable, dimension(:,:,:) :: smallestSeen
   real (kind=double), allocatable, dimension(:,:,:) :: largestSeen
         !   The extremes within each pair's group, tracked in the same pass
         !   that averages them, so the reported disagreement describes the
         !   averaging that actually happened.

   allocate (averaged(numPartials, numPartials, numComponents))
   allocate (smallestSeen(numPartials, numPartials, &
         & numComponents))
   allocate (largestSeen(numPartials, numPartials, &
         & numComponents))

   largestSpread = 0.0_double
   largestValue = 0.0_double

   do energyPoint = 1, numEnergyPoints
      averaged(:,:,:) = 0.0_double
      smallestSeen(:,:,:) = huge(1.0_double)
      largestSeen(:,:,:) = -huge(1.0_double)

      do opIndex = 1, numOps

         ! Final-state partial outer, initial-state partial inner. The
         !   initial-state partial is the leftmost index of both arrays and
         !   so is the one that must run innermost: Fortran stores the
         !   leftmost index fastest, and this keeps both the read and the
         !   write walking memory in order.
         do finPartial = 1, numPartials
            finRotated = permTable(opIndex, finPartial)
            do initPartial = 1, numPartials
               initRotated = permTable(opIndex, initPartial)

               averaged(initPartial, finPartial, :) = &
                     & averaged(initPartial, finPartial, :) &
                     & + values(initRotated, finRotated, :, &
                     & energyPoint)
               smallestSeen(initPartial, finPartial, :) = &
                     & min(smallestSeen(initPartial, &
                     & finPartial, :), &
                     & values(initRotated, finRotated, :, &
                     & energyPoint))
               largestSeen(initPartial, finPartial, :) = &
                     & max(largestSeen(initPartial, &
                     & finPartial, :), &
                     & values(initRotated, finRotated, :, &
                     & energyPoint))
            enddo
         enddo
      enddo

      largestSpread = max(largestSpread, &
            & maxval(largestSeen(:,:,:) &
            & - smallestSeen(:,:,:)))
      largestValue = max(largestValue, &
            & maxval(abs(values(:,:,:,energyPoint))))

      values(:,:,:,energyPoint) = averaged(:,:,:) &
            & / real(numOps, double)
   enddo

   deallocate (averaged)
   deallocate (smallestSeen)
   deallocate (largestSeen)

end subroutine symmetrizePairs

end module O_MathSubs
