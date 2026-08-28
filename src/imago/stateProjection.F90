!! SPDX-License-Identifier: ECL-2.0
!! Copyright (c) 2026 Paul Rulis

!-----------------------------------------------------------------------
! The state projection as a matrix product (PSEUDOCODE 34).
!
! The post-SCF analysis tools -- the density of states, bond order,
!   and the optical properties -- all begin with the same operation:
!   carry every electron state through one valence-by-valence matrix.
!   For the density of states and for bond order that matrix is the
!   overlap S; for the optical properties it is a momentum matrix.
!   The carried states are
!
!!    T(mu,j) = sum over nu of  G(mu,nu) * C(nu,j) ,
!
!   the j-th eigenvector C(:,j) run through the matrix G, formed for
!   every state j at once.  That is a single dense matrix--matrix
!   product, which the level-3 BLAS library evaluates at the speed the
!   hardware allows.  It replaces a hand-rolled loop that, for each
!   pairing of a state with a basis function, walked a whole column of
!   G and summed products -- an arithmetic that moved through memory
!   instead of staying in the processor and so was slow far out of
!   proportion to its operation count (PF9 measured it at thirty-nine
!   minutes on the 1296-atom cell, second only to the secular solve).
!
! Each caller finishes the projection its own way once T is in hand.
!   The density of states and bond order form the Mulliken product
!   Re(conjg(C) .* T) element by element; the optical properties form
!   a different combination of C and T.  So this routine returns T and
!   stops: the finishing step is one line at each caller and is not
!   worth a second shared routine (PSEUDOCODE 34.4).  The routine is
!   written ONCE, here, so that bond order and optics call the same
!   product when their own recasts are taken (TODO PA6, step 5).
!
! This is the serial, single-core recast of DESIGN 9.10's Step 0,
!   taken before any of the work is spread across processors, exactly
!   as PSEUDOCODE 31 recast the density-matrix build before dealing
!   it.  It changes no result by more than the last few digits (the
!   sum is reassociated by the library; PSEUDOCODE 34.7), and it is
!   the serial number the later parallel deal is measured against.
!-----------------------------------------------------------------------
module O_StateProjection

   ! Import the precision kinds shared across the whole code base.
   use O_Kinds

   ! Make sure that no variables are declared accidentally.
   implicit none

   contains

subroutine projectStatesOntoBasis(basisMatrix, eigenVectors, &
      & numBasis, numStates, statesProjected)

   ! Import the precision kinds and the level-3 BLAS interface for
   !   this build.  The overlap (and the bond-order matrix) is
   !   Hermitian in the multi-k complex build and symmetric in the
   !   gamma-point real build, so the specialized Hermitian/symmetric
   !   product zhemm/dsymm is the natural call: it reads only one
   !   triangle of the matrix and so moves half the memory the general
   !   product would, while giving the identical result because the
   !   matrix truly is Hermitian/symmetric.
   use O_Kinds
#ifndef GAMMA
   use zhemmInterface
#else
   use dsymmInterface
#endif

   ! Make sure that no variables are declared accidentally.
   implicit none

   ! Define the passed parameters.
   !   numBasis        the valence dimension: the number of basis
   !                   functions, and the order of the square matrix.
   !   numStates       the number of electron states (columns of C).
   !   basisMatrix     the numBasis-by-numBasis Hermitian (complex) or
   !                   symmetric (real) matrix G the states are carried
   !                   through -- the overlap for the density of states
   !                   and bond order, a momentum matrix for optics.
   !   eigenVectors    the numBasis-by-numStates eigenvector matrix C
   !                   for one k-point and one spin.
   !   statesProjected the numBasis-by-numStates result T = G C.  It is
   !                   the caller's preallocated workspace and is
   !                   written in full (BETA is zero, so its incoming
   !                   contents are not read).
   integer, intent(in) :: numBasis
   integer, intent(in) :: numStates
#ifndef GAMMA
   complex (kind=double), dimension(numBasis,numBasis), &
         & intent(in) :: basisMatrix
   complex (kind=double), dimension(numBasis,numStates), &
         & intent(in) :: eigenVectors
   complex (kind=double), dimension(numBasis,numStates), &
         & intent(out) :: statesProjected
#else
   real (kind=double), dimension(numBasis,numBasis), &
         & intent(in) :: basisMatrix
   real (kind=double), dimension(numBasis,numStates), &
         & intent(in) :: eigenVectors
   real (kind=double), dimension(numBasis,numStates), &
         & intent(out) :: statesProjected
#endif

   ! Form T = G C in one level-3 BLAS product.  The arguments, in
   !   order, are: multiply from the left ('L', so G is the left
   !   factor and is therefore numBasis by numBasis), read the upper
   !   triangle ('U') of the Hermitian/symmetric G, the result has
   !   numBasis rows and numStates columns, the scalar 1 in front of
   !   the product G C, then G and C, the scalar 0 in front of the
   !   incoming (unused) result, and finally the result T.  Every
   !   leading dimension is numBasis because all three matrices are
   !   stored with that many rows.
#ifndef GAMMA
   call zhemm('L', 'U', numBasis, numStates, &
         & (1.0_double, 0.0_double), basisMatrix, numBasis, &
         & eigenVectors, numBasis, (0.0_double, 0.0_double), &
         & statesProjected, numBasis)
#else
   call dsymm('L', 'U', numBasis, numStates, 1.0_double, &
         & basisMatrix, numBasis, eigenVectors, numBasis, &
         & 0.0_double, statesProjected, numBasis)
#endif

end subroutine projectStatesOntoBasis

end module O_StateProjection
