!***********************************************************************
!
!    Copyright (c) 2016, Lawrence Livermore National Security, LLC.
!                        Produced at the Lawrence Livermore National
!                        Laboratory.
!                        Written by Nicolas Schunck, schunck1@llnl.gov
!
!    LLNL-CODE-728299 All rights reserved.
!    LLNL-CODE-573953 All rights reserved.
!
!    Copyright 2017, R. Navarro Perez, N. Schunck, R. Lasseri, C. Zhang,
!                    J. Sarich
!    Copyright 2012, M.V. Stoitsov, N. Schunck, M. Kortelainen, H.A. Nam,
!                    N. Michel, J. Sarich, S. Wild
!    Copyright 2005, M.V. Stoitsov, J. Dobaczewski, W. Nazarewicz, P.Ring
!
!    This file is part of HFBTHO.
!
!    HFBTHO is free software: you can redistribute it and/or modify it
!    under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    HFBTHO is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with HFBTHO. If not, see <http://www.gnu.org/licenses/>.
!
!    OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC
!    LICENSE
!
!    Our Preamble Notice
!
!      A. This notice is required to be provided under our contract
!         with the U.S. Department of Energy (DOE). This work was
!         produced at the Lawrence Livermore National Laboratory under
!         Contract No. DE-AC52-07NA27344 with the DOE.
!      B. Neither the United States Government nor Lawrence Livermore
!         National Security, LLC nor any of their employees, makes any
!         warranty, express or implied, or assumes any liability or
!         responsibility for the accuracy, completeness, or usefulness
!         of any information, apparatus, product, or process disclosed,
!         or represents that its use would not infringe privately-owned
!         rights.
!      C. Also, reference herein to any specific commercial products,
!         process, or services by trade name, trademark, manufacturer
!         or otherwise does not necessarily constitute or imply its
!         endorsement, recommendation, or favoring by the United States
!         Government or Lawrence Livermore National Security, LLC. The
!         views and opinions of authors expressed herein do not
!         necessarily state or reflect those of the United States
!         Government or Lawrence Livermore National Security, LLC, and
!         shall not be used for advertising or product endorsement
!         purposes.
!
!    The precise terms and conditions for copying, distribution and
!    modification are contained in the file COPYING.
!
!***********************************************************************

! ==================================================================== !
!                                                                      !
!                      GAUSS QUADRATURE PACKAGE                        !
!                                                                      !
! ==================================================================== !

!-------------------------------------------------------------------
!> This module provides a set of routines to compute the weights and
!> nodes of several standard quadrature schemes. It also contains
!> several routines specifically related to Hermite polynomials.
!-------------------------------------------------------------------
!>  Subroutines: - gausspoints
!>               - Gaussq(kindi,n,alpha,beta,kpts,endpts,b,t,w)
!>               - Class(kindi,N,ALPHA,BETA,B,A,MUZERO)
!>               - GBTQL2(N,D,E,Z,IERR)
!>               - D_HERM(X,N,HER,DHER,NDIM)
!>               - DEVHER(NOSACT)
!>  Functions: - GBSLVE(SHIFT,N,A,B)
!>             - pr_gamma(x)
!----------------------------------------------------------------------!

Module HFBTHO_gauss

  Use HFBTHO_utilities
  Use HFBTHO

  Implicit None

Contains
  !=======================================================================
  !> The routine determines the points and weights for Gauss quadratures
  !> in the cases of Gauss-Legendre, -Laguerre and -Hermite formulas.
  !> Note that for Gauss-Hermite and Gauss-Laguerre integrations, the
  !> weights are defined in such a way that they contain the weight function
  !> associated with the orthohonal polynomials. For example, in the case
  !> of Gauss-Hermite quadrature, we may write
  !>
  !>  \f{align}{
  !>      \int dx f(x)g(x) & = \int dx e^{-x^2} f(x)e^{x^2} g(x)e^{x^2} \\
  !>                       & = \sum_{i} w_i f(x_i)e^{x_i^2} g(x_i)e^{x_i^2} \\
  !>                       & = \sum_{i} w_i e^{x_i^2} f(x_i)g(x_i)
  !>  \f}
  !>
  !> The weights are thus multiplied by \f$ e^{x_i^2} \f$.
  !=======================================================================
  Subroutine gausspoints
    Implicit None
    Real(pr):: al,be,sparity
    Real(pr), Allocatable :: endpts(:),b(:),t(:),w(:)
    Integer(ipr) :: N,i,j,KINDI,kpts,nparity
    !
    al=0.0_pr; be=0.0_pr; kpts=0
    !
    !--------------------------------------------------------------------
    !------------------>> Gauss-Hermite (positive nodes) <<--------------
    !--------------------------------------------------------------------
    If(Parity) Then
       KINDI=4; N=2*ngh ! Parity conserved
       nparity=ngh; sparity=two
    Else
       KINDI=4; N=ngh   ! Parity not conserved
       nparity=0;   sparity=one
    End If
    Allocate(endpts(2)); Allocate(b(N),t(N),w(N))
    Call Gaussq(KINDI,N,al,be,kpts,endpts,b,t,w)
    If(ierror_flag.Ne.0) Return
    Do i=nparity+1,N
       j=i-nparity
       xh(j)=t(i)
       ! Build in the Gaussian weight function into the weights wh
       wh(j)=sparity*Exp(xh(j)*xh(j)+Log(w(i)))
    End Do
    Deallocate(endpts,b,t,w)
    !--------------------------------------------------------------------
    !---------------------------->> Gauss-Laguerre <<--------------------|
    !--------------------------------------------------------------------
    KINDI=6; N=ngl
    Allocate(endpts(2)); Allocate(b(N),t(N),w(N))
    Call Gaussq(KINDI,N,al,be,kpts,endpts,b,t,w)
    If(ierror_flag.Ne.0) Return
    Do j=1,ngl
       xl(j)=t(j)
       ! Build in the exponential weight function into the weights wl
       wl(j)=Exp(xl(j)+Log(w(j)))
       sxl(j)=Sqrt(xl(j))
    End Do
    Deallocate(endpts,b,t,w)
    !--------------------------------------------------------------------
    !----------------->> Gauss-Legendre (positive nodes) <<--------------
    !--------------------------------------------------------------------
    If(nleg.Gt.0) Then
       KINDI=1; N=2*nleg
       Allocate(endpts(2)); Allocate(b(N),t(N),w(N))
       Call Gaussq(KINDI,N,al,be,kpts,endpts,b,t,w)
       If(ierror_flag.Ne.0) Return
       Do j=1,nleg
          i=nleg+j
          xleg(j)=t(i); wleg(j)=w(i)
       End Do
       Deallocate(endpts,b,t,w)
    End If
    !
  End Subroutine gausspoints
  !=======================================================================
  !> This set of routines computes the nodes t(j) and weights w(j) for
  !> Gaussian-type quadrature rules with pre-assigned nodes. These are
  !> used when one wishes to approximate
  !>
  !>   \f[
  !>       \int_{a}^{b}  f(x) w(x) dx \approx \sum_{j=1}^{n} w_j f(t_j)
  !>   \f]
  !>
  !> (note w(x) and w(j) have no connection with each other). Here w(x)
  !> is one of six possible non-negative weight functions (listed below),
  !> and f(x) is the function to be integrated. Gaussian quadrature is
  !> particularly useful on infinite intervals (with appropriate weight
  !> functions), since then other techniques often fail. Associated with
  !> each weight function w(x) is a set of orthogonal polynomials. The
  !> nodes t(j) are just the zeroes of the proper n-th degree polynomial.
  !>
  !> Inputs (all real numbers are in double precision)
  !>  - kindi ..: an integer between 1 and 6 giving the type of quadrature
  !>       * 1:  Legendre quadrature, w(x) = 1 on \f$ [-1, 1] \f$
  !>       * 2:  Chebyshev quadrature of the first kind
  !>             \f$ w(x) = 1/\sqrt(1 - x^2) \f$ on \f$ [-1, +1] \f$
  !>       * 3:  Chebyshev quadrature of the second kind
  !>             \f$ w(x) = \sqrt(1 - x^2) \f$ on \f$ [-1, 1] \f$
  !>       * 4:  Hermite quadrature, \f$ w(x) = \exp(-x^2) \f$ on
  !>             \f$ ]-\infty, +\infty[ \f$
  !>       * 5:  Jacobi quadrature, \f$ w(x) = (1-x)^{\alpha}(1+x)^{\beta}\f$
  !>             on \f$ [-1, 1] \f$, \f$ \alpha, \beta > -1 \f$.
  !>             Note: kind=2 and 3 are a special case of this.
  !>       * 6:  generalized Laguerre quadrature,
  !>             \f$ w(x) = \exp(-x) x^{\alpha} \f$ on
  !>             \f$ [0, +\infty[ \f$, with \f$ alpha > -1 \f$.
  !>  - n .....: the number of points used for the quadrature rule
  !>  - alpha .: real parameter used only for Gauss-Jacobi and Gauss-
  !>             Laguerre quadrature (otherwise use 0.d0).
  !>  - beta ..: real parameter used only for Gauss-Jacobi quadrature
  !>             (otherwise use 0.d0)
  !>  - kpts ..: (integer) normally 0, unless the left or right end-
  !>             point (or both) of the interval is required to be a
  !>             node (this is called Gauss-Radau or Gauss-Lobatto
  !>             quadrature). Then kpts is the number of fixed
  !>             endpoints (1 or 2).
  !>  - endpts : real array of length 2. Contains the values of any
  !>             fixed endpoints, if kpts = 1 or 2.
  !>  - b .....: real scratch array of length n
  !>
  !> Outputs (both double precision arrays of length n)
  !>  - t .....: the desired nodes.
  !>  - w .....: the desired weights w(j).
  !>
  !> NOTE: Underflow may sometimes occur, but is harmless.
  !>
  !> Adapted from GO library at www.netlib.org
  !=======================================================================
  Subroutine Gaussq(kindi,n,alpha,beta,kpts,endpts,b,t,w)
    Implicit None
    Integer(ipr) :: N,kindi
    Real(pr):: alpha,beta,MUZERO,GAM,T1
    Integer(ipr) :: j1,J2,kpts,ierr
    Real(pr):: b(n),t(n),w(n),endpts(2)
    !--------------------------------------------------------------------
    !--------------------------------------------------------------------
    !
    Call Class(kindi,n,alpha,beta,b,t,muzero)
    !
    If(KPTS.Eq.0) Then
       W=0.0_pr; W(1)=1._pr
       Call GBTQL2(N,T,B,W,IERR)
       W=MUZERO*W*W
       Return
    End If
    If(KPTS.Eq.2) Then
       GAM=GBSLVE(ENDPTS(1),N,T,B)
       T1=((ENDPTS(1)-ENDPTS(2))/(GBSLVE(ENDPTS(2),N,T,B)-GAM))
       B(N-1)=Sqrt(T1)
       T(N)=ENDPTS(1)+GAM*T1
       W=0.0_pr; W(1)=1._pr
       Call GBTQL2(N,T,B,W,IERR)
       W=MUZERO*W*W
       Return
    End If
    T(N)=GBSLVE(ENDPTS(1),N,T,B)*B(N-1)**2+ENDPTS(1)
    W=0.0_pr; W(1)=1._pr
    Call GBTQL2(N,T,B,W,IERR)
    W=MUZERO*W*W
  End Subroutine Gaussq
  !=======================================================================
  !
  !=======================================================================
  Real(pr) Function GBSLVE(SHIFT,N,A,B)
    Implicit None
    Integer(ipr) :: N,NM1,i
    Real(pr) :: ALPHA,SHIFT,A(N),B(N)
    ALPHA=A(1)-SHIFT
    NM1=N-1
    Do I=2,NM1
       ALPHA=A(I)-SHIFT-B(I-1)**2/ALPHA
    End Do
    GBSLVE=1.0_pr/ALPHA
  End Function GBSLVE
  !=======================================================================
  !>  This procedure supplies the coefficients a(j), b(j) of the
  !>  recurrence relation
  !>
  !>   \f[
  !>       b_j p_j(x) = (x-a_j) p_{j-1}(x) - b_{j-1} p_{j-2}(x)
  !>   \f]
  !>
  !>  for the various classical (normalized) orthogonal polynomials,
  !>  and the zero-th moment
  !>
  !>   \f[
  !>       \mu_{0} = \int w(x) dx
  !>   \f]
  !>
  !>  of the given polynomial's weight function w(x). Since the
  !>  polynomials are orthonormalized, the tridiagonal matrix is
  !>  guaranteed to be symmetric.
  !>
  !>  The input parameter alpha is used only for Laguerre and Jacobi
  !>  polynomials, and the parameter beta is used only for Jacobi
  !>  polynomials. The Laguerre and Jacobi polynomials require the Gamma
  !>  function.
  !>
  !>  Adapted from GO library at www.netlib.org
  !=======================================================================
  Subroutine Class(kindi,N,ALPHA,BETA,B,A,MUZERO)
    Implicit None
    Integer(ipr) :: N,kindi,i,NM1
    Real(pr) :: MUZERO,ALPHA,BETA,A(N),B(N)
    Real(pr) :: PI,ABI,DI20,AB,A2B2,FI
    Data PI / 3.1415926535897930_pr /
    NM1=N-1
    Select Case (kindi)
    Case (1) ! Legendre polynomials
       MUZERO=2.0_pr
       Do I=1,NM1
          A(I)=0.0_pr
          ABI=Real(I,Kind=pr)
          B(I)=ABI/Sqrt(4.0_pr*ABI*ABI-1.0_pr)
       End Do
       A(N)=0.0_pr
    Case (2) ! Chebyshev polynomials of the first kind
       MUZERO=PI
       Do I=1,NM1
          A(I)=0.0_pr
          B(I)=0.50_pr
       End Do
       B(1)=Sqrt(0.50_pr)
       A(N)=0.0_pr
    Case (3) ! Chebyshev polynomials of the second kind
       MUZERO=PI/2.0_pr
       Do I=1,NM1
          A(I)=0.0_pr
          B(I)=0.50_pr
       End Do
       A(N)=0.0_pr
    Case (4) ! Hermite polynomials
       MUZERO=Sqrt(PI)
       Do I=1,NM1
          A(I)=0.0_pr
          DI20=I/2.0_pr
          B(I)=Sqrt(DI20)
       End Do
       A(N)=0.0_pr
    Case (5) ! Jacobi polynomials
       AB=ALPHA+BETA
       ABI=2.0_pr+AB
       MUZERO=2.0_pr**(AB+1.0_pr)*pr_gamma(ALPHA+1.0_pr)*pr_gamma(BETA+1.0_pr)/pr_gamma(ABI)
       A(1)=(BETA-ALPHA)/ABI
       B(1)=Sqrt(4.0_pr*(1.0_pr+ALPHA)*(1.0_pr+BETA)/((ABI+1.0_pr)*ABI*ABI))
       A2B2=BETA*BETA-ALPHA*ALPHA
       Do I=2,NM1
          ABI=2.0_pr*I+AB
          A(I)=A2B2/((ABI-2.0_pr)*ABI)
          FI=I
          B(I)=Sqrt(4.0_pr*FI*(FI+ALPHA)*(FI+BETA)*(FI+AB)/((ABI*ABI-1.0_pr)*ABI*ABI))
       End Do
       ABI=2.0_pr*N+AB
       A(N)=A2B2/((ABI-2.0_pr)*ABI)
    Case (6) ! Laguerre polynomials
       MUZERO=pr_gamma(ALPHA+1.0_pr)
       Do I=1,NM1
          FI=I
          A(I)=2.0_pr*FI-1.0_pr+ALPHA
          B(I)=Sqrt(FI*(FI+ALPHA))
       End Do
       A(N)=2.0_pr*N-1.0_pr+ALPHA
    Case default
    End Select
  End Subroutine Class
  !=======================================================================
  !
  !=======================================================================
  Subroutine GBTQL2(N,D,E,Z,IERR)
    Implicit None
    Integer(ipr) :: N,IERR
    Real(pr) :: D(N),E(N),Z(N)
    Integer(ipr) :: I,J,K,L,M,II,MML
    Real(pr) :: MACHEP,P,G,R,S,C,F,B
    !MACHEP=16.0_pr**(-14)
    MACHEP=epsilon(1.0d0)
    IERR=0
    If(N.Eq.1) Return
    E(N)=0.0_pr
    Do L= 1,N
       J=0
       Do
          Do M=L,N
             If(M .Eq. N) Exit
             If(Abs(E(M)) .Le. MACHEP*(Abs(D(M))+Abs(D(M+1)))) Exit
             Continue
          End Do
          P=D(L)
          If(M .Eq. L) Exit
          If(J .Eq. 30) Then
             IERR=L
             Return
          End If
          J=J+1
          G=(D(L+1)-P) / (2.0_pr*E(L))
          R=Sqrt(G*G+1.0_pr)
          G=D(M) - P + E(L)/(G+Sign(R,G))
          S=1.0_pr
          C=1.0_pr
          P=0.0_pr
          MML=M-L
          Do II=1, MML
             I=M-II
             F=S*E(I)
             B=C*E(I)
             If(Abs(F).Ge.Abs(G)) Then
                C=G/F
                R=Sqrt(C*C+1.0_pr)
                E(I+1)=F*R
                S=1.0_pr/R
                C=C*S
             Else
                S=F/G
                R=Sqrt(S*S+1.0_pr)
                E(I+1)=G*R
                C=1.0_pr/R
                S=S*C
             End If
             G=D(I+1)-P
             R=(D(I)-G)*S + 2.0_pr*C*B
             P=S*R
             D(I+1)=G+P
             G=C*R - B
             F=Z(I+1)
             Z(I+1)=S*Z(I) + C*F
             Z(I)=C*Z(I) - S*F
          End Do
          D(L)=D(L)-P
          E(L)=G
          E(M)=0.0_pr
       End Do
    End Do
    Do II=2, N
       I=II-1
       K=I
       P=D(I)
       Do J=II,N
          If(D(J) .Ge. P) Cycle
          K=J
          P=D(J)
       End Do
       If(K .Eq. I) Cycle
       D(K)=D(I)
       D(I)=P
       P=Z(I)
       Z(I)=Z(K)
       Z(K)=P
    End Do
  End Subroutine GBTQL2
  !=======================================================================
  !> D_HERM calculates Hermite polynomials at point X up to N-th
  !> degree. It also calculates  derivatives  DHER(i)  of these
  !> polynomials.
  !>                    DHER(I)=2*I*HER(I-1)
  !=======================================================================
  Subroutine D_HERM(X,N,HER,DHER,NDIM)
    Integer(ipr), Intent(In) :: N, NDIM
    Real(pr), Intent(in) :: X
    Real(pr), Dimension(0:NDIM), Intent(Inout) :: HER,DHER
    Integer(ipr) :: i
    Real(pr) :: F,DF
    !
    HER (1)=1.0_pr
    DHER(1)=0.0_pr
    !
    If(N.Le.0) Return
    !
    HER (2)=X + X
    DHER(2)=2.0_pr
    !
    If((N-1).Le.0) Return
    !
    Do I=2,N-1
       F=X*HER(I)-Real(I-1,Kind=pr)*HER(I-1)
       HER(I+1)=F+F
       DF=Real(I,Kind=pr)*HER(I)
       DHER(I+1)=DF+DF
    End Do
    !
  End Subroutine D_HERM
  !=======================================================================
  !> The routine computes the expansion of a product of two Hermite     !
  !> polynomials up to order NOSACT as a linear combination of Hermite  !
  !> polynomials. The expansion goes from 0 t0 2*NOSACT. Coefficients   !
  !> are stored in the global variable COEF00(0:NOSACT,0:NOSACT,0:ngh). !
  !> The routine also computes the normalization factors of Hermite     !
  !> polynomials up to order 2*NOSACT, which are stored in HERFAC (also !
  !> a global variable.                                                 !
  !=======================================================================
  Subroutine DEVHER(NOSACT)
    Integer(ipr), Intent(IN) :: NOSACT
    !
    Integer(ipr) :: I,IZEROS,NOSCIL,N,M,K,IGAUSS
    Real(pr) :: SQ2,XZER,RESULT
    Real(pr), Dimension(0:NOSACT) :: PHERMI,DHERMI
    Real(pr), Dimension(0:NOSACT,1:ngh) :: HERPLN,DH1PLN,DH2PLN,HERSQ2
    !
    If(.Not.Allocated(COEF00)) Allocate(COEF00(0:NOSACT,0:NOSACT,0:ngh))
    If(.Not.Allocated(HERFAC)) Allocate(HERFAC(0:2*NOSACT))
    !
    HERFAC(0)=Sqrt(Sqrt(4.0_pr*Atan(1.0_pr)))
    Do I=1,2*NOSACT
       HERFAC(I)=HERFAC(I-1)*Sqrt(Real(2*I,Kind=pr))
    End Do
    !
    ! Defining the values of the Hermite polynomials at gauss zeros
    SQ2=Sqrt(2.0_pr)
    Do IZEROS = 1,ngh
       XZER = xh(IZEROS)
       Call D_HERM(XZER,NOSACT,PHERMI,DHERMI,NOSACT)
       Do NOSCIL = 0,NOSACT
          HERPLN(NOSCIL,IZEROS)= PHERMI(NOSCIL)/HERFAC(NOSCIL)
          DH1PLN(NOSCIL,IZEROS)=(DHERMI(NOSCIL)- XZER*PHERMI(NOSCIL))/HERFAC(NOSCIL)
          DH2PLN(NOSCIL,IZEROS)=(XZER*XZER-2*NOSCIL-1)*PHERMI(NOSCIL)/HERFAC(NOSCIL)
       End Do
       XZER = xh(IZEROS) / SQ2
       Call D_HERM(XZER,NOSACT,PHERMI,DHERMI,NOSACT)
       Do NOSCIL = 0,NOSACT
          HERSQ2(NOSCIL,IZEROS)= PHERMI(NOSCIL)/HERFAC(NOSCIL)
       End Do
    End Do
    ! Calculating the expansion coefficients for the polynomials
    Do N=0,NOSACT
       Do M=0,NOSACT
          Do K=M+N,0,-2
             RESULT=0.0_pr
             Do IGAUSS=1,ngh
                RESULT=RESULT+wh(IGAUSS)*HERPLN(N,IGAUSS)*HERPLN(M,IGAUSS)*HERPLN(K,IGAUSS)
             End Do
             COEF00(K,M,N)=RESULT
          End Do
       End Do
    End Do
    ! The code uses dimensioned wave functions and their derivatives.
    ! Therefore, the coefficients are now scaled by appropriate factors
    ! which take into account these dimensions.
    Do N=0,NOSACT
       Do M=0,NOSACT
          Do K=M+N,0,-2
             COEF00(K,M,N)=COEF00(K,M,N)*Sqrt(bz)
          End Do
       End Do
    End Do
    !
  End Subroutine DEVHER
  !=======================================================================
  !>  pr_gamma evaluates \f$ \Gamma(x) \f$ for a real argument.
  !>
  !>  Discussion:
  !>    This function was originally named DGAMMA. However, a number of
  !>    Fortran compilers now include a library function of this name. To
  !>    avoid conflicts, this function was renamed pr_gamma. This routine
  !>    calculates the GAMMA function for a real argument X. Computation
  !>    is based on an algorithm outlined in reference 1 below. The
  !>    program uses rational functions that approximate the GAMMA
  !>    function to at least 20 significant decimal digits. Coefficients
  !>    for the approximation over the interval (1,2) are unpublished.
  !>    Those for the approximation for 12 <= X are from reference 2.
  !>
  !>  Licensing:
  !>    This code is distributed under the GNU LGPL license.
  !>
  !>  Modified:
  !>    18 January 2008
  !>
  !>  Author:
  !>    Original FORTRAN77 version by William Cody, Laura Stoltz.
  !>    FORTRAN90 version by John Burkardt.
  !>
  !>  Reference:
  !>    - 1. William Cody, "An Overview of Software Development for Special
  !>         Functions," in Numerical Analysis Dundee, 1975, Edited by GA
  !>         Watson, Lecture Notes in Mathematics 506, Springer, 1976.
  !>    - 2. John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles
  !>        Mesztenyi, John Rice, Henry Thatcher, Christoph Witzgall,
  !>        "Computer Approximations," Wiley, 1968, LC: QA297.C64.
  !>
  !>  Parameters:
  !>    - Input, real ( kind = 8 ) X, the argument of the function.
  !>    - Output, real ( kind = 8 ) R8_GAMMA, the value of the function.
  !=======================================================================
  Real(pr) Function pr_gamma(x)
    Implicit None
    !
    !  Coefficients for minimax approximation over (12, INF).
    Logical :: parity
    Integer(ipr) :: i,n
    Real(pr), Dimension(7) :: c = (/ -1.910444077728000000000D-03, &
      8.417138778129500000000000D-04,-5.952379913043012000000D-04, &
      7.936507935003502480000000D-04,-2.777777777777681622553D-03, &
      8.333333333333333331554247D-02, 5.708383526100000000000D-03 /)
    Real(pr) :: eps,fact,pi,res,sqrtpi,sum,x,xbig,xden,xinf,xminin,xnum,y,y1,ysq,z
    Real(pr) :: p(8),q(8)
    !
    ! Mathematical constants
    Data sqrtpi /0.9189385332046727417803297D+00/
    Data pi /3.1415926535897932384626434D+00/
    !
    ! Numerator and denominator coefficients for rational minimax
    ! approximation over (1,2).
    Data p / -1.71618513886549492533811D+00,  2.47656508055759199108314D+01, &
             -3.79804256470945635097577D+02,  6.29331155312818442661052D+02, &
              8.66966202790413211295064D+02, -3.14512729688483675254357D+04, &
             -3.61444134186911729807069D+04,  6.64561438202405440627855D+04 /

    Data q / -3.08402300119738975254353D+01,  3.15350626979604161529144D+02, &
             -1.01515636749021914166146D+03, -3.10777167157231109440444D+03, &
              2.25381184209801510330112D+04,  4.75584627752788110767815D+03, &
             -1.34659959864969306392456D+05, -1.15132259675553483497211D+05 /

    parity = .False.; fact = one; n = 0; y = x
    xbig = 171.624D+00
    xminin = Tiny(1.0D0); eps = Epsilon(1.0D0) ; xinf = Huge(1.0D0)
    !
    ! Argument is negative.
    If(y<=zero) Then

       y = - x
       y1 = Aint ( y )
       res = y - y1

       If(res/=zero) Then
          If(y1/=Aint(y1*half)*two) Then
             parity = .True.
          End If
          fact = - pi / Sin(pi*res)
          y = y + one
       Else
          res = xinf
          pr_gamma = res
          Return
       End If

    End If
    !
    ! Argument is positive.
    if (y < eps) Then
       !
       ! Argument < EPS.
       If(xminin <= y) Then
          res = one / y
       Else
          res = xinf
          pr_gamma = res
          Return
       End If

    Else If(y<pp12) Then

       y1 = y

       ! 0.0 < argument < 1.0.
       If(y<one) Then
          z = y
          y = y + one
       Else
          ! 1.0 < argument < 12.0.
          ! Reduce argument if necessary.
          n = Int(y) - 1
          y = y - Real(n,Kind=pr)
          z = y - one
       End If
       !
       ! Evaluate approximation for 1.0 < argument < 2.0.
       xnum = zero
       xden = one
       Do i = 1, 8
          xnum = ( xnum + p(i) ) * z
          xden = xden * z + q(i)
       End Do

       res = xnum / xden + one
       !
       ! Adjust result for case  0.0 < argument < 1.0.
       If(y1<y) Then
          res = res / y1
       Else If(y<y1) Then
          ! Adjust result for case 2.0 < argument < 12.0.
          Do i = 1, n
             res = res * y
             y = y + one
          End Do
       End If

    Else
       !
       !  Evaluate for 12.0 <= argument.
       If(y<=xbig) Then
          ysq = y * y
          sum = c(7)
          Do i = 1, 6
             sum = sum / ysq + c(i)
          End Do
          sum = sum / y - y + sqrtpi
          sum = sum + ( y - half ) * Log(y)
          res = Exp(sum)
       Else
          res = huge ( res )
          pr_gamma = res
          Return
       End If

    End If
    !
    !  Final adjustments and Return.
    If(parity) Then
       res = - res
    End If

    If(fact/=one) Then
       res = fact / res
    End If

    pr_gamma = res

    Return
  End Function pr_gamma  ! Gamma function in double precision
  !=======================================================================
  !
  !=======================================================================
End Module HFBTHO_gauss

