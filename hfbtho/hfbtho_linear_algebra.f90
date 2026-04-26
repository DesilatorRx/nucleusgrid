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
!                      LINEAR ALGEBRA PACKAGE                          !
!                                                                      !
! ==================================================================== !

!-------------------------------------------------------------------
!> This module several routines implementing basic linear algebra
!> and mathematical operations
!----------------------------------------------------------------------
!  Subroutines: - lingd(ma,mx,n,m,a,x,d,Ifl)
!               - csplin(n, x, y, b, c, d)
!               - cseval(n,u,x,y,b,c,d,splf0)
!               - deri(h,n,f1,dunl)
!----------------------------------------------------------------------
Module linear_algebra

  Use HFBTHO_utilities

  Implicit None

Contains
  !=======================================================================
  !> Solves the system of linear equations A*X = B
  !> At the beginning the matrix B is stored in X. During the calculation
  !> it will be overwritten. D is the determinant of A
  !=======================================================================
  Subroutine lingd(ma,mx,n,m,a,x,d,Ifl)
    Integer(ipr) :: ma,mx,n,m,Ifl
    Integer(ipr), Save :: i,j,k,l,k1,n1
    Real(pr) ::  a(ma,m),x(mx,m),d
    Real(pr), Save :: tollim,one,zero,p,q,tol,cp,cq
    Data tollim/1.d-10/,one/1.d0/,zero/0.d0/
    Ifl = 1; p = zero
    Do i=1,n
       q = zero
       Do j=1,n
          q = q + Abs(a(i,j))
       End Do
       If(q.Gt.p) p = q
    End Do
    tol = tollim*p; d   = one
    Do k=1,n
       p = zero
       Do j=k,n
          q = Abs(a(j,k))
          If(q.Lt.p) Cycle
          p = q; i = j
       End Do
       If (p.Lt.tol) Then
          Write (6,200) ('-',j=1,80),tol,i,k,a(i,k),('-',j=1,80)
200     Format (/1x,80a1/' *****  ERROR IN LINGD , TOLERANZ =',e10.4, &
               ' VALUE OF A(',i3,',',i3,') IS ',e10.4/1x,80a1)
          Ifl = -1
          Return
       End If
       cp = one/a(i,k)
       If(i.Ne.k) Then
          d = -d
          Do l=1,m
             cq = x(i,l); x(i,l) = x(k,l); x(k,l) = cq
          End Do
          Do l=k,n
             cq = a(i,l); a(i,l) = a(k,l); a(k,l) = cq
          End Do
       End If
       d = d*a(k,k)
       If(k.Eq.n) Exit
       k1 = k + 1
       Do i=k1,n
          cq=a(i,k)*cp
          Do l=1,m
             x(i,l)=x(i,l)-cq*x(k,l)
          End Do
          Do l=k1,n
             a(i,l)=a(i,l)-cq*a(k,l)
          End Do
       End Do
    End Do
    Do l=1,m
       x(n,l)=x(n,l)*cp
    End Do
    If(n.Eq.1) Return
    n1=n-1
    Do k=1,n1
       cp = one/a(n-k,n-k)
       Do l=1,m
          cq = x(n-k,l)
          Do i=1,k
             cq = cq-a(n-k,n+1-i)*x(n+1-i,l)
          End Do
          x(n-k,l) = cq*cp
       End Do
    End Do
    Return
  End Subroutine lingd
  !=======================================================================
  !> The coefficients \f$ b_i, c_i, d_i, i=1,2,\dots,n \f$ are computed
  !> for a cubic interpolating spline
  !>   \f[
  !>      s(x) = y_i + b_i (x-x_i) + c_i (x-x_i)^2 + d_i (x-x_i)^3
  !>   \f]
  !> for \f$ x_i \leq x \leq x_{i+1} \f$
  !> Input
  !>   - n = the number of data points or knots (n.ge.2)
  !>   - x = the abscissas of the knots in strictly increasing order
  !>   - y = the ordinates of the knots
  !> Output
  !>   - b, c, d  = arrays of spline coefficients as defined above.
  !>   \f[
  !>       y_i = s(x_i), \quad\quad
  !>       b_i = \frac{ds}{dx}(x_i),  \quad\quad
  !>       c_i = \frac{1}{2}\frac{d^{2}s}{dx^{2}}(x_i),  \quad\quad
  !>       d_i = \frac{1}{6}\frac{d^{3}s}{dx^{3}}(x_i)
  !>   \f]
  !> The accompanying function subprogram \ref cseval can be used to
  !> evaluate the spline, its derivative or even its 2nd derivative.
  !=======================================================================
  Subroutine csplin(n, x, y, b, c, d)
    Integer(ipr), Save :: nm1,i,ib
    Integer(ipr) :: n
    Real(pr) :: x(n), y(n), b(n), c(n), d(n)
    Real(pr), Save :: t,zero=0.0d0,two=2.0d0,tr=3.0d0
    ! check input for consistency
    If(n.Lt.2) Stop '-n < 2 in csplin call--'
    nm1 = n-1
    Do i = 1, nm1
       If(x(i).Ge.x(i+1)) Stop 'x not strictly ascending in csplin call'
    End Do
    If (n.Ne.2) Then
       ! set up tridiagonal system
       ! b = diagonal, d = offdiagonal, c = right hand side.
       d(1) = x(2) - x(1); c(2) = (y(2) - y(1))/d(1)
       Do i = 2, nm1
          d(i) = x(i+1) - x(i); b(i) = two*(d(i-1) + d(i))
          c(i+1) = (y(i+1) - y(i))/d(i); c(i) = c(i+1) - c(i)
       End Do
       ! end conditions.  third derivatives at  x(1)  and  x(n)
       ! obtained from divided dIfferences
       b(1) = -d(1); b(n) = -d(n-1); c(1) = zero; c(n) = zero
       If (n.Ne.3) Then
          c(1) =  c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
          c(n) =  c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
          c(1) =  c(1)*d(1)**2/(x(4)-x(1))
          c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
          ! forward elimination
       Else
          Do i = 2, n
             t = d(i-1)/b(i-1); b(i) = b(i) - t*d(i-1); c(i) = c(i) - t*c(i-1)
          End Do
       End If
       ! back substitution
       c(n) = c(n)/b(n)
       Do ib = 1, nm1
          i = n-ib
          c(i) = (c(i) - d(i)*c(i+1))/b(i)
       End Do
       ! compute polynomial coefficients
       b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + two*c(n))
       Do i = 1, nm1
          b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + two*c(i))
          d(i) = (c(i+1) - c(i))/d(i); c(i) = tr*c(i)
       End Do
       c(n) = tr*c(n); d(n) = d(n-1)
       Return
    Else
       b(1) = (y(2)-y(1))/(x(2)-x(1)); c(1) = zero; d(1) = zero
       Return
    End If
  End Subroutine csplin
  !=======================================================================
  !
  !=======================================================================
  Subroutine cseval(n,u,x,y,b,c,d,splf0)
    Integer(ipr) :: n
    Integer(ipr), Save :: i=1,j,k
    Real(pr) :: x(n),y(n),b(n),c(n),d(n),u,splf0
    Real(pr), Save :: dx
    If(i.Ge.n)      i = 1
    If(u.Lt.x(i))   Go To 10
    If(u.Le.x(i+1)) Go To 30
    ! binary search
10  i = 1
    j = n + 1
20  k = (i+j)/2
    If(u.Lt.x(k)) j = k
    If(u.Ge.x(k)) i = k
    If(j.Gt.i+1) Go To 20
    ! evaluate splf0
30 dx = u - x(i)
    splf0 = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
    Return
  End Subroutine cseval
  !=======================================================================
  !> First derivative of 'f1' if the step is 'h'
  !=======================================================================
  Subroutine deri(h,n,f1,dunl)
    Integer(ipr) :: n
    Integer(ipr), Save :: k
    Real(pr) :: h,f1(n),dunl(n)
    Real(pr), Save :: t60,t12
    Real(pr), Save :: t8=8.0d0,t45=45.0d0,t9=9.0d0
    t60 =1.0d0/(h*60.0d0); t12 =1.0d0/(h*12.0d0)
    !
    dunl(1)  =(t8*f1(2)-f1(3)+f1(1))*t12
    dunl(2)  =(t45*(f1(3)-f1(1))-t9*f1(4)+f1(5)-f1(1))*t60
    dunl(3)  =(t45*(f1(4)-f1(2))-t9*(f1(5)-f1(1))+f1(6))*t60
    dunl(n)  =(-t8*f1(n-1)+f1(n)+f1(n-2))*t12
    dunl(n-1)=(t45*(f1(n)-f1(n-2))+t9*f1(n-3)-f1(n)-f1(n-4))*t60
    dunl(n-2)=(t45*(f1(n-1)-f1(n-3))-t9*(f1(n)-f1(n-4))-f1(n-5))*t60
    Do k=4,n-3
       dunl(k) =(t45*(f1(k+1)-f1(k-1))-t9*(f1(k+2)-f1(k-2))+f1(k+3)-f1(k-3))*t60
    End Do
    Return
  End Subroutine deri
  !=======================================================================
  !
  !=======================================================================
End Module linear_algebra

