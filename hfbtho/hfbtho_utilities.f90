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
!                         BASIC UTILITIES PACKAGE                      !
!                                                                      !
! ==================================================================== !

!----------------------------------------------------------------------
!> This module defines basic data types, constants and file units.
!----------------------------------------------------------------------
!  Subroutines: - get_CPU_time(subname,is)
!----------------------------------------------------------------------
Module HFBTHO_utilities

  Implicit None

  Integer, Parameter, Public :: ipr=Kind(1)     ! to set the precision of the DFT solver
  Integer, Parameter, Public :: pr =Kind(1.0d0) ! to set the precision of the DFT solver

  ! I/O
  Integer, Public :: lout = 6, lfile = 7

  ! Global numbers
  Real(pr), Parameter :: zero=0.0_pr,half= 0.5_pr,one=1.0_pr,two  =2.0_pr,three=3.0_pr, &
                         four=4.0_pr,five= 5.0_pr,six=6.0_pr,seven=7.0_pr,eight=8.0_pr, &
                         nine=9.0_pr,ten =10.0_pr
  ! Whole global numbers pp#
  Real(pr), Parameter :: pp12=12.0_pr,pp16=16.0_pr,pp15=15.0_pr,pp20=20.0_pr, &
                         pp24=24.0_pr,pp27=27.0_pr,pp32=32.0_pr,pp64=64.0_pr, &
                         pp40=40.0_pr
  ! Fractional global numbers p#
  Real(pr), Parameter :: p12=one/two,   p13=one/three,  p14=one/four,  p23=two/three, &
                         p43=four/three,p32=three/two,  p34=three/four,p53=five/three,&
                         p18=one/eight, p38=three/eight,p59=five/nine, p52=five/two,  &
                         p54=five/four, p74=seven/four

 Contains
  !=======================================================================
  !
  !=======================================================================
  Subroutine get_CPU_time(subname,is)
    Implicit None
    Integer, Intent(in)   :: is
    Character*(*), Intent(in)  :: subname
    Character(Len=15) :: subprint
    Integer(ipr), Save :: t1,t2,countrate,countmax
    !
    If(is.Eq.0) Then
       Call system_clock(t1,countrate,countmax)
    Else
       Call system_clock(t2,countrate,countmax)
       Write(subprint,'(a15)') subname
       Write(*,'(a,a15,a,F16.6)') '  Time in seconds -> ',subprint,':',(t2-t1)/real(countrate,kind=pr)
    End If
    !
  End Subroutine get_CPU_time
  !=======================================================================
  !
  !=======================================================================
End Module HFBTHO_utilities
