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

#if(READ_FUNCTIONAL==1)
!----------------------------------------------------------------------
!>  This module contains subroutines to read the parameters of an
!>  arbitrary Skyrme-like functional from a file instead of using the
!>  ones that are hardcoded in the UNEDF module (see hfbtho_unedf.f90).
!>  This is particularly useful for making calculations with statistical
!>  ensembles of parameter sets.
!>
!> @author
!>    Rodrigo Navarro Perez
!----------------------------------------------------------------------
!  Subroutines: - read_HFBTHO_Functional
!               - replace_functional
!               - broadcast_functional
!----------------------------------------------------------------------
Module HFBTHO_read_functional
  Use HFBTHO_utilities
  Use HFBTHO
  Implicit None
Contains

  subroutine read_HFBTHO_Functional
    implicit none
    integer(ipr) :: lfunctional=16,i
    open(lfunctional,file='hfbtho_FUNCTIONAL.dat')
    read(lfunctional,*)
    read(lfunctional,*) (functional_vector(i),i=1,n_func_param)
    close(lfunctional)
  end subroutine read_HFBTHO_Functional

  subroutine replace_functional
    implicit none
    RHO_NM   = functional_vector( 1)
    E_NM     = functional_vector( 2)
    K_NM     = functional_vector( 3)
    ASS_NM   = functional_vector( 4)
    LASS_NM  = functional_vector( 5)
    SMASS_NM = functional_vector( 6)
    CrDr(0)  = functional_vector( 7)
    CrDr(1)  = functional_vector( 8)
    CpV0(0)  = functional_vector( 9)
    CpV0(1)  = functional_vector(10)
    CrdJ(0)  = functional_vector(11)
    CrdJ(1)  = functional_vector(12)
  end subroutine replace_functional

#if(USE_MPI==1)
  subroutine broadcast_functional
    implicit none
    include 'mpif.h'
    call mpi_bcast(functional_vector,12,mpi_double_precision,&
         0,mpi_comm_world,ierr_mpi)

  end subroutine broadcast_functional
#endif

end module HFBTHO_read_functional

#endif
