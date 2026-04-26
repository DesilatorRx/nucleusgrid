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
!                        QRPA-pnFAM INTERFACE PACKAGE                  !
!                                                                      !
! ==================================================================== !

!----------------------------------------------------------------------
!> This module provides an interface between HFBTHO and the QRPA-pnFAM
!> code developed by the Chapel Hill group.
!>
!>  @author
!>    Mikka Mustonen, Thomas Shafer
!----------------------------------------------------------------------
!  Subroutines: - save_HFBTHO_solution(plaintext)
!----------------------------------------------------------------------
Module HFBTHO_storage
   implicit none

   ! Version number to check against
   integer, parameter :: VERSION = 9

   contains

   ! ---------------------------------------------------------------------------
   ! Saves the relevant HFBTHO solution information to file for re-use by FAM.
   ! ---------------------------------------------------------------------------
   subroutine save_HFBTHO_solution(plaintext)
      use UNEDF,  only: cr0 => crho, crr => cdrho, ctau, cdrho => crdr, ctj => cj, crdj, use_j2terms
      use HFBTHO, only: nb, id, REqpP, RUqpP, RVqpP, REqpN, RUqpN, RVqpN, nr, nz, nl, ns, npar,    &
                        nghl, wdcori, y_opt, fh, qhla_opt, fi1r_opt, fi1z_opt, fi2d_opt, npr, ala, &
                        alast, del, cpv0, cpv1, bet, ehfb, pwi, fn_T, fp_T, temper, entropy,       &
                        switch_on_temperature, nbx, ka, kd, KqpP, KpwiP, KqpN, KpwiN

      implicit none

      logical, intent(in), optional :: plaintext
      logical :: save_as_text

      integer :: ifh, ierr

      ! Option to save in plain text
      save_as_text = .false.
      if (present(plaintext)) then
         if (plaintext .eqv. .true.) then
            save_as_text = .true.
         end if
      end if

      ! Open the file
      ifh = 77
      if (save_as_text) then
         write(*,'(/a)') ' ### STORING HFB SOLUTION (PLAINTEXT format)'
         write(*,'(a)')  ' Filename: "hfb.solution.txt"'
         open(unit=ifh, file='hfb.solution.txt', status='unknown', iostat=ierr)
      else
         write(*,'(/a)') ' ### STORING HFB SOLUTION (BINARY format)'
         write(*,'(a)')  ' Filename: "hfb.solution"'
         open(unit=ifh, file='hfb.solution', status='unknown', form='unformatted', iostat=ierr)
      end if

      if (ierr /= 0) then
         write(*,'(a,i0)') ' Error saving HFB solution: could not open file to write. Error', ierr
         stop
      end if

      ! ------------------------------------------------------------------------
      ! Store the HFB solution details to 'hfb.solution(.txt if plain text)'
      ! ------------------------------------------------------------------------

      if (.not. save_as_text) then
         ! The solution version number
         write(ifh) VERSION

         ! Basic HFB quantities
         ! N, Z, A where applicable
         write(ifh) npr(:)             ! Particle number
         write(ifh) ala(:)             ! Lambdas ala (quasiparticles measured w.r.t. these, I think)
         write(ifh) alast(:)           ! Lambdas alast (last-bound s.p. energy)
         write(ifh) del(:)             ! Pairing gaps
         write(ifh) pwi                ! Pairing window
         write(ifh) CpV0(:), CpV1(:)   ! Pairing strengths
         write(ifh) bet                ! Total deformation
         write(ifh) ehfb               ! Binding energy

         ! The HFB quasiparticle energies and amplitudes
         write(ifh) nbx
         write(ifh) nb
         write(ifh) id(:)
         write(ifh) REqpP(:)
         write(ifh) RVqpP(:)
         write(ifh) RUqpP(:)
         write(ifh) REqpN(:)
         write(ifh) RVqpN(:)
         write(ifh) RUqpN(:)

         ! Pairing window active q.p. levels
         write(ifh) ka(:,:)
         write(ifh) kd(:,:)
         write(ifh) KqpP(:)
         write(ifh) KpwiP(:)
         write(ifh) KqpN(:)
         write(ifh) KpwiN(:)

         ! Basis state quantum numbers
         write(ifh) nr(:)
         write(ifh) nz(:)
         write(ifh) nl(:)
         write(ifh) ns(:)
         write(ifh) npar(:)

         ! Wave functions and integration data
         write(ifh) nghl                      ! number of integration points
         write(ifh) wdcori(:)                 ! inverse of integration weights
         write(ifh) y_opt(:)                  ! 1/rho in 'fm^(-1)'
         write(ifh) fh(:)                     ! z in 'fm'
         write(ifh) transpose(qhla_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)
         write(ifh) transpose(fi1r_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)
         write(ifh) transpose(fi1z_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)
         write(ifh) transpose(fi2d_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)

         ! Time-even isovector coupling constants from HFB mean field
         write(ifh) cr0(1), crr(1), cdrho(1), ctau(1), ctj(1), crdj(1)
         ! Store .TRUE. => 1, .FALSE. => 0
         write(ifh) merge(1, 0, use_j2terms)

         ! Temperature-dependence of the HFB solution
         ! Store .TRUE. => 1, .FALSE. => 0
         write(ifh) merge(1, 0, switch_on_temperature)
         write(ifh) temper
         write(ifh) entropy
         write(ifh) fp_T(:)
         write(ifh) fn_T(:)
      else
         ! The solution version number
         write(ifh,*) VERSION

         ! Basic HFB quantities
         ! N, Z, A where applicable
         write(ifh,*) npr(:)           ! Requested particle number
         write(ifh,*) ala(:)           ! Lambdas ala (quasiparticles measured w.r.t. these, I think)
         write(ifh,*) alast(:)         ! Lambdas alast (last-bound s.p. energy)
         write(ifh,*) del(:)           ! Pairing gaps
         write(ifh,*) pwi              ! Pairing window
         write(ifh,*) CpV0(:), CpV1(:) ! Pairing strengths
         write(ifh,*) bet              ! Total deformation
         write(ifh,*) ehfb             ! Binding energy

         ! The HFB quasiparticle energies and amplitudes
         write(ifh,*) nbx
         write(ifh,*) nb
         write(ifh,*) id(:)
         write(ifh,*) REqpP(:)
         write(ifh,*) RVqpP(:)
         write(ifh,*) RUqpP(:)
         write(ifh,*) REqpN(:)
         write(ifh,*) RVqpN(:)
         write(ifh,*) RUqpN(:)

         ! Pairing window active q.p. levels
         write(ifh,*) ka(:,:)
         write(ifh,*) kd(:,:)
         write(ifh,*) KqpP(:)
         write(ifh,*) KpwiP(:)
         write(ifh,*) KqpN(:)
         write(ifh,*) KpwiN(:)

         ! Basis state quantum numbers
         write(ifh,*) nr(:)
         write(ifh,*) nz(:)
         write(ifh,*) nl(:)
         write(ifh,*) ns(:)
         write(ifh,*) npar(:)

         ! Wave functions and integration data
        write(ifh,*) nghl                      ! number of integration points
        write(ifh,*) wdcori(:)                 ! inverse of integration weights
        write(ifh,*) y_opt(:)                  ! 1/rho in 'fm^(-1)'
        write(ifh,*) fh(:)                     ! z in 'fm'
        write(ifh,*) transpose(qhla_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)
        write(ifh,*) transpose(fi1r_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)
        write(ifh,*) transpose(fi1z_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)
        write(ifh,*) transpose(fi2d_opt(:,:))  ! 1st index: state, 2nd index: coordinate (before transposing)

        ! Time-even isovector coupling constants from HFB mean field
        write(ifh,*) cr0(1), crr(1), cdrho(1), ctau(1), ctj(1), crdj(1)
        ! Store .TRUE. => 1, .FALSE. => 0
        write(ifh,*) merge(1, 0, use_j2terms)

        ! Temperature-dependence of the HFB solution
        ! Store .TRUE. => 1, .FALSE. => 0
        write(ifh,*) merge(1, 0, switch_on_temperature)
        write(ifh,*) temper
        write(ifh,*) entropy
        write(ifh,*) fp_T(:)
        write(ifh,*) fn_T(:)
      end if

      close(ifh)
      write(*,'(a)') ' Storage completed.'
   end subroutine save_HFBTHO_solution
  !=======================================================================
  !
  !=======================================================================
end module HFBTHO_storage
