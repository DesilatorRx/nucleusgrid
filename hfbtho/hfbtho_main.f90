!***********************************************************************
!
!    Copyright (c) 2012, Lawrence Livermore National Security, LLC.
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

Program hfbthoprog
  Use HFBTHO
  Use HFBTHO_utilities
  Call Main_Program
End Program hfbthoprog
  !=======================================================================
  !
  !=======================================================================
  Subroutine Main_Program
    Use HFBTHO_utilities
    Use HFBTHO
#if(USE_QRPA==1)
    Use HFBTHO_storage
#endif
#if(USE_MPI==1)
    Use HFBTHO_mpi_communication
#endif
#if(DO_MASSTABLE==1 || DRIP_LINES==1)
    Use HFBTHO_large_scale
#endif
#if(READ_FUNCTIONAL==1)
    USE HFBTHO_read_functional
#endif
    USE HFBTHO_localization
    Use HFBTHO_solver
    Implicit None
    Integer(ipr) :: iblocase(2),nkblocase(2,5)
    Integer(ipr) :: i,it,icount,jcount,l,noForce,icalc
    Integer(ipr) :: j, ib, mu, nu
    Logical      :: file_exists
    ! Initialize MPI environment and possibly create subcommunicators
#if(USE_MPI==1)
    call MPI_INIT(ierr_mpi)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr_mpi)
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_taskid, ierr_mpi)
#if(DRIP_LINES==1)
    call Create_MPI_Teams
#endif
#else
    mpi_size = 1
    mpi_taskid = 0
#endif
    !
    write(ID_string,'(1i3.3)') mpi_taskid
    task_error = 0
    !-------------------------------------------------------------------
    ! Read input and namelist data for the requested nucleus
    !-------------------------------------------------------------------
    Call initialize_HFBTHO_NAMELIST
    !-------------------------------------------------------------------
    !  Process 0 reads general input data. In serial mode, that is all
    !  the data needed.
    !-------------------------------------------------------------------
    if(mpi_taskid.eq.0) then
       Call read_HFBTHO_NAMELIST
#if(DO_MASSTABLE==1)
       call read_HFBTHO_MassTable
#endif
#if(DRIP_LINES==1)
       call read_HFBTHO_StableLine
#endif
#if(READ_FUNCTIONAL==1)
       call read_HFBTHO_Functional
#endif
    endif
    !-------------------------------------------------------------------
    ! Broadcast of process-independent input data in parallel mode
    !-------------------------------------------------------------------
#if(USE_MPI==1)
    if(mpi_taskid.eq.0) then
       call Construct_Vectors
       allocate(task_error_gthr(0:mpi_size-1))
    endif
    call broadcast_vectors
    if(nrows.gt.0) then
       if(mpi_taskid.gt.0) call Deconstruct_Vectors
    endif
#if(READ_FUNCTIONAL==1)
       call broadcast_functional
#endif
#endif
    !-------------------------------------------------------------------
    ! Allocation of vectors for mass tables
    !-------------------------------------------------------------------
#if(DO_MASSTABLE==1)
    if(nrows.gt.0) then
#if(USE_MPI==1)
       call allocate_out_vectors
#else
       call allocate_mass_table
#endif
    endif
#endif
    !-------------------------------------------------------------------
    ! Mass table calculation
    !-------------------------------------------------------------------
#if(DO_MASSTABLE==1)
    pi = 4*atan(1.0_pr)
    icalc = 0
    !team leader creates a file for bookkeeping
    if(mpi_taskid.eq.0) then
       open(127,file='TableLog.dat')
    endif
    !loop over elements of the mass table
    do iRow = 0,nRows
       Z_chain = Z_masstable(iRow)
       N_chain = N_masstable(iRow)
       A_chain = Z_chain+N_chain
       beta_deformation = beta_masstable(iRow)
       Q20 = Q20_masstable(iRow)!beta_deformation*sqrt(5/pi)*(A_chain)**(5/3._pr)/100._pr
       write(row_string,'(1a1,1i4.4)') '_',irow
       if(irow.eq.0.and.nrows.gt.0) cycle
       if(mpi_taskid.eq.0) then
          write(127,'(a5,2i5,2f15.8)') row_string,Z_chain,&
               N_chain,Q20,beta_deformation
       endif
       !only do the calculations that correspond to your task id
       if(mod(irow,mpi_size).ne.mpi_taskid) cycle
       proton_number  = Z_chain
       neutron_number = N_chain
       expectation_values(2) = Q20
       basis_deformation = beta_deformation
#endif
    !-------------------------------------------------------------------
    ! Dripline calculation
    !-------------------------------------------------------------------
#if(DRIP_LINES==1)
    beta_step = 1.0_pr/real(number_deformations-1,kind=pr)
    pi = 4*atan(1.0_pr)
    write(team_string,'(1i3.3)') team_color
    !team leader allocates array to recieve energies and creates a file for bookkeeping
    if(team_rank.eq.0) then
       allocate(energy_chain_gthr(0:team_size-1))
       open(127,file='TeamTable'//team_string//'.dat')
    endif
    calc_counter = 0
    !loop over the nuclei in the "stable line"
    do iRow = 1,nRows
       !only do the chains that correspond to your team
       if(mod(irow,number_teams).ne.team_color) cycle
       Minimum_Energy_Prev = 100
       Z_chain = Z_stable_line(iRow)
       N_chain = N_stable_line(iRow)
       !move along the isotopic (or isotonic) chain until the drip line is reached
       do
          A_chain = Z_chain + N_chain
          beta_deformation = -0.5_pr - beta_step
          !loop over the different basis deformations of each nucleus
          do i_deformation = 1,number_deformations
             beta_deformation = beta_deformation + beta_step
             calc_counter = calc_counter + 1
             Q20 = beta_deformation*sqrt(5/pi)*(A_chain)**(5/3._pr)/&
                  100._pr
             write(row_string,'(a1,a3,a1,i6.6)') '_',team_string,'_',calc_counter
             !team leader writes type of calculation for bookkeeping
             if(team_rank.eq.0) then
                if(direction_sl(irow).eq.1) then
                   direction_str = ' isotopic'
                else
                   direction_str = ' isotonic'
                endif
                write(127,'(a11,2i5,2f15.8,a9)') row_string,z_chain,&
                     n_chain,Q20,beta_deformation,direction_str
             endif
             if(number_deformations.eq.1) beta_deformation = 0._pr
             !only calculate what corresponds to each process
             if(mod(i_deformation,team_size).ne.team_rank) cycle
             proton_number  = Z_chain
             neutron_number = N_chain
             expectation_values(2) =  Q20
             basis_deformation = beta_deformation
#endif
       !-------------------------------------------------------------------
       ! Regular HFBTHO execution, whether or not mass table or dripline
       ! mode is activated
       !-------------------------------------------------------------------
       !memo: Namelist /HFBTHO_GENERAL/ number_of_shells,oscillator_length,&
       !                                proton_number,neutron_number,type_of_calculation
       !      Namelist /HFBTHO_INITIAL/ beta2_deformation, beta4_deformation
       !      Namelist /HFBTHO_ITERATIONS/ number_iterations, accuracy
       !      Namelist /HFBTHO_FUNCTIONAL/ functional, add_initial_pairing, type_of_coulomb
       !      Namelist /HFBTHO_CONSTRAINTS/ lambda_values, lambda_active, expectation_values
       !      Namelist /HFBTHO_BLOCKING/ proton_blocking(1:5), neutron_blocking(1:5)
       !      Namelist /HFBTHO_BLOCKING/ switch_to_THO, projection_is_on, gauge_points, delta_N, delta_P
       !      Namelist /HFBTHO_TEMPERATURE/ set_temperature, temperature
       !      Namelist /HFBTHO_FEATURES/ collective_inertia, fission_fragments, pairing_regularization,
       !                                 localiation_functions
       !      Namelist /HFBTHO_NECK/ set_neck_constrain, neck_value
       !      Namelist /HFBTHO_DEBUG/ number_Gauss, number_Laguerre, number_Legendre, &
       !                              force_parity, print_time
       !
       n00_INI                   = number_of_shells       ! number of shells
       b0_INI                    = oscillator_length      ! oscillator length
       q_INI                     = basis_deformation      ! deformation beta_2 of the basis
       npr_INI(1)                = neutron_number         ! N
       npr_INI(2)                = proton_number          ! Z
       kindhfb_INI               = type_of_calculation    ! 1: HFB, -1: HFB+LN
       !
       b2_0                      = beta2_deformation      ! beta2 parameter of the initial WS solution
       b4_0                      = beta4_deformation      ! beta4 parameter of the initial WS solution
       !
       MAX_ITER_INI              = number_iterations      ! max number of iterations
       epsi_INI                  = accuracy               ! convergence of iterations
       inin_INI                  = restart_file           ! restart from file
       !
       skyrme_INI                = TRIM(functional)       ! functional
       Add_Pairing_INI           = add_initial_pairing    ! add pairing starting from file
       icou_INI                  = type_of_coulomb        ! coul: no-(0), dir.only-(1), plus exchange-(2)
       !
       set_pairing               = user_pairing           ! pairing is defined by user if .True.
       V0n_INI                   = vpair_n                ! pairing strength for neutrons
       V0p_INI                   = vpair_p                ! pairing strength for protons
       pwi_INI                   = pairing_cutoff         ! pairing q.p. cutoff
       cpv1_INI                  = pairing_feature        ! Type of pairing: volume, surface, mixed
       !
       nkblocase(1,:)            = neutron_blocking       ! config. of neutron blocked state
       nkblocase(2,:)            = proton_blocking        ! config. of proton blocked state
       !
       iLST_INI                  = switch_to_THO          ! 0:HO, -1:HO->THO, 1:THO
       keypj_INI                 = gauge_points           ! PNP: number of gauge points
       iproj_INI                 = projection_is_on       ! projecting on different nucleus
       npr1pj_INI                = delta_N                ! its neutron number
       npr2pj_INI                = delta_Z                ! its proton number
       !
       switch_on_temperature     = set_temperature        ! switches on temperature mode
       temper                    = temperature            ! value of the temperature
       !
       collective_inertia        = collective_inertia     ! calculate collective mass and zero-point energy
       fission_fragments         = fission_fragments      ! calculate fission fragment characteristics
       pairing_regularization    = pairing_regularization ! activates the regularization of the pairing force
       localization_functions    = localization_functions ! computes localization functions
       !
       neck_constraints          = set_neck_constrain     ! activate the constraint on the neck
       neckRequested             = neck_value             ! set the requested value for the neck
       dfrag_constraints         = set_dfrag_constrain     ! activate the constraint on the fragment distance
       dfragRequested            = dfrag_value             ! set the requested value for the fragment distance
       csi_constraints           = set_csi_constrain     ! activate the constraint on the mass asymmetry
       csiRequested              = csi_value             ! set the requested value for the mass asymmetry

       !
       ngh_INI                   = number_Gauss           ! number of Gauss-Hermite points for z-direction
       ngl_INI                   = number_Laguerre        ! number of Gauss-Laguerre points for rho-direction
       nleg_INI                  = number_Legendre        ! number of Gauss-Legendre points for Coulomb
       basis_HFODD_INI           = compatibility_HFODD    ! flag to enforce same basis as HFODD
       nstate_INI                = number_states          ! total number of states in basis
       Parity_INI                = force_parity           ! reflection symmetry
       IDEBUG_INI                = print_time             ! debug
       DO_FITT_INI               = .False.                ! calculates quantities for reg.optimization
       Print_HFBTHO_Namelist_INI = .False.                ! Print Namelist
       !
       ! Checking consistency of *_INI variables
       Call check_consistency
       If(ierror_flag.Ne.0) Then
          Print *,ierror_info(ierror_flag)
          Write(6,'("Terminating...")')
          Return
       End If
       ! Check if there is at least one constraint
       icount=0
       Do l=1,lambdaMax
          If(lambda_active(l).Gt.0) icount=icount+1
       End Do
       ! If there is at least one constraint, check if any breaks parity
       If(icount.Gt.0) Then
          jcount=0
          Do l=1,lambdaMax,2
             If(lambda_active(l).Gt.0) jcount=jcount+1
          End Do
          If(jcount.Gt.0) Parity_INI=.False.
       Else
          collective_inertia = .False.
       End If
       If(fission_fragments) Parity_INI=.False.
       !
       Call read_UNEDF_NAMELIST(skyrme_INI,noForce)
       ! If functional is used, projection automaticaly switched off
       If(noForce.Eq.0) iproj_INI=0
       ! Read parameters of the energy functionals from a file
#if(READ_FUNCTIONAL==1)
       call replace_functional
#endif
       !---------------------------------------------------------------------------
       ! GROUND STATE BLOCKING WALKER: blocking candidates are predefined
       ! by the parent nucleus and we block them one by one
       !---------------------------------------------------------------------------
       iblocase=0; bloqpdif=zero !blomax=0; blomax will be charged from the previous solution
       Do it=1,2
          If(nkblocase(it,1).Ne.0.And.nkblocase(it,2).Eq.0) Then
             If(it.Eq.1) Then
                iblocase(1)=iblocase(1)+1
                If(iblocase(1).Gt.blomax(1)) iblocase(1)=1
             Else
                If(iblocase(1).Le.1) iblocase(2)=iblocase(2)+1
             End If
             nkblo_INI(it,1)=Sign(iblocase(it),nkblocase(it,1))
             nkblo_INI(it,2)=0
          Else
             ! case of external blocking
             nkblo_INI(it,:)=nkblocase(it,:)
          Endif
       End Do
       !---------------------------------------------------------------------------
       ! MANUAL BLOCKING: manualBlocking=1 in the module (ocasionaly used)
       ! One types which level to be blocked referencing the parent nucleus
       !---------------------------------------------------------------------------
       If(manualBlocking.Ne.0) Then
          Write(*,'(a,5(1pg12.4))') 'Please print the number of the neutron level to block, num='
          Read(*,*) nkblo_INI(1,1)
          Write(*,'(a,5(1pg12.4))') 'Neutron blocked level num=',nkblo_INI(1,1)
          nkblo_INI(1,2)=0
          Write(*,'(a,5(1pg12.4))') 'Please print the number of the proton level to block, num='
          Read(*,*) nkblo_INI(2,1)
          Write(*,'(a,5(1pg12.4))') 'Proton blocked level num=',nkblo_INI(2,1)
          nkblo_INI(2,2)=0
       End If
       !--------------------------------------------------------------------
       ! Run the solver in all cases EVEN/ODDS, FITS/NO-FITS
       !--------------------------------------------------------------------
       Call HFBTHO_DFT_SOLVER
       !--------------------------------------------------------------------
       ! Calculate localization functions
       !--------------------------------------------------------------------
       If (localization_functions) Then
           Call localization()
       End If
       !--------------------------------------------------------------------
       ! Display error messages in case of problems
       !--------------------------------------------------------------------
#if(USE_MPI!=1)
       If (ierror_flag.Ne.0) Then
          Write(*,*)
          Write(*,'(a)') ' ERRORS IN HFBTHO_SOLVER'
          Do i=1,ierror_flag
             Write(*,'(a,i2,2x,a)') ' error_flag=',i,ierror_info(i)
          End Do
          Write(*,*)
       Else
          Write(*,*)
          Write(*,'(a)') ' HFBTHO_SOLVER ended without errors'
          Write(*,*)
       End If
#else
       If (ierror_flag.Ne.0) Then
          task_error = 1
          Write(*,*)
          Write(*,'(a33,i4)') ' ERRORS IN HFBTHO_SOLVER, process',mpi_taskid
          Write(*,'(2(a3,i4),a18,f6.1)')' Z=',Z_chain,' N=',N_chain, ' basis deformation', beta_deformation
          Do i=1,ierror_flag
             Write(*,'(a8,i4,a12,i2,2x,a)') ' process',mpi_taskid,', error_flag=',i,ierror_info(i)
          End Do
          Write(*,*)
          ierror_flag = 0
       End If
#endif
    !-------------------------------------------------------------------
    ! Mass table calculation
    !-------------------------------------------------------------------
#if(DO_MASSTABLE==1)
       if(nrows.gt.0) then
          write(*,'(1a27,1i3,1a4,1i3,a)') 'task '//id_string//' finished row '//row_string
#if(USE_MPI==1)
          call fill_out_vectors(icalc)
#else
          call fill_mass_table(icalc)
#endif
          icalc = icalc + 1
       endif
       Close(lfile) ! close the output
    enddo
    if(nrows.gt.0) then
#if(USE_MPI==1)
       call gather_results
#endif
       call print_mass_table
    endif
    if(team_rank.eq.0) then
       close(127)
    endif
#endif
    !-------------------------------------------------------------------
    ! Dripline calculation
    !-------------------------------------------------------------------
#if(DRIP_LINES==1)
             !Calculations without Lipkin-Nogami
             if(kindhfb_INI.gt.0) then
                Energy_chain = ehfb
             !Calculations with Lipkin-Nogami
             else
                Energy_chain = etot
             endif
             Close(lfile) ! close the output
          ! end loop over deformations
          enddo
          call find_minimum_energy
          separation_2N = Minimum_Energy_Prev - Minimum_Energy
          Minimum_Energy_Prev = Minimum_Energy
          if(separation_2N.lt.0._pr) then
             ! Drip line has been reached
             exit
          else
             ! Drip line has not been reached
             if(direction_sl(irow).eq.1) then
                if(N_chain.ge.310) exit
                N_chain = N_chain + 2
             else
                if(Z_chain.ge.120) exit
                Z_chain = Z_chain + 2
             endif
          endif
       ! end isotopic (or isotonic) chain
       enddo
       !team leader announces isotopic (or isotonic) chain finished
       if(team_rank.eq.0) then
          if(direction_sl(irow).eq.1) then
             write(*,128) 'team',team_color,' finished isotopic chain Z=',Z_chain,' at N=',N_chain
          else
             write(*,128) 'team',team_color,' finished isotonic chain N=',N_chain,' at Z=',Z_chain
          endif
       endif
    !end loop over nuclei inside the "stable line
    enddo
    !team leader closes bookkeeping file
    if(team_rank.eq.0) then
       close(127)
    endif
128 format(a4,i3,a27,i4,a6,i4)

#endif

#if(USE_MPI==1)
    call mpi_gather(task_error,1,mpi_integer,&
               task_error_gthr,1,mpi_integer,0,MPI_COMM_WORLD,&
               ierr_mpi)
    if(mpi_taskid.eq.0) then
       if(sum(task_error_gthr).eq.0) then
          Write(*,*)
          Write(*,'(a)') ' HFBTHO_SOLVER ended without errors'
          Write(*,*)
       endif
    endif
#endif
    !--------------------------------------------------------------------
    ! For single-nucleus calculation (so far) records the HFB solution
    ! in a format readable by the pnFAM code of Mustonen & Shafer
    !--------------------------------------------------------------------
#if(DO_MASSTABLE==0)
#if(USE_QRPA==1)
If (lout.Lt.lfile) Close(lfile) ! close the output
    call save_HFBTHO_solution  ! for pnFAM
#endif
#endif
    !
  End Subroutine Main_Program
