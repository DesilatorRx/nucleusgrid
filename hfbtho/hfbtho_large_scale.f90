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

#if(DRIP_LINES==1 || DO_MASSTABLE==1)

! ==================================================================== !
!                                                                      !
!                  LARGE-SCALE CALCULATION PACKAGE                     !
!                                                                      !
! ==================================================================== !
!----------------------------------------------------------------------
!>  This module contains various utility functions to perform large
!>  scale calculations with HFBTHO, be it full nuclear chart from
!>  dripline to dripline or a subset of nuclei with different
!>  deformations. Calculations of this type require USE_MPI=1.
!>  Refer to hfbtho_main.f90 for further information.
!>
!>  @author
!>  Rodrigo Navarro Perez
!----------------------------------------------------------------------
!  Subroutines: - read_HFBTHO_MassTable
!               - allocate_mass_table
!               - fill_mass_table
!               - print_mass_table
!----------------------------------------------------------------------
Module HFBTHO_large_scale
  Use HFBTHO_utilities
  Use HFBTHO
  Implicit None
Contains

  !=======================================================================
  !>  Opens and reads the 'hfbtho_MASSTABLE.dat' file that indicates
  !>  the nuclei (with different deformations) to be calculated with
  !>  HFBTHO
  !>    - The first integer read from the file gives the number of
  !>      calculations to be made (called nrows along the source code)
  !>    - Each row contains the number of protons (Z), number of neutrons
  !>      (N), requested value for the quadrupole moment (Q20) and basis
  !>      deformation (beta)
  !>    - If nrows is 0 (or negative) the rest of the file will be
  !>      ignored and a regular HFBTHO calculation will be made with the
  !>      parameters given in the 'hfbtho_NAMELIST.dat' file.
  !>    - If nrows is greater than zero the first nrows calculations will
  !>      be made with the parameters given in each row. Other input
  !>      parameters will be take from the 'hfbtho_NAMELIST.dat' file.
  !>  nrows should not be larger than the actual number of rows on the
  !>  file.
  !=======================================================================
  subroutine read_HFBTHO_MassTable
    implicit none
    integer(ipr) :: lmasstable=15,i,j
    open(lmasstable,file='hfbtho_MASSTABLE.dat')
    read(lmasstable,*)
    read(lmasstable,*) nRows
    read(lmasstable,*)
    allocate(   Z_masstable(0:nRows))
    allocate(   N_masstable(0:nRows))
    allocate( Q20_masstable(0:nRows))
    allocate(beta_masstable(0:nRows))
       Z_masstable(0) =  proton_number
       N_masstable(0) = neutron_number
     Q20_masstable(0) = expectation_values(2)
    beta_masstable(0) = basis_deformation
    do i=1,nRows
       read(lmasstable,*) Z_masstable(i),N_masstable(i),Q20_masstable(i)&
            ,beta_masstable(i)
    enddo
    close(lmasstable)
  end subroutine read_HFBTHO_MassTable
  !=======================================================================
  !>  Allocates the arrays necessary to construct the output mass table
  !=======================================================================
  subroutine allocate_mass_table
    implicit none
    allocate(ierrors_out(1:nrows))
    allocate(Z_out(1:nrows))
    allocate(N_out(1:nrows))
    allocate(Q20_in(1:nrows))
    allocate(beta_in(1:nrows))
    allocate(E_HFB_out(1:nrows))
    allocate(Q20Z_out(1:nrows))
    allocate(Q20N_out(1:nrows))
    allocate(Q20T_out(1:nrows))
  end subroutine allocate_mass_table
  !=======================================================================
  !>  The results from the HFBTHO calculation are stored in arrays.
  !>  This subroutine is used in the case that the mass table calculated
  !>  without MPI parallelism
  !=======================================================================
  subroutine fill_mass_table(icalc)
    implicit none
    integer, intent(in) :: icalc !< Calculations performed so far
    ierrors_out(icalc+1) = ierror_flag
    Z_out(icalc+1) = Z_masstable(iRow)
    N_out(icalc+1) = N_masstable(iRow)
    Q20_in(icalc+1) = Q20_masstable(iRow)
    beta_in(icalc+1) = beta_masstable(iRow)
    E_HFB_out(icalc+1) = ehfb
    Q20Z_out(icalc+1) = qmoment(2,1)
    Q20N_out(icalc+1) = qmoment(2,2)
    Q20T_out(icalc+1) = qmoment(2,3)
  end subroutine fill_mass_table
  !=======================================================================
  !>  The output mass table is written to the 'MassTableOut.dat' file
  !>  and on screen
  !=======================================================================
  subroutine print_mass_table
    implicit none
    integer :: i
!!#if(DO_MASSTABLE==1)
    if(mpi_taskid.eq.0) then
!!#endif
       open(100,file='MassTableOut.dat')
       do i = 1,nrows
          write(  *,'(4i5,6f15.8)') i,ierrors_out(i),Z_out(i),N_out(i),Q20_in(i),beta_in(i),E_HFB_out(i),Q20Z_out(i),Q20N_out(i),Q20T_out(i)
          write(100,'(4i5,6f15.8)') i,ierrors_out(i),Z_out(i),N_out(i),Q20_in(i),beta_in(i),E_HFB_out(i),Q20Z_out(i),Q20N_out(i),Q20T_out(i)
       enddo
       close(100)
!!#if(DO_MASSTABLE==1)
    endif
!!#endif
  end subroutine print_mass_table
  !=======================================================================
  !>  This subroutine reads the location of the valley of stability as
  !>  well as the direction (proton or neutron) in which dripline
  !>  calculations should proceed.
  !=======================================================================
  subroutine read_HFBTHO_StableLine
    implicit none
    integer(ipr) :: ldripline=15,i,j
    open(ldripline,file='hfbtho_STABLELINE.dat')
    read(ldripline,*)
    read(ldripline,*) nRows
    read(ldripline,*)
    allocate(Z_stable_line(0:nRows))
    allocate(N_stable_line(0:nRows))
    allocate(direction_sl(0:nRows))
    Z_stable_line(0) =  proton_number
    N_stable_line(0) = neutron_number
    do i=1,nRows
       read(ldripline,*) Z_stable_line(i),N_stable_line(i),direction_sl(i)
    enddo
    close(ldripline)
  end subroutine read_HFBTHO_StableLine

end module HFBTHO_large_scale

#endif



#if(USE_MPI==1)

!----------------------------------------------------------------------
!>  This module contain utility functions and subroutines to
!>  exchange information among MPI processes involved in a mass table
!>  or dripline calculation.
!>
!> @author
!> Rodrigo Navarro Perez
!----------------------------------------------------------------------
!  Subroutines: - Create_MPI_teams
!               - Construct_Vectors
!               - broadcast_vectors
!               - allocate_mpi_vectors
!               - allocate_out_vectors
!               - Deconstruct_Vectors
!               - gather_results
!               - fill_out_vectors(icalc)
!----------------------------------------------------------------------
Module HFBTHO_mpi_communication
  Use HFBTHO_utilities
  Use HFBTHO
  Implicit None
  Include 'mpif.h'
Contains

  !=======================================================================
  !>  Creates new communicators for teams of MPI processes to calculate
  !>  each nucleus with different deformations. The processes are
  !>  distributed as evenly as possible and each team will have at
  !>  most as many members as the number of deformations that will be used
  !>  (having more members would be inefficient)
  !>
  !>  A group and communicator are created with one member of each team
  !>  that will act as team leader and communicate with a 'super-leader'
  !>  to distribute the nucleus to be calculated by each group
  !=======================================================================
  subroutine Create_MPI_Teams
    implicit none
    integer :: i,GROUP_leaders,GROUP_world
    integer, allocatable  :: leaders_ranks(:)
    integer :: leaders_size, leaders_rank
    !Determine number of teams (each team must have at most 11 members)
    number_teams = (mpi_size-1)/number_deformations + 1
    team_color = mod(mpi_taskid,number_teams)
    !Create Team communicators by splitting the world communicator
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,team_color,mpi_taskid,COMM_team,&
         ierr_mpi)
    !Get size team and rank within the team
    call MPI_COMM_SIZE(COMM_team, team_size, ierr_mpi)
    call MPI_COMM_RANK(COMM_team, team_rank, ierr_mpi)

    ! !The processes ranked first in the world group will be team leaders
    ! allocate(leaders_ranks(0:number_teams-1))
    ! do i = 0,number_teams-1
    !    leaders_ranks(i) = i
    ! enddo
    ! ! Get a handle on the world group
    ! call MPI_COMM_GROUP(MPI_COMM_WORLD,GROUP_world,ierr_mpi)
    ! ! Put leaders into the group and create leaders communicator
    ! call MPI_GROUP_INCL(GROUP_world,number_teams,leaders_ranks,&
    !      GROUP_leaders,ierr_mpi)
    ! call MPI_COMM_CREATE(MPI_COMM_WORLD,GROUP_leaders,COMM_leaders,&
    !      ierr_mpi)

    ! !For leaders only,
    ! !get number of leathers (must be equal to number of teams) and
    ! !rank within the leaders (must be equal to the rank within the world)
    ! if(COMM_leaders.ne.MPI_COMM_NULL) then
    !    call MPI_COMM_SIZE(COMM_leaders, leaders_size, ierr_mpi)
    !    call MPI_COMM_RANK(COMM_leaders, leaders_rank, ierr_mpi)
    ! else
    !    leaders_size = -1
    !    leaders_rank = -1
    ! endif
  end subroutine Create_MPI_Teams
  !=======================================================================
  !>  After the master process has read the 'hfbtho_NAMELIST.dat' and
  !>  'hfbtho_MASSTABLE.dat' files, it constructs (allocates and fills out)
  !>  a set of vectors to be broadcast to the other processes.
  !>  - vector_int_mpi contains all the integers to broadcasted
  !>  - vector_real_mpi contains all the reals to broadcasted
  !>  - vector_logicals_mpi contains all the logicals to broadcasted
  !>
  !>  Vector_sizes contains 4 integers, the size of the three previous
  !>  arrays and the number of rows on the mass table.
  !=======================================================================
  subroutine Construct_Vectors
    implicit none
    integer :: i
    allocate(vector_sizes(1:4))
    vector_sizes(1) = 27 + 2*lambdamax +  n_int_masstable_in*nRows
    vector_sizes(2) = 13 +   lambdamax + n_real_masstable_in*nRows
    vector_sizes(3) = 12
    vector_sizes(4) =  nRows
    allocate( vector_int_mpi(1:vector_sizes(1)))
    allocate( vector_real_mpi(1:vector_sizes(2)))
    allocate( vector_log_mpi(1:vector_sizes(3)))
    ! Inputs of type 'integer'
    vector_int_mpi( 1) = number_of_shells
    vector_int_mpi( 2) = proton_number
    vector_int_mpi( 3) = neutron_number
    vector_int_mpi( 4) = type_of_calculation
    vector_int_mpi( 5) = number_iterations
    vector_int_mpi( 6) = type_of_coulomb
    vector_int_mpi( 7) = restart_file
    vector_int_mpi( 8) = projection_is_on
    vector_int_mpi( 9) = gauge_points
    vector_int_mpi(10) = delta_Z
    vector_int_mpi(11) = delta_N
    vector_int_mpi(12) = switch_to_THO
    vector_int_mpi(13) = number_Gauss
    vector_int_mpi(14) = number_Laguerre
    vector_int_mpi(15) = number_Legendre
    vector_int_mpi(16) = number_states
    vector_int_mpi(17) = print_time
    do i = 1,5
       vector_int_mpi(17+i) =  proton_blocking(i)
       vector_int_mpi(22+i) = neutron_blocking(i)
    enddo
    do i = 1,lambdamax
       vector_int_mpi(27+i) = lambda_values(i)
       vector_int_mpi(27+lambdamax+i) = lambda_active(i)
    enddo
#if(DO_MASSTABLE==1)
    do i = 1,nRows
       vector_int_mpi(27+2*lambdamax+n_int_masstable_in*(i-1)+1) = Z_masstable(i)
       vector_int_mpi(27+2*lambdamax+n_int_masstable_in*(i-1)+2) = N_masstable(i)
    enddo
#endif
#if(DRIP_LINES==1)
    do i = 1,nRows
       vector_int_mpi(27+2*lambdamax+n_int_masstable_in*(i-1)+1) = Z_stable_line(i)
       vector_int_mpi(27+2*lambdamax+n_int_masstable_in*(i-1)+2) = N_stable_line(i)
       vector_int_mpi(27+2*lambdamax+n_int_masstable_in*(i-1)+3) = direction_sl(i)
    enddo
#endif
    ! Inputs of type 'real'
    vector_real_mpi( 1) = oscillator_length
    vector_real_mpi( 2) = basis_deformation
    vector_real_mpi( 3) = beta2_deformation
    vector_real_mpi( 4) = beta4_deformation
    vector_real_mpi( 5) = accuracy
    vector_real_mpi( 6) = temperature
    vector_real_mpi( 7) = vpair_n
    vector_real_mpi( 8) = vpair_p
    vector_real_mpi( 9) = pairing_cutoff
    vector_real_mpi(10) = pairing_feature
    vector_real_mpi(11) = neck_value
    vector_real_mpi(12) = dfrag_value
    vector_real_mpi(13) = csi_value
    do i = 1,lambdamax
       vector_real_mpi(11+i) = expectation_values(i)
    enddo
#if(DO_MASSTABLE==1)
    do i = 1,nRows
       vector_real_mpi(11+lambdamax+n_real_masstable_in*(i-1)+1) = Q20_masstable(i)
       vector_real_mpi(11+lambdamax+n_real_masstable_in*(i-1)+2) = beta_masstable(i)
    enddo
#endif
    ! Inputs of type 'logical'
    vector_log_mpi( 1) = add_initial_pairing
    vector_log_mpi( 2) = set_temperature
    vector_log_mpi( 3) = collective_inertia
    vector_log_mpi( 4) = fission_fragments
    vector_log_mpi( 5) = pairing_regularization
    vector_log_mpi( 6) = localization_functions
    vector_log_mpi( 7) = set_neck_constrain
    vector_log_mpi( 8) = compatibility_HFODD
    vector_log_mpi( 9) = force_parity
    vector_log_mpi(10) = user_pairing
    vector_log_mpi(11) = set_dfrag_constrain
    vector_log_mpi(12) = set_csi_constrain
    
  end subroutine Construct_Vectors
  !=======================================================================
  !>  The master process broadcasts first the vector_sizes arrays, the
  !>  other processes allocate the other arrays (int, reals and logicals)
  !>  then the master process broadcasts integers, reals, logicals and
  !>  a single string
  !=======================================================================
  subroutine broadcast_vectors
    implicit none
    if(mpi_taskid.gt.0) allocate(vector_sizes(1:4))
    call mpi_bcast(vector_sizes,4,mpi_integer,0,mpi_comm_world,ierr_mpi)
    if(mpi_taskid.gt.0) then
       nRows = vector_sizes(4)
       call allocate_mpi_vectors
    endif
    call mpi_bcast( vector_int_mpi,vector_sizes(1),mpi_integer,0,mpi_comm_world,ierr_mpi)
    call mpi_bcast(vector_real_mpi,vector_sizes(2),mpi_double_precision,0,mpi_comm_world,ierr_mpi)
    call mpi_bcast( vector_log_mpi,vector_sizes(3),mpi_logical,0,mpi_comm_world,ierr_mpi)
    call mpi_bcast(functional,30,mpi_character,0,mpi_comm_world,ierr_mpi)
  end subroutine broadcast_vectors
  !=======================================================================
  !>  Non-master processes allocate the arrays that will be received
  !>  from the master's broadcast.
  !=======================================================================
  subroutine allocate_mpi_vectors
    implicit none
    allocate( vector_int_mpi(1:vector_sizes(1)))
    allocate(vector_real_mpi(1:vector_sizes(2)))
    allocate( vector_log_mpi(1:vector_sizes(3)))
#if(DO_MASSTABLE==1)
    allocate(   Z_masstable(0:nRows))
    allocate(   N_masstable(0:nRows))
    allocate( Q20_masstable(0:nRows))
    allocate(beta_masstable(0:nRows))
#endif
#if(DRIP_LINES==1)
    allocate(Z_stable_line(0:nRows))
    allocate(N_stable_line(0:nRows))
    allocate(direction_sl(0:nRows))
#endif
  end subroutine allocate_mpi_vectors
  !=======================================================================
  !>  After receiving the master's broadcast, the non-master processes
  !>  'deconstruct' the arrays by setting the parameters to be used
  !>  in the HFBTHO calculations
  !=======================================================================
  subroutine Deconstruct_Vectors
    implicit none
    integer :: i
    ! Inputs of type 'integer'
    number_of_shells    = vector_int_mpi( 1)
    proton_number       = vector_int_mpi( 2)
    neutron_number      = vector_int_mpi( 3)
    type_of_calculation = vector_int_mpi( 4)
    number_iterations   = vector_int_mpi( 5)
    type_of_coulomb     = vector_int_mpi( 6)
    restart_file        = vector_int_mpi( 7)
    projection_is_on    = vector_int_mpi( 8)
    gauge_points        = vector_int_mpi( 9)
    delta_Z             = vector_int_mpi(10)
    delta_N             = vector_int_mpi(11)
    switch_to_THO       = vector_int_mpi(12)
    number_Gauss        = vector_int_mpi(13)
    number_Laguerre     = vector_int_mpi(14)
    number_Legendre     = vector_int_mpi(15)
    number_states       = vector_int_mpi(16)
    print_time          = vector_int_mpi(17)
    do i = 1,5
       proton_blocking(i)  = vector_int_mpi(17+i)
       neutron_blocking(i) = vector_int_mpi(22+i)
    enddo
    do i = 1,lambdamax
       lambda_values(i) = vector_int_mpi(27+i)
       lambda_active(i) = vector_int_mpi(27+lambdamax+i)
    enddo
#if(DO_MASSTABLE==1)
    do i = 1,nRows
       Z_masstable(i) = vector_int_mpi(27+2*lambdamax+n_int_masstable_in*(i-1)+1)
       N_masstable(i) = vector_int_mpi(27+2*lambdamax+n_int_masstable_in*(i-1)+2)
    enddo
#endif
#if(DRIP_LINES==1)
    do i = 1,nRows
       Z_stable_line(i) = vector_int_mpi(27+2*lambdamax+n_int_masstable_in*(i-1)+1)
       N_stable_line(i) = vector_int_mpi(27+2*lambdamax+n_int_masstable_in*(i-1)+2)
       direction_sl(i)  = vector_int_mpi(27+2*lambdamax+n_int_masstable_in*(i-1)+3)
    enddo
#endif
    ! Inputs of type 'real'
    oscillator_length = vector_real_mpi( 1)
    basis_deformation = vector_real_mpi( 2)
    beta2_deformation = vector_real_mpi( 3)
    beta4_deformation = vector_real_mpi( 4)
    accuracy          = vector_real_mpi( 5)
    temperature       = vector_real_mpi( 6)
    vpair_n           = vector_real_mpi( 7)
    vpair_p           = vector_real_mpi( 8)
    pairing_cutoff    = vector_real_mpi( 9)
    pairing_feature   = vector_real_mpi(10)
    neck_value        = vector_real_mpi(11)
    dfrag_value       = vector_real_mpi(12)
    csi_value         = vector_real_mpi(13)
    
    
    do i = 1,lambdamax
       expectation_values(i) =  vector_real_mpi(12+i)
    enddo
#if(DO_MASSTABLE==1)
    do i = 1,nRows
        Q20_masstable(i) = vector_real_mpi(13+lambdamax+n_real_masstable_in*(i-1)+1)
       beta_masstable(i) = vector_real_mpi(13+lambdamax+n_real_masstable_in*(i-1)+2)
    enddo
#endif
    ! Inputs of type 'logical'
    add_initial_pairing    = vector_log_mpi( 1)
    set_temperature        = vector_log_mpi( 2)
    collective_inertia     = vector_log_mpi( 3)
    fission_fragments      = vector_log_mpi( 4)
    pairing_regularization = vector_log_mpi( 5)
    localization_functions = vector_log_mpi( 6)
    set_neck_constrain     = vector_log_mpi( 7)
    compatibility_HFODD    = vector_log_mpi( 8)
    force_parity           = vector_log_mpi( 9)
    user_pairing           = vector_log_mpi(10)
    set_dfrag_constrain    = vector_log_mpi(11)
    set_csi_constrain      = vector_log_mpi(12)
  end subroutine Deconstruct_Vectors
  !=======================================================================
  !>  Every process allocates the vectors that will be gathered by
  !>  the master.
  !=======================================================================
  subroutine allocate_out_vectors
    implicit none
    deallocate(vector_sizes,vector_int_mpi,vector_real_mpi,vector_log_mpi)
    nrows_task = nrows/mpi_size
    if(mpi_taskid.gt.0.and.mpi_taskid.le.mod(nrows,mpi_size)) then
       nrows_task = nrows_task + 1
    endif
    allocate(vector_sizes(1:2))
    vector_sizes(1) = n_int_masstable_out*nrows_task
    vector_sizes(2) = n_real_masstable_out*nrows_task
    allocate( vector_int_mpi(1:vector_sizes(1)))
    allocate(vector_real_mpi(1:vector_sizes(2)))
  end subroutine allocate_out_vectors
  !=======================================================================
  !>  Get and broadcast the value of the minimum energy for a given nucleus
  !=======================================================================
  subroutine find_minimum_energy()
    implicit none
    !team leader gathers energies and looks for the minimum
    call mpi_gather(Energy_Chain,1,mpi_double_precision,&
                    Energy_Chain_gthr,1,mpi_double_precision,0,COMM_team,&
                    ierr_mpi)
    if(team_rank.eq.0) then
       Minimum_Energy = minval(Energy_Chain_gthr,1)
    endif
    !team leader broadcasts the minimum
    call mpi_bcast(Minimum_Energy,1,mpi_double_precision,0,COMM_team,ierr_mpi)
  end subroutine find_minimum_energy
  !=======================================================================
  !>  After each HFBTHO calculation, each process saves the results in
  !>  vectors that will be gathered later by the master.
  !=======================================================================
  subroutine fill_out_vectors(icalc)
    implicit none
    integer, intent(in) :: icalc
       vector_int_mpi(n_int_masstable_out*icalc+1) = ierror_flag
       vector_int_mpi(n_int_masstable_out*icalc+2) = Z_masstable(iRow)
       vector_int_mpi(n_int_masstable_out*icalc+3) = N_masstable(iRow)
       vector_real_mpi(n_real_masstable_out*icalc+1) = Q20_masstable(iRow)
       vector_real_mpi(n_real_masstable_out*icalc+2) = beta_masstable(iRow)
       if(kindhfb_INI.gt.0) then
          vector_real_mpi(n_real_masstable_out*icalc+3) = ehfb
       else
          vector_real_mpi(n_real_masstable_out*icalc+3) = etot
       endif
       vector_real_mpi(n_real_masstable_out*icalc+4) = qmoment(2,1)
       vector_real_mpi(n_real_masstable_out*icalc+5) = qmoment(2,2)
       vector_real_mpi(n_real_masstable_out*icalc+6) = qmoment(2,3)
  end subroutine fill_out_vectors
  !=======================================================================
  !>  The master process gathers the result from all calculations
  !>  and fills out the resulting mass table including number of protons,
  !>  number of neutrons, input quadrupole moment, basis deformation, HFB
  !>  energy, output quadrupole moment from protons, neutrons and total
  !=======================================================================
  subroutine gather_results
    implicit none
    integer :: i,j
    if(mpi_taskid.eq.0) then
       allocate(vector_sizes_gthr(1:2*mpi_size))
       allocate(vector_sizes_int_gthr(0:mpi_size-1))
       allocate(vector_sizes_real_gthr(0:mpi_size-1))
       allocate(vector_disp_int_gthr(0:mpi_size-1))
       allocate(vector_disp_real_gthr(0:mpi_size-1))
    endif
    call mpi_gather(vector_sizes,2,mpi_integer,vector_sizes_gthr,2,&
         mpi_integer,0,mpi_comm_world,ierr_mpi)

    if(mpi_taskid.eq.0) then
       disp_int = 0
       disp_real = 0
       do i = 0,mpi_size-1
          vector_sizes_int_gthr(i)  = vector_sizes_gthr(2*i+1)
          vector_sizes_real_gthr(i) = vector_sizes_gthr(2*i+2)
          if(i.eq.0) cycle
          vector_disp_int_gthr(i)  = disp_int
          vector_disp_real_gthr(i) = disp_real
          disp_int = disp_int + vector_sizes_gthr(2*i+1)
          disp_real = disp_real + vector_sizes_gthr(2*i+2)
       enddo
       vector_disp_int_gthr(0)  = disp_int
       vector_disp_real_gthr(0) = disp_real
       allocate( vector_int_gthr(1:sum( vector_sizes_int_gthr)))
       allocate(vector_real_gthr(1:sum(vector_sizes_real_gthr)))
    endif

    call mpi_gatherv(vector_int_mpi,vector_sizes(1),mpi_integer,&
         vector_int_gthr,vector_sizes_int_gthr,vector_disp_int_gthr,&
         mpi_integer,0,mpi_comm_world,ierr_mpi)

    call mpi_gatherv(vector_real_mpi,vector_sizes(2),&
         mpi_double_precision,vector_real_gthr,vector_sizes_real_gthr,&
         vector_disp_real_gthr,mpi_double_precision,0,mpi_comm_world,&
         ierr_mpi)

    if(mpi_taskid.eq.0) then
       allocate(ierrors_out(1:nrows))
       allocate(Z_out(1:nrows))
       allocate(N_out(1:nrows))
       allocate(Q20_in(1:nrows))
       allocate(beta_in(1:nrows))
       allocate(E_HFB_out(1:nrows))
       allocate(Q20Z_out(1:nrows))
       allocate(Q20N_out(1:nrows))
       allocate(Q20T_out(1:nrows))
       do i = 1,mpi_size
          do j = 1,vector_sizes_int_gthr(mod(i,mpi_size))/n_int_masstable_out
             ierrors_out(mpi_size*(j-1)+i) = vector_int_gthr(vector_disp_int_gthr(mod(i,mpi_size))+n_int_masstable_out*(j-1)+1)
                   Z_out(mpi_size*(j-1)+i) = vector_int_gthr(vector_disp_int_gthr(mod(i,mpi_size))+n_int_masstable_out*(j-1)+2)
                   N_out(mpi_size*(j-1)+i) = vector_int_gthr(vector_disp_int_gthr(mod(i,mpi_size))+n_int_masstable_out*(j-1)+3)
          enddo
          do j = 1,vector_sizes_real_gthr(mod(i,mpi_size))/n_real_masstable_out
                Q20_in(mpi_size*(j-1)+i) = vector_real_gthr(vector_disp_real_gthr(mod(i,mpi_size))+n_real_masstable_out*(j-1)+1)
               beta_in(mpi_size*(j-1)+i) = vector_real_gthr(vector_disp_real_gthr(mod(i,mpi_size))+n_real_masstable_out*(j-1)+2)
             E_HFB_out(mpi_size*(j-1)+i) = vector_real_gthr(vector_disp_real_gthr(mod(i,mpi_size))+n_real_masstable_out*(j-1)+3)
              Q20Z_out(mpi_size*(j-1)+i) = vector_real_gthr(vector_disp_real_gthr(mod(i,mpi_size))+n_real_masstable_out*(j-1)+4)
              Q20N_out(mpi_size*(j-1)+i) = vector_real_gthr(vector_disp_real_gthr(mod(i,mpi_size))+n_real_masstable_out*(j-1)+5)
              Q20T_out(mpi_size*(j-1)+i) = vector_real_gthr(vector_disp_real_gthr(mod(i,mpi_size))+n_real_masstable_out*(j-1)+6)
          enddo
       enddo
    endif

  end subroutine gather_results

end module HFBTHO_mpi_communication
#endif
