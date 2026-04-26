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
!                      INPUT/OUTPUT PACKAGE                            !
!                                                                      !
! ==================================================================== !

!-------------------------------------------------------------------
!> This module contains all routines dealing with input (reading
!> data from disk) and output (writing data on disk).
!>
!> @author
!> Nicolas Schunck
!----------------------------------------------------------------------
!  Subroutines: - FileLabels()
!               - inout()
!               - write_version()
!               - read_data_old()
!               - read_data()
!               - write_data()
!               - write_data_old()
!               - blosort()
!  Functions: - version_number()
!             - check_file()
!----------------------------------------------------------------------
Module HFBTHO_io

  Use HFBTHO_utilities
  Use HFBTHO
  Use hfbtho_gogny
  Use HFBTHO_collective

  Implicit None

  Integer(ipr), PUBLIC, SAVE :: VERSION_DATA = 3 !< Version number of the binary file. VERSION_DATA=2 corresponds to version
                                                 !< 2.00d of HFBTHO; VERSION_DATA=1 to version 1.66 (not implemented)
  Character(Len=50), PUBLIC, SAVE :: welfile !< Name of the binary file

Contains
  !=======================================================================
  !> Subroutine that defines a filename, e.g., filelabel='s070_040'
  !=======================================================================
  Subroutine FileLabels(NPRI,ININL,FILELABEL)
    Implicit None
    Integer(ipr) :: it,ininabs,ininl,nprt(2),NPRI(2)
    Character(1)  :: sinin
    Character(3)  :: snpr(2)
    Character(8)  :: filelabel
    !
    ininabs=iabs(ininl)
    If(ininabs.Eq.4.or.ininabs.Eq.400) sinin='t'
    If(ininabs.Eq.3.or.ininabs.Eq.300) sinin='o'
    If(ininabs.Eq.2.or.ininabs.Eq.200) sinin='p'
    If(ininabs.Eq.1.or.ininabs.Eq.100) sinin='s'
    !
    nprt=npri
    Do it=itmin,itmax
       If(npri(it).Ne.2*(npri(it)/2)) nprt(it)=nprt(it)+iparenti(it)  !iparent=-/+ means particles/holes
       write(snpr(it),'(i3.3)') nprt(it)
       If(nprt(it).Lt.10  ) Then
          Write(snpr(it),'(a2,i1)') '00',nprt(it)
       Else
          If(nprt(it).Lt.100 ) Then
             Write(snpr(it),'(a1,i2)') '0',nprt(it)
          Else
             Write(snpr(it),'(i3)') nprt(it)
          End If
       End If
    End Do
    !
    Write(filelabel,'(a1,a3,a1,a3)')  sinin,snpr(1),'_',snpr(2)
    !
#if(DO_MASSTABLE==1)
    If(iLST1.Le.0) Write(welfile,'(a8,a5,a4)')  FILELABEL,row_string,'.hel'
    If(iLST1.Gt.0) Write(welfile,'(a8,a5,a4)')  FILELABEL,row_string,'.tel'
#elif(DRIP_LINES==1)
    If(iLST1.Le.0) Write(welfile,'(a8,a11,a4)')  FILELABEL,row_string,'.hel'
    If(iLST1.Gt.0) Write(welfile,'(a8,a11,a4)')  FILELABEL,row_string,'.tel'
#else
    If(iLST1.Le.0) Write(welfile,'("hfbtho_output.hel")')
    If(iLST1.Gt.0) Write(welfile,'("hfbtho_output.tel")')
#endif
    !
  End Subroutine FileLabels
  !=======================================================================
  !> This subroutine is the central interface to all input/output
  !> operations in HFBTHO. Depending on the value of its input parameter
  !> it will either read a binary file to get input data, or it will write
  !> the current data in a binary file. In the reading phase, the code
  !> tries to read a version number. If this read is successful, it is
  !> assumed the binary file is in version VERSION_DATA
  !=======================================================================
  Subroutine inout(is,iexit)
    Implicit None
    Integer(ipr), INTENT(IN) :: is !> Integer specifying if the file is read (=1) or written (=2)
    Integer(ipr), INTENT(INOUT) :: iexit !> Integer giving the exit status of the routine (0: OK, >0: not OK)
    Character(50) :: action
    Character(8) :: filelabel
    Integer(ipr) :: VERSION_read
    ! label organization
    Call FileLabels(NPR,ININ,FILELABEL)
    If(ierror_flag.Ne.0) Return
    !---------------------------------------------------------------------
    ! Read data to start the calculation
    !---------------------------------------------------------------------
    If(is.Eq.1) Then
       action = 'Read'
       ! Start from scratch
       If(inin.Gt.0) Then
          iexit=1
       Else
          ! Check the file is valid
          iexit = check_file(welfile, action)
          If(iexit==0) Then
             ! Read the version number of the file
             VERSION_read = version_number(welfile)
             ! If version is current, try to read the data
             If(VERSION_read==3) Then
                Call read_data(iexit)
             Else
                ! Close and reopen the file
                iexit = check_file(welfile, action)
                ! Try to read the file assuming it is in the old format
                Call read_data_old(iexit)
             End If
          Else
             iexit=1
          End If
       End If
    End If
    !---------------------------------------------------------------------
    ! Write data on disk
    !---------------------------------------------------------------------
    If(is.Eq.2.And.iasswrong(3).Eq.0) Then
       action = 'Write'
       ! Check the file is valid
       iexit = check_file(welfile, action)
       If(iexit==0) Then
          Call write_version()
          Call write_data()
       Else
          iexit=1
       End If
    End If
  End Subroutine inout
  !=======================================================================
  !> Subroutine that just writes the version number on disk
  !=======================================================================
  Subroutine write_version()
    Implicit None
    Write(lwou) VERSION_DATA
  End Subroutine write_version
  !=======================================================================
  !> Function that reads and returns the version number from the input
  !> file. Returns -1 in case the version cannot be read.
  !=======================================================================
  Integer(ipr) Function version_number(filename)
    Implicit None
    Character(50), INTENT(IN) :: filename
    Integer(ipr) :: ierr,version_read,iexit
    version_read = 1
    Read(lwin,IOSTAT=ierr) version_read
    If(ierr.NE.0) Then
       version_read = -1
    End If
    version_number = version_read
  End Function version_number
  !=======================================================================
  !> Function that checks the status of the file. If the file exists and
  !> can be opened, it is opened and the functions returns an exit status
  !> of 0. If the file should be read but does not exist, the function
  !> returns 1; if the file should be written but does not exist, the
  !> function opens a new file and returns 0.
  !=======================================================================
  Integer(ipr) Function check_file(filename, action)
    Implicit None
    Character(50), INTENT(IN) :: filename !> Name of the file to check
    Character(50), INTENT(IN) :: action !> Equal to 'Read' or 'Write', defines what to do if the file does not exist
    !
    Logical :: file_exists,file_opened
    Integer(ipr) :: ierr,iexit
    iexit = 0
    file_exists=.False.; inquire(file=filename, exist=file_exists); ierr=0
    If(Trim(action) == 'Read') Then
       If(file_exists) Then
          file_opened=.False.; inquire(unit=lwin, opened=file_opened)
          If(file_opened) Then
             Close(lwin)
          End If
          Open(lwin,file=filename,status='old',form='unformatted',IOSTAT=ierr)
          If(ierr.NE.0) Then
             iexit = 1
          End If
       Else
          iexit = 1
       End If
    End If
    If(Trim(action) == 'Write') Then
       If(file_exists) Then
          file_opened=.False.; inquire(unit=lwou, opened=file_opened)
          If(file_opened) Then
             Close(lwou)
          End If
          Open(lwou,file=filename,status='old',form='unformatted',IOSTAT=ierr)
          If(ierr.NE.0) Then
             iexit = 1
          End If
       Else
          Open(lwou,file=filename,status='new',form='unformatted',IOSTAT=ierr)
          iexit = 0
       End If
    End If
    check_file = iexit
  End Function check_file
  !=======================================================================
  !> Subroutine that reads the data for 'old' binary files corresponding
  !> to HFBTHO version 200d
  !=======================================================================
  Subroutine read_data_old(iexit)
    Implicit None
    Integer(ipr), INTENT(INOUT) :: iexit
    Integer(ipr)  :: is,iw,n1,n2,nd,ib,bloall1,lambdaMax1,ierr,counterLine
    Real(pr)      :: tz1(2),b01,bz1,bp1,beta1,v0r(2),v1r(2),pwir
    Integer(ipr)  :: npr1,npr11,ngh1,ngl1,n001,nt1,ibro,i
    Integer(ipr)  :: ntx1,nb1,nhhdim1,NLANSA0,NLANSA1,NZA2NRA,NZA1,NLA1
    Integer(ipr)  :: ID1(nbx)
    !---------------------------------------------------------------------
    ! Read data
    !---------------------------------------------------------------------
    Do iw=lout,lfile
       Write(iw,*)
       Write(iw,*) ' Reading from wel_file: ',welfile
       Write(iw,*)
    End Do
    counterLine = 0; iexit = 0
    Read(lwin,Err=100,End=100) npr11,npr1,ngh1,ngl1,n001,nb1,nt1
    counterLine = counterLine+1
    If(Abs(n001).Ne.Abs(n00).And.nb1.Ne.nb) go to 100
    Read(lwin,Err=100,End=100) b01,bz1,bp1,beta1,si,etot,rms,bet,xmix,v0r,v1r,pwir, &
                               del,ept,ala,ala2,alast,tz1,varmas,varmasNZ,pjmassNZ, &
                               ass,skass
    brin=zero; bbroyden='L'; !si=one;
    counterLine = counterLine+1
    Read(lwin,Err=100,End=100) ntx1,nb1,nhhdim1
    counterLine = counterLine+1
    Read(lwin,Err=100,End=100) lambdaMax1
    counterLine = counterLine+1
    Read(lwin,Err=100,End=100) multLag
    counterLine = counterLine+1
    Read(lwin,Err=100,End=100) id1
    counterLine = counterLine+1
    Read(lwin,Err=100,End=100) brin
    counterLine = counterLine+1
    !
    ! Add small pairing de=de+0.1 in the no-LN case to prevent pairing collapse
    If(kindhfb.Eq.1.And.Add_Pairing) Then
       ibro=0
       Do ib=1,NB
          ND=ID1(ib)
          I=ibro
          Do N1=1,ND
             Do N2=1,N1
                I=I+1
                brin(i+nhhdim2)=brin(i+nhhdim2)+0.10_pr
                brin(i+nhhdim3)=brin(i+nhhdim3)+0.10_pr
             End Do !N2
          End Do !N1
          ibro=i
       End Do !IB
    End If
    Do ib=1,NB
       ND=ID1(ib)
       Do N1=1,ND
          Read(lwin,ERR=100,End=100) NLANSA0,NLANSA1,NZA2NRA,NZA1,NLA1
       End Do
    End Do
    counterLine = counterLine+1
    ! blocking
    Read(lwin,ERR=100,End=100)  bloall1
    counterLine = counterLine+1
    Read(lwin,ERR=100,End=100)  bloblo,blo123,blok1k2,blomax,bloqpdif
    counterLine = counterLine+1
    If(bloall1.Ne.bloall) go to 100
    !tel
    If(iLST.Gt.0) Then
       Read(lwin,ERR=100,End=100) decay,rmm3,cmm3,amm3,bmm3,itass,iqqmax
       If(Allocated(fdsx)) Deallocate(fdsx,fdsy,fdsy1,fdsy2,fdsy3,  &
            fspb0,fspc0,fspd0,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,  &
            fspb3,fspc3,fspd3)
       Allocate(fdsx(iqqmax),fdsy(iqqmax),fdsy1(iqqmax),  &
            fdsy2(iqqmax),fdsy3(iqqmax),fspb0(iqqmax),fspc0(iqqmax),  &
            fspd0(iqqmax),fspb1(iqqmax),fspc1(iqqmax),fspd1(iqqmax),  &
            fspb2(iqqmax),fspc2(iqqmax),fspd2(iqqmax),fspb3(iqqmax),  &
            fspc3(iqqmax),fspd3(iqqmax))
       Read(lwin,ERR=100,End=100) fdsx,fdsy,fdsy1,fdsy2,fdsy3,fspb0,fspc0,fspd0  &
            ,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,fspb3,fspc3,fspd3
    End If
    Close(lwin)
    Return
    !
100 Continue
    iexit=1
     Do iw=lout,lfile
        Write(iw,'(1x,a,a,a)')
        Write(iw,'(1x,a,a,a)')   ' The file ',welfile,' is corrupted!'
        Write(iw,'(1x,a,i2,a)')  ' Problem occurs at line ',counterLine,'        '
        Write(iw,'(1x,a,a,a)')   ' STARTING FROM SCRATCH WITH ININ=IABS(ININ)!'
        Write(iw,'(1x,a,a,a)')
     End Do
  End Subroutine read_data_old
  !=======================================================================
  !> This subroutine reads a binary file containing the results of a
  !> HFBTHO calculation. In the new format (VERSION_DATA=3), the binary
  !> file is structured by keywords and contains the HF and pairing field
  !> on the Gauss quadrature mesh. This new format allows restarts even
  !> when (i) the basis is different (as long as the quadrature mesh is
  !> the same) (ii) the number of constraints is different (iii) metadata
  !> such as number of protons, neutrons, characteristics of the EDF, etc.
  !> are different. The price to pay for this flexibility is an increased
  !> size of the binary file.
  !=======================================================================
  Subroutine read_data(iexit)
    Implicit None
    Integer(ipr), INTENT(INOUT) :: iexit
    Character(Len=8) :: key
    Logical :: different_basis, different_mesh
    Integer(ipr) :: counterLine,iw,ib,nd,n1,icons,lambda,jcons,lambda_r,bloall_r
    Integer(ipr) :: switch_to_THO_r, projection_is_on_r
    Logical :: collective_inertia_r, fission_fragments_r, pairing_regularization_r, localization_functions_r, &
               set_temperature_r, set_neck_constrain_r, set_dfrag_constrain_r, set_csi_constrain_r
    Integer(ipr) :: Z_r,N_r,n00_r,nb_r,nt_r,ngh_r,ngl_r,nleg_r
    Integer(ipr) :: nr_r,nz_r,nl_r,ns_r,numberCons_r,lambdaMax_r
    Integer(ipr) :: NDCOMP_r,NOCOMP_r
    Real(pr) :: b0_r,bz_r,bp_r,neckLag_r,neckRequested_r,pwi_r,tz_r,varmas_r,varmasNZ_r
    Real(pr) :: dfragLag_r,dfragRequested_r,csiLag_r,csiRequested_r
    ! Arrays
    Integer(ipr), Allocatable :: ID_r(:),multLambda_r(:)
    Real(pr), Allocatable :: multRequested_r(:),multLag_r(:)
    Real(pr), Allocatable :: xh_r(:),xl_r(:),wh_r(:),wl_r(:)
    !---------------------------------------------------------------------
    ! Read data
    !---------------------------------------------------------------------
    Do iw=lout,lfile
       Write(iw,*)
       Write(iw,*) ' Reading from wel_file: ',welfile
       Write(iw,*)
    End Do
    counterLine = 0; iexit = 0
    ! Loop over all keywords in the file
    Do
       Read(lwin,Err=100,End=99) key
       ! Metadata: N, Z, force, optional flags
       If(Trim(key)=='Metadata') Then
          Read(lwin,Err=100,End=100) Z_r,N_r
          Read(lwin,Err=100,End=100) collective_inertia_r, fission_fragments_r, pairing_regularization_r, localization_functions_r
          Read(lwin,Err=100,End=100) switch_to_THO_r, projection_is_on_r, set_temperature_r
          Read(lwin,Err=100,End=100) set_neck_constrain_r, set_dfrag_constrain_r, set_csi_constrain_r
       End If
       ! Basis: deformations, number of quanta, integration mesh
       If(Trim(key)=='HO-Basis') Then
          Read(lwin,Err=100,End=100) b0_r,bz_r,bp_r
          Read(lwin,Err=100,End=100) n00_r,nb_r,nt_r,ngh_r,ngl_r,nleg_r
          If(ngh_r.Ne.ngh.Or.ngl_r.Ne.ngl) Then
             Write(6,'("Error in read_data - Inconsistent quadrature mesh!")')
             Write(6,'("The code will start from scratch")')
             iexit=1
          End If
          Allocate(xh_r(1:ngh_r),xl_r(1:ngl_r),wh_r(1:ngh_r),wl_r(1:ngl_r))
          Read(lwin,Err=100,End=100) xh_r,xl_r,wh_r,wl_r
       End If
       ! Quantum numbers and \Omega-blocks information
       If(Trim(key)=='QuantNum') Then
          Allocate(ID_r(nb_r))
          Read(lwin,Err=100,End=100) ID_r
          Do ib=1,nb_r
             nd=ID_r(ib)
             Do n1=1,nd
                Read(lwin,Err=100,End=100) nr_r,nz_r,nl_r,ns_r
             End Do
          End Do
          ! Check if the file has a different basis
          different_basis = nb.Ne.nb_r .Or. b0_r.Ne.b0 .Or. bz_r.Ne.bz .Or. bp_r.Ne.bp &
                                       .Or. n00_r.Ne.n00 .Or. nt_r.Ne.nt
       End If
       ! Various
       If(Trim(key)=='Various.') Then
          Read(lwin,Err=100,End=100) si,etot,rms,bet,xmix
          Read(lwin,Err=100,End=100) pwi_r,del,ept,ala,ala2,alast
          Read(lwin,Err=100,End=100) tz_r,varmas_r,varmasNZ_r,pjmassNZ,ass,skass
          siold=si
       End If
       ! Constraints: requested values, Lagrange parameters
       If(Trim(key)=='Constrai') Then
          Read(lwin,Err=100,End=100) numberCons_r,lambdaMax_r
          Allocate(multLambda_r(1:numberCons_r),multRequested_r(0:lambdaMax_r),multLag_r(1:lambdaMax_r))
          Read(lwin,Err=100,End=100) multLambda_r
          Read(lwin,Err=100,End=100) multRequested_r
          Read(lwin,Err=100,End=100) multLag_r
          If(set_neck_constrain_r) Then
             Read(lwin,Err=100,End=100) neckRequested_r
             Read(lwin,Err=100,End=100) neckLag_r
          End If
          If(set_dfrag_constrain_r) Then
             Read(lwin,Err=100,End=100) dfragRequested_r
             Read(lwin,Err=100,End=100) dfragLag_r
          End If
          If(set_csi_constrain_r) Then
             Read(lwin,Err=100,End=100) csiRequested_r
             Read(lwin,Err=100,End=100) csiLag_r
          End If
          
          
          Read(lwin,Err=100,End=100) qfield
          ! Reset values of Lagrange parameters based on the values read on disk
          Do icons=1,numberCons
             lambda=multLambda(icons)
             Do jcons=1,numberCons_r
                lambda_r=multLambda_r(jcons)
                If(lambda==lambda_r) Then
                   multLag(lambda) = multLag_r(lambda_r)
                End If
             End Do
          End Do
          If(set_neck_constrain .And. set_neck_constrain_r) Then
             neckLag = neckLag_r
          End If
          If(set_dfrag_constrain .And. set_dfrag_constrain_r) Then
             dfragLag = dfragLag_r
          End If
          If(set_csi_constrain .And. set_csi_constrain_r) Then
             csiLag = csiLag_r
          End If
          
       End If
       ! Density matrix and pairing tensor in coordinate space
       If(Trim(key)=='Densits.') Then
          Read(lwin,Err=100,End=100) ro     ! 2*rho
          Read(lwin,Err=100,End=100) aka    ! Kappa
       End If
       ! HF and pairing field in coordinate space
       If(Trim(key)=='FieldsN.') Then
          Read(lwin,Err=100,End=100) vn     ! RHO_ij
          Read(lwin,Err=100,End=100) vhbn   ! TAU_ij
          Read(lwin,Err=100,End=100) vrn    ! NABLAr RHO__ij
          Read(lwin,Err=100,End=100) vzn    ! NABLAz RHO__ij
          Read(lwin,Err=100,End=100) vdn    ! DELTA RHO_ij
          Read(lwin,Err=100,End=100) vsn    ! NABLA . J__ij
          Read(lwin,Err=100,End=100) vSFIZn ! JFIZ_ij
          Read(lwin,Err=100,End=100) vSZFIn ! JZFI_ij
          Read(lwin,Err=100,End=100) vSFIRn ! JFIR_ij
          Read(lwin,Err=100,End=100) vSRFIn ! JRFI_ij
          Read(lwin,Err=100,End=100) dvn    ! \Delta_ij
          End If
       If(Trim(key)=='FieldsP.') Then
          Read(lwin,Err=100,End=100) vp
          Read(lwin,Err=100,End=100) vhbp
          Read(lwin,Err=100,End=100) vrp
          Read(lwin,Err=100,End=100) vzp
          Read(lwin,Err=100,End=100) vdp
          Read(lwin,Err=100,End=100) vsp
          Read(lwin,Err=100,End=100) vSFIZp
          Read(lwin,Err=100,End=100) vSZFIp
          Read(lwin,Err=100,End=100) vSFIRp
          Read(lwin,Err=100,End=100) vSRFIp
          Read(lwin,Err=100,End=100) dvp
       End If
       ! Blocking
       If(Trim(key)=='Blocking') Then
          Do ib=1,2
             Call blosort(ib,blomax(ib))
          End Do
          Read(lwin,Err=100,End=100) bloall_r
          Read(lwin,Err=100,End=100) bloblo,blo123,blok1k2,blomax,bloqpdif
       End If
       ! Temperature
       If(Trim(key)=='Temperat') Then
          Read(lwin,Err=100,End=100) temper,entropy
          Read(lwin,Err=100,End=100) fp_T
          Read(lwin,Err=100,End=100) fn_T
       End If
       ! THO
       If(Trim(key)=='THObasis') Then
          If(Allocated(fdsx)) Then
             Read(lwin,Err=100,End=100) decay,rmm3,cmm3,amm3,bmm3,itass,iqqmax
             Read(lwin,Err=100,End=100) fdsx,fdsy,fdsy1,fdsy2,fdsy3,fspb0,fspc0,fspd0,  &
                                        fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,fspb3,fspc3,fspd3
          End If
       End If
       ! Pairing regularization
       If(Trim(key)=='Regular.') Then
          Read(lwin,Err=100,End=100) MEFFn
          Read(lwin,Err=100,End=100) MEFFp
          Read(lwin,Err=100,End=100) geff_inv
       End If
       ! Gogny force
       If(Trim(key)=='GognyVNN') Then
          ! Finite-range basis is implemented in configuration space and the
          ! basis on file must be the same as the current one for smooth restart
          If(different_basis) iexit=1
          Read(lwin,Err=100,End=100) NumVz,NumVr
          Read(lwin,Err=100,End=100) vrGogny
          Read(lwin,Err=100,End=100) vzGogny
       End If
       ! Collective inertia
       If(Trim(key)=='CollMass') Then
          Read(lwin,Err=100,End=100) NDCOMP_r,NOCOMP_r
          Read(lwin,Err=100,End=100) SK1
          Read(lwin,Err=100,End=100) SK2
          Read(lwin,Err=100,End=100) SK3
          Read(lwin,Err=100,End=100) ATDMAS
          Read(lwin,Err=100,End=100) GCMMAS
          Read(lwin,Err=100,End=100) E0_ATD
          Read(lwin,Err=100,End=100) E0_GCM
       End If
    End Do
    !
 99 Continue
    Close(lwin)
    Return
    !
100 Continue
    iexit=1
    Close(lwin)
    Do iw=lout,lfile
       Write(iw,'(1x,a,a,a)')
       Write(iw,'(1x,a,a,a)')   ' The file ',welfile,' is corrupted!'
       Write(iw,'(1x,a,a8,a)')  ' Problem occurs for key ',key,'        '
       Write(iw,'(1x,a,a,a)')   ' STARTING FROM SCRATCH WITH ININ=IABS(ININ)!'
       Write(iw,'(1x,a,a,a)')
    End Do
    !
  End Subroutine read_data
  !=======================================================================
  !> This subroutine writes a binary file containing the results of a
  !> HFBTHO calculation. In the new format (VERSION_DATA=3), the binary
  !> file is structured by keywords and contains the HF and pairing field
  !> on the Gauss quadrature mesh.
  !=======================================================================
  Subroutine write_data()
    Implicit None
    Integer(ipr) :: ib,nd,n1,ibasis,nla,nra,nza,nsa

    ! Metadata: N, Z, force, optional flags
    Write(lwou) 'Metadata'
    Write(lwou) npr(2),npr(1)
    Write(lwou) collective_inertia, fission_fragments, pairing_regularization, localization_functions
    Write(lwou) switch_to_THO, projection_is_on, set_temperature, set_neck_constrain, &
         set_dfrag_constrain, set_csi_constrain
    ! Basis: deformations, number of quanta, integration mesh
    Write(lwou) 'HO-Basis'
    Write(lwou) b0,bz,bp
    Write(lwou) n00,nb,nt,ngh,ngl,nleg
    Write(lwou) xh,xl,wh,wl
    ! Quantum numbers and \Omega-blocks information
    Write(lwou) 'QuantNum'
    Write(lwou) ID
    ibasis=0
    Do ib=1,nb
       nd=ID(ib)
       Do n1=1,nd
          ibasis=ibasis+1
          nla=NL(ibasis); nra=NR(ibasis); nza=NZ(ibasis); nsa=NS(ibasis)
          Write(lwou) nra,nza,nla,nsa
       End Do
    End Do
    ! Pairing
    Write(lwou) 'Various.'
    Write(lwou) si,etot,rms,bet,xmix
    Write(lwou) pwi,del,ept,ala,ala2,alast
    Write(lwou) tz,varmas,varmasNZ,pjmassNZ,ass,skass
    ! Constraints: requested values, Lagrange parameters
    Write(lwou) 'Constrai'
    Write(lwou) numberCons,lambdaMax
    Write(lwou) multLambda
    Write(lwou) multRequested
    Write(lwou) multLag
    If(set_neck_constrain) Then
       Write(lwou) neckRequested
       Write(lwou) neckLag
    End If
    If(set_dfrag_constrain) Then
       Write(lwou) dfragRequested
       Write(lwou) dfragLag
    End If
    If(set_csi_constrain) Then
       Write(lwou) csiRequested
       Write(lwou) csiLag
    End If
    Write(lwou) qfield
    ! Density matrix and pairing tensor in coordinate space
    Write(lwou) 'Densits.'
    Write(lwou) ro     ! Rho
    Write(lwou) aka    ! Kappa
    ! HF and pairing field in coordinate space
    Write(lwou) 'FieldsN.'
    Write(lwou) vn     ! RHO_ij
    Write(lwou) vhbn   ! TAU_ij
    Write(lwou) vrn    ! NABLAr RHO__ij
    Write(lwou) vzn    ! NABLAz RHO__ij
    Write(lwou) vdn    ! DELTA RHO_ij
    Write(lwou) vsn    ! NABLA . J__ij
    Write(lwou) vSFIZn ! JFIZ_ij
    Write(lwou) vSZFIn ! JZFI_ij
    Write(lwou) vSFIRn ! JFIR_ij
    Write(lwou) vSRFIn ! JRFI_ij
    Write(lwou) dvn    ! \Delta_ij
    Write(lwou) 'FieldsP.'
    Write(lwou) vp
    Write(lwou) vhbp
    Write(lwou) vrp
    Write(lwou) vzp
    Write(lwou) vdp
    Write(lwou) vsp
    Write(lwou) vSFIZp
    Write(lwou) vSZFIp
    Write(lwou) vSFIRp
    Write(lwou) vSRFIp
    Write(lwou) dvp
    ! Blocking
    Write(lwou) 'Blocking'
    Do ib=1,2
       Call blosort(ib,blomax(ib))
    End Do
    Write(lwou) bloall
    Write(lwou) bloblo,blo123,blok1k2,blomax,bloqpdif
    ! Temperature
    If(set_temperature) Then
       Write(lwou) 'Temperat'
       Write(lwou) temper,entropy
       Write(lwou) fp_T
       Write(lwou) fn_T
    End If
    ! THO
    If(switch_to_THO.Ne.0) Then
       Write(lwou) 'THObasis'
       If(iLST.Gt.0) Then
          If(Allocated(fdsx)) Then
             Write(lwou) decay,rmm3,cmm3,amm3,bmm3,itass,iqqmax
             Write(lwou) fdsx,fdsy,fdsy1,fdsy2,fdsy3,fspb0,fspc0,fspd0  &
                  ,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,fspb3,fspc3,fspd3
          End If
       End If
    End If
    ! Pairing regularization
    If(pairing_regularization) Then
       Write(lwou) 'Regular.'
       Write(lwou) MEFFn
       Write(lwou) MEFFp
       Write(lwou) geff_inv
    End If
    ! Gogny force
    If(finite_range) Then
       Write(lwou) 'GognyVNN'
       Write(lwou) NumVz,NumVr
       Write(lwou) vrGogny
       Write(lwou) vzGogny
    End If
    ! Collective inertia
    If(collective_inertia) Then
       Write(lwou) 'CollMass'
       Write(lwou) NDCOMP,NOCOMP
       Write(lwou) SK1
       Write(lwou) SK2
       Write(lwou) SK3
       Write(lwou) ATDMAS
       Write(lwou) GCMMAS
       Write(lwou) E0_ATD
       Write(lwou) E0_GCM
    End If
    Close(lwou)
    !
  End Subroutine write_data
  !=======================================================================
  !> For debugging only: this subroutine writes the binary file using the
  !> old convention of HFBTHO version 200d
  !=======================================================================
  Subroutine write_data_old()
    Implicit None
    Integer(ipr)  :: iw,N1,ND,ib,ibasis
    Integer(ipr)  :: npr1,npr11,NLANSA1,NLA,NRA,NZA,NSA
    !
    npr11=npr(1); npr1=npr(2)
    Write(lwou) npr11,npr1,ngh,ngl,n00,nb,nt
    Write(lwou) b0,bz,bp,beta0,si,etot,rms,bet,xmix,CpV0,CpV1,pwi,  &
                del,ept,ala,ala2,alast,tz,varmas,varmasNZ,pjmassNZ, &
                ass,skass
    Write(lwou) ntx,nb,nhhdim
    Write(lwou) lambdaMax
    Write(lwou) multLag
    Write(lwou) id
    Write(lwou) brin
    ibasis=0
    Do ib=1,NB
       ND=ID(ib)
       Do N1=1,ND
          ibasis=ibasis+1
          NLA=NL(ibasis); NRA=NR(ibasis); NZA=NZ(ibasis); NSA=NS(ibasis); NLANSA1=(-1)**(NZA+NLA)
          Write(lwou) 2*NLA+NSA,NLANSA1,NZA+2*NRA+NLA,NZA,NLA
       End Do
    End Do
    !---------------------------------------------------------------------
    ! blocking: sort blocking candidates first
    !---------------------------------------------------------------------
    Do ib=1,2
       Call blosort(ib,blomax(ib))
    End Do
    Write(lwou) bloall
    Write(lwou) bloblo,blo123,blok1k2,blomax,bloqpdif
    !tel
    If(iLST.Gt.0) Then
       If(Allocated(fdsx)) Then
          Write(lwou) decay,rmm3,cmm3,amm3,bmm3,itass,iqqmax
          Write(lwou) fdsx,fdsy,fdsy1,fdsy2,fdsy3,fspb0,fspc0,fspd0  &
               ,fspb1,fspc1,fspd1,fspb2,fspc2,fspd2,fspb3,fspc3,fspd3
       End If
    End If
    Close(lwou)
    Do iw=lout,lfile
       Write(iw,'(a,a,a)')
       Write(iw,'(a,a,a)') '  Writing to wel_file: ',welfile
       Write(iw,'(a,a,a)') ' __________________________________  '
       Write(iw,'(a,a,a)') '  The tape ',welfile,' recorded:     '
       Write(iw,'(a,a,a)') '  nucname,npr,ngh,ngl,n00,nb,nt      '
       Write(iw,'(a,a,a)') '  b0,beta0,si,etot,rms,bet,xmix      '
       Write(iw,'(a,a,a)') '  pairing:     CpV0,CpV1,pwi         '
       Write(iw,'(a,a,a)') '  delta:       del,ept               '
       Write(iw,'(a,a,a)') '  lambda:      ala,ala2,alast,tz     '
       Write(iw,'(a,a,a)') '  asymptotic:  varmas,ass,skass      '
       Write(iw,'(a,a,a)') '  ntx,nb,nhhdim,id,N_rz,n_r,n_z      '
       Write(iw,'(a,a,a)') '  Omega2,Sigma2,Parity,Lambda        '
       Write(iw,'(a,a,a)') '  matrices(inbro):    hh,de          '
       Write(iw,'(a,a,a)') '  *all blocking candidates           '
       If(Allocated(fdsx)) Write(iw,'(a,a,a)') '  *all THO arrays                    '
       Write(iw,'(a,a,a)') ' __________________________________  '
       Write(iw,'(a,a,a)')
    End Do
    !
  End Subroutine write_data_old
  !=======================================================================
  !> Sorting blocking candidates before writing them to disk
  !=======================================================================
  Subroutine blosort(it,n)
    Implicit None
    Integer(ipr) :: it,ip,n,i,k,j
    Real(pr) :: p
    Do i=1,n
       k=i; p=bloqpdif(i,it)
       If (i.Lt.n) Then
          Do j=i+1,n
             If (bloqpdif(j,it).Lt.p) Then
                k=j; p=bloqpdif(j,it)
             End If
          End Do
          If (k.Ne.i) Then
             bloqpdif(k,it)=bloqpdif(i,it); bloqpdif(i,it)=p
             ip=bloblo(k,it);  bloblo(k,it)=bloblo(i,it);  bloblo(i,it)=ip
             ip=blo123(k,it);  blo123(k,it)=blo123(i,it);  blo123(i,it)=ip
             ip=blok1k2(k,it); blok1k2(k,it)=blok1k2(i,it); blok1k2(i,it)=ip
          End If
       End If
    End Do
  End Subroutine blosort
  !=======================================================================
  !> Generate the values of a scalar field on an irregular Cartesian mesh
  !> given its values on a cylindrical mesh
  !=======================================================================
  Subroutine polar_to_cartesian()
    Implicit None
  End Subroutine polar_to_cartesian
  !=======================================================================
  !
  !=======================================================================
End Module HFBTHO_io
