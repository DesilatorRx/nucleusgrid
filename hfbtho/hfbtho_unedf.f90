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
!                             UNEDF PACKAGE                            !
!                                                                      !
! ==================================================================== !

!----------------------------------------------------------------------
!> This module provides the main \theCode DFT solver. It includes
!> routines for the calculation and diagonalization of the HFB matrix;
!> the definition of densities on the quadrature mesh and in configuration
!> space; the self-consistent loop the calculation of expectation values of
!> observables. It also includes a package of a few routines to perform
!> particle number projection.
!>
!>  @author
!>    Markus Kortelainen, Mario Stoitsov, Nicolas Schunck
!----------------------------------------------------------------------
!  Subroutines: - f01234(lpr)
!               - thofun(key,r,f,f1,f2,f3,fj,lpr,units)
!               - densitr(it,xr,yr,yrP,msw)
!               - gaupolr(z,x)
!----------------------------------------------------------------------
Module UNEDF
  Use HFBTHO_utilities
  Implicit None
  !
  Character(16), Private :: Version='17'
  !
  ! Version History
  !--------------------------------------------------------------------------------------
  ! ver#17:(Mario)   use_TMR_pairing=0/1 standard/TMR pairing added
  !                  to Namelist. Using:
  !                  CpV0(0)=G,    CpV0(1)=a
  !                  CpV1(0)=vfacn,CpV1(1)=vfacp
  ! ver#16:(Mario)   #ifndef hide_dme preprocessing directive included
  ! ver#15:(Markus)  Added parameter CExPar, used in Coul. excange term.
  !                  Also, all the channels included in direct Hartree
  ! ver#14:(Markus)  Added function Vexternal for the external field,
  !                  and use_j2terms to switch off tensor terms.
  !                  Direct Hartree set to zero.
  ! Ver#13:(Mario)   Added ac2,ac3,acoord
  ! ver#12:(Mario)   hartree term temprorary dropped. rDr NNN terms taken
  !                  with a factor of 1/2
  ! ver#11:(Mario)   Gaussian approximation to the Hartree term added,
  ! [3/10/2010]      hatree_rc removed. NB! Function HartreeV is an
  !                  elemental function with possible array arguments
  ! ver#10: (Markus) Added e2charg (e^2 for Coulomb) to the public variables
  ! ver#9: (Mario)   Hartree 'CHrho' calculated in INM with rc='hatree_rc'
  ! [2/2/2010]       is subtracted from Crho(0)at DMEorder >= 0.
  !                  CHrho added to the public list, 'hatree_rc' added
  !                  to interaction parameters and the namelist.
  !                  In the case DMEorder=-1 (standard Skyrme)
  !                  both, 'CHrho' and 'hatree_rc', do not play.
  !                  New function HartreeV(u) defines Hatree energy as
  !                  E(Hartree)=(1/2)*Int[rho_0(r)*V(|r-r'|)*rho_0(r')]
  !                  HartreeV(u) is zero for u=<'hatree_rc'
  ! ver#8: (Markus)  Hartree DME terms dropped out.
  ! ver#7: (Markus)  Added switch to turn off the 3N terms.
  !        (Mario)   Added Abs to density and gradient dependent LDA
  !                  Public :: DMEorder,DMElda,use_DME3N_terms
  ! ver#6: (Mario)   Skyrme transformation added.
  ! ver#5: (Mario)   Print_Namelist=T/F added to the namelist
  ! ver#4: (Markus)  Added natural units to the module. Used only for printing.
  ! ver#3: (Mario)   Uamplitudes(0:3,0:7) in normal order
  !
  ! t for Uamplitudes(t,*)
  ! 0 -> 0,0
  ! 1 -> 1,1
  ! 2 -> 0,1
  ! 3 -> 1,0
  ! n for Uamplitudes(*,n)
  ! 0 -> U
  ! 1 -> dU/dRHO_0
  ! 2 -> dU/dRHO_1
  ! 3 -> d2U/(dRHO_0*dRHO_0)
  ! 4 -> d2U/(dRHO_1*dRHO_1)
  ! 5 -> d2U/(dRHO_0*dRHO_1)
  ! 6 -> dU/d(TAU_0)
  ! 7 -> dU/d(Delta RHO_0)
  !
  ! TESTED MATTHEMETICA<=>BIRUC & SCOTT; MATTHEMETICA<=>Module UNEDF (energy amplitudes only)
  ! ver#2: (Mario) Pairing included
  !  - set_functional_parameters(fname,lpr)
  !  - pairing incorporated into CpV0(0:1),CpV1(0:1)
  !    as public variables also serving two public amplitudes
  !     Urhorhopr(0:1,0)=CpV0(0:1)+CpV1(0:1)*rho(0)
  !     Urhorhopr(0:1,1)=CpV1(0:1)
  !     so, they can be used with appropriate values by the DME solver
  !  -need improvement later,
  !      currently HFBTHO uses CpV0(0:1), CpV0(0:1)  as before
  !      just substituting V0,V1 in pn-representation
  !      CpV0*(1-CpV1/0.16*rho_0)and this defines
  !      the default values in the module CpV0=V0,CpV1=1/2)
  !  -NAMELIST and input/output modified. RESERVED NAMES ARE:
  !      -namelist forbiden:
  !          'UNRDF'  - best UNEDF
  !          'SKYRME' - best SKYRME
  !      -namelist inforced but not for C-parameters (use_INM=F)
  !       or NM-parameters (use_INM=T) defined by the solver
  !          'FITS'
  !      -namelist inforced (one can overwrite all):
  !          'ANY OTHER NAME'
  !       i.e., the solver defines C-/NM- only using 'FITS'
  ! ver#1: (Mario) Complete rewrite consistent with HFBTHO
  !  -CB-LDA added
  !  -INM added
  !  -HFBTHO BENCHMARK: LN, ZR(110) prolate solution with SLY4,
  !   mixed pairing and tensor terms. Agreement with previouse
  !   implemetation to the last significant digit in the cases:
  !      - Standard Skyrme
  !      - LO+LDA
  !      - LO+CB-LDA
  !      - (NrNr=0,rDj=0), (rDr=0,jDr=0), 0.5(NrNr=-rDr,jDr=-rDj)
  !   -use_j2terms removed, i.e., in the SKYRME case CJ=0 removes all
  !    tensor terms, while in DME tensor terms are always present
  ! ver#0: (Marcus) Basic coding from scratch
  !   -DME(u) consistent with Mathematica numbers
  !   -including small 'u' approximation
  !--------------------------------------------------------------------------------------
  !
  ! === PUBLIC VARIABLES ===
  !
  ! Use pointers to prevent conflicts with UNEDF public variabes
  ! Example: Use UNEDF, pr=>my_pr, ipr=>my_ipr, Crho=>my_Crho ...
  !
  !--------------------------------------------------------------------------------------
  !
  ! === PUBLIC VARIABLES ===
  !
  Logical, Public :: use_charge_density, use_cm_cor,use_DME3N_terms,   &
                     use_j2terms,use_full_cm_cor,use_INM,use_Namelist, &
                     Print_Namelist,finite_range,hb0_charge_dependent
  Integer(ipr), Public :: DMEorder,DMElda,use_TMR_pairing
  Real(pr), Public, Dimension(0:3,0:7) :: Urhorho,Urhotau,UrhoDrho,Unablarho  ! ph DME amplitudes
  Real(pr), Public, Dimension(0:3,0:7) :: UJnablarho,UrhonablaJ,UJJ
  Real(pr), Public, Dimension(0:3,0:7) :: Urhorhopr                           ! pp amplitudes
  Real(pr), Public, Dimension(0:1) :: UEnonstdr,UFnonstdr,URnonstdr           ! Other amplitudes
  Real(pr), Public :: hbzero,sigma,e2charg,CExPar                             ! hbr^2/2m, DD sigma, e^2 charge, coul.exch.
  Real(pr), Public :: hbzeron,hbzerop                                         ! hbr^2/2m_n,hbr^2/2m_p
  Real(pr), Public, Dimension(0:1) :: Crho,Cdrho,Ctau,CrDr,CrdJ,CJ,CpV0,CpV1  ! basic coupling constants
  Real(pr), Public :: E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM,P_NM,KA_NM
  Real(pr), Public :: CHrho                                                   ! Crho(0) from the Hartree term in NM
  Real(pr), Public :: mpi,gA,fpi,c1,c3,c4,cd,ce,LambdaX
  Real(pr), PUBLIC :: t0s,t0a,drs,dra,ts,ta,t3alp,t3al0,t3alm,t324,alp,alm,wla0, &
                      wla1,TA7,TA8,TB7,TB8,tv1,tv2,tv3,tv4,tv5,tv6,ts1,ts2,t4o3
  !Gogny parameters
  real(pr), Public, allocatable, dimension(:) :: mu_g, W_g, B_g, H_g, M_g
  integer(ipr) :: n_g=2
  !
  ! === PRIVATE VARIABLES ===
  !
  Real(pr), Private, Dimension(0:1) :: nuCrho,nuCdrho,nuCtau,nuCrDr  ! basic coupling constants in natural units
  Real(pr), Private, Dimension(0:1) :: nuCrdJ,nuCJ,nuCpV0,nuCpV1     !
  Real(pr), Private :: t0,t1,t2,t3,x0,x1,x2,x3,b4,b4p,te,to
  Real(pr), Private :: nuLambda,nufpi                                ! parameters associated to natural units
  Real(pr), Private, Dimension(0:1) :: Cnrho,CJdr                    ! hidden and always zero
  Integer(ipr), Private :: i_cut                                     ! dmeorder: -1=Standard Skyrme, 0=LO, 1=NLO, 2=N2LO
  Real(pr), Private :: Pi,eps                                        ! dmelda: 0=Kf-LDA, 1=CB-LDA
  Real(pr), Private :: kfconst,CK                                    ! (3Pi^2/2)^(1/3)
  Real(pr), Parameter, Private :: mevfm=197.30_pr;
  Real(pr), Private :: rho(0:1),tau(0:1),nrho2(0:1),lrho(0:1)
  Real(pr), Private :: mpi2,fpi2,fpi4,gA2,gA4,gA6,CHartree
  Real(pr), Private :: arhorho,brhorho,arhodrho,brhodrho,arhotau,brhotau,ajj,bjj,adrdr,bdrdr
  Real(pr), Private :: darhorho,dbrhorho,darhodrho,dbrhodrho,darhotau,dbrhotau,dajj,dbjj,dadrdr,dbdrdr
  Real(pr), Private :: ddarhodrho,ddbrhodrho,ddarhotau,ddbrhotau,ddarhorho,ddbrhorho
  Real(pr), Private :: hrho0rho0,hrho1rho1,hdr0dr0,hdr1dr1,hrho0Drho0,hrho1Drho0, &
                       hrho1Drho1,hrho0tau0,hrho1tau0,hrho1tau1,hJ0dr0,hrho0DJ0,hJ1dr1,hrho1DJ1, &
                       hJ0dr1,hrho1DJ0,hJ1dr0,hJ0J0,hJ0J1,hJ1J1
  Real(pr), Private :: dhrho0rho0,dhrho1rho1,dhdr0dr0,dhdr1dr1,dhrho0Drho0, &
                       dhrho1Drho0,dhrho1Drho1,dhrho0tau0,dhrho1tau0,dhrho1tau1,dhJ0dr0,dhrho0DJ0, &
                       dhJ1dr1,dhrho1DJ1,dhJ0dr1,dhrho1DJ0,dhJ1dr0,dhJ0J0,dhJ0J1,dhJ1J1
  Real(pr), Private :: ddhrho0rho0,ddhrho1rho1,ddhrho0Drho0,ddhrho1Drho0, &
                       ddhrho1Drho1,ddhrho0tau0,ddhrho1tau0,ddhrho1tau1
  Real(pr), Private, Dimension(3,3,33) :: ctr0r0,ctr1r1,ctdr0dr0,ctdr1dr1, & ! coefficients for 3N part
                                          ctr0Dr0,ctr1Dr0,ctr1Dr1,ctr0t0,ctr1t0,ctr1t1,ctJ0dr0,ctr0dJ0,ctJ1dr1, &
                                          ctr1dJ1,ctJ0dr1,ctr1dJ0,ctJ1dr0,ctJ0J0,ctJ0J1,ctJ1J1
  Real(pr), Private :: u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12
  Real(pr), Private :: ual,lual,atu,asqu,asqu4
  Real(pr), Private :: ac2,ac3,acoord
  Parameter(acoord=0.50_pr,ac2=4.0_pr*(acoord**2-acoord+0.50_pr),ac3=2.0_pr*(acoord**2-acoord+0.50_pr))
  Character (30) :: FunctionalName
  !
  Real(pr), Private :: A1_1,A1_2,A1_3,A1_4,A1_5,b1_1,b1_2,b1_3,b1_4,b1_5
  Real(pr), Private :: A3_1,A3_2,A3_3,A3_4,A3_5,b3_1,b3_2,b3_3,b3_4,b3_5
  Real(pr), Private :: h0mpi6,h0mpi6c1,h0mpi6c3,h0mpi6NM,h0mpi6c1NM,h0mpi6c3NM
  !
  Namelist /UNEDF_NAMELIST/ FunctionalName, DMEorder, DMElda, use_INM, hbzero, use_TMR_pairing,   &
                            Crho, Cdrho, Ctau, CrDr, CrdJ, CJ, sigma, CpV0, CpV1, e2charg,        &
                            E_NM, K_NM, SMASS_NM, RHO_NM, ASS_NM, LASS_NM, VMASS_NM,              &
                            mpi, gA, fpi, c1, c3, c4, cd, ce, LambdaX,                            &
                            use_cm_cor, use_charge_density, use_DME3N_terms, use_j2terms, CExPar, &
                            Print_Namelist
 Contains
  !
  !=======================================================================
  !
  !=======================================================================
  Subroutine calculate_U_parameters(rho0_in,rho1_in,tau0_in,tau1_in,laprho0,laprho1,nablarho0s,nablarho1s)
    Implicit None
    Real(pr), Intent(in) :: rho0_in,rho1_in,tau0_in,tau1_in
    Real(pr), Intent(in), Optional :: &
         nablarho0s,nablarho1s,laprho0,laprho1
    Integer(ipr) :: t,i,j,k,l
    Real(pr) :: u,du,ddu,dtu,dlu
    Real(pr) :: ph,aux,daux,ddaux
    Real(pr) :: y,dy,ddy,marc,dmarc,ddmarc,mlog,dmlog,ddmlog
    Real(pr) :: ucut,ucut3n
    !
    ucut=0.1_pr; ucut3n=0.6_pr
    !
    rho(0)=rho0_in; rho(1)=rho1_in;
    tau(0)=tau0_in; tau(1)=tau1_in;
    !
    lrho=0.0_pr; nrho2=0.0_pr;
    If (Present(laprho0)) lrho(0)=laprho0
    If (Present(laprho1)) lrho(1)=laprho1
    If (Present(nablarho0s)) nrho2(0)=nablarho0s
    If (Present(nablarho1s)) nrho2(1)=nablarho1s
    !
    arhorho=0.0_pr; darhorho=0.0_pr; ddarhorho=0.0_pr
    brhorho=0.0_pr; dbrhorho=0.0_pr; ddbrhorho=0.0_pr
    arhodrho=0.0_pr; darhodrho=0.0_pr; ddarhodrho=0.0_pr
    brhodrho=0.0_pr; dbrhodrho=0.0_pr; ddbrhodrho=0.0_pr
    arhotau=0.0_pr; darhotau=0.0_pr; ddarhotau=0.0_pr
    brhotau=0.0_pr; dbrhotau=0.0_pr; ddbrhotau=0.0_pr
    adrdr=0.0_pr; dadrdr=0.0_pr
    bdrdr=0.0_pr; dbdrdr=0.0_pr
    ajj=0.0_pr; dajj=0.0_pr
    bjj=0.0_pr; dbjj=0.0_pr
    !
    hrho0rho0=0.0_pr; hrho1rho1=0.0_pr; hdr0dr0=0.0_pr; hdr1dr1=0.0_pr
    hrho0Drho0=0.0_pr; hrho1Drho0=0.0_pr; hrho1Drho1=0.0_pr
    hrho0tau0=0.0_pr; hrho1tau0=0.0_pr; hrho1tau1=0.0_pr
    hJ0dr0=0.0_pr; hrho0DJ0=0.0_pr; hJ1dr1=0.0_pr; hrho1DJ1=0.0_pr
    hJ0dr1=0.0_pr; hrho1DJ0=0.0_pr; hJ1dr0=0.0_pr
    hJ0J0=0.0_pr; hJ0J1=0.0_pr; hJ1J1=0.0_pr
    dhrho0rho0=0.0_pr; dhrho1rho1=0.0_pr; dhdr0dr0=0.0_pr; dhdr1dr1=0.0_pr
    dhrho0Drho0=0.0_pr; dhrho1Drho0=0.0_pr; dhrho1Drho1=0.0_pr
    dhrho0tau0=0.0_pr; dhrho1tau0=0.0_pr; dhrho1tau1=0.0_pr
    dhJ0dr0=0.0_pr; dhrho0DJ0=0.0_pr; dhJ1dr1=0.0_pr; dhrho1DJ1=0.0_pr
    dhJ0dr1=0.0_pr; dhrho1DJ0=0.0_pr; dhJ1dr0=0.0_pr
    dhJ0J0=0.0_pr; dhJ0J1=0.0_pr; dhJ1J1=0.0_pr
    ddhrho0rho0=0.0_pr; ddhrho1rho1=0.0_pr
    ddhrho0Drho0=0.0_pr; ddhrho1Drho0=0.0_pr; ddhrho1Drho1=0.0_pr
    ddhrho0tau0=0.0_pr; ddhrho1tau0=0.0_pr; ddhrho1tau1=0.0_pr
    !
    u=0.0_pr; du=0.0_pr; ddu=0.0_pr; dtu=0.0_pr; dlu=0.0_pr
    !
    Urhorho=0.0_pr   ; Urhotau=0.0_pr
    UrhoDrho=0.0_pr  ; Unablarho=0.0_pr
    UJnablarho=0.0_pr; UrhonablaJ=0.0_pr
    Urhorhopr=0.0_pr ; UJJ=0.0_pr
    UEnonstdr=0.0_pr ; UFnonstdr=0.0_pr ; URnonstdr=0.0_pr
    !
    ! Notations for Uamplitudes(0:3,0:7)
    ! t for Uamplitudes(t,*)
    ! 0 -> 0,0
    ! 1 -> 1,1
    ! 2 -> 0,1
    ! 3 -> 1,0
    ! n for Uamplitudes(*,n)
    ! 0 -> U
    ! 1 -> dU/dRHO_0
    ! 2 -> dU/dRHO_1
    ! 3 -> d2U/(dRHO_0*dRHO_0)
    ! 4 -> d2U/(dRHO_1*dRHO_1)
    ! 5 -> d2U/(dRHO_0*dRHO_1)
    ! 6 -> dU/d(TAU_0)
    ! 7 -> dU/d(Delta RHO_0)
    !
    !! 2N terms
    Do t=0,1
       ph=1.0_pr
       If(t.Eq.1) ph=-1.0_pr
       Urhorho(t,0)=Crho(t)+Cdrho(t)*rho(0)**sigma &
            +0.50_pr*(arhorho+ph*brhorho)*mevfm
       Urhotau(t,0)=Ctau(t)+0.50_pr*(arhotau+ph*brhotau)*mevfm        !! These two determine the
       UrhoDrho(t,0)=Crdr(t)+ac2*0.50_pr*(arhoDrho+ph*brhoDrho)*mevfm !! effective mass (when recoupling to p-n)??
       UJJ(t,0)=CJ(t)+0.50_pr*(ajj+ph*bjj)*mevfm
       Unablarho(t,0)=Cnrho(t)+0.50_pr*(adrdr+ph*bdrdr)*mevfm
       UrhonablaJ(t,0)=Crdj(t)
       UJnablarho(t,0)=Cjdr(t)

       Urhorho(t,1)=sigma*Cdrho(t)*(rho(0)**sigma)/(rho(0)+eps) &
            +0.50_pr*(darhorho+ph*dbrhorho)*du*mevfm
       Urhotau(t,1)=0.50_pr*(darhotau+ph*dbrhotau)*du*mevfm
       UrhoDrho(t,1)=ac2*0.50_pr*(darhoDrho+ph*dbrhoDrho)*du*mevfm
       UJJ(t,1)=0.50_pr*(dajj+ph*dbjj)*du*mevfm
       Unablarho(t,1)=0.50_pr*(dadrdr+ph*dbdrdr)*du*mevfm

       Urhorho(t,6)=0.50_pr*(darhorho+ph*dbrhorho)*dtu*mevfm
       Urhotau(t,6)=0.50_pr*(darhotau+ph*dbrhotau)*dtu*mevfm
       UrhoDrho(t,6)=ac2*0.50_pr*(darhoDrho+ph*dbrhoDrho)*dtu*mevfm
       UJJ(t,6)=0.50_pr*(dajj+ph*dbjj)*dtu*mevfm
       Unablarho(t,6)=0.50_pr*(dadrdr+ph*dbdrdr)*dtu*mevfm

       Urhorho(t,7)=0.50_pr*(darhorho+ph*dbrhorho)*dlu*mevfm
       Urhotau(t,7)=0.50_pr*(darhotau+ph*dbrhotau)*dlu*mevfm
       UrhoDrho(t,7)=ac2*0.50_pr*(darhoDrho+ph*dbrhoDrho)*dlu*mevfm
       UJJ(t,7)=0.50_pr*(dajj+ph*dbjj)*dlu*mevfm
       Unablarho(t,7)=0.50_pr*(dadrdr+ph*dbdrdr)*dlu*mevfm

       Urhorho(t,3)=sigma*(sigma-1.0_pr)*Cdrho(t)*(rho(0)**sigma)/(rho(0)**2+eps) &
            +0.50_pr*(darhorho+ph*dbrhorho)*ddu*mevfm &
            +0.50_pr*(ddarhorho+ph*ddbrhorho)*du*du*mevfm
       Urhotau(t,3)=0.50_pr*(darhotau+ph*dbrhotau)*ddu*mevfm &
            +0.50_pr*(ddarhotau+ph*ddbrhotau)*du*du*mevfm
       UrhoDrho(t,3)=ac2*0.50_pr*(darhoDrho+ph*dbrhoDrho)*ddu*mevfm &
            +ac2*0.50_pr*(ddarhoDrho+ph*ddbrhoDrho)*du*du*mevfm

    End Do
    Urhorhopr(0,0) = CpV0(0)*(1.0_pr-CpV1(0)*rho(0)/0.16_pr)          &
         +CpV0(1)*(1.0_pr-CpV1(1)*rho(0)/0.16_pr)
    Urhorhopr(1,0) = CpV0(0)*(1.0_pr-CpV1(0)*rho(0)/0.16_pr)          &
         +CpV0(1)*(1.0_pr-CpV1(1)*rho(0)/0.16_pr)
    Urhorhopr(2,0) = (CpV0(0)*(1.0_pr-CpV1(0)*rho(0)/0.16_pr)         &
         -CpV0(1)*(1.0_pr-CpV1(1)*rho(0)/0.16_pr))*2.0_pr
    Urhorhopr(0,1) = (-CpV0(0)*CpV1(0)-CpV0(1)*CpV1(1))/0.16_pr
    Urhorhopr(1,1) = (-CpV0(0)*CpV1(0)-CpV0(1)*CpV1(1))/0.16_pr
    Urhorhopr(2,1) = 2.0_pr*(-CpV0(0)*CpV1(0)+CpV0(1)*CpV1(1))/0.16_pr
    Urhorhopr=Urhorhopr/16.0_pr
    UEnonstdr=0.0_pr; UFnonstdr=0.0_pr; URnonstdr=0.0_pr
    !
    If (.Not.use_j2terms) Then
       UJJ=0.0_pr
    End If
    !
  End Subroutine calculate_U_parameters
  !=======================================================================
  !
  !=======================================================================
  Subroutine read_UNEDF_NAMELIST(fname,noForces)
    Use HFBTHO_utilities, Only: lout
    !--------------------------------------------------------------------------------
    ! RESERVED NAMES ARE:
    !  -namelist forbiden:
    !          'UNEDF'  - best UNEDF
    !          'SKYRME' - best SKYRME
    !  -namelist inforced but not for C-parameters (use_INM=F)
    !   or NM-parameters (use_INM=T) defined by the solver
    !          'FITS'
    !  -namelist inforced (one can overwrite all):
    !          'ANY OTHER NAME'
    ! i.e., the DME solver defines C-/NM- only using 'FITS'
    !--------------------------------------------------------------------------------
    Implicit None
    Character (30), Intent(inout) :: fname
    Character (30) :: inforcedname
    Logical        :: regularization
    Integer(ipr)   :: ios,lnamelist=16,noForces
    !
    ! parameters
    eps     = Spacing(1.0_pr)
    Pi      = 4.0_pr*Atan(1.0_pr)
    kfconst =(1.50_pr*Pi**2)**(1.0_pr/3.0_pr)    ! (3Pi^2/2)^(1/3)
    CK      = 3.0_pr/5.0_pr*kfconst**2
    !
    use_Namelist=.True.
    Do
       !---------------------------------------------------------------------
       ! Some default values for all cases
       !---------------------------------------------------------------------
       Print_Namelist=.False.
       FunctionalName="Bla-Bla"
       ! kind of the functional
       use_INM                = .False.
       use_DME3N_terms        = .False.
       use_charge_density     = .False.
       regularization         = .False.
       use_cm_cor             = .False.
       use_full_cm_cor        = .False.
       use_j2terms            = .False.
       hb0_charge_dependent   = .False.
       finite_range           = .False.
       use_TMR_pairing        =  0
       DMEorder               = -1
       DMElda                 =  0
       ! Coupling constants: ph channel
       Crho(0)  = -727.0933239596374733_pr; Crho(1)  =  474.8709969984467989_pr
       CDrho(0) =  612.1037411660222460_pr; CDrho(1) = -705.7204872069220301_pr
       Ctau(0)  =   33.8846741217252401_pr; Ctau(1)  =   32.4047409594248919_pr
       CrDr(0)  =  -76.9962031249999939_pr; CrDr(1)  =   15.6571351249999999_pr
       CrdJ(0)  =  -92.2500000000000000_pr; CrdJ(1)  =  -30.7500000000000000_pr
       CJ(0)    =   17.2096115000000012_pr; CJ(1)    =   64.5758124999999978_pr
       Cnrho    =    0.0000000000000000_pr; CJdr     =    0.0000000000000000_pr
       ! Coupling constants: pp channel
       CpV0     = -258.2000000000000000_pr; CpV1     =    0.5000000000000000_pr
       ! Various
       sigma    =    0.3062227576210547_pr;
       hbzero   =   20.7355300000000007_pr;
       e2charg  =    1.4399784085965135_pr ; CExPar = 1.0_pr
       ! DME
       mpi=   138.03_pr/197.3_pr; fpi = 92.4_pr/197.3_pr; gA = 1.29_pr
       c1 =    -0.81_pr/1000.0_pr*197.3_pr
       c3 =    -3.40_pr/1000.0_pr*197.3_pr
       c4 =     3.40_pr/1000.0_pr*197.3_pr
       cd = -2062.00_pr/1000.0_pr
       ce =  -625.00_pr/1000.0_pr
       ! Natural units
       LambdaX  = 700.0_pr/197.3_pr; nuLambda = 700.0_pr; nufpi = 93.0_pr
       ! Nuclear matter
       E_NM     = -15.972149141444596410_pr; RHO_NM   =  0.159538756711733398_pr
       K_NM     = 229.900964482603626493_pr; SMASS_NM =  1.439546988976078357_pr
       ASS_NM   =  32.004302815052007247_pr; LASS_NM  = 45.961751480461579433_pr
       VMASS_NM =   1.249838547196253424_pr
       !---------------------------------------------------------------------
       ! Select the functional: start with interaction
       !---------------------------------------------------------------------
       noForces=0 ! No forces to start with
       Call skforce(fname,noForces)
       !
       If (noForces.Eq.1) Then
           inforcedname='FORCE'
           use_Namelist=.False.
       Else
           FUNCTIONAL: Select Case (Trim(fname))
           Case ('FITS')
              inforcedname='FITS'
              use_Namelist=.False.
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              P_NM     =    0.000000000000000_pr
              VMASS_NM =    1.249838000000000_pr
           Case ('UNE0')
              inforcedname='UNE0'
              use_Namelist=.False.
              ! kind of the functional
              use_INM    = .True.
              use_cm_cor = .True.
              ! Surface coefficients
              CrDr(0)  =  -55.260600000000000_pr
              CrDr(1)  =  -55.622600000000000_pr
              CpV0(0)  = -170.374000000000000_pr
              CpV0(1)  = -199.202000000000000_pr
              CrdJ(0)  =  -79.530800000000000_pr
              CrdJ(1)  =   45.630200000000000_pr
              CJ(0)    =    0.000000000000000_pr
              CJ(1)    =    0.000000000000000_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.160526000000000_pr
              E_NM     =  -16.055900000000000_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  230.000000000000000_pr
              ASS_NM   =   30.542900000000000_pr
              LASS_NM  =   45.080400000000000_pr
              SMASS_NM =    0.900000000000000_pr
              VMASS_NM =    1.249838000000000_pr
           Case ('UNE1')
              inforcedname='UNE1'
              use_Namelist=.False.
              ! kind of the functional
              use_INM  = .True.
              ! Surface coefficients
              CrDr(0)  =  -45.135131022237300_pr
              CrDr(1)  = -145.382167908057000_pr
              CpV0(0)  = -186.065399575124000_pr
              CpV0(1)  = -206.579593890106000_pr
              CrdJ(0)  =  -74.026333176459900_pr
              CrdJ(1)  =  -35.658261114791700_pr
              CJ(0)    =    0.000000000000000_pr
              CJ(1)    =    0.000000000000000_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.158706769332587_pr
              E_NM     =  -15.800000000000000_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  220.000000000000000_pr
              ASS_NM   =   28.986789057772100_pr
              LASS_NM  =   40.004790480413600_pr
              SMASS_NM =    0.992423332283364_pr
              VMASS_NM =    1.249838574232270_pr
           Case ('UNE2')
              inforcedname='UNE2'
              use_Namelist=.False.
              ! kind of the functional
              use_INM     = .True.
              use_j2terms = .True.
              ! Surface coefficients
              CrDr(0)  =  -46.831409147060600_pr
              CrDr(1)  = -113.163790795259000_pr
              CpV0(0)  = -208.889001962571000_pr
              CpV0(1)  = -230.329984038628000_pr
              CrdJ(0)  =  -64.308862415783800_pr
              CrdJ(1)  =  -38.650194685135500_pr
              CJ(0)    =  -54.433363597372100_pr
              CJ(1)    =  -65.903031044593800_pr
              CExPar   =    1.000000000000000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.156310622197074_pr
              E_NM     =  -15.800000000000000_pr
              P_NM     =    0.000000000000000_pr
              K_NM     =  239.929568022437000_pr
              ASS_NM   =   29.131006470773700_pr
              LASS_NM  =   40.000000000000000_pr
              SMASS_NM =    1.073763804147980_pr
              VMASS_NM =    1.249838574232270_pr
           Case ('N0LO')
              inforcedname='N0LO'
              use_Namelist=.False.
              ! kind of the functional
              use_INM         = .True.
              use_j2terms     = .False.
              use_DME3N_terms = .False.
              DMEorder        = 0
              ! Surface coefficients
              CrDr(0)  =  -67.437_pr
              CrDr(1)  =   21.551_pr
              CpV0(0)  = -241.203_pr
              CpV0(1)  = -252.818_pr
              CrdJ(0)  =  -95.451_pr
              CrdJ(1)  =  -65.906_pr
              CJ(0)    =    0.000_pr
              CJ(1)    =    0.000_pr
              CExPar   =    1.000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.1595_pr
              E_NM     =  -15.9700_pr
              P_NM     =    0.0000_pr
              K_NM     =  229.9000_pr
              ASS_NM   =   32.0000_pr
              LASS_NM  =   45.9600_pr
              SMASS_NM =    1.4400_pr
              VMASS_NM =    1.2500_pr
           Case ('N1LO')
              inforcedname='N1LO'
              use_Namelist=.False.
              ! kind of the functional
              use_INM         = .True.
              use_j2terms     = .False.
              use_DME3N_terms = .False.
              DMEorder        = 1
              ! Surface coefficients
              CrDr(0)  =  -63.996_pr
              CrDr(1)  =    -9.276_pr
              CpV0(0)  = -241.484_pr
              CpV0(1)  = -252.222
              CrdJ(0)  =  -95.463_pr
              CrdJ(1)  =  -60.800_pr
              CJ(0)    =    0.000_pr
              CJ(1)    =    0.000_pr
              CExPar   =    1.000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.1595_pr
              E_NM     =  -15.9700_pr
              P_NM     =    0.0000_pr
              K_NM     =  229.9000_pr
              ASS_NM   =   32.0000_pr
              LASS_NM  =   45.9600_pr
              SMASS_NM =    1.4400_pr
              VMASS_NM =    1.2500_pr
           Case ('N2LO')
              inforcedname='N2LO'
              use_Namelist=.False.
              ! kind of the functional
              use_INM         = .True.
              use_j2terms     = .False.
              use_DME3N_terms = .True.
              DMEorder        = 2
              ! Surface coefficients
              CrDr(0)  = -197.132_pr
              CrDr(1)  =  -12.503_pr
              CpV0(0)  = -272.164_pr
              CpV0(1)  = -193.188_pr
              CrdJ(0)  = -193.188_pr
              CrdJ(1)  =   37.790_pr
              CJ(0)    =    0.000_pr
              CJ(1)    =    0.000_pr
              CExPar   =    1.000_pr
              ! Various
              Cnrho    =    0.000000000000000_pr
              CJdr     =    0.000000000000000_pr
              hbzero   =   20.735530000000000_pr
              e2charg  =    1.439978408596513_pr
              ! Associated INM parameters
              RHO_NM   =    0.1595_pr
              E_NM     =  -15.9700_pr
              P_NM     =    0.0000_pr
              K_NM     =  229.9000_pr
              ASS_NM   =   32.0000_pr
              LASS_NM  =   45.9600_pr
              SMASS_NM =    1.4400_pr
              VMASS_NM =    1.2500_pr
           Case default
              inforcedname=fname
              use_Namelist=.True.
           End Select FUNCTIONAL
       End If
       !---------------------------------------------------------------------
       ! Exit loop condition
       !---------------------------------------------------------------------
       If(.Not.use_Namelist) Exit
       !---------------------------------------------------------------------
       ! Read namelists
       !---------------------------------------------------------------------
       Open(lnamelist,file='UNEDF_NAMELIST.DAT',DELIM='APOSTROPHE') ! 'QUOTE'
       Read(UNIT=lnamelist,NML=UNEDF_NAMELIST,iostat=ios)
       If(ios.Ne.0) Then
          ! WRong entry within UNEDF_NAMELIST.DAT file
          Write(*,'(1X,/,A)') 'ATTENTION: WRONG INPUT!'
          Write(*,*) 'THE INPUT DATA WITH LABEL FUNCTIONALNAME=''',Trim(INFORCEDNAME),''''
          Write(*,*) 'INSIDE THE UNEDF_NAMELIST.DAT FILE IS WRONG.'
          Write(*,*) 'PLESE CORECT AND TRY AGAIN!'
          Stop 'PROGRAM STOP IN read_UNEDF_NAMELIST'
       End If
       Close(lnamelist)
       If(Trim(FunctionalName).Eq.Trim(inforcedname)) Exit
    End Do
    !---------------------------------------------------------------------
    ! See what the namelists modified
    !---------------------------------------------------------------------
    INFORCED_FUNCTIONAL: Select Case (Trim(inforcedname))
    Case ("FORCE")
       FunctionalName='FORCE'
    Case ("UNE0")
       FunctionalName='UNE0'
    Case ("UNE1")
       FunctionalName='UNE1'
    Case ("UNE2")
       FunctionalName='UNE2'
    Case ("N0LO")
       FunctionalName='N0LO'
    Case ("N1LO")
       FunctionalName='N1LO'
    Case ("N2LO")
       FunctionalName='N2LO'
    Case ("FITS")
       FunctionalName='FITS'
    Case default
       ! Missing entry within hfbtho_NAMELIST.dat file
       If(Trim(FunctionalName).Ne.Trim(inforcedname)) Then
          Write(*,'(1X,/,A)') 'ATTENTION: MISSING INPUT!'
          Write(*,*) 'THE INPUT DATA WITH LABEL FUNCTIONALNAME=''',Trim(INFORCEDNAME),''''
          Write(*,*) 'IS MISSING INSIDE THE UNEDF_NAMELIST.DAT FILE.'
          Write(*,*) 'PLEASE CORECT AND TRY AGAIN!'
          Stop 'PROGRAM STOP IN SET_FUNCTIONAL_PARAMETERS'
       End If
    End Select INFORCED_FUNCTIONAL
    !
  End Subroutine read_UNEDF_NAMELIST
  !=======================================================================
  !> Set up parameters of the Gogny force
  !=======================================================================
  Subroutine gogny_force(fname)
    Implicit None
    Character (30), Intent(inout) :: fname
    !
    INTERACTION: Select Case (Trim(fname))
    !---------------------------------------------------------------------
    ! PRC 21, 1568 (1980)
    !---------------------------------------------------------------------
    Case ('D1')
       n_g = 2
       if(.not.allocated(mu_g)) &
            allocate(mu_g(1:n_g),W_g(1:n_g),B_g(1:n_g),H_g(1:n_g),M_g(1:n_g))
       mu_g= [ 0.7_pr, 1.2_pr ]
       W_g = [ -402.40_pr, -21.30_pr ]
       B_g = [ -100.00_pr, -11.77_pr ]
       H_g = [ -496.20_pr,  37.27_pr ]
       M_g = [  -23.56_pr, -68.81_pr ]
    !---------------------------------------------------------------------
    ! CPC 63, 365 (1991)
    !---------------------------------------------------------------------
    Case ('D1S')
       n_g = 2
       if(.not.allocated(mu_g)) &
            allocate(mu_g(1:n_g),W_g(1:n_g),B_g(1:n_g),H_g(1:n_g),M_g(1:n_g))
       mu_g= [ 0.7_pr, 1.2_pr ]
       W_g = [ -1720.30_pr,  103.64_pr ]
       B_g = [  1300.00_pr, -163.48_pr ]
       H_g = [ -1813.53_pr,  162.81_pr ]
       M_g = [  1397.60_pr, -223.93_pr ]
    !---------------------------------------------------------------------
    ! PRC 21, 1568 (1980)
    !---------------------------------------------------------------------
    Case ('D1p')
       n_g = 2
       if(.not.allocated(mu_g)) &
            allocate(mu_g(1:n_g),W_g(1:n_g),B_g(1:n_g),H_g(1:n_g),M_g(1:n_g))
       mu_g= [ 0.7_pr, 1.2_pr ]
       W_g = [ -402.40_pr, -21.30_pr ]
       B_g = [ -100.00_pr, -11.77_pr ]
       H_g = [ -496.20_pr,  37.27_pr ]
       M_g = [  -23.56_pr, -68.81_pr ]
    !---------------------------------------------------------------------
    !
    !---------------------------------------------------------------------
    Case ('D1N')
       n_g = 2
       if(.not.allocated(mu_g)) &
            allocate(mu_g(1:n_g),W_g(1:n_g),B_g(1:n_g),H_g(1:n_g),M_g(1:n_g))
       mu_g= [ 0.8_pr, 1.2_pr ]
       W_g = [ -2047.61_pr,  293.02_pr ]
       B_g = [  1700.00_pr, -300.78_pr ]
       H_g = [ -2414.93_pr,  414.59_pr ]
       M_g = [  1519.35_pr, -316.84_pr ]
    !---------------------------------------------------------------------
    ! Default
    !---------------------------------------------------------------------
    Case default
        !Write(6,'("No Gogny interaction defined in routine skforce()")')
    End Select INTERACTION
    !
    Return
  End Subroutine gogny_force
  !=======================================================================
  !> Set up Pairing & Skyrme force parameters and their combinations
  !=======================================================================
  Subroutine skforce(fname,noForces)
    Implicit None
    Integer(ipr) :: noForces
    Real(pr) :: A,wls,TA7,TA8
    Real(pr) :: zero,one,two,three,four,five,six,seven,eight,nine
    Real(pr) :: half,pp16,pp24
    Character (30), Intent(inout) :: fname
    !
    zero = 0.0_pr; one = 1.0_pr; two = 2.0_pr; three = 3.0_pr; four = 4.0_pr
    five = 5.0_pr; six = 6.0_pr; seven = 7.0_pr; eight = 8.0_pr; nine = 9.0_pr
    half = 0.5_pr; pp16 = 16.0_pr; pp24 = 24.0_pr
    !
    ! Default for all forces if not modified
    hbzero = 1.0d0/0.04823_pr ! DMSHB0=1/hbzero
    sigma = one
    t0 = zero; x0 = zero
    t1 = zero; x1 = zero
    t2 = zero; x2 = one
    t3 = zero; x3 = one
    wls= zero; b4 = wls/two; b4p=wls/two
    te = zero; to = zero
    CExPar=1.0_pr
    !
    noForces=0 ! No forces at all
    !
    INTERACTION: Select Case (Trim(fname))
    !---------------------------------------------------------------------
    ! SIII force, Beiner et al., NPA238 (1975) 29
    !---------------------------------------------------------------------
    Case ('SIII')
        ! ph-Force
        noForces=1
        use_cm_cor = .True.
        hbzero = 20.73533_pr
        t0 = -.1128750d+04; x0 = +0.4500000_pr
        t1 = +.3950000d+03; x1 = +0.0000000_pr
        t2 = -.9500000d+02; x2 = +0.0000000_pr
        t3 = +.1400000d+05; x3 = +1.0000000_pr
        wls= +.1200000d+03; sigma = one
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -265.2500_pr, -340.0625_pr /)
    !---------------------------------------------------------------------
    ! SKM* forces
    !---------------------------------------------------------------------
    Case ('SKM*')
        ! ph-Force
        noForces=1
        use_cm_cor = .True.
        hbzero = 20.73_pr
        t0 = -.2645000d+04; x0 = +.0900000_pr
        t1 = +.4100000d+03; x1 = +.0000000_pr
        t2 = -.1350000d+03; x2 = +.0000000_pr
        t3 = +.1559500d+05; x3 = +.0000000_pr
        wls= +.1300000d+03; sigma = one/six
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -265.2500_pr, -340.0625_pr /)
    !---------------------------------------------------------------------
    ! SKP force, Dobaczewski et al., NPA422 (1984) 103
    !---------------------------------------------------------------------
    Case ('SKP')
        ! ph-Force
        noForces=1
        use_cm_cor  = .True.
        use_j2terms = .True.
        hbzero = 20.730_pr
        t0 =-0.2931696d+04; x0 = 0.2921515_pr
        t1 = 0.3206182d+03; x1 = 0.6531765_pr
        t2 =-0.3374091d+03; x2 =-0.5373230_pr
        t3 = 0.1870896d+05; x3 = 0.1810269_pr
        wls= 0.1000000d+03; sigma=one/six
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -265.2500_pr, -340.0625_pr /)
    !---------------------------------------------------------------------
    ! SLY4 force
    !---------------------------------------------------------------------
    Case ('SLY4')
        ! ph-Force
        noForces=1
        use_cm_cor = .True.
        hbzero = 20.735530_pr
        t0 =-0.2488913d+04; x0 = 0.8340000_pr
        t1 = 0.4868180d+03; x1 =-0.3440000_pr
        t2 =-0.5463950d+03; x2 =-1.0000000_pr
        t3 = 0.1377700d+05; x3 = 1.3540000_pr
        wls= 0.1230000d+03; sigma=one/six
        b4 = wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -325.2500_pr, -340.0625_pr /) ! HFB
    !---------------------------------------------------------------------
    ! SLY5 force
    !---------------------------------------------------------------------
    Case ('SLY5')
        ! ph-Force
        noForces=1
        use_cm_cor  = .True.
        use_j2terms = .True.
        hbzero = 20.73553_pr
        t0 =-0.2483450d+04; x0 = 0.7760000_pr
        t1 = 0.4842300d+03; x1 =-0.3170000_pr
        t2 =-0.5566900d+03; x2 =-1.0000000_pr
        t3 = 0.1375700d+05; x3 = 1.2630000_pr
        wls= 0.1250000d+03; sigma=one/six
        b4 = wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -291.5000_pr, -297.7402_pr /) ! HFB
    !---------------------------------------------------------------------
    ! SLY6 forces
    !---------------------------------------------------------------------
    Case ('SLY6')
        ! ph-Force
        noForces=1
        use_cm_cor      = .True.
        use_full_cm_cor = .True.
        hbzero = 20.73553_pr
        t0 =-0.2479500d+04; x0 = 0.8250000_pr
        t1 = 0.4621800d+03; x1 =-0.4650000_pr
        t2 =-0.4486100d+03; x2 =-1.0000000_pr
        t3 = 0.1367300d+05; x3 = 1.3550000_pr
        wls= 0.1220000d+03; sigma=one/six
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -291.5000_pr, -297.7402_pr /) ! HFB
    !---------------------------------------------------------------------
    ! SLY6 forces
    !---------------------------------------------------------------------
    Case ('SLY7')
        ! ph-Force
        noForces=1
        use_cm_cor      = .True.
        use_j2terms     = .True.
        use_full_cm_cor = .True.
        hbzero = 20.73553_pr
        t0 =-0.2480800d+04; x0 = 0.8480000_pr
        t1 = 0.4612900d+03; x1 =-0.4920000_pr
        t2 =-0.4339300d+03; x2 =-1.0000000_pr
        t3 = 0.1366900d+05; x3 = 1.3930000_pr
        wls= 0.1250000d+03; sigma=one/six
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -291.5000_pr, -297.7402_pr /) ! HFB
    !---------------------------------------------------------------------
    ! SKI3 force, P.G.-Reinhard et al. Nucl. Phys. A584 (1995) 467-488
    !---------------------------------------------------------------------
    Case ('SKI3')
        ! ph-Force
        noForces=1
        use_cm_cor      = .True.
        use_full_cm_cor = .True.
        hbzero = 20.7525d0
        t0 =-0.176288d+04; x0 = 0.30830_pr
        t1 = 0.561608d+03; x1 =-1.17220_pr
        t2 =-0.227090d+03; x2 =-1.09070_pr
        t3 = 0.810620d+04; x3 = 1.29260_pr
        sigma=one/four
        b4 = 94.254_pr; b4p=zero
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -357.2324_pr, -388.5625_pr /)
    !---------------------------------------------------------------------
    ! SKO forces
    !---------------------------------------------------------------------
    Case ('SKO')
        ! ph-Force
        noForces=1
        use_cm_cor      = .True.
        use_full_cm_cor = .True.
        hbzero = 20.735530_pr
        t0 =-0.21036530d+04; x0 = -0.2107010_pr
        t1 = 0.30335200d+03; x1 = -2.8107520_pr
        t2 = 0.79167400d+03; x2 = -1.4615950_pr
        t3 = 0.13553252d+05; x3 = -0.4298810_pr
        wls= 0.12300000d+03; sigma=one/four
        b4 = 0.17657800d+03; b4p=-0.1987490d+03
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -259.0391_pr, -274.8433_pr /)
    !---------------------------------------------------------------------
    ! SKX forces, A.Brown; Phys.Rev. C58 (1998) 220
    !---------------------------------------------------------------------
    Case ('SKX')
        ! ph-Force
        noForces=1
        use_cm_cor  = .True.
        use_j2terms = .True.
        hbzero = 20.73_pr
        t0 = -1445.300_pr; x0 = 0.340_pr
        t1 =   246.900_pr; x1 = 0.580_pr
        t2 =  -131.800_pr; x2 = 0.127_pr
        t3 = 12103.900_pr; x3 = 0.030_pr
        sigma=one/two
        b4 = 0.0743d+03; b4p=zero
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -259.0391_pr, -274.8433_pr /)
    !---------------------------------------------------------------------
    ! HFB9 forces
    !---------------------------------------------------------------------
    Case ('HFB9')
        ! ph-Force
        noForces=1
        use_cm_cor  = .True.
        use_j2terms = .True.
        hbzero = 20.73553_pr
        t0 =-0.20439180d+04; x0 = 0.5149210_pr
        t1 = 0.41159870d+03; x1 =-0.9537990_pr
        t2 =-0.19418860d+03; x2 =-0.3322490_pr
        t3 = 0.12497170d+05; x3 = 0.8994350_pr
        wls= 0.14990000d+03; sigma=one/four
        b4=wls/two; b4p=wls/two
        ! pp-Forces
        CpV1=0.50_pr
        CpV0=(/ -263.5000_pr, -274.9668_pr /)
    !---------------------------------------------------------------------
    ! PRC 21, 1568 (1980)
    !---------------------------------------------------------------------
    Case ('D1')
        ! ph-Force
        noForces=1
        use_cm_cor   = .True.
        use_j2terms  = .False.
        finite_range = .True.
        hbzero = 20.73667552957479_pr
        t3 = 6.0_pr*1350.00_pr; x3 = one;
        wls= 115.000_pr;
        sigma=one/three
        b4=wls/two; b4p=wls/two
        ! No delta-pairing here. The pairing will be in the finite range part
        CpV1= zero
        CpV0= zero
        ! W, B, H, M parameters of the finite-range two-body
        Call gogny_force(fname)
    !---------------------------------------------------------------------
    ! CPC 63, 365 (1991)
    !---------------------------------------------------------------------
    Case ('D1S')
        ! ph-Force
        noForces=1
        use_cm_cor   = .True.
        use_j2terms  = .False.
        finite_range = .True.
        hbzero = 20.73667552957479_pr
        t3 = 6.0_pr*1390.600_pr; x3 = one;
        wls= 130.000_pr;
        sigma=one/three
        b4=wls/two; b4p=wls/two
        ! No delta-pairing here. The pairing will be in the finite range part
        CpV1= zero
        CpV0= zero
        ! W, B, H, M parameters of the finite-range two-body
        Call gogny_force(fname)
    !---------------------------------------------------------------------
    ! PRC 21, 1568 (1980)
    !---------------------------------------------------------------------
    Case ('D1p')
        ! ph-Force
        noForces=1
        use_cm_cor   = .True.
        use_j2terms  = .False.
        finite_range = .True.
        hbzero = 20.73667552957479_pr
        t3 = 6.0_pr*1350.00_pr; x3 = one;
        wls= 130.000_pr;
        sigma=one/three
        b4=wls/two; b4p=wls/two
        ! No delta-pairing here. The pairing will be in the finite range part
        CpV1= zero
        CpV0= zero
        ! W, B, H, M parameters of the finite-range two-body
        Call gogny_force(fname)
    !---------------------------------------------------------------------
    !
    !---------------------------------------------------------------------
    Case ('D1N')
        ! ph-Force
        noForces=1
        use_cm_cor   = .True.
        use_j2terms  = .False.
        finite_range = .True.
        hbzero = 20.73667552957479_pr
        t3 = 6.0_pr*1609.50_pr; x3 = one;
        wls= 115.000_pr;
        sigma=one/three
        b4=wls/two; b4p=wls/two
        ! No delta-pairing here. The pairing will be in the finite range part
        CpV1= zero
        CpV0= zero
        ! W, B, H, M parameters of the finite-range two-body
        Call gogny_force(fname)
    !---------------------------------------------------------------------
    !T0X0 A very simple skyrme functional
    !---------------------------------------------------------------------
    Case ('T0X0')
        noForces=1
        use_cm_cor  = .True.
        t0 =  0.0_pr; x0 = 0.0_pr;
!        t0 =  -1128.75_pr; x0 = 0.45_pr;
        hbzero = 20.73667622931579_pr
        ! pp-Forces
        CpV1= zero
        CpV0= zero !No pairing here.
    !---------------------------------------------------------------------
    ! Default
    !---------------------------------------------------------------------
    Case default
        Write(6,'("No Skyrme interaction defined in routine skforce()")')
    End Select INTERACTION
    !
    If (noForces.Eq.1) Then
        ! obtain coupling constants
        Call C_from_t()
        ! Frequent combinations entering the energy
        tv1   =  t0*(one+half*x0)*half;    tv2 = t0*(x0+half)*half
        tv3   =  t3*(one+half*x3)/12.0_pr; tv4 = t3*(x3+half)/12.0_pr
        tv5   = (t1*(one+half*x1)+t2*(one+half*x2))/four
        tv6   = (t2*(half+x2)-t1*(half+x1))/four
        ts1   = (t2*(one+half*x2)-three*t1*(one+half*x1))/pp16
        ts2   = (t1*(half+x1)*three+t2*(half+x2))/pp16
        t4o3  =  four/three; t324 = t3/pp24
        ! Frequent combinations entering the potential
        t0s   =  t0*(one-x0)*half; t0a = t0*(one+x0*half)
        drs   = (t2*(one+x2)-t1*(one-x1))*three/pp16
        dra   = (t2*(one+half*x2)-three*t1*(one+half*x1))/eight
        ts    = (t1*(one-x1) + three*t2*(one+x2))/eight
        ta    = (t1*(one+half*x1) + t2*(one+half*x2))/four
        t3alp = t3*(two+sigma)*(two+x3)/pp24
        t3al0 = t3*(x3+half)/six; t3alm = t3*sigma*(one+two*x3)/pp24
        alp   = one + sigma; alm = sigma - one
        wla0  = CrdJ(0)+CrdJ(1); wla1  = CrdJ(0)-CrdJ(1);
        TA7   = zero; TA8 = zero
        If(use_j2terms) Then
           TA7=(T1*(ONE-X1)-T2*(ONE+X2))/eight + five*to/four
           TA8=-(T1*X1+T2*X2)/four             + five*(te+to)/four
        End If
        TB7 = TA7; TB8 = TA8*half
    End If
    !
    Return
  End Subroutine skforce
  !=======================================================================
  !> Define functional parameters
  !=======================================================================
  Subroutine set_functional_parameters(fname,lpr)
    Implicit None
    Logical, Intent(in) :: lpr
    Character (30), Intent(inout) :: fname
    Logical :: regularization
    Integer(ipr), Parameter :: lin=15
    !
    ! parameters
    FunctionalName=fname
    eps=Spacing(1.0_pr)
    Pi=4.0_pr*Atan(1.0_pr)
    kfconst=(1.50_pr*Pi**2)**(1.0_pr/3.0_pr)    ! (3Pi^2/2)^(1/3)
    CK=3.0_pr/5.0_pr*kfconst**2
    nuLambda=700.0_pr ; nufpi = 93.0_pr
    !
    Call Make_Parameter_Free_Useful_Combinations()
    !
    ! exact Hartree CHrho from INM
    CHrho=0.0_pr; !!!!If (dmeorder.eq.3) Call CHrho_from_NM()
    !
    If(use_INM) Then
       Call calculate_C_from_NM(E_NM,K_NM,SMASS_NM,RHO_NM,ASS_NM,LASS_NM,VMASS_NM)
    Else
       Crho(0)=Crho(0)+CHrho
    End If
    If(.Not.finite_range) Call calculate_NM_properties()
    !
    Crho(0)=Crho(0)-CHrho
    !
    Call calculate_natural_units()
    !
    ! Print output
    !If(lpr) Then
    !   Call print_functional_parameters()
    !End If
    !
  End Subroutine set_functional_parameters
  !=======================================================================
  !
  !=======================================================================
  Subroutine print_functional_parameters()
    Use HFBTHO_utilities, Only: lout,lfile
    Implicit None
    Integer(ipr) :: iw
    !
    Do iw=lout,lfile
       Write(iw,'("  ---------------------------------------")')
       Write(iw,'("           UNEDF Module Version: ",a)') Trim(Version)
       Write(iw,'("         M.Kortelainen & M.Stoitsov ")')
       Write(iw,'("  ---------------------------------------")')
       !
       Write(iw,*)
       Write(iw,'(2x,a," functional")') Trim(FunctionalName)
       Write(iw,'("  ----------------------------------------")')
       Write(iw,'("  Crho(0)= ",g26.18,"; Crho(1)= ",g26.18)') Crho
       Write(iw,'("  CDrho(0)=",g26.18,"; CDrho(1)=",g26.18)') CDrho
       Write(iw,'("  Ctau(0)= ",g26.18,"; Ctau(1)= ",g26.18)') Ctau
       Write(iw,'("  CrDr(0)= ",g26.18,"; CrDr(1)= ",g26.18)') Crdr
       Write(iw,'("  CrdJ(0)= ",g26.18,"; CrdJ(1)= ",g26.18)') CrdJ
       Write(iw,'("  CJ(0)=   ",g26.18,"; CJ(1)=   ",g26.18)') CJ
       Write(iw,'("  CpV0(0)= ",g26.18,"; CpV0(1)= ",g26.18)') CpV0
       Write(iw,'("  CpV1(0)= ",g26.18,"; CpV1(1)= ",g26.18)') CpV1
       Write(iw,'("  sigma=   ",g26.18,"; hbzero=  ",g26.18)') sigma,hbzero
       Write(iw,'("  e^2 chrg=",g26.18,"; CExPar=  ",g26.18)') e2charg,CExPar
       Write(iw,'("  c.m. correction: ",L1,", chr. density in direct Coul: ",L1)') use_cm_cor,use_charge_density
       Write(iw,'("  use tensor terms: ",L1)') use_j2terms
       ! Finite-range force (Gogny force)
       If(finite_range) Then
          Write(iw,*)
          Write(iw,'("  Finite-range potential")')
          Write(iw,'("  ----------------------------------------")')
          Write(iw,'("  mu=",2f26.18)') mu_g
          Write(iw,'("  W=",2f26.18)') W_g
          Write(iw,'("  B=",2f26.18)') B_g
          Write(iw,'("  H=",2f26.18)') H_g
          Write(iw,'("  M=",2f26.18)') M_g
       End If
       ! Natural units
       Write(iw,*)
       Write(iw,'("  Coupling constants in natural units")')
       Write(iw,'("  ----------------------------------------")')
       If(.Not.finite_range) Then
          Write(iw,'("  Crho_nu(0)= ",g26.18,"; Crho_nu(1)= ",g26.18)') nuCrho
       End If
       Write(iw,'("  CDrho_nu(0)=",g26.18,"; CDrho_nu(1)=",g26.18)') nuCDrho
       If(.Not.finite_range) Then
          Write(iw,'("  Ctau_nu(0)= ",g26.18,"; Ctau_nu(1)= ",g26.18)') nuCtau
          Write(iw,'("  CrDr_nu(0)= ",g26.18,"; CrDr_nu(1)= ",g26.18)') nuCrdr
       End If
       Write(iw,'("  CrdJ_nu(0)= ",g26.18,"; CrdJ_nu(1)= ",g26.18)') nuCrdJ
       If(.Not.finite_range) Then
          Write(iw,'("  CJ_nu(0)=   ",g26.18,"; CJ_nu(1)=   ",g26.18)') nuCJ
          Write(iw,'("  CpV0_nu(0)= ",g26.18,"; CpV0_nu(1)= ",g26.18)') nuCpV0
          Write(iw,'("  CpV1_nu(0)= ",g26.18,"; CpV1_nu(1)= ",g26.18)') nuCpV1
       End If
       Write(iw,'("  fpi_nu=     ",g26.18,"; Lambda_nu=  ",g26.18)') nufpi,nuLambda
       ! DME
       If(dmeorder.Ge.0) Then
          Write(iw,*)
          Write(iw,'("  DME parameters")')
          Write(iw,'("  ----------------------------------------")')
          Write(iw,'("       gA=",f12.6," mpi [1/fm]=",f12.6," fpi [1/fm]=",f12.6)') gA,mpi,fpi
          Write(iw,'("  c1 [fm]=",f12.6,"    c3 [fm]=",f12.6,"    c4 [fm]=",f12.6)') c1,c3,c4
          Write(iw,'("       cd=",f12.6,"         ce=",f12.6," LamX[1/fm]=",f12.6)') cd,ce,LambdaX
          Write(iw,'("  ->CHrho=",f12.6)') CHrho
          If(dmeorder.Ge.2) Write(iw,'("  use 3N terms: ",L1)') use_DME3N_terms
       End If
       ! Nuclear matter (not implemented for finite-range potential yet)
       If(.Not.finite_range) Then
          Write(iw,*)
          Write(iw,'("  Nuclear matter properties")')
          Write(iw,'("  ----------------------------------------")')
          Write(iw,'("  E_NM=    ",g26.18,"; K_NM=     ",g26.18)') E_NM,K_NM
          Write(iw,'("  P_NM=    ",g26.18,"; RHO_NM=   ",g26.18)') P_NM,RHO_NM
          Write(iw,'("  ASS_NM=  ",g26.18,"; LASS_NM=  ",g26.18)') ASS_NM,LASS_NM
          Write(iw,'("  SMASS_NM=",g26.18,"; VMASS_NM= ",g26.18)') SMASS_NM,VMASS_NM
          !
          Call t_from_C()
          !
       End If
       ! (t,x) parametrization of the Skyrme functional
          Write(iw,*)
       Write(iw,'("  Associated (t,x)-coupling constants")')
          Write(iw,'("  ----------------------------------------")')
       If(.Not.finite_range) Then
          Write(iw,'("  t0=    ",g26.18,"; x0=     ",g26.18)') t0,x0
          Write(iw,'("  t1=    ",g26.18,"; x1=     ",g26.18)') t1,x1
          Write(iw,'("  t2=    ",g26.18,"; x2=     ",g26.18)') t2,x2
       End If
       Write(iw,'("  t3=    ",g26.18,"; x3=     ",g26.18)') t3,x3
       Write(iw,'("  b4=    ",g26.18,"; b4p=    ",g26.18)') b4,b4p
       Write(iw,'("  te=    ",g26.18,"; to=     ",g26.18)') te,to
       Write(iw,'("  sigma= ",g26.18,"; hbzero= ",g26.18)') sigma,hbzero
       !
       If(Print_Namelist) Then
          Write(iw,*)
          SELECTED_FUNCTIONAL: Select Case (Trim(FunctionalName))
          Case ("UNEDF","SKYRME")
                Write(iw,'("NAMELIST CONTENT (cannot be modified for functional names UNEDF,SKYRME)")')
                Write(iw,'("-----------------------------------------------------------------------")')
          Case ("FITS")
                Write(iw,'("NAMELIST CONTENT (advanced usage: modify all but not C-, NM-, and more...)")')
                Write(iw,'("--------------------------------------------------------------------------")')
          Case default
                Write(iw,'("NAMELIST CONTENT (copy/past to UNEDF_NAMELIST.DAT and modify)")')
                Write(iw,'("-------------------------------------------------------------")')
          End Select SELECTED_FUNCTIONAL
          Write(*,'(" !NB: FUNCTIONALNAME should be always in quotations")')
          Write(*,UNEDF_NAMELIST)
       End If
    End Do
  End Subroutine print_functional_parameters
  !=======================================================================
  !> Calculates coupling constants in natural units
  !=======================================================================
  Subroutine calculate_natural_units
    Implicit None
    nuCrho = Crho*(nufpi**2)/(mevfm**3)
    nuCdrho = Cdrho*(nufpi**2)*((nuLambda*nufpi*nufpi)**sigma)/(mevfm**(3.0_pr*(1.0_pr+sigma)))
    nuCtau = Ctau*((nufpi*nuLambda)**2)/(mevfm**5)
    nuCrDr = CrDr*((nufpi*nuLambda)**2)/(mevfm**5)
    nuCrdJ = CrdJ*((nufpi*nuLambda)**2)/(mevfm**5)
    nuCJ = CJ*((nufpi*nuLambda)**2)/(mevfm**5)
    nuCpV0 = CpV0*(nufpi**2)/(mevfm**3)
    nuCpV1 = CpV1*(nufpi**4)*nuLambda/(mevfm**6)
  End Subroutine calculate_natural_units
  !=======================================================================
  !> Calculates volume C-constants (and sigma) form NM properties
  !>
  !> Input: E,K,SMASS,RHO,ASS,LASS,VMASS,sigma_NM(optional)
  !>
  !> Output: Crho(0),Crho(1),Cdrho(0),Cdrho(1),Ctau(0),Ctau(0),sigma(optional)
  !>
  !> Options:
  !>  - When sigma_NM exists then 'sigma'=sigma_NM
  !>  - When sigma_NM does not exist then 'sigma' is defined from NM
  !=======================================================================
  Subroutine calculate_C_from_NM(E,K,SMASS,RHO,ASS,LASS,VMASS,sigma_NM)
    Implicit None
    Real(pr), Intent(in) :: E,K,SMASS,RHO,ASS,LASS,VMASS
    Real(pr), Intent(in), Optional :: sigma_NM
    Real(pr) :: aRho0Rho0,daRho0Rho0,ddaRho0Rho0,aRho1Rho1,daRho1Rho1,ddaRho1Rho1
    Real(pr) :: aRho0Tau0,daRho0Tau0,ddaRho0Tau0,aRho1Tau1,daRho1Tau1,ddaRho1Tau1
    Real(pr) :: u,tauc,rho2
    Real(pr),Parameter :: c13=1.0_pr/3.0_pr,c23=2.0_pr/3.0_pr
    !
    tauc=CK*RHO**c23; u=(kfconst/mpi)*RHO**c13; rho2=rho**2
    !
    Call calculate_U_parameters(RHO,RHO,tauc*RHO,tauc*RHO,0.0_pr,0.0_pr)
    !
    aRho0Rho0=0.50_pr*(aRhoRho+bRhoRho)*mevfm
    aRho1Rho1=0.50_pr*(aRhoRho-bRhoRho)*mevfm
    aRho0Tau0=0.50_pr*(aRhoTau+bRhoTau)*mevfm
    aRho1Tau1=0.50_pr*(aRhoTau-bRhoTau)*mevfm
    daRho0Rho0=0.50_pr*(daRhoRho+dbRhoRho)*mevfm
    daRho1Rho1=0.50_pr*(daRhoRho-dbRhoRho)*mevfm
    daRho0Tau0=0.50_pr*(daRhoTau+dbRhoTau)*mevfm
    daRho1Tau1=0.50_pr*(daRhoTau-dbRhoTau)*mevfm
    ddaRho0Rho0=0.50_pr*(ddaRhoRho+ddbRhoRho)*mevfm
    ddaRho1Rho1=0.50_pr*(ddaRhoRho-ddbRhoRho)*mevfm
    ddaRho0Tau0=0.50_pr*(ddaRhoTau+ddbRhoTau)*mevfm
    ddaRho1Tau1=0.50_pr*(ddaRhoTau-ddbRhoTau)*mevfm
    !
    ! set/calculate sigma
    If (Present(sigma_NM)) Then
        sigma=sigma_NM
    Else
        sigma=((1.0_pr/3.0_pr)*(-K+tauc*hbzero*(-3.0_pr+4.0_pr*SMASS)-9.0_pr*E+9.0_pr*RHO2*hRho0Rho0 &
             +21.0_pr*tauc*RHO2*hRho0Tau0+u*RHO*(daRho0Rho0+5.0_pr*tauc*daRho0Tau0 &
             +7.0_pr*RHO*dhRho0Rho0+11.0_pr*tauc*RHO*dhRho0Tau0+u*ddaRho0Rho0 &
             +u*tauc*ddaRho0Tau0+u*RHO*ddhRho0Rho0+u*tauc*RHO*ddhRho0Tau0))) &
             /(tauc*hbzero*(-3.0_pr+2.0_pr*SMASS)+3.0_pr*E+3.0_pr*RHO2*hRho0Rho0 &
             +3.0_pr*tauc*RHO2*hRho0Tau0+u*RHO*(daRho0Rho0+tauc*daRho0Tau0 &
             + RHO*dhRho0Rho0+tauc*RHO*dhRho0Tau0))
    End If
    !
    Crho(0)=(c13*(tauc*hbzero*(-3.0_pr+(2.0_pr-3.0_pr*sigma)*SMASS) &
        +3.0_pr*(1.0_pr+sigma)*E-3.0_pr*sigma*RHO*aRho0Rho0 &
        +3.0_pr*(1.0_pr-sigma)*RHO2*hRho0Rho0+3.0_pr*tauc*RHO2*hRho0Tau0 &
        +u*RHO*(daRho0Rho0+tauc*daRho0Tau0+RHO*dhRho0Rho0 &
        +tauc*RHO*dhRho0Tau0)))/(sigma*RHO)
    Cdrho(0)=(c13*RHO**(-1.0_pr-sigma)*(tauc*hbzero*(3.0_pr-2.0_pr*SMASS)&
        -3.0_pr*E-3.0_pr*RHO**2*hRho0Rho0-3.0_pr*tauc*RHO2*hRho0Tau0&
        -u*RHO*(daRho0Rho0+tauc*daRho0Tau0+RHO*dhRho0Rho0 &
        +tauc*RHO*dhRho0Tau0)))/sigma
    Ctau(0)=(hbzero*(SMASS-1.0_pr)-RHO*(aRho0Tau0+RHO*hRho0Tau0))/RHO
    !
    Crho(1)=(27.0_pr*ASS*(1.0_pr+sigma)-9.0_pr*LASS &
        +5.0_pr*tauc*hbzero*(5.0_pr-6.0_pr*VMASS+3.0_pr*sigma*(-4.0_pr+3.0_pr*VMASS)) &
        +20.0_pr*tauc*(2.0_pr-3.0_pr*sigma)*RHO*aRho0Tau0 &
        +RHO*(-27.0_pr*sigma*aRho1Rho1+5.0_pr*tauc*(11.0_pr-12.0_pr*sigma)*RHO*hRho0Tau0 &
        -27.0_pr*(-1.0_pr+sigma)*RHO*hRho1Rho1+9.0_pr*tauc*(5.0_pr-3.0_pr*sigma)*RHO*hRho1Tau0 &
        +45.0_pr*tauc*RHO*hRho1Tau1+40.0_pr*tauc*Ctau(0)-60.0_pr*tauc*sigma*Ctau(0) &
        +5.0_pr*u*tauc*daRho0Tau0+9.0_pr*u*daRho1Rho1+15.0_pr*u*tauc*daRho1Tau1 &
        +5.0_pr*u*tauc*RHO*dhRho0Tau0+9.0_pr*u*RHO*dhRho1Rho1+9.0_pr*u*tauc*RHO*dhRho1Tau0 &
        +15.0_pr*u*tauc*RHO*dhRho1Tau1))/(27.0_pr*sigma*RHO)
    Cdrho(1)=-(RHO**(-1.0_pr-sigma)*(27.0_pr*ASS-9.0_pr*LASS &
        +5.0_pr*tauc*hbzero*(5.0_pr-6.0_pr*VMASS)+40.0_pr*tauc*RHO*aRho0Tau0 &
        +55.0_pr*tauc*RHO2*hRho0Tau0+27.0_pr*RHO**2*hRho1Rho1+45.0_pr*tauc*RHO2*hRho1Tau0 &
        +45.0_pr*tauc*RHO2*hRho1Tau1+40.0_pr*tauc*RHO*Ctau(0) +5.0_pr*u*tauc*RHO*daRho0Tau0 &
        +9.0_pr*u*RHO*daRho1Rho1+15.0_pr*u*tauc*RHO*daRho1Tau1 &
        +5.0_pr*u*tauc*RHO2*dhRho0Tau0+9.0_pr*u*RHO2*dhRho1Rho1 &
        +9.0_pr*u*tauc*RHO2*dhRho1Tau0 +15.0_pr*u*tauc*RHO2*dhRho1Tau1))/(27.0_pr*sigma)
    Ctau(1)=(hbzero-hbzero*VMASS+RHO*(aRho0Tau0-aRho1Tau1+RHO*hRho0Tau0-RHO*hRho1Tau1+Ctau(0)))/RHO
    !
  End Subroutine calculate_C_from_NM
  !=======================================================================
  !> Calculates INM properties
  !=======================================================================
  Subroutine calculate_NM_properties()
    Implicit None
    Real(pr) :: aRho0Rho0,daRho0Rho0,ddaRho0Rho0,aRho1Rho1,daRho1Rho1,ddaRho1Rho1
    Real(pr) :: aRho0Tau0,daRho0Tau0,ddaRho0Tau0,aRho1Tau1,daRho1Tau1,ddaRho1Tau1
    Real(pr) :: u,tauc,rho_NM2
    Real(pr), Parameter :: c13=1.0_pr/3.0_pr,c23=2.0_pr/3.0_pr
    !
    RHO_NM=find_NM_RHOC()
    !
    aRho0Rho0=0.50_pr*(aRhoRho+bRhoRho)*mevfm
    aRho1Rho1=0.50_pr*(aRhoRho-bRhoRho)*mevfm
    aRho0Tau0=0.50_pr*(aRhoTau+bRhoTau)*mevfm
    aRho1Tau1=0.50_pr*(aRhoTau-bRhoTau)*mevfm
    daRho0Rho0=0.50_pr*(daRhoRho+dbRhoRho)*mevfm
    daRho1Rho1=0.50_pr*(daRhoRho-dbRhoRho)*mevfm
    daRho0Tau0=0.50_pr*(daRhoTau+dbRhoTau)*mevfm
    daRho1Tau1=0.50_pr*(daRhoTau-dbRhoTau)*mevfm
    ddaRho0Rho0=0.50_pr*(ddaRhoRho+ddbRhoRho)*mevfm
    ddaRho1Rho1=0.50_pr*(ddaRhoRho-ddbRhoRho)*mevfm
    ddaRho0Tau0=0.50_pr*(ddaRhoTau+ddbRhoTau)*mevfm
    ddaRho1Tau1=0.50_pr*(ddaRhoTau-ddbRhoTau)*mevfm
    tauc=CK*RHO_NM**c23; u=(kfconst/mpi)*RHO_NM**c13; rho_NM2=rho_NM**2
    !
    ! Symmetric Nuclear Matter
    E_NM=tauc*hbzero+RHO_NM*(aRho0Rho0+RHO_NM*hRho0Rho0+Crho(0)+RHO_NM**sigma*Cdrho(0)) &
      +tauc*RHO_NM*(aRho0Tau0+RHO_NM*hRho0Tau0+Ctau(0))
    P_NM=c13*RHO_NM**2*((2.0_pr*tauc*hbzero)/RHO_NM+3.0_pr*aRho0Rho0+5.0_pr*tauc*aRho0Tau0 &
      +6.0_pr*RHO_NM*hRho0Rho0+8.0_pr*tauc*RHO_NM*hRho0Tau0+3.0_pr*Crho(0) &
      +3.0_pr*(1+sigma)*RHO_NM**sigma*Cdrho(0)+5.0_pr*tauc*Ctau(0)+u*daRho0Rho0 &
      +u*tauc*daRho0Tau0+u*RHO_NM*dhRho0Rho0+u*tauc*RHO_NM*dhRho0Tau0)
    SMASS_NM=1.0_pr+(RHO_NM*(aRho0Tau0+RHO_NM*hRho0Tau0+Ctau(0)))/hbzero
    K_NM=9.0_pr*sigma*(1+sigma)*RHO_NM**(1+sigma)*Cdrho(0) &
      +(-2.0_pr*tauc*hbzero+10.0_pr*tauc*RHO_NM*aRho0Tau0+18.0_pr*RHO_NM2*hRho0Rho0 &
      +40.0_pr*tauc*RHO_NM**2*hRho0Tau0+4.0_pr*u*RHO_NM*daRho0Rho0 &
      +RHO_NM*(10.0_pr*tauc*Ctau(0)+u*(8.0_pr*tauc*daRho0Tau0+u*ddaRho0Rho0 &
      +(10.0_pr*RHO_NM*dhRho0Rho0+14.0_pr*tauc*RHO_NM*dhRho0Tau0 &
      +(u*tauc*ddaRho0Tau0+u*RHO_NM*ddhRho0Rho0+u*tauc*RHO_NM*ddhRho0Tau0)))))
    !
    ! Asymmetric Nuclear Matter
    ASS_NM=RHO_NM2*hRho1Rho1+RHO_NM*(aRho1Rho1+Crho(1)+RHO_NM**sigma*Cdrho(1)) &
       +(tauc*(5.0_pr*hbzero+RHO_NM*(5.0_pr*aRho0Tau0+15.0_pr*aRho1Tau1+5.0_pr*RHO_NM*hRho0Tau0 &
       +9.0_pr*RHO_NM*hRho1Tau0+5.0_pr*(3.0_pr*RHO_NM*hRho1Tau1+Ctau(0)+3.0_pr*Ctau(1)))))/9.0_pr
    VMASS_NM=(hbzero+RHO_NM*(aRho0Tau0-aRho1Tau1+RHO_NM*hRho0Tau0-RHO_NM*hRho1Tau1+Ctau(0)-Ctau(1)))/hbzero
    LASS_NM=6.0_pr*RHO_NM2*hRho1Rho1+3.0_pr*RHO_NM*(aRho1Rho1+Crho(1)+(1.0_pr+sigma)*RHO_NM**sigma*Cdrho(1)) &
       +u*RHO_NM*daRho1Rho1 +u*RHO_NM2*dhRho1Rho1 &
       +(tauc*(10.0_pr*hbzero+8.0_pr*RHO_NM2*(5.0_pr*hRho0Tau0+9.0_pr*hRho1Tau0+15.0_pr*hRho1Tau1) &
       +25.0_pr*RHO_NM*(aRho0Tau0+3.0_pr*aRho1Tau1+Ctau(0)+3*Ctau(1)) &
       +5.0_pr*u*RHO_NM*(daRho0Tau0+3.0_pr*daRho1Tau1) &
       +u*RHO_NM2*(5.0_pr*dhRho0Tau0+9.0_pr*dhRho1Tau0+15.0_pr*dhRho1Tau1)))/9.0_pr
    KA_NM=18.0_pr*RHO_NM2*hRho1Rho1+9.0_pr*sigma*(1.0_pr+sigma)*RHO_NM**(1.0_pr+sigma)*Cdrho(1) &
       +4.0_pr*u*RHO_NM*daRho1Rho1 +10.0_pr*u*RHO_NM2*dhRho1Rho1 &
       + u**2*RHO_NM*ddaRho1Rho1+u**2*RHO_NM2*ddhRho1Rho1 &
       +(tauc*(-10.0_pr*hbzero+40.0_pr*RHO_NM2*(5.0_pr*hRho0Tau0+9.0_pr*hRho1Tau0+15.0_pr*hRho1Tau1) &
       +50.0_pr*RHO_NM*(aRho0Tau0+3.0_pr*aRho1Tau1+Ctau(0)+3*Ctau(1)) &
       +40.0_pr*u*RHO_NM*(daRho0Tau0+3.0_pr*daRho1Tau1) &
       +14.0_pr*u*RHO_NM2*(5.0_pr*dhRho0Tau0+9.0_pr*dhRho1Tau0 &
       +15.0_pr*dhRho1Tau1)+5.0_pr*u**2*RHO_NM*(ddaRho0Tau0 &
       +3.0_pr*ddaRho1Tau1)+u**2*RHO_NM2*(5.0_pr*ddhRho0Tau0+9*ddhRho1Tau0+15*ddhRho1Tau1)))/9.
     !
  End Subroutine calculate_NM_properties
  !=======================================================================
  !> Find the INM saturation density RHO_NM using the Secant Method
  !=======================================================================
  Real(pr) Function find_NM_RHOC()
    Implicit None
    Integer(pr) :: iter
    Real(pr) :: aRho0Rho0,daRho0Rho0,ddaRho0Rho0,aRho1Rho1,daRho1Rho1,ddaRho1Rho1
    Real(pr) :: aRho0Tau0,daRho0Tau0,ddaRho0Tau0,aRho1Tau1,daRho1Tau1,ddaRho1Tau1
    Real(pr) :: kfconstmpi,u,tauc
    Real(pr) :: rhom0,rhom,rhom2,w,w0,step,energy
    Real(pr),Parameter :: c13=1.0_pr/3.0_pr,c23=2.0_pr/3.0_pr
    !
    kfconstmpi=kfconst/mpi; step=-0.0010_pr; iter=0
    ! initial point
    rhom=0.170_pr; tauc=CK*rhom**c23; u=kfconstmpi*rhom**c13; rhom2=rhom**2
    !
    Call calculate_U_parameters(rhom,rhom,tauc*rhom,tauc*rhom,0.0_pr,0.0_pr)
    !
    aRho0Rho0=0.50_pr*(aRhoRho+bRhoRho)*mevfm; daRho0Rho0=0.50_pr*(daRhoRho+dbRhoRho)*mevfm
    aRho0Tau0=0.50_pr*(aRhoTau+bRhoTau)*mevfm; daRho0Tau0=0.50_pr*(daRhoTau+dbRhoTau)*mevfm
    w0=c13*rhom2*((2.0_pr*tauc*hbzero)/rhom+3.0_pr*aRho0Rho0+5.0_pr*tauc*aRho0Tau0 &
      +6.0_pr*rhom*hRho0Rho0+8.0_pr*tauc*rhom*hRho0Tau0+3.0_pr*Crho(0) &
      +3.0_pr*(1.0_pr+sigma)*rhom**sigma*Cdrho(0)+5.0_pr*tauc*Ctau(0)+u*daRho0Rho0 &
      +u*tauc*daRho0Tau0+u*rhom*dhRho0Rho0+u*tauc*rhom*dhRho0Tau0)
    rhom0=rhom; rhom=rhom+step
    !
    ! secant method
    Do While(Abs(step).Ge.eps*100.0_pr)
       iter=iter+1
       tauc=CK*rhom**c23; u=kfconstmpi*rhom**c13; rhom2=rhom**2
       !
       Call calculate_U_parameters(rhom,rhom,tauc*rhom,tauc*rhom,0.0_pr,0.0_pr)
       !
       aRho0Rho0=0.50_pr*(aRhoRho+bRhoRho)*mevfm; daRho0Rho0=0.50_pr*(daRhoRho+dbRhoRho)*mevfm
       aRho0Tau0=0.50_pr*(aRhoTau+bRhoTau)*mevfm; daRho0Tau0=0.50_pr*(daRhoTau+dbRhoTau)*mevfm
       w=c13*rhom2*((2.0_pr*tauc*hbzero)/rhom+3.0_pr*aRho0Rho0+5.0_pr*tauc*aRho0Tau0 &
         +6.0_pr*rhom*hRho0Rho0+8.0_pr*tauc*rhom*hRho0Tau0+3.0_pr*Crho(0) &
         +3.0_pr*(1.0_pr+sigma)*rhom**sigma*Cdrho(0)+5.0_pr*tauc*Ctau(0)+u*daRho0Rho0 &
         +u*tauc*daRho0Tau0+u*rhom*dhRho0Rho0+u*tauc*rhom*dhRho0Tau0)
       step=-w*(rhom-rhom0)/(w-w0)
       rhom0=rhom; w0=w; rhom=rhom+step
       If(iter.Gt.100) Stop 'STOP(In find_NM_RHOC)'
       !energy=tauc*hbzero+rhom*(aRho0Rho0+rhom*hRho0Rho0+Crho(0)+rhom**sigma*Cdrho(0)) &
       ! +tauc*rhom*(aRho0Tau0+rhom*hRho0Tau0+Ctau(0))
       !Write(6,'(a,15(1pg12.4))') ' rhom0,rhom,step,e,w=',rhom0,rhom,step,energy,w
    End Do
    find_NM_RHOC=rhom
  End Function find_NM_RHOC
  !=======================================================================
  !
  !=======================================================================
  Subroutine C_from_t()
    !--------------------------------------------------------------------------------
    ! C- from (t,x)-
    !--------------------------------------------------------------------------------
    Implicit None
    Crho(0)  =   3.0_pr/8.0_pr  * t0
    Cdrho(0) =  (1.0_pr/16.0_pr)* t3
    Crho(1)  = -(1.0_pr/4.0_pr) * t0*(0.50_pr+x0)
    Cdrho(1) = -(1.0_pr/24.0_pr)* t3*(0.50_pr+x3)
    Ctau(0)  =  (3.0_pr/16.0_pr)* t1+(1.0_pr/4.0_pr)*t2*(5.0_pr/4.0_pr+x2)
    Ctau(1)  = -(1.0_pr/8.0_pr) * t1*(0.5+x1)+(1.0_pr/8.0_pr)*t2*(0.50_pr+x2)
    CrDr(0)  =  (1.0_pr/16.0_pr)* t2*(5.0_pr/4.0_pr+x2)-(9.0_pr/64.0_pr)*t1
    CrDr(1)  =  (3.0_pr/32.0_pr)* t1*(0.5+x1)+(1.0_pr/32.0_pr)*t2*(0.50_pr+x2)
    CJ(0)    = -(1.0_pr/16.0_pr)*(t1*(2.0_pr*x1-1.0_pr)+t2*(2.0_pr*x2+1)-5*te-15*to)
    CJ(1)    = -(1.0_pr/16.0_pr)*(t2 -t1 + 5.0_pr*te -5.0_pr*to )
    CrdJ(0)  = -b4-(0.50_pr)*b4p
    CrdJ(1)  = -0.50_pr*b4p
  End Subroutine C_from_t
  !=======================================================================
  !
  !=======================================================================
  Subroutine t_from_C()
    !--------------------------------------------------------------------------------
    ! (t,x)- from C-
    !--------------------------------------------------------------------------------
    Implicit None
    t0  =  (8.0_pr/3)*Crho(0)
    t1  =  4.0_pr/3.0_pr*(Ctau(0)-4.0_pr*CrDr(0))
    t2  =  4.0_pr/3.0_pr*(3.0_pr*Ctau(0)-6.0_pr*Ctau(1)+4.0_pr*CrDr(0)-8.0_pr*CrDr(1))
    t3  =  16.0_pr*Cdrho(0)
    x0  = -0.50_pr*(3.0_pr*Crho(1)/Crho(0)+1.0_pr)
    x1  =  2.0_pr*(-Ctau(0)-3.0_pr*Ctau(1)+4.0_pr*CrDr(0)+12.0_pr*CrDr(1))/t1/3.0_pr
    x2  = -2.0_pr*(3.0_pr*Ctau(0)-15.0_pr*Ctau(1)+4.0_pr*CrDr(0)-20.0_pr*CrDr(1))/t2/3.0_pr
    x3  = -0.50_pr*(3.0_pr*Cdrho(1)/Cdrho(0)+1.0_pr)
    b4  =  CrdJ(1)-CrdJ(0)
    b4p = -2.0_pr*CrdJ(1)
    te  = (4.0_pr/15.0_pr)*(3.0_pr*CJ(0)-9.0_pr*CJ(1)-4.0_pr*CrDr(0)+12.0_pr*CrDr(1)-2.0_pr*Ctau(0)+6.0_pr*Ctau(1))
    to  = (4.0_pr/15.0_pr)*(3.0_pr*CJ(0)+3.0_pr*CJ(1)+4.0_pr*CrDr(0)+4.0_pr*CrDr(1))
  End Subroutine t_from_C
  !=======================================================================
  !
  !=======================================================================
  Subroutine CHrho_from_NM()
    !--------------------------------------------------------------------------------
    ! CHrho from NM, E_NM(Hartree)=CHrho*RHO_NM
    !--------------------------------------------------------------------------------
    Implicit None
    Real(pr) :: z3=1.50_pr
    !
    !!CHrho= &
    !!+h0mpi6c3NM*(A3_1/b3_1**z3+A3_2/b3_2**z3+A3_3/b3_3**z3+A3_4/b3_4**z3+A3_5/b3_5**z3) &
    !!+h0mpi6c1NM*(A1_1/b1_1**z3+A1_2/b1_2**z3+A1_3/b1_3**z3+A1_4/b1_4**z3+A1_5/b1_5**z3)
    CHrho = 0.0_pr
    !
  End Subroutine CHrho_from_NM
  !=======================================================================
  !
  !=======================================================================
  Elemental Function HartreeV00(u)
    !--------------------------------------------------------------------------------
    ! HartreeV(u), E(Hartree)=(1/2)*Int[rho_0(r)*V(|r-r'|)*rho_0(r')]
    !--------------------------------------------------------------------------------
    Implicit None
    Real(pr), Intent(in) :: u
    Real(pr)             :: x2,HartreeV00
    !
    !!x2=(u*mpi)**2
    !
    !!HartreeV=h0mpi6c1*(Exp(-x2*b1_1)*A1_1+Exp(-x2*b1_2)*A1_2+Exp(-x2*b1_3)*A1_3+Exp(-x2*b1_4)*A1_4+Exp(-x2*b1_5)*A1_5)+&
    !!h0mpi6c3*(Exp(-x2*b3_1)*A3_1+Exp(-x2*b3_2)*A3_2+Exp(-x2*b3_3)*A3_3+Exp(-x2*b3_4)*A3_4+Exp(-x2*b3_5)*A3_5)
    !
    HartreeV00=0.0_pr
    !
  End Function HartreeV00
  !=======================================================================
  !
  !=======================================================================
  Elemental Function HartreeV01(u)
    !--------------------------------------------------------------------------------
    ! HartreeV(u), E(Hartree)=(1/2)*Int[rho_0(r)*V(|r-r'|)*rho_1(r')]
    !--------------------------------------------------------------------------------
    Implicit None
    Real(pr), Intent(in) :: u
    Real(pr)             :: x2,HartreeV01
    !
    HartreeV01=0.0_pr
    !
  End Function HartreeV01
  !=======================================================================
  !
  !=======================================================================
  Elemental Function HartreeV11(u)
    !--------------------------------------------------------------------------------
    ! HartreeV(u), E(Hartree)=(1/2)*Int[rho_1(r)*V(|r-r'|)*rho_1(r')]
    !--------------------------------------------------------------------------------
    Implicit None
    Real(pr), Intent(in) :: u
    Real(pr)             :: x2,HartreeV11
    !
    HartreeV11=0.0_pr
    !
  End Function HartreeV11
  !=======================================================================
  !
  !=======================================================================
  Elemental Function ThetaFunction2(u)
    !--------------------------------------------------------------------------------
    ! ThetaFunction2(u)=0 or 1  when x2<2  or x2>2
    !--------------------------------------------------------------------------------
    Implicit None
    Real(pr), Intent(in) :: u
    Real(pr)             :: x2,ThetaFunction2
    !
    x2=(u*mpi)
    !
    ThetaFunction2=0.0_pr
    If(x2.Gt.2.0_pr) ThetaFunction2=1.0_pr
    !
  End Function ThetaFunction2
  !=======================================================================
  !
  !=======================================================================
  Subroutine Make_Parameter_Free_Useful_Combinations()
    !--------------------------------------------------------------------------------
    ! Make Useful combinations
    !--------------------------------------------------------------------------------
    Implicit None
    !
    If(dmeorder.Ge.0) Then
       !
       mpi2=mpi**2
       gA2=gA**2; gA4=gA2**2; gA6=gA2**3;
       fpi2=fpi**2; fpi4=fpi2**2;
       CHartree =mevfm*(3.0_pr*gA2)/(32.0_pr*fpi4*Pi**2)
       h0mpi6=197.30_pr*(mpi**6)*(3.0_pr*gA*gA)/(32.0_pr*fpi**4*Pi**2)
       h0mpi6c1=h0mpi6*c1;         h0mpi6c3=h0mpi6*c3
       !
       h0mpi6NM=197.30_pr*(3.0_pr*(mpi**3)*gA2)/(64.0_pr*fpi**4*Sqrt(Pi))
       h0mpi6c1NM=h0mpi6NM*c1;     h0mpi6c3NM=h0mpi6NM*c3
       !
       A3_1=42.7132145164590_pr;   A3_2=0.670441422115440_pr; A3_3=0.0525713896514650_pr;
       A3_4=0.0012545731701320_pr; A3_5=5.81008627207380_pr
       b3_1=3.0809379008590_pr;    b3_2=0.905186811964580_pr; b3_3=0.474514509597610_pr;
       b3_4=0.228138177966090_pr;  b3_5=1.66931540698090_pr;
       !
       A1_1=2.5000830618386_pr;    A1_2=0.619542286897850_pr; A1_3=0.169682589033730_pr;
       A1_4=0.0276112113725470_pr; A1_5=0.00108164458809540_pr
       b1_1=1.75854210706510_pr;   b1_2=0.88882524524657_pr;  b1_3=0.46377235143756_pr;
       b1_4=0.247665887704790_pr;  b1_5=0.132222413002680_pr
       !
    End If
    !
  End Subroutine Make_Parameter_Free_Useful_Combinations
  !=======================================================================
  !
  !=======================================================================
  Elemental Function Vexternal(t,x,y,z)
    !
    Implicit None
    Integer(ipr), Intent(in) :: t  !! isospin index: 0=isoscalar, 1=isovector
    Real(pr), Intent(in) :: x,y,z  !! position in cartesian basis
    Real(pr) :: Vexternal
    !
    Vexternal = 0.0_pr
    !
  End Function Vexternal
  !=======================================================================
  !
  !=======================================================================
End Module UNEDF

