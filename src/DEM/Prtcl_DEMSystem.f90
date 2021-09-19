module Prtcl_DEMSystem
  use MPI
  use Prtcl_decomp_2d
#ifdef CFDDEM
  use decomp_2d,only: nrank
  use Prtcl_DumpPrtcl
#endif
  use Prtcl_TypeDef
  use Prtcl_Timer
  use Prtcl_LogInfo
  use Prtcl_Parameters
  use Prtcl_Property
  use Prtcl_Geometry
  use Prtcl_Variables
  use Prtcl_Comm
  use Prtcl_CL_and_CF
  use Prtcl_ContactSearch
  use Prtcl_Integration
  use Prtcl_ContactSearchPW
  use Prtcl_IOAndVisu
  implicit none
  private
    
  !// DEMSystem class 
  type DEMSystem
    integer :: iterNumber   = 0  ! iteration number 
        
    !// timers
    type(timer):: m_total_timer
    type(timer):: m_pre_iter_timer
    type(timer):: m_comm_cs_timer
    type(timer):: m_CSCF_PP_timer
    type(timer):: m_CSCF_PW_timer
    type(timer):: m_Acceleration_timer
    type(timer):: m_integration_timer
    type(timer):: m_write_prtcl_timer
    type(timer):: m_comm_exchange_timer
  contains
    procedure:: Initialize => DEMS_Initialize
    
    ! iterating simulation for n time steps
    procedure:: iterate     => DEMS_iterate
    
    ! performing pre-iterations 
    procedure:: preIteration    => DEMS_preIteration
    
    !calculating acceleration of particles - linear and angular
    procedure:: clc_Acceleration    => DEMS_clc_Acceleration
  end type DEMSystem
  type(DEMSystem),public::DEM
    
contains

!********************************************************************************
!   Initializing DEMSystem object with particles which are inserted from a 
!   predefined plane 
!********************************************************************************
  subroutine  DEMS_Initialize(this,chDEMPrm)
    implicit none
    class(DEMSystem)::this
    character(*),intent(in)::chDEMPrm
    integer:: ierror
    character(256):: ch
    
    !// Initializing main log info and visu
    this%IterNumber=DEM_Opt%ifirst-1
    ch= 'mkdir -p '//DEM_opt%Res_Dir//' '//DEM_opt%RestartDir//' 2> /dev/null'
    if(nrank==0) call system(trim(adjustl(ch)))
    call DEMLogInfo%InitLog(DEM_opt%Res_Dir, DEM_opt%RunName,DEM_opt%LF_file_lvl, DEM_opt%LF_cmdw_lvl)
    if(nrank==0) call Write_DEM_Opt_to_Log()

    ! Step1: Physical property
    call DEMProperty%InitPrtclProperty(chDEMPrm)
    call DEMProperty%InitWallProperty(chDEMPrm)
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step1: Physical properties of particels and walls are set.",1)
      call DEMLogInfo%OutInfo("Physical properties contains "// trim( num2str(DEM_opt%numPrtcl_Type ) ) // &
                              " particle types and "//trim( num2str(DEM_opt%numWall_type ) )// " wall types.",2)
    endif

    ! Step2: set the geometry
    call DEMGeometry%MakeGeometry(chDEMPrm)
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step2: Geometry is set", 1 )
      call DEMLogInfo%OutInfo("Geometry Contains "//trim(num2str(DEMGeometry%num_pWall))//" Plane walls.", 2)
    endif

    ! Step3: initilize all the particle variables
    call GPrtcl_list%AllocateAllVar()
#ifdef CFDDEM
    if(.not.DEM_Opt%RestartFlag) then
      if(DEM_Opt%numPrtcl>0) call DEM_IO%ReadInitialCoord()
      if(nrank==0) then
        call DEMLogInfo%OutInfo("Step3: Initial Particle coordinates are READING into DEMSystem ...", 1 )
        call DEMLogInfo%OutInfo("Number of particles avaiable in the system:"//trim(num2str(DEM_opt%numPrtcl)),2)
      endif
      DEM_opt%np_InDomain = DEM_opt%numPrtcl
      if(DEM_Opt%numPrtclFix>0) call DEM_IO%ReadFixedCoord()
    else
      if(DEM_Opt%numPrtcl>0) then
        call DEM_IO%Read_Restart()
      else
        DEM_opt%np_InDomain=0
      endif
      if(nrank==0) then
        call DEMLogInfo%OutInfo("Step3: Particles are READING from the Resarting file ...", 1 )
        call DEMLogInfo%OutInfo("Number of particles avaiable in domain:"//trim(num2str(DEM_opt%np_InDomain)),2)
      endif
      if(DEM_Opt%numPrtclFix>0) call DEM_IO%ReadFixedRestart()
    endif
#else
    if(.not.DEM_Opt%RestartFlag) then
      call GPrtcl_list%MakingAllPrtcl(chDEMPrm)
      if(nrank==0) then
        call DEMLogInfo%OutInfo("Step3: Particles are MAKING into DEMSystem ...", 1 )
        call DEMLogInfo%OutInfo("Number of particles avaiable in the system:"//trim(num2str(DEM_opt%numPrtcl)),2)
      endif
      DEM_opt%np_InDomain = DEM_opt%numPrtcl
    else
      call DEM_IO%Read_Restart()
      if(nrank==0) then
        call DEMLogInfo%OutInfo("Step3: Particles are READING from the Resarting file ...", 1 )
        call DEMLogInfo%OutInfo("Number of particles avaiable in domain:"//trim(num2str(DEM_opt%np_InDomain)),2)
      endif
    endif
    if(DEM_Opt%numPrtclFix>0) call DEM_IO%ReadFixedCoord()
#endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    ! Step4: Initializing visu
    call DEM_IO%Init_visu(chDEMPrm)
#ifdef CFDDEM
    call DEM_PDump%Initilize(chDEMPrm)
#endif

    ! Step5: initialize the inter-processors communication
    call DEM_Comm%InitComm()
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step4: Initializing the inter-processors communication . . . ", 1 )
    endif

    ! Step6: Initializing contact list and contact force 
    call  GPPW_CntctList%InitContactList()
    if(DEM_Opt%RestartFlag) then
      if(DEM_opt%np_InDomain>0)call DEM_IO%RestartCL()
    endif
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step5: Initializing contact list and contact force models . . . ", 1 )
      if( DEM_opt%CF_Type == CFT_LSD ) then
        ch = "linear spring-dashpot with limited tangential displacement"
      elseif( DEM_opt%CF_Type == CFT_nLin ) then
        ch = "non-linear visco-elastic model with limited tangential displacement"
      endif
      call DEMLogInfo%OutInfo("Contact force model is "//trim(ch), 2 )

      if(DEM_opt%PI_Method==PIM_FE) then
        ch = "Forward Euler             "
      elseif(DEM_opt%PI_Method==PIM_AB2) then
        ch = "Adams Bashforth: 2nd Order"
      elseif(DEM_opt%PI_Method==PIM_AB3) then
        ch = "Adams Bashforth: 3nd Order"
      endif
      call DEMLogInfo%OutInfo("Linear   movement Integration scheme is : "//trim(ch),2)

      if(DEM_opt%PRI_Method==PIM_FE) then
        ch = "Forward Euler             "
      elseif(DEM_opt%PRI_Method==PIM_AB2) then
        ch = "Adams Bashforth: 2nd Order"
      elseif(DEM_opt%PRI_Method==PIM_AB3) then
        ch = "Adams Bashforth: 3nd Order"
      endif
      call DEMLogInfo%OutInfo("Rotating movement Integration scheme is : "//trim(ch),2)
    endif

    ! Step7: Initializing contact search method
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step6: Initializing contact search method . . . ", 1)
      call DEMLogInfo%OutInfo("Particle-Particle contact search intialization...",2)
      call DEMLogInfo%OutInfo("Particle-Wall contact search intialization...",2)
    endif
    call DEMContactSearch%InitContactSearch()
    call DEMContactSearchPW%InitContactSearchPW()
    
    ! Step8: timers for recording the execution time of different parts of program
    if(nrank==0) call DEMLogInfo%OutInfo("Step7: Initializing timers . . . ", 1 )
    call this%m_total_timer%reset()
    call this%m_pre_iter_timer%reset()
    call this%m_comm_cs_timer%reset()
    call this%m_CSCF_PP_timer%reset()
    call this%m_CSCF_PW_timer%reset()
    call this%m_Acceleration_timer%reset()
    call this%m_integration_timer%reset()
    call this%m_comm_exchange_timer%reset()
    
  end subroutine DEMS_Initialize

  !********************************************************************************
  !   iterating over time 
  !   calls all the required methods to do numIter iterations in the DEM system
  !********************************************************************************
  subroutine DEMS_iterate(this,itime)
    implicit none
    class(DEMSystem) this
    integer,intent(in)::itime

    ! locals 
    character(256):: chLine
    integer:: Consv_Cont(2),Consv_Cont1(2),ierror,npwcs(4)

    ! body
    call this%m_total_timer%start()

    ! pre-iteration adjustments 
    call this%m_pre_iter_timer%start()
    call this%preIteration()
    call this%m_pre_iter_timer%finish()

    ! inter-processor commucation for contact search ( ghost particle )
    call this%m_comm_cs_timer%start()
    call DEM_Comm%Comm_For_Cntct()
    call this%m_comm_cs_timer%finish()

    ! finding contacts among particels, and then calculating contact forces
    call this%m_CSCF_PP_timer%start()
    call DEMContactSearch%FindContacts()
    call this%m_CSCF_PP_timer%finish()

    ! finding contacts between particles and walls, and then calculating contact forces
    call this%m_CSCF_PW_timer%start()
    call DEMContactSearchPW%FindContactsPW()
    call this%m_CSCF_PW_timer%finish() 

    ! calculating linear and angular accelerations
    call this%m_Acceleration_timer%start()
    call this%clc_Acceleration()
    call this%m_Acceleration_timer%finish()

    ! correcting position and velocities 
    call this%m_integration_timer%start()
    call Prtcl_Integrate()
    call this%m_integration_timer%finish()
   
    ! inter-processor commucation for exchange
    call this%m_comm_exchange_timer%start()
    call DEM_Comm%Comm_For_Exchange()
    call GPPW_CntctList%RemvReleased()
    call this%m_comm_exchange_timer%finish()
    this%iterNumber = this%iterNumber + 1

    ! writing results to the output file and Restart file
    call this%m_write_prtcl_timer%start()
    call MPI_ALLREDUCE(GPrtcl_list%nlocal, DEM_Opt%np_InDomain, 1, int_type,MPI_SUM,MPI_COMM_WORLD,ierror)
#ifndef CFDDEM
    if( mod(this%IterNumber,DEM_opt%SaveVisu)== 0)   call DEM_IO%dump_visu(itime)
#else
    if( mod(this%IterNumber,DEM_opt%SaveVisu)== 0)   call DEM_IO%dump_visu(itime/icouple)
    if( mod(this%IterNumber,DEM_PDump%WriteCacheFreq)== 0)  call DEM_PDump%WriteCache(itime)
#endif

    if( mod(this%IterNumber,DEM_opt%BackupFreq)== 0 .or. itime==DEM_opt%ilast) then
      call DEM_IO%Write_Restart(itime)
#ifdef CFDDEM
      call DEM_IO%WriteFixedRestart(itime)
      call DEM_PDump%PrtclVarDump(itime)
#endif
      call DEM_IO%Delete_Prev_Restart(itime)
    endif
    call this%m_write_prtcl_timer%finish()
    call this%m_total_timer%finish()
        
    ! output to log file and terminal/command window
    IF((this%IterNumber==DEM_Opt%ifirst .or. mod(this%IterNumber,DEM_opt%Cmd_LFile_Freq)==0) ) THEN
      Consv_Cont1 =  DEMContactSearch%get_numContact()
      call MPI_REDUCE(Consv_Cont1, Consv_Cont,       2,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
      call MPI_REDUCE(GPPW_CntctList%numCntcts,npwcs,4,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
      if(nrank/=0) return
    
      ! command window and log file output
      call DEMLogInfo%OutInfo("DEM performed "//trim(num2str(this%IterNumber))//" iterations up to here!",1)
      call DEMLogInfo%OutInfo("Execution time [tot, last, ave] [sec]: "//trim(num2str(this%m_total_timer%tot_time))//", "// &
      trim(num2str(this%m_total_timer%last_time ))//", "//trim(num2str(this%m_total_timer%average())),2)

      call DEMLogInfo%OutInfo("PreItertion time [tot, ave]        : "//trim(num2str(this%m_pre_iter_timer%tot_time))//", "// &
      trim(num2str(this%m_pre_iter_timer%average())),3)

      call DEMLogInfo%OutInfo("Comm_For_Contact [tot, ave]        : "//trim(num2str(this%m_comm_cs_timer%tot_time))//", "// &
      trim(num2str(this%m_comm_cs_timer%average())), 3)

      call DEMLogInfo%OutInfo("CS and CF P-P time [tot, ave]      : "//trim(num2str(this%m_CSCF_PP_timer%tot_time))//", "// &
      trim(num2str(this%m_CSCF_PP_timer%average())), 3)

      call DEMLogInfo%OutInfo("CS and CF P-W time [tot, ave]      : "//trim(num2str(this%m_CSCF_PW_timer%tot_time))//", "// &
      trim(num2str(this%m_CSCF_PW_timer%average())), 3)

      call DEMLogInfo%OutInfo("Acceleration time [tot, ave]       : "//trim(num2str(this%m_Acceleration_timer%tot_time))//", "// &
      trim(num2str(this%m_Acceleration_timer%average())), 3)

      call DEMLogInfo%OutInfo("Integration time [tot, ave]        : "//trim(num2str(this%m_integration_timer%tot_time))//", "// &
      trim(num2str(this%m_integration_timer%average())), 3)

      call DEMLogInfo%OutInfo("Comm_For_Exchange [tot, ave]       : "//trim(num2str(this%m_comm_exchange_timer%tot_time))//", "// &
      trim(num2str(this%m_comm_exchange_timer%average())), 3)

      call DEMLogInfo%OutInfo("Write to file time [tot, ave]      : "//trim(num2str(this%m_write_prtcl_timer%tot_time))//", "// &
      trim(num2str(this%m_write_prtcl_timer%average())), 3)
     
      chLine = "Particle number in  domain:  "//trim(num2str(DEM_Opt%np_InDomain))
      call DEMLogInfo%OutInfo(chLine, 2)        

      call DEMLogInfo%OutInfo("Contact information", 2)
      chLine = "No. consrvtv. contacts, same level | cross level: "// trim(num2str(Consv_Cont(1)))//" | "//trim(num2str(Consv_Cont(2)))
      call DEMLogInfo%OutInfo(chLine, 3)
      chLine = "No. exact contacts P-P | P-GP | P-FP | P-W : "//trim(num2str(npwcs(1)))//" | "//trim(num2str(npwcs(2)))//" | "//trim(num2str(npwcs(3)))//" | "//trim(num2str(npwcs(4)))
      call DEMLogInfo%OutInfo(chLine, 3)
    ENDIF

  end subroutine DEMS_iterate
    
  !**********************************************************************
  ! DEMS_preIteration
  !**********************************************************************
  subroutine DEMS_preIteration(this)
    implicit none
    class(DEMSystem)::this
        
    ! update wall neighbor list if necessary
    call DEMContactSearchPW%UpdateNearPrtclsPW(this%iterNumber)
    GPrtcl_cntctForce =zero_r3
    GPrtcl_torque= zero_r3
    call GPPW_CntctList%PreIteration()
  end subroutine DEMS_preIteration

  !**********************************************************************
  ! calculating acceleration of particles (linear and angular) 
  !**********************************************************************
#ifdef CFDDEM
  subroutine DEMS_clc_Acceleration(this)
    implicit none
    class(DEMSystem) this

    ! locals
    integer:: i,itype
    type(real3):: Fpforce
    real(RK)::Mass, Inertia,MassInFluid,MassTot,PrtclDensity

    do i=1,GPrtcl_list%nlocal
      itype   = GPrtcl_pType(i)
      Mass    = DEMProperty%Prtcl_PureProp(itype)%Mass
      MassTot    = DEMProperty%Prtcl_PureProp(itype)%MassOfFluid*half +Mass
      MassInFluid= DEMProperty%Prtcl_PureProp(itype)%MassInFluid
      Inertia    = DEMProperty%Prtcl_PureProp(itype)%Inertia

      Fpforce = 1.50_RK*GPrtcl_FpForce(i)-0.50_RK*GPrtcl_FpForce_old(i)
      GPrtcl_linAcc(1,i) = ((GPrtcl_cntctForce(i)+Fpforce)/MassTot)+ MassInFluid/MassTot*DEM_opt%gravity
      GPrtcl_rotAcc(1,i) = ((GPrtcl_torque(i)/Inertia))
    enddo
  end subroutine DEMS_clc_Acceleration

#else
  subroutine DEMS_clc_Acceleration(this)
    implicit none
    class(DEMSystem) this

    ! locals
    integer:: i,itype
    real(RK)::Mass, Inertia
    
    do i=1,GPrtcl_list%nlocal
      itype   = GPrtcl_pType(i)
      Mass    = DEMProperty%Prtcl_PureProp(itype)%Mass
      Inertia = DEMProperty%Prtcl_PureProp(itype)%Inertia
      GPrtcl_linAcc(1,i) = ((GPrtcl_cntctForce(i))/Mass)+DEM_opt%gravity
      GPrtcl_rotAcc(1,i) = ((GPrtcl_torque(i)/Inertia))
    enddo
  end subroutine DEMS_clc_Acceleration
#endif

  !**********************************************************************
  ! Write_DEM_Opt_to_Log
  !**********************************************************************
  subroutine Write_DEM_Opt_to_Log()
    implicit none

    ! locals
    character(64):: RunName, Res_Dir,RestartDir
    integer::numPrtcl, numPrtclFix, CS_Method, CF_Type, CT_Model, PI_Method, PRI_Method
    integer::numPrtcl_Type, numWall_type, CntctList_Size, CS_numlvls, Base_wall_id
    integer::Wall_max_update_iter, SaveVisuDEM,BackupFreqDEM, Cmd_LFile_Freq, LF_file_lvl, LF_cmdw_lvl
    real(RK)::dtDEM,Wall_neighbor_ratio,Prtcl_cs_ratio
    type(real3):: gravity, minpoint, maxpoint
    logical,dimension(3)::IsPeriodic
    integer:: GeometrySource,ifirstDEM,ilastDEM
    character(64):: Geom_Dir 
    logical::RestartFlag
    NAMELIST /DEMOptions/ RestartFlag,numPrtcl, numPrtclFix, dtDEM, gravity, minpoint,maxpoint, CS_Method,   &
                          CF_Type, CT_Model, PI_Method, PRI_Method, numPrtcl_Type, numWall_type, CS_numlvls, &
                          CntctList_Size, Base_wall_id, Wall_max_update_iter, Wall_neighbor_ratio,RunName,   &
                          Res_Dir,RestartDir,BackupFreqDEM, SaveVisuDEM, Cmd_LFile_Freq, LF_file_lvl,        &
                          LF_cmdw_lvl, GeometrySource, Geom_Dir,ifirstDEM,ilastDEM,Prtcl_cs_ratio,IsPeriodic
#ifdef CFDDEM
    NAMELIST/CFDDEMCoupling/icouple,UpdateDEMflag,is_clc_Lift,is_clc_Basset,is_clc_Basset_fixed,is_clc_ViscousForce,&
                            is_clc_PressureGradient,is_clc_ViscousForce,is_clc_FluidAcc,FluidAccCoe,SaffmanConst,RatioSR
    NAMELIST/BassetOptions/ mWinBasset, mTailBasset, BassetAccuracy, BassetTailType 
#endif

    RestartFlag = DEM_opt%RestartFlag 
    numPrtcl    = DEM_opt%numPrtcl   
    numPrtclFix = DEM_opt%numPrtclFix  
    dtDEM       = DEM_opt%dt       
    ifirstDEM   = DEM_opt%ifirst   
    ilastDEM    = DEM_opt%ilast    
    gravity     = DEM_opt%gravity  
    minpoint    = DEM_opt%SimDomain_min 
    maxpoint    = DEM_opt%SimDomain_max 
    IsPeriodic  = DEM_opt%IsPeriodic 
           
    Prtcl_cs_ratio=  DEM_opt%Prtcl_cs_ratio 
    CS_Method  =  DEM_opt%CS_Method 
    CF_Type    = DEM_opt%CF_Type   
    CT_Model   = DEM_opt%CT_Model  
    PI_Method  = DEM_opt%PI_Method  
    PRI_Method = DEM_opt%PRI_Method
           
    numPrtcl_Type = DEM_opt%numPrtcl_Type
    numWall_type  = DEM_opt%numWall_type 
          
    CntctList_Size =  DEM_opt%CntctList_Size 
    CS_numlvls     = DEM_opt%CS_numlvls    
           
    Base_wall_id  =  DEM_opt%Base_wall_id 
    Wall_max_update_iter = DEM_opt%Wall_max_update_iter
    Wall_neighbor_ratio  = DEM_opt%Wall_neighbor_ratio 
           
    RunName = DEM_opt%RunName 
    Res_Dir = DEM_opt%Res_Dir 
    RestartDir    = DEM_opt%RestartDir 
    BackupFreqDEM = DEM_opt%BackupFreq 
    SaveVisuDEM   = DEM_opt%SaveVisu 
    Cmd_LFile_Freq= DEM_opt%Cmd_LFile_Freq 
    LF_file_lvl   = DEM_opt%LF_file_lvl 
    LF_cmdw_lvl   = DEM_opt%LF_cmdw_lvl 
           
    GeometrySource=  DEM_opt%GeometrySource 
    Geom_Dir  = DEM_opt%Geom_Dir 
    write(DEMLogInfo%nUnit, nml=DEMOptions)
#ifdef CFDDEM
    write(DEMLogInfo%nUnit, nml=CFDDEMCoupling)
    write(DEMLogInfo%nUnit, nml=BassetOptions)
#endif
   
  end subroutine Write_DEM_Opt_to_Log

end module Prtcl_DEMSystem
