module Prtcl_DumpPrtcl
  use MPI
  use Prtcl_TypeDef
  use Prtcl_LogInfo
  use Prtcl_Variables
  use Prtcl_Parameters
  use Prtcl_CL_and_CF
  use decomp_2d,only: nrank
  implicit none
  private
  
  integer,parameter:: Prtcl_Dump_Flag = 93
  integer,parameter:: RKP=4
  character(64)::DumpPrtclDir
  logical,save:: DumpPrtclFlag
  logical,save:: ResetDumpFlag
  integer,save:: nDumpPrtclSize,mDumpPrtclSize,WriteCacheFreq

  type DumpPrtclVar
    integer:: itime
    integer:: id
    integer:: pType
    integer:: CntctFlag
    real(RKP):: Pos(3)
    real(RKP):: linVel(3)
    real(RKP):: RotVel(3)
  end type DumpPrtclVar
  type(DumpPrtclVar),dimension(:),allocatable:: DumpPrtclVec
  
  type Prtcl_Dump
    integer:: WriteCacheFreq=100000
  contains
    procedure:: Initilize
    procedure:: WriteCache   
    procedure:: PrtclVarDump   
  end type Prtcl_Dump
  type(Prtcl_Dump),public:: DEM_PDump

contains

  !**********************************************************************
  ! Initilize
  !**********************************************************************
  subroutine Initilize(this,chFile)
    implicit none
    class(Prtcl_Dump)::this
    character(*),intent(in)::chFile

    ! locals
    character(256):: ch
    real(RK)::yDump=zero
    integer:: pid,nUnitFile,ierror,disp0
    type(DumpPrtclVar):: DumpPrtclVar_t
    namelist/DumpPrtclOptions/DumpPrtclFlag,ResetDumpFlag,yDump,mDumpPrtclSize,DumpPrtclDir,WriteCacheFreq

    nUnitFile = GetFileUnit()
    open(unit=nUnitFile, file=chFile, status='old', form='formatted', IOSTAT=ierror)
    if(ierror /= 0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"InitDumpPrtclVar","Cannot open file:"//trim(adjustl(chFile)))
    read(nUnitFile, nml=DumpPrtclOptions)
    close(nUnitFile,IOSTAT=ierror)

    if(.not. DumpPrtclFlag) return
    if(ResetDumpFlag) then
      do pid=1,GPrtcl_list%nlocal
        if(GPrtcl_PosR(pid)%y >= yDump) then
          GPrtcl_usrMark(pid)=Prtcl_Dump_Flag
        else
          GPrtcl_usrMark(pid)=1
        endif
      enddo
    endif
    DEM_PDump%WriteCacheFreq =  WriteCacheFreq

    nDumpPrtclSize=0
    if(mDumpPrtclSize<10000 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"InitDumpPrtclVar","So small mDumpPrtclSize:"//trim(num2str(mDumpPrtclSize)) )
    allocate(DumpPrtclVec(mDumpPrtclSize),Stat=ierror )
    if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"InitDumpPrtclVar","Allocation failed")

    write(ch,"(A)") 'mkdir -p '//trim(adjustl(DumpPrtclDir))//' 2> /dev/null'
    if(nrank==0) then 
      call system(trim(adjustl(ch)))
      write(DEMLogInfo%nUnit, nml=DumpPrtclOptions)
    endif

  end subroutine Initilize

  !**********************************************************************
  ! WriteCache
  !**********************************************************************
  subroutine WriteCache(this,itime)
    implicit none
    class(Prtcl_Dump)::this
    integer,intent(in)::itime

    ! locals
    integer::pid,nlocal

    if(.not. DumpPrtclFlag) return
    nlocal=GPrtcl_list%nlocal
    if(nlocal>mDumpPrtclSize) call DEMLogInfo%CheckForError(ErrT_Abort,"WriteCache","So small mDumpPrtclSize !!!")
    do pid=1,nlocal
      if(GPrtcl_usrMark(pid)/=Prtcl_Dump_Flag) cycle
      nDumpPrtclSize=nDumpPrtclSize+1

      DumpPrtclVec(nDumpPrtclSize)%itime     = itime
      DumpPrtclVec(nDumpPrtclSize)%id        = GPrtcl_id(pid)
      DumpPrtclVec(nDumpPrtclSize)%pType     = GPrtcl_pType(pid)
      DumpPrtclVec(nDumpPrtclSize)%CntctFlag = GPPW_CntctList%IsCntct(pid)
      DumpPrtclVec(nDumpPrtclSize)%Pos(1)    = real(GPrtcl_PosR(pid)%x, RKP)           ! real01
      DumpPrtclVec(nDumpPrtclSize)%Pos(2)    = real(GPrtcl_PosR(pid)%y, RKP)           ! real02
      DumpPrtclVec(nDumpPrtclSize)%Pos(3)    = real(GPrtcl_PosR(pid)%z, RKP)           ! real03
      DumpPrtclVec(nDumpPrtclSize)%linVel(1) = real(GPrtcl_linVel(1,pid)%x, RKP)       ! real04
      DumpPrtclVec(nDumpPrtclSize)%linVel(2) = real(GPrtcl_linVel(1,pid)%y, RKP)       ! real05
      DumpPrtclVec(nDumpPrtclSize)%linVel(3) = real(GPrtcl_linVel(1,pid)%z, RKP)       ! real06
      DumpPrtclVec(nDumpPrtclSize)%RotVel(1)     = real(GPrtcl_RotVel(1,pid)%x, RKP)   ! real07
      DumpPrtclVec(nDumpPrtclSize)%RotVel(2)     = real(GPrtcl_RotVel(1,pid)%y, RKP)   ! real08
      DumpPrtclVec(nDumpPrtclSize)%RotVel(3)     = real(GPrtcl_RotVel(1,pid)%z, RKP)   ! real09

      if(nDumpPrtclSize==mDumpPrtclSize) call this%PrtclVarDump(itime)
    enddo

  end subroutine WriteCache

  !**********************************************************************
  ! PrtclVarDump
  !**********************************************************************
  subroutine PrtclVarDump(this,itime)
    implicit none
    class(Prtcl_Dump)::this
    integer,intent(in)::itime

    ! locals
    integer::ierror,nUnit
    character(256)::chFile

    if(.not.DumpPrtclFlag)return
    if(nDumpPrtclSize==0) return
    write(chFile,'(A,I5.5,A,I10.10)')trim(DumpPrtclDir)//'rank',nrank,'_',itime

    nUnit=GetFileUnit()
    open(unit=nUnit,file=trim(chFile),status='replace',form='unformatted',access='stream',IOSTAT=ierror)
    IF(ierror/=0) THEN
      call DEMLogInfo%CheckForError(ErrT_Abort,"PrtclVarDump","Cannot open file: "//trim(chFile))
    ELSE
      write(nUnit)DumpPrtclVec(1:nDumpPrtclSize)
    ENDIF
    close(nUnit,IOSTAT=ierror)

    nDumpPrtclSize=0
  end subroutine PrtclVarDump

end module Prtcl_DumpPrtcl
