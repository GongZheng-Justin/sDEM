!********************************************************************!
!*    file name  : Prtcl_Variables.f90                              *!
!*    module name: Prtcl_Variables                                  *!  
!*                                                                  *!
!*    purpose:                                                      *! 
!*      1) All datas required to represent spherical particles      *!
!*      2) Initialize all the particle variables                    *!
!*                                                                  *!
!*  Author: Zheng Gong           Date: 23:Feb:2020                  *!
!*                                                                  *!
!********************************************************************!

module Prtcl_Variables
  use MPI
  use Prtcl_TypeDef
  use Prtcl_decomp_2d
  use Prtcl_LogInfo
  use Prtcl_Parameters
  use Prtcl_Property
  implicit none
  private

  integer,dimension(:),allocatable,public:: GPrtcl_id
  integer,dimension(:),allocatable,public:: GPrtcl_pType
  integer,dimension(:),allocatable,public:: GPrtcl_usrMark
  type(real4),dimension(:),allocatable,public::   GPrtcl_PosR
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_linVel
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_linAcc
  type(real3),dimension(:),allocatable,public::   GPrtcl_theta
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_rotVel
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_rotAcc
  type(real3),dimension(:),allocatable,public::   GPrtcl_cntctForce
  type(real3),dimension(:),allocatable,public::   GPrtcl_torque

  integer,dimension(:),allocatable,public::       GPFix_id
  integer,dimension(:),allocatable,public::       GPFix_pType
  type(real4),dimension(:),allocatable,public::   GPFix_PosR

  type VarList
    integer:: nlocal     ! number of particles in local processor
    integer:: mlocal     ! the possible maxium # of particles in local processor
    integer:: mlocalFix  ! number of fixed particles in local processor

    integer:: nGhost_CS
    integer:: mGhost_CS
    integer:: nGhostFix_CS

    integer:: tsize  ! size for translational veloctity and acceleration
    integer:: rsize  ! size for rotational veloctity and acceleration
  contains
    procedure:: AllocateAllVar     => GL_AllocateAllVar
    procedure:: MakingAllPrtcl     => GL_MakingAllPrtcl
    procedure:: ReallocatePrtclVar => GL_ReallocatePrtclVar
    procedure:: copy               => GL_copy
  end type VarList
  type(VarList),public::GPrtcl_list 

  ! Ghost particle variables
  integer,dimension(:),allocatable,public:: GhostP_id
  integer,dimension(:),allocatable,public:: GhostP_pType
  type(real4),dimension(:),allocatable,public:: GhostP_PosR
  type(real3),dimension(:),allocatable,public:: GhostP_linVel
  type(real3),dimension(:),allocatable,public:: GhostP_rotVel
    
contains

  !**********************************************************************
  ! GL_AllocateAllVar
  !**********************************************************************
  subroutine GL_AllocateAllVar(this)
    implicit none
    class(VarList)::this

    ! locals
    type(real3)::SimLen
    integer:: mlocal,numPrtcl,mlocalFix,numPrtclFix
    integer:: iErr01,iErr02,iErr03,iErr04,iErr05,iErr06,iErr07
    integer:: iErr08,iErr09,iErr10,iErr11,iErr12,iErr13,iErrSum
    real(RK)::xst,xed,yst,yed,zst,zed,vol_tot,vol_local

    numPrtcl    = DEM_opt%numPrtcl
    numPrtclFix = DEM_opt%numPrtclFix

    ! step0: determine initial misze
    xst=DEM_decomp%xSt; xed=DEM_decomp%xEd
    yst=DEM_decomp%ySt; yed=DEM_decomp%yEd
    zst=DEM_decomp%zSt; zed=DEM_decomp%zEd

    ! step1: allocating memory for particles
    if(DEM_opt%PI_Method==PIM_FE) then
      this%tsize=1
    elseif (DEM_opt%PI_Method==PIM_AB2) then
      this%tsize=2
    elseif(DEM_opt%PI_Method==PIM_AB3) then
      this%tsize=3
    endif
    if(DEM_opt%PRI_Method==PIM_FE) then
      this%rsize=1
    elseif (DEM_opt%PRI_Method==PIM_AB2) then
      this%rsize=2
    elseif(DEM_opt%PRI_Method==PIM_AB3) then
      this%rsize=3
    endif 

    SimLen = DEM_Opt%SimDomain_max - DEM_Opt%SimDomain_min
    vol_tot= SimLen%x * SimLen%y * SimLen%z
    
    vol_local =(xed-xst)*(yed-yst)*(zed-zst)
    mlocal = int(vol_local/vol_tot*real(numPrtcl,kind=RK))
    mlocal = int(1.5_RK*mlocal)
    mlocal = min(mlocal, numPrtcl)
    mlocal = max(mlocal, 10)
        
    allocate(GPrtcl_id(mlocal),                Stat=iErr01)
    allocate(GPrtcl_pType(mlocal),             Stat=iErr02)
    allocate(GPrtcl_usrMark(mlocal),           Stat=iErr03)
    allocate(GPrtcl_PosR(mlocal),              Stat=iErr04)
    allocate(GPrtcl_linVel(this%tsize,mlocal), Stat=iErr05)
    allocate(GPrtcl_linAcc(this%tsize,mlocal), Stat=iErr06)
    allocate(GPrtcl_theta(mlocal),             Stat=iErr07)
    allocate(GPrtcl_rotVel(this%rsize,mlocal), Stat=iErr08)
    allocate(GPrtcl_rotAcc(this%rsize,mlocal), Stat=iErr09)
    allocate(GPrtcl_cntctForce(mlocal),        Stat=iErr10)
    allocate(GPrtcl_torque(mlocal),            Stat=iErr11)
    iErrSum=abs(iErr01)+abs(iErr02)+abs(iErr03)+abs(iErr04)+abs(iErr05)+abs(iErr06)+&
            abs(iErr07)+abs(iErr08)+abs(iErr09)+abs(iErr10)+abs(iErr11)
    if(iErrSum/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"GL_AllocateAllVar: ","Allocation failed 1")

    mlocalFix = int(vol_local/vol_tot*real(numPrtclFix,kind=RK))
    mlocalFix = int(1.5_RK*mlocalFix)
    mlocalFix = max(mlocalFix, 10)
    mlocalFix = min(mlocalFix, numPrtclFix)
    allocate(GPFix_id(mlocalFix),              Stat=iErr01)
    allocate(GPFix_pType(mlocalFix),           Stat=iErr02)
    allocate(GPFix_PosR(mlocalFix),            Stat=iErr03)
    iErrSum=abs(iErr01)+abs(iErr02)+abs(iErr03)
    if(iErrSum/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"GL_AllocateAllVar: ","Allocation failed 2")

    GPrtcl_id = 0
    GPrtcl_pType = 1
    GPrtcl_usrMark = 1
    GPrtcl_PosR = zero_r4
    GPrtcl_linVel  = zero_r3
    GPrtcl_linAcc  = zero_r3
    GPrtcl_theta   = zero_r3
    GPrtcl_rotVel  = zero_r3
    GPrtcl_rotAcc  = zero_r3
    GPrtcl_cntctForce = zero_r3
    GPrtcl_torque     = zero_r3
    GPFix_id = 0
    GPFix_pType = 1
    GPFix_PosR = zero_r4

    ! step2: Initialize this%mlocal,this%nlocal
    this%nlocal    = 0
    this%mlocal    = mlocal
    this%mlocalFix = mlocalFix

  end subroutine GL_AllocateAllVar

  !**********************************************************************
  ! GL_MakingAllPrtcl
  !**********************************************************************
  subroutine GL_MakingAllPrtcl(this, chFile)
    implicit none
    class(VarList)::this
    character(*),intent(in)::chFile
        
    ! locals         
    type(real3):: MkPrtclMinpoint_real3, MkPrtclMaxpoint_real3
    real(RK),dimension(3):: MkPrtclMinpoint, MkPrtclMaxpoint
    real(RK):: Distance_Ratio,VelMag
    character,dimension(3):: Fill_order
    NAMELIST /ParticleMakingOption/  MkPrtclMinpoint,MkPrtclMaxpoint, Fill_order, Distance_Ratio, VelMag
        
    integer::i,j,k,code,bin_id,numPrtcl,nUnitFile,myistat,nlocal,nlocal_sum
    type(real3):: cntr,l1_vec,l2_vec,l3_vec,lmin_p,lmax_p,dx,pos1,pos2
    integer,dimension(:),allocatable:: sum_bin
    real(RK)::xst,xed,yst,yed,zst,zed
    real(RK)::sum_divided, maxRad,rvelt,randt
    real(RK),dimension(3):: vel_dir

    nUnitFile = d_read_DEMParam_unit    
    open(unit=nUnitFile, file=chFile,status='old',form='formatted',IOSTAT=myistat)
    if(myistat/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl", "Cannot open file: "//trim(chFile))
    read(nUnitFile, nml=ParticleMakingOption)
    write(d_DEMLogInfo_unit, nml=ParticleMakingOption)
    close(nUnitFile,IOSTAT=myistat)
    MkPrtclMinpoint_real3 = MkPrtclMinpoint
    MkPrtclMaxpoint_real3 = MkPrtclMaxpoint
    numPrtcl = DEM_opt%numPrtcl
        
    xst=DEM_decomp%xSt; xed=DEM_decomp%xEd
    yst=DEM_decomp%ySt; yed=DEM_decomp%yEd
    zst=DEM_decomp%zSt; zed=DEM_decomp%zEd

    ! step3: Assign Prtcl_Type, Prtcl_id and Prtcl_PosR
    call Fill_Vectors(Fill_order, l1_vec, l2_vec, l3_vec)
    call system_clock(count=code); code=0
    call random_seed(size = j)
    call random_seed(put = code+63946*(/ (i - 1, i = 1, j) /))
    maxRad = maxval( DEMProperty%Prtcl_PureProp%Radius )
    lmin_p = MkPrtclMinpoint_real3 + Distance_Ratio* maxRad*one_r3
    lmax_p = MkPrtclMaxpoint_real3 - Distance_Ratio* maxRad*one_r3
    dx = Distance_Ratio *two*maxRad* one_r3
        
    allocate(sum_bin(DEM_opt%numPrtcl_Type))
    sum_bin(1)=DEMProperty%nPrtcl_in_Bin(1)
    do j=2, DEM_opt%numPrtcl_Type
      sum_bin(j)=sum_bin(j-1)+DEMProperty%nPrtcl_in_Bin(j)
    enddo

    cntr = lmin_p       ! start point
    nlocal=0
    do i=1,numPrtcl
      do j=1,DEM_opt%numPrtcl_Type
        if(i<=sum_bin(j)) then
          bin_id = j; exit
        endif
      enddo

      if(cntr%x>=xst.and.cntr%x<xed.and.cntr%y>=yst.and.cntr%y<yed.and.cntr%z>=zst.and.cntr%z<zed) then
        if(nlocal>=this%mlocal) call this%ReallocatePrtclVar(nlocal)
        nlocal=nlocal+1
        GPrtcl_id(nlocal) = i
        GPrtcl_pType(nlocal)= bin_id
        GPrtcl_PosR(nlocal) = cntr
        GPrtcl_PosR(nlocal)%w = DEMProperty%Prtcl_PureProp(bin_id)%Radius
      endif

      cntr = cntr + l1_vec * dx
      if((l1_vec.dot.cntr)>=(l1_vec.dot.lmax_p))then
        cntr = (lmin_p*l1_vec) + ((cntr+dx)*l2_vec)+(cntr*l3_vec)
        if((l2_vec.dot.cntr)>=(l2_vec.dot.lmax_p))then
          cntr = (cntr*l1_vec)+(lmin_p*l2_vec)+((cntr+dx)*l3_vec)
          if((l3_vec.dot.cntr) >= (l3_vec.dot.lmax_p) .and. nrank==0) then
            call DEMLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: ","Not enough space for positioning" ) 
          endif
        endif
      endif
    enddo

    call MPI_REDUCE(nlocal,nlocal_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,code)
    if(nlocal_sum/= numPrtcl .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: "," nlocal_sum/= numPrtcl " )
    endif
    this%nlocal = nlocal
    deallocate(sum_bin)

    ! randomize the GPrtcl_PosR
    do k=1,3
      do i=1,this%nlocal
        call random_number(randt)
        j=int(randt*this%nlocal)+1
        pos1=GPrtcl_PosR(i)
        pos2=GPrtcl_PosR(j)
        GPrtcl_PosR(i)=pos2
        GPrtcl_PosR(j)=pos1
      enddo
    enddo
        
    ! step3: assign the remaining variables
    do i=1,this%nlocal
      call random_number(vel_dir)
      vel_dir=two*vel_dir-one
      rvelt=sqrt(vel_dir(1)**2+vel_dir(2)**2+vel_dir(3)**2)
      if(rvelt>1.0E-10_RK) then
        vel_dir(1)=vel_dir(1)/rvelt
        vel_dir(2)=vel_dir(2)/rvelt
        vel_dir(3)=vel_dir(3)/rvelt
      else
        vel_dir(1)= -one*l3_vec%x
        vel_dir(2)= -one*l3_vec%y
        vel_dir(3)= -one*l3_vec%z
      endif
      GPrtcl_linVel(1,i)%x=VelMag*vel_dir(1)
      GPrtcl_linVel(1,i)%y=VelMag*vel_dir(2)
      GPrtcl_linVel(1,i)%z=VelMag*vel_dir(3)
    enddo

  end subroutine GL_MakingAllPrtcl

  !**********************************************************************
  ! Fill order Vector 
  !**********************************************************************
  subroutine Fill_Vectors( fill_order , l1 , l2, l3 )
    implicit none
    character,dimension(3) :: fill_order
    type(real3),intent(out):: l1, l2, l3
    
    if(fill_order(1)=="x" .or. fill_order(1)=="X") then
      l1 = real3(one,zero,zero)
    elseif(fill_order(1)=="y" .or. fill_order(1)=="Y") then
      l1 = real3(zero,one,zero)
    elseif(fill_order(1)=="z" .or. fill_order(1)=="Z") then
      l1 = real3(zero,zero,one)
    else
      call DEMLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: ","Fill_order wrong: 1 ")
    endif
            
    if(fill_order(2)=="x".or.fill_order(2)=="X") then
      l2 = real3(one,zero,zero)
    elseif(fill_order(2)=="y".or.fill_order(2)=="Y") then
      l2 = real3(zero,one,zero)
    elseif(fill_order(2)=="z".or.fill_order(2)=="Z") then
      l2 = real3(zero,zero,one)
    else
      call DEMLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: ","Fill_order wrong: 2 ")
    endif            
            
    if(fill_order(3)=="x".or.fill_order(3)=="X") then
      l3 = real3(one,zero,zero)
    elseif(fill_order(3)=="y".or.fill_order(3)=="Y") then
      l3 = real3(zero,one,zero) 
    elseif(fill_order(3)=="z".or.fill_order(3)=="Z") then
      l3 = real3(zero,zero,one)
    else
      call DEMLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: ","Fill_order wrong: 3 ")
    endif
  end subroutine Fill_Vectors

  !**********************************************************************
  ! copy i2 to i1
  !**********************************************************************
  subroutine GL_copy(this,i1,i2)
    implicit none
    class(VarList)::this
    integer,intent(in)::i1,i2

    if(i1==i2) return
    GPrtcl_id(i1)       = GPrtcl_id(i2)
    GPrtcl_pType(i1)    = GPrtcl_pType(i2)
    GPrtcl_usrMark(i1)  = GPrtcl_usrMark(i2)
    GPrtcl_PosR(i1)     = GPrtcl_PosR(i2)
    GPrtcl_linVel(:,i1) = GPrtcl_linVel(:,i2)
    GPrtcl_linAcc(:,i1) = GPrtcl_linAcc(:,i2)
    GPrtcl_theta(i1)    = GPrtcl_theta(i2)
    GPrtcl_rotVel(:,i1) = GPrtcl_rotVel(:,i2)
    GPrtcl_rotAcc(:,i1) = GPrtcl_rotAcc(:,i2)
  end subroutine GL_copy

  !**********************************************************************
  ! Reallocate particle varaibles 
  !**********************************************************************
  subroutine GL_ReallocatePrtclVar(this,np_new)
    implicit none
    class(VarList)::this
    integer,intent(in)::np_new

    ! locals
    integer::i
    integer:: sizep,sizen
    integer,dimension(:),allocatable:: IntVec
    type(real3),dimension(:),allocatable::Real3Vec
    type(real4),dimension(:),allocatable::Real4Vec
    type(real3),dimension(:,:),allocatable::Real3Arr

    sizep= this%mlocal
    sizen= int(1.2_RK*real(sizep,kind=RK))
    sizen= max(sizen, np_new+1)
    sizen= min(sizen,DEM_Opt%numPrtcl)
    this%mlocal=sizen

    ! ======= integer vector part =======
    call move_alloc(GPrtcl_id, IntVec)
    allocate(GPrtcl_id(sizen))
    GPrtcl_id(1:sizep)=IntVec
    GPrtcl_id(sizep+1:sizen)=0

    call move_alloc(GPrtcl_pType, IntVec)
    allocate(GPrtcl_pType(sizen))
    GPrtcl_pType(1:sizep)=IntVec
    GPrtcl_pType(sizep+1:sizen)=1

    call move_alloc(GPrtcl_usrMark, IntVec)
    allocate(GPrtcl_usrMark(sizen))
    GPrtcl_usrMark(1:sizep)=IntVec
    GPrtcl_usrMark(sizep+1:sizen)=1
    deallocate(IntVec)

    ! ======= real3 vercor part =======
    call move_alloc(GPrtcl_theta,Real3Vec)
    allocate(GPrtcl_theta(sizen))
    GPrtcl_theta(1:sizep)=Real3Vec
    GPrtcl_theta(sizep+1:sizen)=zero_r3
    deallocate(Real3Vec)
 
    deallocate(GPrtcl_cntctForce)
    allocate(GPrtcl_cntctForce(sizen))

    deallocate(GPrtcl_torque)
    allocate(GPrtcl_torque(sizen))
    
    ! ======= real3 matrix part =======
    call move_alloc(GPrtcl_linVel,Real3Arr)
    allocate(GPrtcl_linVel(this%tsize,sizen))
    GPrtcl_linVel(1:this%tsize,1:sizep)=Real3Arr
    GPrtcl_linVel(1:this%tsize,sizep+1:sizen)=zero_r3

    call move_alloc(GPrtcl_linAcc,Real3Arr)
    allocate(GPrtcl_linAcc(this%tsize,sizen))
    GPrtcl_linAcc(1:this%tsize,1:sizep)=Real3Arr
    GPrtcl_linAcc(1:this%tsize,sizep+1:sizen)=zero_r3

    call move_alloc(GPrtcl_rotVel,Real3Arr)
    allocate(GPrtcl_rotVel(this%rsize,sizen))
    GPrtcl_rotVel(1:this%rsize,1:sizep)=Real3Arr
    GPrtcl_rotVel(1:this%rsize,sizep+1:sizen)=zero_r3

    call move_alloc(GPrtcl_rotAcc,Real3Arr)
    allocate(GPrtcl_rotAcc(this%rsize,sizen))
    GPrtcl_rotAcc(1:this%rsize,1:sizep)=Real3Arr
    GPrtcl_rotAcc(1:this%rsize,sizep+1:sizen)=zero_r3  
    deallocate(Real3Arr) 

    ! ======= real4 vercor part =======
    call move_alloc(GPrtcl_PosR,Real4Vec)
    allocate(GPrtcl_PosR(sizen))
    GPrtcl_PosR(1:sizep)=Real4Vec
    GPrtcl_PosR(sizep+1:sizen)=zero_r4
    deallocate(Real4Vec)

    call DEMLogInfo%CheckForError(ErrT_Pass," ReallocatePrtclVar"," Need to reallocate particle variables")
    call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    call DEMLogInfo%OutInfo("Previous matirx length is :"//trim(num2str(sizep)),3)
    call DEMLogInfo%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),3)
  end subroutine GL_ReallocatePrtclVar

end module Prtcl_Variables
