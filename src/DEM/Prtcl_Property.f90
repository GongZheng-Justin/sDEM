module Prtcl_Property
  use Prtcl_TypeDef
  use Prtcl_LogInfo
  use Prtcl_Parameters
  use Prtcl_decomp_2d,only: nrank
  implicit none
  private
  
  real(RK),parameter:: d_Density       = 2500.0_RK
  real(RK),parameter:: d_Radius        = 0.005_RK  
  real(RK),parameter:: d_PoissonRatio  = 0.25_RK
  real(RK),parameter:: d_YoungsModulus = 1.0E6_RK
    
  real(RK),parameter:: d_FrictionCoe      = 0.3_RK
  real(RK),parameter:: d_FrictionCoe_Roll = 0.1_RK
  real(RK),parameter:: d_RestitutionCoe_n = 0.75_RK
  real(RK),parameter:: d_RestitutionCoe_t = one    
  real(RK),parameter:: d_rel_vel          = one    
    
  type PureProperty
    real(RK):: Radius   = d_Radius
    real(RK):: Density  = d_Density
    real(RK):: YoungsModulus =d_YoungsModulus
    real(RK):: PoissonRatio  = d_PoissonRatio
        
    real(RK):: ShearModulus 
    real(RK):: Mass
    real(RK):: Volume
    real(RK):: Inertia
  end type PureProperty

  type,public:: BinaryProperty
    real(RK):: YoungsModEff ! effective Young's modulus
    real(RK):: ShearModEff  ! effective shear modulus
    real(RK):: RadEff       ! effective Radiusius
    real(RK):: MassEff      ! effective Mass
        
    real(RK):: RestitutionCoe_n = d_RestitutionCoe_n 
    real(RK):: RestitutionCoe_t = d_RestitutionCoe_t        
    real(RK):: FrictionCoe      = d_FrictionCoe
    real(RK):: FrictionCoe_Roll = d_FrictionCoe_Roll
    real(RK):: rel_vel_kn       = d_rel_vel      ! relative particle velocity for calculating spring stiffness

    real(RK):: LinearSpringStiffness_n
    real(RK):: LinearSpringStiffness_t        
    real(RK):: dmp_n = zero
    real(RK):: dmp_t = zero
  end type BinaryProperty
    
  type PhysicalProperty
    integer,allocatable,dimension(:) :: nPrtcl_in_Bin
    integer,allocatable,dimension(:) :: CS_Hrchl_level ! level for NBS-Munjiza-Hierarchy Contact Search Method  
    type(pureProperty),allocatable,dimension(:) :: Prtcl_PureProp
    type(pureProperty),allocatable,dimension(:) :: Wall_PureProp
    type(BinaryProperty),allocatable,dimension(:,:):: Prtcl_BnryProp    
    type(BinaryProperty),allocatable,dimension(:,:):: PrtclWall_BnryProp         
  contains
    procedure:: InitPrtclProperty
    procedure:: InitWallProperty
    procedure,private:: clc_BnryPrtcl_Prop  => PhP_clc_BnryPrtcl_Prop
  end type PhysicalProperty

  type(PhysicalProperty),public::DEMProperty

contains

  !*******************************************************
  ! initializing the size distribution with property 
  !*******************************************************
  subroutine InitPrtclProperty( this,  chFile )
    implicit none
    class(PhysicalProperty)::this
    character(*),intent(in)::chFile
        
    ! locals
    real(RK),dimension(:),allocatable:: Bin_Divided, Density, Diameter,YoungsModulus_P, PoissonRatio_P
    real(RK)::FrictionCoe_PP, FrictionCoe_Roll_PP, RestitutionCoe_n_PP, RestitutionCoe_t_PP
    namelist/ParticlePhysicalProperty/Bin_Divided, Density, Diameter,YoungsModulus_P, PoissonRatio_P, &
                                      FrictionCoe_PP,FrictionCoe_Roll_PP,RestitutionCoe_n_PP,RestitutionCoe_t_PP
    integer:: i,j,code,nPType,nUnitFile,myistat, sum_prtcl,bin_pnum,prdiff,bin_id
    real(RK):: sum_divided,rtemp
    type(PureProperty):: pari, parj 
        
    nPType  = DEM_opt%numPrtcl_Type
    allocate( Bin_Divided(nPType))
    allocate( Density(nPType))
    allocate( Diameter(nPType))
    allocate( YoungsModulus_P(nPType))
    allocate( PoissonRatio_P(nPType))
        
    allocate( this%CS_Hrchl_level(nPType))
    allocate( this%nPrtcl_in_Bin(nPType))
    allocate( this%Prtcl_PureProp(nPType))
    allocate( this%Prtcl_BnryProp(nPType,nPType)) 
        
    nUnitFile = d_read_DEMParam_unit    
    open(unit=nUnitFile, file=chFile, status='old', form='formatted', IOSTAT=myistat)
    if(myistat /= 0 ) call DEMLogInfo%CheckForError(ErrT_Abort,"InitPrtclProperty","Cannot open file:"//trim(chFile))
    read(nUnitFile, nml=ParticlePhysicalProperty)
    write(d_DEMLogInfo_unit, nml=ParticlePhysicalProperty)
    close(nUnitFile,IOSTAT=myistat)
        
    ! calculate this%nPrtcl_in_Bin
    sum_divided=zero
    do i=1,nPType
      sum_divided = sum_divided + Bin_Divided(i)
    enddo
    sum_prtcl = 0
    do i=1,nPType
      bin_pnum = int(DEM_opt%numPrtcl*Bin_Divided(i)/sum_divided)
      this%nPrtcl_in_Bin(i)= bin_pnum
      sum_prtcl = sum_prtcl + bin_pnum
    enddo
    prdiff = DEM_opt%numPrtcl - sum_prtcl
    if( prdiff > 0 ) then
      call system_clock(count=code)
      call random_seed(size = j)
      call random_seed(put = code+63946*(/(i-1,i=1,j)/))
      do i=1, prdiff 
        call random_number(rtemp)
        bin_id = int(rtemp*nPType) + 1
        this%nPrtcl_in_Bin(bin_id) = this%nPrtcl_in_Bin(bin_id)  + 1
        enddo
     endif
        
     ! calculate particle properties
     do i = 1, nPType
       this%Prtcl_PureProp(i)%Density = Density(i)
       this%Prtcl_PureProp(i)%Radius = half * Diameter(i)
       this%Prtcl_PureProp(i)%YoungsModulus = YoungsModulus_P(i)
       this%Prtcl_PureProp(i)%PoissonRatio = PoissonRatio_P(i)

       this%Prtcl_PureProp(i)%Volume = four/three*Pi*(this%Prtcl_PureProp(i)%Radius)**3
       this%Prtcl_PureProp(i)%Mass = Density(i)*this%Prtcl_PureProp(i)%Volume
       this%Prtcl_PureProp(i)%Inertia = two/five*this%Prtcl_PureProp(i)%Mass*(this%Prtcl_PureProp(i)%Radius**2)
       this%Prtcl_PureProp(i)%ShearModulus = YoungsModulus_P(i)/(two*(one+PoissonRatio_P(i)))
     end do
        
     do i=1, nPType
       do j=1,nPType
        pari = this%Prtcl_PureProp(i)
        parj = this%Prtcl_PureProp(i)
        this%Prtcl_BnryProp(i,j)= this%clc_BnryPrtcl_Prop(pari,parj,FrictionCoe_PP,FrictionCoe_Roll_PP, &
                                                          RestitutionCoe_n_PP,RestitutionCoe_t_PP,.false.) 
      end do
    end do        
        
  end subroutine InitPrtclProperty

  !*****************************************************************************    
  ! setting the physical property of wall 
  !*****************************************************************************
  subroutine InitWallProperty(this, chFile )
    implicit none
    class(PhysicalProperty)::this
    character(*),intent(in)::chFile
        
    ! locals
    real(RK),dimension(:),allocatable:: YoungsModulus_W, PoissonRatio_W
    real(RK)::FrictionCoe_PW, FrictionCoe_Roll_PW, RestitutionCoe_n_PW, RestitutionCoe_t_PW
    namelist /WallPhysicalProperty/ YoungsModulus_W, PoissonRatio_W, FrictionCoe_PW,&
                                    FrictionCoe_Roll_PW,RestitutionCoe_n_PW,RestitutionCoe_t_PW        
    integer:: i, j, nPType,nWType,nUnitFile,myistat
    type(PureProperty) pari, wall
        
    nPType = DEM_opt%numPrtcl_Type
    nWType = DEM_opt%numWall_type
    allocate(YoungsModulus_W(nWType) )
    allocate(PoissonRatio_W(nWType) )
        
    allocate( this%Wall_PureProp(nWType) )
    allocate( this%PrtclWall_BnryProp(nPType, nWType) )
 
    nUnitFile = d_read_DEMParam_unit    
    open(unit=nUnitFile, file=chFile, status='old', form='formatted', IOSTAT=myistat)
    if(myistat /= 0 .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"InitWallProperty" ,"Cannot open file:"//trim(chFile))
    endif
    read(nUnitFile, nml=WallPhysicalProperty)
    write(d_DEMLogInfo_unit, nml=WallPhysicalProperty)
    close(nUnitFile)        
        
    do i=1,nWType
       this%Wall_PureProp(i)%Density = 1.0E50_RK
       this%Wall_PureProp(i)%Radius  = 1.0E50_RK
       this%Wall_PureProp(i)%YoungsModulus = YoungsModulus_W(i)
       this%Wall_PureProp(i)%PoissonRatio = PoissonRatio_W(i)
           
       this%Wall_PureProp(i)%Mass = 1.0E50_RK
       this%Wall_PureProp(i)%Volume  = 1.0E50_RK
       this%Wall_PureProp(i)%Inertia = 1.0E50_RK
       this%Wall_PureProp(i)%ShearModulus = YoungsModulus_W(i)/(two*(one+PoissonRatio_W(i)))
    end do

    do i = 1, nPType
      do j = 1, nWType
        pari = this%Prtcl_PureProp(i)
        wall = this%Wall_PureProp(j)
        this%PrtclWall_BnryProp(i,j)=this%clc_BnryPrtcl_Prop(pari,wall,FrictionCoe_PW,FrictionCoe_Roll_PW, &
                                                             RestitutionCoe_n_PW, RestitutionCoe_t_PW, .true. ) 
      end do
    end do       
  end subroutine InitWallProperty

  !*****************************************************************************    
  ! Calculating the binary contact properties
  !*****************************************************************************
  function PhP_clc_BnryPrtcl_Prop(this,pari,parj,FrictionCoe,FrictionCoe_Roll,RestitutionCoe_n, &
                                  RestitutionCoe_t,iswall) result( Bnry )
    implicit none
    class(PhysicalProperty) this
    class(PureProperty),intent(in):: pari, parj
    real(RK),intent(in):: FrictionCoe, FrictionCoe_Roll, RestitutionCoe_n, RestitutionCoe_t
    logical,intent(in) :: iswall
    type(BinaryProperty):: Bnry
    real(RK):: kappa, pri,prj
    
    if(.not.iswall) then
      Bnry%RadEff = (pari%Radius * parj%Radius)/(pari%Radius + parj%Radius)
      Bnry%MassEff = (pari%Mass * parj%Mass )/( pari%Mass+parj%Mass )
    else
      Bnry%RadEff = pari%Radius
      Bnry%MassEff = pari%Mass         
    endif
    pri = pari%PoissonRatio ! Possion ratio i
    prj = parj%PoissonRatio ! Possion ratio j
        
    Bnry%YoungsModEff=one/((one-pri**2)/pari%YoungsModulus + (one-prj**2)/parj%YoungsModulus)   ! 2.47, p31
    Bnry%ShearModEff= one/((two- pri)/pari%ShearModulus + (two-prj)/parj%ShearModulus)          ! 2.71, p43
    Bnry%LinearSpringStiffness_n = 1.2024_RK*(sqrt(Bnry%MassEff) * (Bnry%YoungsModEff**2)*Bnry%RadEff* &
                                                Bnry%rel_vel_kn )**(two/five) ! 2.46, p31
            
    kappa = ((one-pri)/pari%ShearModulus + (one-prj)/parj%ShearModulus) / &
            ((one-half*pri)/pari%ShearModulus +(one-half*prj)/parj%ShearModulus )
    Bnry%LinearSpringStiffness_t = kappa*Bnry%LinearSpringStiffness_n   ! 2.52, p34
    
    Bnry%RestitutionCoe_n = RestitutionCoe_n
    Bnry%RestitutionCoe_t = RestitutionCoe_t
    Bnry%FrictionCoe = FrictionCoe
    Bnry%FrictionCoe_Roll = FrictionCoe_Roll
    Bnry%dmp_n = Damping_normal(Bnry) 
    Bnry%dmp_t = Damping_tangential(Bnry)
  end function

  !*****************************************************************************    
  ! Damping_normal
  !*****************************************************************************
  function Damping_normal( Bnry ) result (res)
    implicit none  
    type(BinaryProperty):: bnry
    real(RK):: res, K_hertz
        
    select case( DEM_opt%CF_Type )
    case(CFT_LSD_nl, CFT_LSD_l)
      res = -two*log(Bnry%RestitutionCoe_n)*sqrt(Bnry%MassEff*Bnry%LinearSpringStiffness_n)/sqrt((log(Bnry%RestitutionCoe_n))**2 + pi**two)  ! 2.51, p34
    case(CFT_nLin_nl, CFT_nLin_l )
      K_hertz = four/three*Bnry%YoungsModEff*sqrt(Bnry%RadEff)    ! 2.64, p40
      res=-2.2664_RK*log(Bnry%RestitutionCoe_n)*sqrt(Bnry%MassEff*K_hertz)/sqrt(log(Bnry%RestitutionCoe_n)**2 + 10.1354_RK) ! 2.63, p40
    end select
    
  end function Damping_normal

  !*****************************************************************************    
  ! Damping_tangential
  !*****************************************************************************
  function Damping_tangential( Bnry ) result (res)
    implicit none
    type(BinaryProperty) bnry
    real(RK) res
        
    res = zero
    if(nrank==0)print*, "No equation is considered for tangential damping yet"
  end function Damping_tangential

end module Prtcl_Property
