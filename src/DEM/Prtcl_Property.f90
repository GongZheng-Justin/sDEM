module Prtcl_Property
  use Prtcl_TypeDef
  use Prtcl_LogInfo
  use Prtcl_Parameters
#ifdef CFDDEM
  use decomp_2d,only: nrank
  use m_Parameters,only: FluidDensity
#else
  use Prtcl_decomp_2d,only: nrank
#endif
  implicit none
  private
    
  type PureProperty
    real(RK):: Radius = 0.005_RK
    real(RK):: Density= 2500.0_RK
    real(RK):: PoissonRatio= 0.25_RK
    real(RK):: YoungsModulus=5.0E6_RK
    real(RK):: Mass
    real(RK):: Volume
    real(RK):: Inertia
#ifdef CFDDEM
    real(RK):: MassInFluid
    real(RK):: MassOfFluid
#endif
  end type PureProperty

  type,public:: BinaryProperty
    real(RK):: RadEff       ! effective Radiusius
    real(RK):: MassEff      ! effective Mass

    real(RK):: StiffnessCoe_n
    real(RK):: StiffnessCoe_t        
    real(RK):: DampingCoe_n = zero
    real(RK):: DampingCoe_t = zero
    real(RK):: RestitutionCoe_n = 0.95_RK  ! Normal Resitution Coefficient
    real(RK):: RestitutionCoe_t = 0.95_RK  ! Tangential Resitution Coefficient       
    real(RK):: FrictionCoe_s = 0.80_RK     ! Coefficient of static  friction
    real(RK):: FrictionCoe_k = 0.15_RK     ! Coefficient of kinetic friction
    real(RK):: FrictionCoe_Roll = 0.10_RK  ! Coefficient of rolling friction
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
  end type PhysicalProperty

  type(PhysicalProperty),public::DEMProperty

contains

  !*******************************************************
  ! initializing the size distribution with property 
  !*******************************************************
  subroutine InitPrtclProperty(this,chFile)
    implicit none
    class(PhysicalProperty)::this
    character(*),intent(in)::chFile
        
    ! locals
    real(RK),dimension(:),allocatable:: Bin_Divided,Density,Diameter,YoungsModulus_P,PoissonRatio_P
    real(RK)::FrictionCoe_s_PP,FrictionCoe_k_PP,FrictionCoe_Roll_PP,RestitutionCoe_n_PP,RestitutionCoe_t_PP
    namelist/ParticlePhysicalProperty/Bin_Divided, Density, Diameter,YoungsModulus_P,PoissonRatio_P, &
             FrictionCoe_s_PP,FrictionCoe_k_PP,FrictionCoe_Roll_PP,RestitutionCoe_n_PP,RestitutionCoe_t_PP
    integer:: i,j,code,nPType,nUnitFile,myistat, sum_prtcl,bin_pnum,prdiff,bin_id
    real(RK):: sum_divided,rtemp
    type(PureProperty)::pari,parj
    type(BinaryProperty)::Bnry
        
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
        
    nUnitFile = getFileUnit()    
    open(unit=nUnitFile, file=chFile, status='old', form='formatted', IOSTAT=myistat)
    if(myistat /= 0 ) call DEMLogInfo%CheckForError(ErrT_Abort,"InitPrtclProperty","Cannot open file:"//trim(chFile))
    read(nUnitFile, nml=ParticlePhysicalProperty)
    write(DEMLogInfo%nUnit, nml=ParticlePhysicalProperty)
    close(nUnitFile,IOSTAT=myistat)
        
    ! calculate this%nPrtcl_in_Bin
    call system_clock(count=code) !code=0
    call random_seed(size = j)
    call random_seed(put = code+63946*(/(i-1,i=1,j)/))
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
#ifdef CFDDEM
      this%Prtcl_PureProp(i)%MassInFluid= (Density(i)-FluidDensity)*this%Prtcl_PureProp(i)%Volume
      this%Prtcl_PureProp(i)%MassOfFluid= FluidDensity*this%Prtcl_PureProp(i)%Volume
#endif
     end do
        
     do i=1,nPType
       do j=1,nPType
        pari = this%Prtcl_PureProp(i)
        parj = this%Prtcl_PureProp(i)
        Bnry= clc_BnryPrtcl_Prop(pari,parj,FrictionCoe_s_PP,FrictionCoe_k_PP,FrictionCoe_Roll_PP,RestitutionCoe_n_PP,RestitutionCoe_t_PP,.false.)
        this%Prtcl_BnryProp(i,j)= Bnry  
      enddo
    enddo
  end subroutine InitPrtclProperty

  !*****************************************************************************    
  ! setting the physical property of wall 
  !*****************************************************************************
  subroutine InitWallProperty(this,chFile )
    implicit none
    class(PhysicalProperty)::this
    character(*),intent(in)::chFile
        
    ! locals
    real(RK),dimension(:),allocatable:: YoungsModulus_W,PoissonRatio_W
    real(RK)::FrictionCoe_s_PW,FrictionCoe_k_PW,FrictionCoe_Roll_PW,RestitutionCoe_n_PW,RestitutionCoe_t_PW
    namelist /WallPhysicalProperty/YoungsModulus_W,PoissonRatio_W,FrictionCoe_s_PW,FrictionCoe_k_PW, &
              FrictionCoe_Roll_PW,RestitutionCoe_n_PW,RestitutionCoe_t_PW     
    integer:: i, j, nPType,nWType,nUnitFile,myistat
    type(PureProperty) pari, wall
        
    nPType = DEM_opt%numPrtcl_Type
    nWType = DEM_opt%numWall_type
    allocate(YoungsModulus_W(nWType) )
    allocate(PoissonRatio_W(nWType) )
    allocate( this%Wall_PureProp(nWType) )
    allocate( this%PrtclWall_BnryProp(nPType, nWType) )
 
    nUnitFile = GetFileUnit()    
    open(unit=nUnitFile, file=chFile, status='old', form='formatted', IOSTAT=myistat)
    if(myistat /= 0 .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"InitWallProperty: " ,"Cannot open file:"//trim(chFile))
    endif
    read(nUnitFile, nml=WallPhysicalProperty)
    close(nUnitFile)        
        
    do i=1,nWType
       this%Wall_PureProp(i)%Density = 1.0E50_RK
       this%Wall_PureProp(i)%Radius  = 1.0E50_RK
       this%Wall_PureProp(i)%YoungsModulus= YoungsModulus_W(i)
       this%Wall_PureProp(i)%PoissonRatio = PoissonRatio_W(i)
           
       this%Wall_PureProp(i)%Mass    = 1.0E50_RK
       this%Wall_PureProp(i)%Volume  = 1.0E50_RK
       this%Wall_PureProp(i)%Inertia = 1.0E50_RK
    enddo

    do i = 1, nPType
      do j = 1, nWType
        pari = this%Prtcl_PureProp(i)
        wall = this%Wall_PureProp(j)
        this%PrtclWall_BnryProp(i,j)=clc_BnryPrtcl_Prop(pari,wall,FrictionCoe_s_PW,FrictionCoe_k_PW,FrictionCoe_Roll_PW, &
                                                             RestitutionCoe_n_PW,RestitutionCoe_t_PW,.true.)
      enddo
    enddo       
  end subroutine InitWallProperty

  !*****************************************************************************    
  ! Calculating the binary contact properties
  !*****************************************************************************
  function clc_BnryPrtcl_Prop(pari,parj,FrictionCoe_s,FrictionCoe_k,FrictionCoe_Roll, &
                              RestitutionCoe_n,RestitutionCoe_t,iswall) result(Bnry)
    implicit none
    class(PureProperty),intent(in):: pari, parj
    real(RK),intent(in):: FrictionCoe_s,FrictionCoe_k,FrictionCoe_Roll,RestitutionCoe_n,RestitutionCoe_t
    logical,intent(in) :: iswall

    type(BinaryProperty)::Bnry
    real(RK)::kappa,pri,prj,ShearModi,ShearModj,YoungsModEff,ShearModEff,Eta,K_hertz
    
    if(.not.iswall) then
      Bnry%RadEff = (pari%Radius * parj%Radius)/(pari%Radius + parj%Radius)
      Bnry%MassEff= (pari%Mass * parj%Mass )/( pari%Mass+parj%Mass )
    else
      Bnry%RadEff = pari%Radius
      Bnry%MassEff= pari%Mass         
    endif
    Bnry%RestitutionCoe_n = RestitutionCoe_n
    Bnry%RestitutionCoe_t = RestitutionCoe_t
    Bnry%FrictionCoe_s = FrictionCoe_s
    Bnry%FrictionCoe_k = FrictionCoe_k
    Bnry%FrictionCoe_Roll = FrictionCoe_Roll

    pri=pari%PoissonRatio
    prj=parj%PoissonRatio
    ShearModi=pari%YoungsModulus/(two*(one+pri))
    ShearModj=parj%YoungsModulus/(two*(one+prj))
    YoungsModEff=one/((one-pri*pri)/pari%YoungsModulus + (one-prj*prj)/parj%YoungsModulus) ! 2.47, p31
    ShearModEff =one/((two- pri)/ShearModi + (two-prj)/ShearModj)                          ! 2.71, p43
    kappa = ((one-pri)/ShearModi + (one-prj)/ShearModj) / &
            ((one-half*pri)/ShearModi +(one-half*prj)/ShearModj )
    kappa=abs(kappa)

    Eta=log(RestitutionCoe_n); Eta=Eta*Eta
    if(DEM_opt%CF_Type == CFT_LSD) then                                             
      Bnry%StiffnessCoe_n = 1.2024_RK*(sqrt(Bnry%MassEff)*(YoungsModEff**2)*Bnry%RadEff)**(0.4_RK)        ! 2.46, p31
      Bnry%StiffnessCoe_t = kappa*Bnry%StiffnessCoe_n                                                     ! 2.52, p34
      Bnry%DampingCoe_n=-two*log(RestitutionCoe_n)*sqrt(Bnry%MassEff*Bnry%StiffnessCoe_n)/sqrt(PI*PI+Eta) ! 2.44, p31
      Bnry%DampingCoe_t= zero ! No equation is considered for tangential damping yet

    elseif(DEM_opt%CF_Type == CFT_nLin) then
      K_hertz =four/three*YoungsModEff*sqrt(Bnry%RadEff)               ! 2.62, P39
      Bnry%StiffnessCoe_n= K_hertz
      Bnry%StiffnessCoe_t= sixteen/three*ShearModEff*sqrt(Bnry%RadEff) ! 2.72, P44
      Bnry%DampingCoe_n  = -2.2664_RK*log(RestitutionCoe_n)*sqrt(Bnry%MassEff*K_hertz)/sqrt(Eta+10.1354_RK) ! 2.63, p40
      Bnry%DampingCoe_t  = zero ! No equation is considered for tangential damping yet
    endif
  end function

end module Prtcl_Property
