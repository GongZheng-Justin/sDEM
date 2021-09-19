module Prtcl_Decomp_2d
  use MPI
  use Prtcl_TypeDef
  use Prtcl_Parameters
#ifdef CFDDEM
  use decomp_2d,only: nrank,nproc,ystart,yend
  use m_Parameters,only: p_row,p_col
  use m_MeshAndMetries,only : dx,dz
#endif
  implicit none
  private

#ifndef CFDDEM
  integer,save,public:: nrank  ! local MPI rank 
  integer,save,public:: nproc  ! total number of processors
#endif

  integer,save,public:: int_type,real_type,real3_type,real4_type
  integer,save,public:: int_byte,real_byte,real3_byte,real4_byte

  TYPE Prtcl_DECOMP_INFO

    ! define neighboring blocks
    ! second dimension 8 neighbour processors:
    !        1:4, 4 edge neighbours; 5:8, 4 cornor neighbours; 0, current processor(if any)
    integer::Prtcl_Pencil
    integer::prow
    integer::pcol
    integer::coord1
    integer::coord2
    integer,dimension(4):: ProcNgh
    real(RK)::xSt,ySt,zSt !min domain of the current and neighbor processors 
    real(RK)::xEd,yEd,zEd !max domain of the current and neighbor processors
  contains
    procedure:: Init_DECOMP => PDI_Init_DECOMP
  end type Prtcl_DECOMP_INFO

  TYPE(Prtcl_DECOMP_INFO), public :: DEM_decomp

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine PDI_Init_DECOMP(this, chFile)
    implicit none
    class(Prtcl_DECOMP_INFO)::this
    character(*),intent(in)::chFile

    !locals
    character::pencil
    logical::read_partition
    type(real3)::pMin, Pmax,SimLen
    integer::row,col,nUnitFile, myistat, coord1,coord2,idTop,idBottom,idLeft,idRight
    NAMELIST /PrtclDomainDecomp/row,col,pencil,read_partition
#ifndef CFDDEM
    real(RK)::xpart,ypart,zpart
#endif

    call Init_Prctl_MPI_TYPE()

    nUnitFile = GetFileUnit()    
    open( unit=nUnitFile, file=chFile,status='old',form='formatted',IOSTAT=myistat )
    if(myistat /= 0 .and. nrank==0) then
      print*, "  PDI_Init_DECOMP: Cannot open file"//trim(adjustl(chFile))
      STOP
    endif
    read(nUnitFile, nml=PrtclDomainDecomp)
    close(nUnitFile)

#ifndef CFDDEM
    if(row*col /= nproc .and. nrank==0) then
      STOP "  PDI_Init_DECOMP: Invalid 2D processor grid- nproc/= row*col"
    endif
#else
    row= p_row
    col= p_col
#endif
    this%prow= row
    this%pcol= col

    if(pencil == "x" .or. pencil == "X") then
      this%Prtcl_Pencil = x_axis
    elseif(pencil == "y" .or. pencil == "Y") then
      this%Prtcl_Pencil = y_axis
    elseif(pencil == "z" .or. pencil == "Z") then
      this%Prtcl_Pencil = z_axis
    else
      if(nrank==0) STOP "  PDI_Init_DECOMP: Invalid pencil sign"
    endif

! ------------------------- neigbor information begins------------------------
! 
! 2D domain decomposition method in DEM_MPI(from left to right:  x-pencil, y-pencil, z-pencil)
! 
!   If we have 6 processors = 3 row * 2 col
! 
!     the arrangement of the subdomains(nrank) is as follow:
!       y               x               x
!       |        4 5    |        4 5    |        4 5   
!       |        2 3    |        2 3    |        2 3
!       |_ _ _z  0 1    |_ _ _z  0 1    |_ _ _y  0 1
!
!     the arrangement of the coord1 is as follow:
!       y               x               x
!       |        2 2    |        2 2    |        2 2
!       |        1 1    |        1 1    |        1 1
!       |_ _ _z  0 0    |_ _ _z  0 0    |_ _ _y  0 0
! 
!     the arrangement of the coord2 is as follow:
!       y               x               x
!       |        0 1    |        0 1    |        0 1
!       |        0 1    |        0 1    |        0 1
!       |_ _ _ z 0 1    |_ _ _z  0 1    |_ _ _y  0 1
! 
!     neighbor index:
! 
!       y               x               x
!       |        6 3 5  |        6 3 5  |        6 3 5
!       |        2 0 1  |        2 0 1  |        2 0 1
!       |_ _ _ z 7 4 8  |_ _ _z  7 4 8  |_ _ _y  7 4 8
! 
!        Here 0 means the center subdomain, and 1-8 stands for the relative location of the eight neighbors

    
    coord1 = int ( nrank / col)
    coord2 = mod ( nrank,  col)
    this%coord1 = coord1
    this%coord2 = coord2

    ! Firstly, all the boundaries are assumed to be periodic
    idTop     = mod(coord1+1,    row)  ! top
    idBottom  = mod(coord1+row-1,row)  ! bottom 
    idLeft    = mod(coord2+col-1,col)  ! left
    idRight   = mod(coord2+1,    col)  ! right
    this%ProcNgh(1) = coord1   * col + idRight 
    this%ProcNgh(2) = coord1   * col + idLeft  
    this%ProcNgh(3) = idTop    * col + coord2  
    this%ProcNgh(4) = idBottom * col + coord2  
   ! this%ProcNgh(5) = idTop    * col + idRight
   ! this%ProcNgh(6) = idTop    * col + idLeft
   ! this%ProcNgh(7) = idBottom * col + idLeft
   ! this%ProcNgh(8) = idBottom * col + idRight

   ! Secondly, modify the edge neighbour ids
    IF(coord1==0) THEN
      if((.not.DEM_Opt%IsPeriodic(2)) .and. this%Prtcl_Pencil == x_axis) this%ProcNgh(4)=MPI_PROC_NULL
      if((.not.DEM_Opt%IsPeriodic(1)) .and. this%Prtcl_Pencil == y_axis) this%ProcNgh(4)=MPI_PROC_NULL
      if((.not.DEM_Opt%IsPeriodic(1)) .and. this%Prtcl_Pencil == z_axis) this%ProcNgh(4)=MPI_PROC_NULL
    ENDIF
    IF(coord1==row-1) THEN
      if((.not.DEM_Opt%IsPeriodic(2)) .and. this%Prtcl_Pencil == x_axis) this%ProcNgh(3)=MPI_PROC_NULL
      if((.not.DEM_Opt%IsPeriodic(1)) .and. this%Prtcl_Pencil == y_axis) this%ProcNgh(3)=MPI_PROC_NULL
      if((.not.DEM_Opt%IsPeriodic(1)) .and. this%Prtcl_Pencil == z_axis) this%ProcNgh(3)=MPI_PROC_NULL
    ENDIF
    IF(coord2==0) THEN
      if((.not.DEM_Opt%IsPeriodic(3)) .and. this%Prtcl_Pencil == x_axis) this%ProcNgh(2)=MPI_PROC_NULL  
      if((.not.DEM_Opt%IsPeriodic(3)) .and. this%Prtcl_Pencil == y_axis) this%ProcNgh(2)=MPI_PROC_NULL
      if((.not.DEM_Opt%IsPeriodic(2)) .and. this%Prtcl_Pencil == z_axis) this%ProcNgh(2)=MPI_PROC_NULL
    ENDIF
    IF(coord2==col-1) THEN
      if((.not.DEM_Opt%IsPeriodic(3)) .and. this%Prtcl_Pencil == x_axis) this%ProcNgh(1)=MPI_PROC_NULL 
      if((.not.DEM_Opt%IsPeriodic(3)) .and. this%Prtcl_Pencil == y_axis) this%ProcNgh(1)=MPI_PROC_NULL  
      if((.not.DEM_Opt%IsPeriodic(2)) .and. this%Prtcl_Pencil == z_axis) this%ProcNgh(1)=MPI_PROC_NULL
    ENDIF

    pMin = DEM_Opt%SimDomain_min
    pMax = DEM_Opt%SimDomain_max
    SimLen = pMax - pMin
#ifndef CFDDEM
    IF(read_partition) THEN

    ELSE
      if(this%Prtcl_Pencil ==x_axis) then
        ypart = SimLen%y/real(row,kind=RK)
        zpart = SimLen%z/real(col,kind=RK)
        this%xSt= pMin%x
        this%xEd= pMax%x
        this%ySt= pMin%y + ypart*real(coord1,  kind=RK)
        this%yEd= pMin%y + ypart*real(coord1+1,kind=RK)
        this%zSt= pMin%z + zpart*real(coord2,  kind=RK)
        this%zEd= pMin%z + zpart*real(coord2+1,kind=RK)
      elseif(this%Prtcl_Pencil ==y_axis) then
        xpart = SimLen%x/real(row,kind=RK)
        zpart = SimLen%z/real(col,kind=RK)
        this%xSt= pMin%x + xpart*real(coord1,  kind=RK)
        this%xEd= pMin%x + xpart*real(coord1+1,kind=RK)
        this%ySt= pMin%y
        this%yEd= pMax%y
        this%zSt= pMin%z + zpart*real(coord2,  kind=RK)
        this%zEd= pMin%z + zpart*real(coord2+1,kind=RK)
      elseif(this%Prtcl_Pencil ==z_axis) then
        xpart = SimLen%x/real(row,kind=RK)
        ypart = SimLen%y/real(col,kind=RK)
        this%xSt= pMin%x + xpart*real(coord1,  kind=RK)
        this%xEd= pMin%x + xpart*real(coord1+1,kind=RK)
        this%ySt= pMin%y + ypart*real(coord2,  kind=RK)
        this%yEd= pMin%y + ypart*real(coord2+1,kind=RK)
        this%zSt= pMin%z
        this%zEd= pMax%z
      endif
    ENDIF
#else
    if(this%Prtcl_Pencil /=y_axis .and. nrank==0) then
      STOP "In the CFDDEM model, only y-pencil is avaiable now"
    endif
    this%xSt= real(ystart(1)-1, kind=RK)*dx
    this%xEd= real(yend(1),     kind=RK)*dx
    this%ySt= pMin%y
    this%yEd= pMax%y
    this%zSt= real(ystart(3)-1, kind=RK)*dz
    this%zEd= real(yend(3),     kind=RK)*dz
#endif

  end subroutine PDI_Init_DECOMP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Init_Prctl_MPI_TYPE()
    implicit none
    integer::ierror
    integer,dimension(4)::disp,blocklen,blocktype
  
    ! integer
    int_type = MPI_INTEGER
    call MPI_TYPE_SIZE(int_type,   int_byte,   ierror)

    ! real
    if(RK==4) then
      real_type = MPI_REAL
    else
      real_type = MPI_DOUBLE_PRECISION
    endif
    call MPI_TYPE_SIZE(real_type,  real_byte,  ierror)

    ! real3 type
    blocklen(1:3)=1
    blocktype(1:3)=real_type
    disp(1)=0
    disp(2)=disp(1)+real_byte
    disp(3)=disp(2)+real_byte
    call MPI_TYPE_STRUCT(3,blocklen(1:3),disp(1:3),blocktype(1:3),real3_type,ierror)
    call MPI_TYPE_COMMIT(real3_type,ierror)
    call MPI_TYPE_SIZE(real3_type, real3_byte, ierror)

    ! real4 type
    blocklen(1:4)=1
    blocktype(1:4)=real_type
    disp(1)=0
    disp(2)=disp(1)+real_byte
    disp(3)=disp(2)+real_byte
    disp(4)=disp(3)+real_byte
    call MPI_TYPE_STRUCT(4,blocklen(1:4),disp(1:4),blocktype(1:4),real4_type,ierror)
    call MPI_TYPE_COMMIT(real4_type,ierror)
    call MPI_TYPE_SIZE(real4_type, real4_byte, ierror)
  end subroutine Init_Prctl_MPI_TYPE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module Prtcl_Decomp_2d
