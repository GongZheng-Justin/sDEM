module Prtcl_IOAndVisu
  use MPI
  use Prtcl_TypeDef
  use Prtcl_LogInfo
  use Prtcl_Parameters
  use Prtcl_Variables
  use Prtcl_CL_and_CF
  use Prtcl_Decomp_2d
  use Prtcl_Comm
  use Prtcl_Property
  implicit none
  private

  integer,parameter:: IK = 4
  integer,parameter:: END_FLAG = -999
  integer,save:: Prev_BackUp_itime  = 1234567
  logical::saveXDMFOnce, save_ID,       save_Diameter, save_Type,   save_UsrMark,    save_LinVel
  logical::save_LinAcc,  save_Theta,    save_RotVel,   save_RotAcc, save_CntctForce, save_Torque

  type::part_io_size_vec
    integer,dimension(1)::  sizes
    integer,dimension(1)::  subsizes
    integer,dimension(1)::  starts
  end type part_io_size_vec
  type part_io_size_mat
    integer,dimension(2)::  sizes
    integer,dimension(2)::  subsizes
    integer,dimension(2)::  starts
  end type part_io_size_mat

  type:: Prtcl_IO_Visu
  contains
    procedure:: Init_visu     =>  PIO_Init_visu
    procedure:: Final_visu    =>  PIO_Final_visu
    procedure:: Dump_visu     =>  PIO_Dump_visu
    procedure:: Read_Restart  =>  PIO_Read_Restart
    procedure:: ReadFixedCoord=>  PIO_ReadFixedCoord
    procedure:: RestartCL     =>  PIO_RestartCL
    procedure:: Write_Restart =>  PIO_Write_Restart
    procedure:: Delete_Prev_Restart =>  PIO_Delete_Prev_Restart
    procedure,private:: Write_XDMF  =>  PIO_Write_XDMF
  end type Prtcl_IO_Visu
  type(Prtcl_IO_Visu),public:: DEM_IO

  ! useful interfaces
  interface Prtcl_dump
    module procedure Prtcl_dump_int_one,     Prtcl_dump_int_vector
    module procedure Prtcl_dump_real_one,    Prtcl_dump_real_vector 
    module procedure Prtcl_dump_real3_vector,Prtcl_dump_real3_matrix
  end interface Prtcl_dump

contains

  !**********************************************************************
  ! PIO_Init_visu
  !**********************************************************************
  subroutine PIO_Init_visu(this, chFile)
    implicit none 
    class(Prtcl_IO_Visu)::this
    character(*),intent(in)::chFile

    ! locals
    integer::nUnitFile,myistat,indent,nflds,ifld
    NAMELIST /PrtclVisuOption/ saveXDMFOnce,save_ID,save_Diameter,save_Type,save_UsrMark,save_LinVel,   &
                               save_LinAcc,save_Theta,save_RotVel,save_RotAcc,save_CntctForce,save_Torque
    character(256)::XdmfFile

    nUnitFile = d_read_DEMParam_unit    
    open(unit=nUnitFile, file=chFile,status='old',form='formatted',IOSTAT=myistat)
    if(myistat/=0)call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Init_visu", "Cannot open file: "//trim(chFile))
    read(nUnitFile, nml=PrtclVisuOption)
    write(d_DEMLogInfo_unit, nml=PrtclVisuOption)
    close(nUnitFile,IOSTAT=myistat)

    ! initialize the XDMF/XDF file
    if(nrank/=0) return
    nUnitFile = d_write_xdmf_unit
    XdmfFile = trim(DEM_opt%Res_Dir)//"PartVisuFor"//trim(DEM_opt%RunName)//".xmf"
    open(unit=nUnitFile, file=XdmfFile,status='replace',form='formatted',IOSTAT=myistat)
    if(myistat /= 0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Init_visu","Cannot open file:  "//trim(XdmfFile))
    ! XDMF/XMF Title
    write(nUnitFile,'(A)') '<?xml version="1.0" ?>'
    write(nUnitFile,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(nUnitFile,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(nUnitFile,'(A)') '<Domain>'

    ! Time series
    indent =  4
    nflds = (DEM_Opt%ilast - DEM_Opt%ifirst +1)/DEM_Opt%SaveVisu  + 1
    write(nUnitFile,'(A)')repeat(' ',indent)//'<Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
    indent = indent + 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'<Time TimeType="List">'
    indent = indent + 4
    write(nUnitFile,'(A,I6,A)')repeat(' ',indent)//'<DataItem Format="XML" NumberType="Int" Dimensions="',nflds,'">' 
    write(nUnitFile,'(A)',advance='no') repeat(' ',indent)
    do ifld = 1,nflds
      write(nUnitFile,'(I9)',advance='no')  (ifld-1)*DEM_Opt%SaveVisu + DEM_Opt%ifirst-1
    enddo
    write(nUnitFile,'(A)')repeat(' ',indent)//'</DataItem>'
    indent = indent - 4
    write(nUnitFile,fmt='(A)')repeat(' ',indent)//'</Time>'
    close(nUnitFile,IOSTAT=myistat)
    if( .not. saveXDMFOnce) return

    do ifld = 1,nflds
      call this%Write_XDMF((ifld-1)*DEM_Opt%SaveVisu + DEM_Opt%ifirst-1)
    enddo

    ! XDMF/XMF Tail
    open(unit=nUnitFile, file=XdmfFile,status='old',position='append',form='formatted',IOSTAT=myistat)
    write(nUnitFile,'(A)') '    </Grid>'
    write(nUnitFile,'(A)') '</Domain>'
    write(nUnitFile,'(A)') '</Xdmf>'
    close(nUnitFile,IOSTAT=myistat)

  end subroutine PIO_Init_visu

  !**********************************************************************
  ! PIO_Delete_Prev_Restart
  !**********************************************************************
  subroutine PIO_Delete_Prev_Restart(this,itime)
    implicit none
    class(Prtcl_IO_Visu)::this
    integer:: itime

    ! lcoals
    character(24)::ch
    character(256)::chFile

    if(nrank/=0) return
    write(ch,'(I10.10)')Prev_BackUp_itime   
    chFile = trim(DEM_opt%RestartDir)//"RestartFor"//trim(DEM_opt%RunName)//trim(adjustl(ch))
    call system("rm "//trim(adjustl(chFile))//" 2> /dev/null")
    Prev_BackUp_itime = itime
  end subroutine PIO_Delete_Prev_Restart

  !**********************************************************************
  ! PIO_RestartCL
  !**********************************************************************
  subroutine PIO_RestartCL(this)
    implicit none
    class(Prtcl_IO_Visu)::this

    ! locals
    character(24)::ch
    character(256)::chFile
    integer:: CntctVec(20),intvec(2),i,j,ncv,nlink_total1,int_t
    integer:: itime,fh,code,nlocal,np,tsize,rsize,nreal3,pbyte
    real(RK)::xst,xed,yst,yed,zst,zed,realt
    integer(kind=MPI_OFFSET_KIND)::disp,disp_pos,disp_CL,disp_TanStart,dispt,disp_Tan
    type(real3)::real3t,TanDelVel(20)

    itime = DEM_Opt%ifirst - 1

    ! Begin to write Restart file
    write(ch,'(I10.10)')itime
    chFile = trim(DEM_opt%RestartDir)//"RestartFor"//trim(DEM_opt%RunName)//trim(adjustl(ch))
    call MPI_FILE_OPEN(MPI_COMM_WORLD, chFile, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, code)
    if(code /= 0 .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_RestartCL","Cannot open file: "//trim(chFile))
    endif
    disp = 0_MPI_OFFSET_KIND;

    np = DEM_Opt%np_InDomain
    disp=disp+int_byte  ! firstly skip the np_InDomain
    call MPI_FILE_READ_AT(fh,disp,nlink_total1, 1, int_type, MPI_STATUS_IGNORE, code);disp=disp+int_byte;

    tsize=GPrtcl_list%tsize; rsize=GPrtcl_list%rsize;
    nlocal=0; nreal3 = 2*(1+tsize+rsize)
    pbyte=int_byte*3 + nreal3*real3_byte  ! corresponding to the subroutine 'PIO_Write_Restart'
    xst=DEM_decomp%xSt; xed=DEM_decomp%xEd
    yst=DEM_decomp%ySt; yed=DEM_decomp%yEd
    zst=DEM_decomp%zSt; zed=DEM_decomp%zEd

    disp_pos     = disp + int_byte*3*np
    disp_CL      = disp + pbyte     *np
    disp_TanStart= disp_CL + (2*nlink_total1 + np)*int_byte
    DO i=1,np
      call MPI_FILE_READ_AT(fh,disp_pos,real3t,1,real3_type,MPI_STATUS_IGNORE,code);disp_pos=disp_pos +real3_byte
      IF(real3t%x>=xst .and. real3t%x<xed .and. real3t%y>=yst .and. real3t%y<yed .and. real3t%z>=zst .and. real3t%z<zed) THEN

        nlocal = nlocal + 1
        ! clc ncv(number of particles/walls which have overlap with this particle)
        dispt= disp_CL; ncv=0
        do
          call MPI_FILE_READ_AT(fh,dispt,int_t,1,int_type,MPI_STATUS_IGNORE,code)
          if( int_t==END_FLAG ) then
            dispt = dispt+int_byte; exit
          endif
          ncv=ncv+1; dispt = dispt+int_byte*2
        enddo
        
        do j=1,ncv
          call MPI_FILE_READ_AT(fh,disp_CL,intvec,2,int_type,MPI_STATUS_IGNORE,code);disp_CL=disp_CL+int_byte*2
          disp_Tan = disp_TanStart + real3_byte*(intvec(2)-1)
          call MPI_FILE_READ_AT(fh,disp_Tan,real3t, 1, real3_type, MPI_STATUS_IGNORE, code)
          CntctVec(j) = intvec(1)
          TanDelVel(j)= real3t
        enddo
        disp_CL=disp_CL+int_byte   ! END_FLAG
        if(ncv>0) call GPPW_CntctList%Add_RestartCntctlink(nlocal,ncv,CntctVec,TanDelVel)
      ELSE
        do
          call MPI_FILE_READ_AT(fh,disp_CL,int_t,1,int_type,MPI_STATUS_IGNORE,code)
          if( int_t==END_FLAG ) then
            disp_CL=disp_CL+int_byte; exit
          endif
          ncv=ncv+1; disp_CL=disp_CL+int_byte*2
        enddo
      ENDIF
    ENDDO
    call MPI_FILE_CLOSE(fh,code)

  end subroutine PIO_RestartCL

  !**********************************************************************
  ! PIO_ReadFixedCoord
  !**********************************************************************
  subroutine PIO_ReadFixedCoord(this)
    implicit none
    class(Prtcl_IO_visu)::this
 
    ! locals
    type(real3)::real3t
    character(256)::chFile
    real(RK)::xst,xed,yst,yed,zst,zed,radius,diam
    integer(kind=MPI_OFFSET_KIND):: byte_total,disp_Pos,disp_Diam,disp_Type,disp,filesize
    integer::code,fh,numPrtclFix,i,nfix,nfix_sum,ptype,nfixNew,bgn_ind,color,key,PrtclFix_WORLD,prank
    integer,allocatable,dimension(:):: nP_in_bin, nP_in_bin_reduce,IntVec
    real(RK),allocatable,dimension(:):: RealVec
    type(real3),allocatable,dimension(:):: real3Vec
    type(real4),allocatable,dimension(:):: real4Vec
    type(part_io_size_vec)::pvsize

    chFile = trim(DEM_opt%RestartDir)//"FixedSpheresCoord.dat"
    call MPI_FILE_OPEN(MPI_COMM_WORLD, chFile, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, code)
    if(code /= 0 .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixedCoord","Cannot open file: "//trim(chFile))
    endif

    numPrtclFix = DEM_opt%numPrtclFix
    xst=DEM_decomp%xSt; xed=DEM_decomp%xEd
    yst=DEM_decomp%ySt; yed=DEM_decomp%yEd
    zst=DEM_decomp%zSt; zed=DEM_decomp%zEd

    ! the data storage sequence in file "FixedSpheresCoord.dat" is as follow:
    ! Position(real3 type), Diameter(real type), Prtcl_Type(integer type)
    byte_total= (real3_byte+real_byte+int_byte)*numPrtclFix
    call MPI_FILE_GET_SIZE(fh,filesize,code)
    if(filesize /= byte_total .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixedCoord","file byte wrong")
    endif

    nfix=0
    disp_Pos = 0_MPI_OFFSET_KIND
    disp_Diam= disp_Pos  +numPrtclFix*real3_byte
    disp_Type= disp_Diam +numPrtclFix*real_byte
    allocate(nP_in_bin(DEM_opt%numPrtcl_Type));        nP_in_bin=0
    allocate(nP_in_bin_reduce(DEM_opt%numPrtcl_Type)); nP_in_bin_reduce=0
    DO i=1, numPrtclFix
      call MPI_FILE_READ_AT(fh,disp_Diam,    diam,   1, real_type,  MPI_STATUS_IGNORE, code)
      call MPI_FILE_READ_AT(fh,disp_Type,    ptype,  1, int_type,   MPI_STATUS_IGNORE, code)
      radius = DEMProperty%Prtcl_PureProp(ptype)%Radius
      if( abs(two*radius/diam -one)>1.0E-6 .and. nrank==0 )then
        call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixedCoord","Diameter not coordinate")
      endif

      call MPI_FILE_READ_AT(fh,disp_pos,  real3t, 1, real3_type, MPI_STATUS_IGNORE, code)
      IF(real3t%x>=xst .and. real3t%x<xed .and. real3t%y>=yst .and. real3t%y<yed .and. real3t%z>=zst .and. real3t%z<zed) THEN

        if(nfix>=GPrtcl_list%mlocalFix)  then
          nfixNew= int(1.2_RK*real(nfix,kind=RK))
          nfixNew= max(nfixNew, nfix+1)
          nfixNew= min(nfixNew,DEM_Opt%numPrtclFix)
          GPrtcl_list%mlocalFix= nfixNew

          call move_alloc(GPFix_id,IntVec)
          allocate(GPFix_id(nfixNew))
          GPFix_id(1:nfix)=IntVec
          call move_alloc(GPFix_pType,IntVec)
          allocate(GPFix_pType(nfixNew))
          GPFix_pType(1:nfix)=IntVec
          deallocate(IntVec)
          call move_alloc(GPFIx_PosR,real4Vec)
          allocate(GPFix_PosR(nfixNew))
          GPFix_PosR(1:nfix)=real4Vec
          deallocate(real4Vec)
        endif
        nfix=nfix+1
         
        GPFix_id(nfix)     = i  + DEM_opt%numPrtcl         ! NOTE HERE
        GPFix_pType(nfix)  = PType
        GPFix_PosR(nfix)   = real3t
        GPFix_PosR(nfix)%w = radius
        nP_in_bin(PType)   = nP_in_bin(PType)+1
      ENDIF
      disp_Pos     = disp_Pos     + real3_byte
      disp_Diam    = disp_Diam    + real_byte
      disp_Type    = disp_Type    + int_byte
    ENDDO
    call MPI_ALLREDUCE(nP_in_bin, nP_in_bin_reduce,DEM_opt%numPrtcl_Type,int_type,MPI_SUM,MPI_COMM_WORLD,code)
    DEMProperty%nPrtcl_in_Bin= DEMProperty%nPrtcl_in_Bin+ nP_in_bin_reduce

    if(nfix>0) then
      call move_alloc(GPFix_id,IntVec)
      allocate(GPFix_id(nfix))
      GPFix_id=IntVec(1:nfix)
      call move_alloc(GPFix_pType,IntVec)
      allocate(GPFix_pType(nfix))
      GPFix_pType=IntVec(1:nfix)
      deallocate(IntVec)
      call move_alloc(GPFIx_PosR,real4Vec)
      allocate(GPFix_PosR(nfix))
      GPFix_PosR=real4Vec(1:nfix)
      deallocate(real4Vec)
    else
      deallocate(GPFix_id,GPFix_pType,GPFix_PosR)
    endif
    call MPI_FILE_CLOSE(fh,code)

    GPrtcl_list%mlocalFix = nfix
    call MPI_REDUCE(nfix,nfix_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,code)
    if(nfix_sum/= numPrtclFix .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixedCoord"," nfix_sum/= numPrtclFix " )
    endif
  end subroutine PIO_ReadFixedCoord

  !**********************************************************************
  ! PIO_Read_Restart
  !**********************************************************************
  subroutine PIO_Read_Restart(this)
    implicit none
    class(Prtcl_IO_Visu)::this

    ! locals
    character(24)::ch
    character(256)::chFile
    integer::itime,fh,code,nlocal,np,i,j,k,int_t,itype,tsize,rsize,nlocal_sum,real3VecLen
    real(RK)::xst,xed,yst,yed,zst,zed,realt
    integer(kind=MPI_OFFSET_KIND)::disp,disp_id,disp_type,disp_mark,disp_pos,disp_linvel,disp_linacc
    integer(kind=MPI_OFFSET_KIND)::disp_theta,disp_rotvel,disp_rotacc
    type(real3)::real3t
    type(real3),allocatable,dimension(:)::real3Vec
    
    itime = DEM_Opt%ifirst - 1  

    ! Begin to write Restart file
    write(ch,'(I10.10)')itime
    chFile = trim(DEM_opt%RestartDir)//"RestartFor"//trim(DEM_opt%RunName)//trim(adjustl(ch))
    call MPI_FILE_OPEN(MPI_COMM_WORLD, chFile, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, code)
    if(code /= 0 .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart","Cannot open file: "//trim(chFile))
    endif
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_READ_AT(fh,disp,np, 1, int_type, MPI_STATUS_IGNORE, code);disp=disp+int_byte;
    if(np>DEM_Opt%numPrtcl .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart"," np_InDomain > numPrtcl " )
    endif
    DEM_Opt%np_InDomain = np
    disp=disp+int_byte   ! skip nlink_total1

    tsize=GPrtcl_list%tsize; rsize=GPrtcl_list%rsize;
    real3VecLen= max(tsize,rsize)
    allocate( real3Vec(real3VecLen) )
    xst=DEM_decomp%xSt; xed=DEM_decomp%xEd
    yst=DEM_decomp%ySt; yed=DEM_decomp%yEd
    zst=DEM_decomp%zSt; zed=DEM_decomp%zEd

    nlocal=0;
    disp_id    = disp
    disp_type  = disp_id    + int_byte  *np
    disp_mark  = disp_type  + int_byte  *np
    disp_pos   = disp_mark  + int_byte  *np
    disp_linvel= disp_pos   + real3_byte*np
    disp_linacc= disp_linvel+ real3_byte*np*tsize
    disp_theta = disp_linacc+ real3_byte*np*tsize
    disp_rotvel= disp_theta + real3_byte*np
    disp_rotacc= disp_rotvel+ real3_byte*np*rsize
    DO i=1,np
      call MPI_FILE_READ_AT(fh,disp_pos,real3t, 1, real3_type, MPI_STATUS_IGNORE, code)
      IF(real3t%x>=xst .and. real3t%x<xed .and. real3t%y>=yst .and. real3t%y<yed .and. real3t%z>=zst .and. real3t%z<zed) THEN

        if(nlocal>=GPrtcl_list%mlocal)  call GPrtcl_list%ReallocatePrtclVar(nlocal)
        nlocal=nlocal+1

        ! id
        call MPI_FILE_READ_AT(fh,disp_id,  int_t, 1, int_type, MPI_STATUS_IGNORE, code) 
        GPrtcl_id(nlocal)=int_t

        ! pType
        call MPI_FILE_READ_AT(fh,disp_type,itype, 1, int_type, MPI_STATUS_IGNORE, code) 
        GPrtcl_pType(nlocal)=itype

        ! Usr_Mark
        call MPI_FILE_READ_AT(fh,disp_mark,int_t, 1, int_type, MPI_STATUS_IGNORE, code) 
        GPrtcl_UsrMark(nlocal)=int_t

        ! PosR
        GPrtcl_PosR(nlocal)= real3t
        GPrtcl_PosR(nlocal)%w=DEMProperty%Prtcl_PureProp(itype)%Radius

        ! LinVec
        call MPI_FILE_READ_AT(fh,disp_linVel, real3Vec, tsize, real3_type, MPI_STATUS_IGNORE, code) 
        GPrtcl_LinVel(1:tsize,nlocal) = real3Vec(1:tsize)

        ! LinAcc
        call MPI_FILE_READ_AT(fh,disp_linAcc, real3Vec, tsize, real3_type, MPI_STATUS_IGNORE, code) 
        GPrtcl_LinAcc(1:tsize,nlocal) = real3Vec(1:tsize)

        ! theta
        call MPI_FILE_READ_AT(fh,disp_theta, real3t, 1, real3_type, MPI_STATUS_IGNORE, code) 
        GPrtcl_theta(nlocal) = real3t       

        ! RotVel
        call MPI_FILE_READ_AT(fh,disp_rotVel, real3Vec, rsize, real3_type, MPI_STATUS_IGNORE, code) 
        GPrtcl_RotVel(1:rsize,nlocal) = real3Vec(1:rsize)

        ! rotAcc
        call MPI_FILE_READ_AT(fh,disp_rotAcc, real3Vec, rsize, real3_type, MPI_STATUS_IGNORE, code) 
        GPrtcl_RotAcc(1:rsize,nlocal) = real3Vec(1:rsize)

      ENDIF

      disp_id    = disp_id    + int_byte
      disp_type  = disp_type  + int_byte
      disp_mark  = disp_mark  + int_byte
      disp_pos   = disp_pos   + real3_byte
      disp_linvel= disp_linvel+ real3_byte*tsize
      disp_linacc= disp_linacc+ real3_byte*tsize
      disp_theta = disp_theta + real3_byte
      disp_rotvel= disp_rotvel+ real3_byte*rsize
      disp_rotacc= disp_rotacc+ real3_byte*rsize
    ENDDO

    GPrtcl_list%nlocal = nlocal
    call MPI_REDUCE(nlocal,nlocal_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,code)
    if(nlocal_sum/= np .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart"," nlocal_sum/= np_InDomain " )
    endif

    call MPI_FILE_CLOSE(fh,code)
  end subroutine PIO_Read_Restart

  !**********************************************************************
  ! PIO_Write_Restart
  !**********************************************************************
  subroutine PIO_Write_Restart(this,itime)
    implicit none 
    class(Prtcl_IO_Visu)::this
    integer,intent(in)::itime
   
    ! locals 
    character(24)::ch
    character(256)::chFile
    integer :: i,j,k,nlocal,bgn_ind,write_num,prank
    integer :: nlinkSum1,nlink_total1,nlink_ind1,nlinkSum2,nlink_ind2
    integer :: code,fh,tsize,rsize,color,key,Prtcl_WORLD1,Prtcl_WORLD2,ncv,prev,now,idv1,idv2
    integer(kind=MPI_OFFSET_KIND)::disp,bgn_byte,filesize
    type(real3),allocatable,dimension(:)::real3Vec
    integer,allocatable,dimension(:)::IntVec
    type(real3)::TanDel
    integer,dimension(40)::CntctVec
    type(part_io_size_vec)::pvsize
    type(part_io_size_mat)::pmsize

    ! Calculate the bgn_ind and nlink_ind
    nlocal    = GPrtcl_list%nlocal
    bgn_ind   = clc_bgn_ind(nlocal)

    nlinkSum1 = GPPW_CntctList%Get_numCntcts()
    nlinkSum2 = GPPW_CntctList%Get_numTanDel()
    nlink_ind1= clc_bgn_ind(nlinkSum1)
    nlink_ind2= clc_bgn_ind(nlinkSum2)
    call MPI_ALLREDUCE(nlinkSum1, nlink_total1, 1, int_type,MPI_SUM, MPI_COMM_WORLD,code )
    call GPPW_CntctList%Prepare_Restart(nlink_ind2)

    ! Create the Prtcl_GROUP
    color = 1; key=nrank
    if(nlocal<=0) color=2
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,Prtcl_WORLD1,code)
    if(color==2) return

    ! Begin to write Restart file
    write(ch,'(I10.10)')itime
    chFile = trim(DEM_opt%RestartDir)//"RestartFor"//trim(DEM_opt%RunName)//trim(adjustl(ch))
    call MPI_FILE_OPEN(Prtcl_WORLD1, chFile, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, code)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,code)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND

    ! Write DEM_Opt%np_InDomain,nlink_total1, in the begining of the Restart file
    call MPI_COMM_RANK(PRTCL_WORLD1, prank, code)
    if(prank==0) then
       write_num  =  1
    else
       write_num  =  0
    endif
    call Prtcl_dump(fh,disp, DEM_Opt%np_InDomain, write_num)
    call Prtcl_dump(fh,disp, nlink_total1,        write_num)

    ! Beigin to write
    pvsize%sizes(1)     = DEM_Opt%np_InDomain
    pvsize%subsizes(1)  = nlocal
    pvsize%starts(1)    = bgn_ind
    call Prtcl_dump(fh,disp, GPrtcl_id(1:nlocal),      pvsize)
    call Prtcl_dump(fh,disp, GPrtcl_pType(1:nlocal),   pvsize)
    call Prtcl_dump(fh,disp, GPrtcl_UsrMark(1:nlocal), pvsize)
    allocate(real3Vec(nlocal))
    do i=1,nlocal
      real3Vec(i)=GPrtcl_PosR(i)
    enddo
    call Prtcl_dump(fh,disp, real3Vec(1:nlocal),  pvsize)
    deallocate(real3Vec)

    tsize=GPrtcl_list%tsize;       rsize=GPrtcl_list%rsize;
    pmsize%sizes(1)     = tsize;   pmsize%sizes(2)     = DEM_Opt%np_InDomain
    pmsize%subsizes(1)  = tsize;   pmsize%subsizes(2)  = nlocal
    pmsize%starts(1)    = 0 ;      pmsize%starts(2)    = bgn_ind  
    call Prtcl_dump(fh,disp, GPrtcl_LinVel(1:tsize,1:nlocal),  pmsize)
    call Prtcl_dump(fh,disp, GPrtcl_LinAcc(1:tsize,1:nlocal),  pmsize)

    pmsize%sizes(1)     = rsize;   pmsize%sizes(2)     = DEM_Opt%np_InDomain
    pmsize%subsizes(1)  = rsize;   pmsize%subsizes(2)  = nlocal
    pmsize%starts(1)    = 0 ;      pmsize%starts(2)    = bgn_ind  
    call Prtcl_dump(fh,disp, GPrtcl_theta(1:nlocal),  pvsize) 
    call Prtcl_dump(fh,disp, GPrtcl_RotVel(1:rsize,1:nlocal),  pmsize)
    call Prtcl_dump(fh,disp, GPrtcl_RotAcc(1:rsize,1:nlocal),  pmsize)

    pvsize%sizes(1)     = 2*nlink_total1 + DEM_Opt%np_InDomain
    pvsize%subsizes(1)  = 2*nlinkSum1    + nlocal
    pvsize%starts(1)    = 2*nlink_ind1   + bgn_ind
    allocate(IntVec(pvsize%subsizes(1)) )
    idv1=0; idv2=0
    DO i=1,nlocal
      call GPPW_CntctList%Gather_Cntctlink_Restart(i,CntctVec,ncv)
      if(ncv>0) then
        idv2=idv1+2*ncv;
        IntVec(idv1+1:idv2)=CntctVec(1:2*ncv)
        idv1=idv2  
      endif      
      idv1=idv1+1
      IntVec(idv1)=END_FLAG
    ENDDO
    call Prtcl_dump(fh,disp, IntVec,  pvsize)
    deallocate(IntVec)
    call MPI_FILE_CLOSE(fh, code)

    ! Begin to write Contact List file
    color = 1; key=nrank
    if(nlinkSum2<=0) color=2
    call MPI_COMM_SPLIT(Prtcl_WORLD1,color,key,Prtcl_WORLD2,code)
    call MPI_COMM_FREE( PRTCL_WORLD1, code)
    if(color==2) return

    prev=1; j=0; k=0
    disp = disp + nlink_ind2 * real3_byte
    allocate(real3Vec(nlocal))
    call MPI_FILE_OPEN(Prtcl_WORLD2, chFile, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, code)
    DO
      call GPPW_CntctList%GetNext_TanDel(TanDel,prev,now);  prev=now+1
      j=j+1; k=k+1; real3Vec(j)=TanDel
      if(k==nlinkSum2) then
        call MPI_FILE_WRITE_AT(fh,disp, real3Vec, j, real3_type, MPI_STATUS_IGNORE, code)
        disp = disp +j*real3_byte; exit
      endif
      if(j==nlocal) then
        call MPI_FILE_WRITE_AT(fh,disp, real3Vec, j, real3_type, MPI_STATUS_IGNORE, code)
        disp = disp +j*real3_byte; j=0 
      endif
    ENDDO
    deallocate(real3Vec)

    call MPI_BARRIER( PRTCL_WORLD2, code)
    call MPI_FILE_CLOSE(fh, code)
    call MPI_COMM_FREE( PRTCL_WORLD2, code)

  end subroutine PIO_Write_Restart

  !**********************************************************************
  ! Purpose:
  !   Create a xdmf/xmf file in order to view the simulation results
  !     by Paraview directly
  ! 
  ! Original Author: 
  !	  Pedro Costa
  ! 
  ! Modified by:
  !	  Zheng Gong
  ! 
  ! Original Source file is downloaded from ( April 2020 ):
  !   https://github.com/p-costa/gen_xdmf_particles
  !              
  !**********************************************************************
  subroutine PIO_Write_XDMF(this,itime) 
    implicit none 
    class(Prtcl_IO_Visu)::this
    integer,intent(in)::itime

    !locals
    integer(kind=MPI_OFFSET_KIND)::disp
    integer:: indent,nUnitFile,myistat,np,dims,iprec
    character(256)::chFile

    if(nrank/=0) return 
    np=DEM_Opt%np_InDomain
    nUnitFile = d_write_xdmf_unit
    chFile = trim(DEM_opt%Res_Dir)//"PartVisuFor"//trim(DEM_opt%RunName)//".xmf"
    open(unit=nUnitFile, file=chFile,status='old',position='append',form='formatted',IOSTAT=myistat)
    if(myistat/=0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Write_XDMF","Cannot open file: "//trim(chFile))
    endif

    indent = 8; disp = 0_MPI_OFFSET_KIND
    chFile = "PartVisuFor"//trim(DEM_opt%RunName)
    dims=3; iprec=RK
    write(nUnitFile,'(A,I10.10,A)')repeat(' ',indent)//'<Grid Name="T',itime,'" GridType="Uniform">'
    indent = indent + 4
    write(nUnitFile,'(A,I9,A)')repeat(' ',indent)//'<Topology TopologyType="Polyvertex" NodesPerElement="',np,'"/>'
    write(nUnitFile,'(A)')repeat(' ',indent)//'<Geometry GeometryType="'//"XYZ"//'">'
    indent = indent + 4
    write(nUnitFile,'(A,I1,A,I2,I9,A,I15,A)')repeat(' ',indent)// '<DataItem Format="Binary"' // &
          ' DataType="Float" Precision="',iprec,'" Endian="Native"' // &
          ' Dimensions="',dims,np,'" Seek="',disp,'">'
    disp = disp+np*dims*iprec
    indent = indent + 4
    write(nUnitFile,'(A,I10.10)')repeat(' ',indent)//trim(chFile),itime
    indent = indent - 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'</DataItem>'
    indent = indent - 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'</Geometry>'

    IF(save_ID) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"ID","Scalar","Int",disp)
    ENDIF
    IF(save_Diameter) THEN
      dims=1; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"Diameter","Scalar","Float",disp)
    ENDIF
    IF(save_Type) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"Type","Scalar","Int",disp)
    ENDIF
    IF(save_UsrMark) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"UsrMark","Scalar","Int",disp)
    ENDIF
    IF(save_LinVel) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"LinVel","Vector","Float",disp)
    ENDIF
    IF(save_LinAcc) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"LinAcc","Vector","Float",disp)
    ENDIF
    IF(save_Theta) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"Theta","Vector","Float",disp)
    ENDIF
    IF(save_RotVel) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"RotVel","Vector","Float",disp)
    ENDIF
    IF(save_RotAcc) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"RotAcc","Vector","Float",disp)
    ENDIF
    IF(save_CntctForce) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"CntctForce","Vector","Float",disp)
    ENDIF
    IF(save_Torque) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"Torque","Vector","Float",disp)
    ENDIF
    
    write(nUnitFile,'(A)')'        </Grid>'
    close(nUnitFile,IOSTAT=myistat)
 
  end subroutine PIO_Write_XDMF

  !**********************************************************************
  ! Write_XDMF_One
  !**********************************************************************
  subroutine Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,chName,chAttribute,chDataType,disp)
    implicit none
    integer,intent(in)::nUnitFile,dims,iprec,np,itime
    character(*),intent(in)::chFile,chName,chAttribute,chDataType
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
    
    ! locals
    integer:: indent
    indent = 12

    write(nUnitFile,'(A)')repeat(' ',indent)//'<Attribute Type="'//trim(chAttribute)//'" Center="Node" Name="'//trim(chName)//'">'
    indent = indent + 4
    write(nUnitFile,'(3A,I1,A,I2,I9,A,I15,A)')repeat(' ',indent)// '<DataItem Format="Binary"' // &
          ' DataType="',trim(chDataType),'" Precision="',iprec,'" Endian="Native"' // &
          ' Dimensions="',dims,np,'" Seek="',disp,'">'
    disp = disp+np*dims*iprec
    indent = indent + 4
    write(nUnitFile,'(A,I10.10)')repeat(' ',indent)//trim(chFile),itime
    indent = indent - 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'</DataItem>'
    indent = indent - 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'</Attribute>'

  end subroutine Write_XDMF_One

  !**********************************************************************
  ! PIO_Final_visu
  !**********************************************************************
  subroutine PIO_Final_visu(this)
    implicit none 
    class(Prtcl_IO_Visu)::this

    ! locals
    integer::nUnitFile,myistat
    character(256)::XdmfFile

    if(nrank/=0 .or. saveXDMFOnce) return
    nUnitFile = d_write_xdmf_unit
    XdmfFile = trim(DEM_opt%Res_Dir)//"PartVisuFor"//trim(DEM_opt%RunName)//".xmf"
    open(unit=nUnitFile, file=XdmfFile,status='old',position='append',form='formatted',IOSTAT=myistat)
    if(myistat /= 0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Final_visu","Cannot open file: "//trim(XdmfFile))
    endif
    ! XDMF/XMF Tail
    write(nUnitFile,'(A)') '    </Grid>'
    write(nUnitFile,'(A)') '</Domain>'
    write(nUnitFile,'(A)') '</Xdmf>'
    close(nUnitFile)

  end subroutine PIO_Final_visu

  !**********************************************************************
  ! PIO_Dump_visu
  !**********************************************************************
  subroutine PIO_Dump_visu(this, itime)
    implicit none
    class(Prtcl_IO_Visu)::this
    integer,intent(in)::itime

    ! locals
    character(24)::ch
    character(256)::chFile
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer :: code,fh,i,color,key,Prtcl_WORLD,nlocal,bgn_ind
    type(part_io_size_vec)::pvsize
    type(real3),allocatable,dimension(:)::real3Vec
    integer,allocatable,dimension(:)::intVec
    real(RK),allocatable,dimension(:)::realVec

    ! write xdmf file first
    if(.not.saveXDMFOnce) call this%Write_XDMF(itime)

    ! update the bgn_ind
    nlocal = GPrtcl_list%nlocal
    bgn_ind=clc_bgn_ind(nlocal)

    ! create the Prtcl_GROUP
    color = 1; key=nrank
    if(nlocal<=0) color=2
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,Prtcl_WORLD,code)
    if(color==2) return
    
    ! begin to dump
    write(ch,'(I10.10)')itime
    chFile = trim(DEM_opt%Res_Dir)//"PartVisuFor"//trim(DEM_opt%RunName)//trim(adjustl(ch))
    call MPI_FILE_OPEN(Prtcl_WORLD, chFile, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, code)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,code)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    pvsize%sizes(1)     = DEM_Opt%np_InDomain
    pvsize%subsizes(1)  = nlocal
    pvsize%starts(1)    = bgn_ind
    if(nlocal<=0) return

    allocate(real3Vec(nlocal))
    do i=1,nlocal
      real3Vec(i)=GPrtcl_PosR(i)
    enddo
    call Prtcl_dump(fh,disp, real3Vec(1:nlocal),  pvsize)
    deallocate(real3Vec)
    if(save_ID) call Prtcl_dump(fh,disp, GPrtcl_id(1:nlocal),  pvsize)
    if(save_Diameter) then
      allocate(realVec(nlocal))
      do i=1,nlocal
        realVec(i)= two*GPrtcl_PosR(i)%w
      enddo   
      call Prtcl_dump(fh,disp, realVec(1:nlocal),  pvsize)
      deallocate(realVec)
    endif
    if(save_Type)     call Prtcl_dump(fh,disp, GPrtcl_pType(1:nlocal),     pvsize)

    if(save_UsrMark)  then
      allocate(intvec(nlocal))
      do i=1,nlocal
        intvec(i)=GPrtcl_UsrMark(i)
      enddo
      call Prtcl_dump(fh,disp, intvec,   pvsize)
      deallocate(intvec)
    endif

    if(save_LinVel)     call Prtcl_dump(fh,disp, GPrtcl_LinVel(1,1:nlocal),  pvsize)
    if(save_LinAcc)     call Prtcl_dump(fh,disp, GPrtcl_LinAcc(1,1:nlocal),  pvsize)
    if(save_Theta)      call Prtcl_dump(fh,disp, GPrtcl_Theta(1:nlocal),     pvsize)
    if(save_RotVel)     call Prtcl_dump(fh,disp, GPrtcl_RotVel(1,1:nlocal),  pvsize)
    if(save_RotAcc)     call Prtcl_dump(fh,disp, GPrtcl_RotAcc(1,1:nlocal),  pvsize)
    if(save_CntctForce) call Prtcl_dump(fh,disp, GPrtcl_CntctForce(1:nlocal),pvsize)
    if(save_Torque)     call Prtcl_dump(fh,disp, GPrtcl_Torque(1:nlocal),    pvsize)
    call MPI_FILE_CLOSE(fh, code) 
    call MPI_COMM_FREE( PRTCL_WORLD, code)

  end subroutine PIO_Dump_visu

  !**********************************************************************
  ! clc_bgn_ind
  !**********************************************************************
  function clc_bgn_ind(nlocal) result(bgn_ind)
    implicit none
    integer,intent(in)::nlocal
    integer::bgn_ind

    ! locals
    integer::code
    integer,dimension(MPI_STATUS_SIZE) :: status1
  
    bgn_ind=0
    if(nproc<=1) return
    IF(nrank==0)THEN
      call MPI_SEND(bgn_ind+nlocal, 1, int_type, nrank+1, 0, MPI_COMM_WORLD, code)
    ELSEIF(nrank /= nproc-1) THEN
      call MPI_RECV(bgn_ind,        1, int_type, nrank-1, 0, MPI_COMM_WORLD, status1,  code)
      call MPI_SEND(bgn_ind+nlocal, 1, int_type, nrank+1, 0, MPI_COMM_WORLD, code)
    ELSE
      call MPI_RECV(bgn_ind,        1, int_type, nrank-1, 0, MPI_COMM_WORLD, status1,  code)
    ENDIF 
  end function clc_bgn_ind

  !**********************************************************************
  ! Prtcl_dump_int_one
  !**********************************************************************
  subroutine Prtcl_dump_int_one(fh,disp,var,write_num)
    implicit none
    integer,intent(in)::fh  
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
    integer,intent(in)::var
    integer,intent(in)::write_num    

    ! locals
    integer :: code
  
    call MPI_FILE_SET_VIEW(fh, disp, int_type,  int_type,'native',MPI_INFO_NULL,code)
    call MPI_FILE_WRITE_ALL(fh, var, write_num, int_type, MPI_STATUS_IGNORE, code)
    disp = disp + int_byte
  end subroutine Prtcl_dump_int_one

  !**********************************************************************
  ! Prtcl_dump_int_vector
  !**********************************************************************
  subroutine Prtcl_dump_int_vector(fh,disp,var,pvsize)
    implicit none
    integer,intent(in)::fh
    type(part_io_size_vec),intent(in)::pvsize  
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
    integer,dimension(1:pvsize%subsizes(1)),intent(in)::var

    ! locals
    integer:: code,newtype
    integer,dimension(1) :: sizes, subsizes, starts

    ! calculate sizes, subsizes and starts
    sizes    = pvsize%sizes
    subsizes = pvsize%subsizes
    starts   = pvsize%starts

    ! write the particle revelant integer vector
    call MPI_TYPE_CREATE_SUBARRAY(1, sizes, subsizes, starts, MPI_ORDER_FORTRAN, int_type, newtype, code)
    call MPI_TYPE_COMMIT(newtype,code)
    call MPI_FILE_SET_VIEW(fh,disp,int_type, newtype,'native',MPI_INFO_NULL,code)
    call MPI_FILE_WRITE_ALL(fh, var, subsizes(1),int_type, MPI_STATUS_IGNORE, code)
    call MPI_TYPE_FREE(newtype,code)
    disp = disp + sizes(1) * int_byte
  end subroutine Prtcl_dump_int_vector

  !**********************************************************************
  ! Prtcl_dump_real_one
  !**********************************************************************
  subroutine Prtcl_dump_real_one(fh,disp,var,write_num)
    implicit none
    integer,intent(in)::fh  
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
    real(RK),intent(in)::var
    integer,intent(in)::write_num
  
    ! lcoals
    integer :: code

    call MPI_FILE_SET_VIEW(fh, disp,real_type, real_type,'native',MPI_INFO_NULL,code)
    call MPI_FILE_WRITE_ALL(fh,var, write_num, real_type, MPI_STATUS_IGNORE, code)
    disp = disp + real_byte
  end subroutine Prtcl_dump_real_one

  !**********************************************************************
  ! Prtcl_dump_real_vector
  !**********************************************************************
  subroutine Prtcl_dump_real_vector(fh,disp,var,pvsize)
    implicit none
    integer,intent(in)::fh
    type(part_io_size_vec),intent(in)::pvsize  
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
    real(RK),dimension(1:pvsize%subsizes(1)),intent(in)::var

    ! locals
    integer:: code,newtype
    integer,dimension(1) :: sizes, subsizes, starts

    ! calculate sizes, subsizes and starts
    sizes    = pvsize%sizes
    subsizes = pvsize%subsizes
    starts   = pvsize%starts

    ! write the particle revelant real vector
    call MPI_TYPE_CREATE_SUBARRAY(1, sizes, subsizes, starts, MPI_ORDER_FORTRAN, real_type, newtype, code)
    call MPI_TYPE_COMMIT(newtype,code)
    call MPI_FILE_SET_VIEW(fh,disp,real_type, newtype,'native',MPI_INFO_NULL,code)
    call MPI_FILE_WRITE_ALL(fh, var, subsizes(1),real_type, MPI_STATUS_IGNORE, code)
    call MPI_TYPE_FREE(newtype,code)
    disp = disp + sizes(1) * real_byte
  end subroutine Prtcl_dump_real_vector

  !**********************************************************************
  ! Prtcl_dump_real3_vector
  !**********************************************************************
  subroutine Prtcl_dump_real3_vector(fh,disp,var,pvsize)
    implicit none
    integer,intent(in)::fh
    type(part_io_size_vec),intent(in)::pvsize  
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
    type(real3),dimension(1:pvsize%subsizes(1)),intent(in)::var

    ! locals
    integer:: code,newtype
    integer,dimension(1) :: sizes, subsizes, starts

    ! calculate sizes, subsizes and starts
    sizes    = pvsize%sizes
    subsizes = pvsize%subsizes
    starts   = pvsize%starts

    ! write the particle revelant real3 vector
    call MPI_TYPE_CREATE_SUBARRAY(1, sizes, subsizes, starts, MPI_ORDER_FORTRAN, real3_type, newtype, code)
    call MPI_TYPE_COMMIT(newtype,code)
    call MPI_FILE_SET_VIEW(fh,disp,real3_type, newtype,'native',MPI_INFO_NULL,code)
    call MPI_FILE_WRITE_ALL(fh, var, subsizes(1),real3_type, MPI_STATUS_IGNORE, code)
    call MPI_TYPE_FREE(newtype,code)
    disp = disp + sizes(1) * real3_byte
  end subroutine Prtcl_dump_real3_vector

  !**********************************************************************
  ! Prtcl_dump_real3_matrix
  !**********************************************************************
  subroutine Prtcl_dump_real3_matrix(fh,disp,var,pmsize)
    implicit none
    integer,intent(in)::fh
    type(part_io_size_mat),intent(in)::pmsize    
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
    type(real3),dimension(1:pmsize%subsizes(1),1:pmsize%subsizes(2)),intent(in)::var

    ! locals
    integer :: code,newtype
    integer, dimension(2) :: sizes, subsizes, starts

    ! calculate sizes, subsizes and starts
    sizes     = pmsize%sizes
    subsizes  = pmsize%subsizes
    starts    = pmsize%starts

    ! write the particle relevant real matrix
    call MPI_TYPE_CREATE_SUBARRAY(2, sizes, subsizes, starts, MPI_ORDER_FORTRAN, real3_type, newtype, code)
    call MPI_TYPE_COMMIT(newtype,code)
    call MPI_FILE_SET_VIEW(fh,disp,real3_type, newtype,'native',MPI_INFO_NULL,code)
    call MPI_FILE_WRITE_ALL(fh, var, subsizes(1)*subsizes(2),real3_type, MPI_STATUS_IGNORE, code)
    call MPI_TYPE_FREE(newtype,code)
    disp = disp + sizes(1) * sizes(2) * real3_byte
  end subroutine Prtcl_dump_real3_matrix

end module Prtcl_IOAndVisu
