module Prtcl_CL_and_CF
  use Prtcl_TypeDef
  use Prtcl_logInfo
  use Prtcl_decomp_2d,only: nrank
  use Prtcl_Property
  use Prtcl_Geometry
  use Prtcl_Variables
  use Prtcl_Parameters
  implicit none
  private
    
  integer,dimension(:),allocatable:: Bucket
  integer,dimension(:),allocatable:: id_j
  integer,dimension(:),allocatable:: next
  integer,dimension(:),allocatable:: CntctStatus
  type(real3),dimension(:),allocatable:: tang_del

  integer,dimension(:),allocatable:: id_i
  integer,dimension(:),allocatable:: Head_Cp  ! cp: counterpart
  integer,dimension(:),allocatable:: next_Cp
    
  type ContactList
    integer:: mBucket
    integer:: numCntcts(4)
    integer:: max_numCntcts
    integer:: nextInsert
  contains
    procedure:: InitContactList => CL_InitContactList
    procedure:: reallocateCL    => CL_reallocateCL
    procedure:: AddContactPP    => CL_AddContactPP
    procedure:: AddContactPPG   => CL_AddContactPPG
    procedure:: AddContactPPFix => CL_AddContactPPFix
    procedure:: AddContactPW    => CL_AddContactPW
    procedure:: Find_Insert
    procedure:: PreIteration    => CL_PreIteration
    procedure:: RemvReleased    => CL_RemvReleased
    procedure:: copy            => CL_copy
    procedure:: IsCntct         => CL_IsCntct
    procedure:: getPrtcl_nlink  => CL_getPrtcl_nlink
    procedure:: Gather_Cntctlink=> CL_Gather_Cntctlink
    procedure:: Add_Cntctlink   => CL_Add_Cntctlink
    procedure:: printCL         => CL_printCL
  
    procedure:: Get_numCntcts   => CL_Get_numCntcts
    procedure:: Get_numTanDel   => CL_Get_numTanDel
    procedure:: GetNext_TanDel  => CL_GetNext_TanDel
    procedure:: Prepare_Restart => CL_Prepare_Restart
    procedure:: Add_RestartCntctlink     => CL_Add_RestartCntctlink 
    procedure:: Gather_Cntctlink_Restart => CL_Gather_Cntctlink_Restart
  end type ContactList
  type(ContactList),public :: GPPW_CntctList 
    
contains
    
  !**********************************************************************
  ! Initializing the contact list for ContactList class
  !**********************************************************************
  subroutine CL_InitContactList(this)
    implicit none
    class(ContactList) :: this
    integer:: max_numCntcts,i,iErrSum
    integer:: iErr1,iErr2,iErr3,iErr4,iErr5,iErr6,iErr7,iErr8
    
    this%numCntcts = 0
    max_numCntcts = DEM_opt%cntctList_Size * GPrtcl_list%mlocal

    this%Max_numCntcts = max_numCntcts
    this%mBucket = GPrtcl_list%mlocal
        
    ! initializing the linked list to stores id pairs
    allocate(Bucket(this%mBucket),    Stat=iErr1)
    allocate(id_j(max_numCntcts),     Stat=iErr2)
    allocate(CntctStatus(max_numCntcts), Stat=iErr3)
    allocate(Next(max_numCntcts),     Stat=iErr4)
    allocate(Tang_del(max_numCntcts), Stat=iErr5)

    allocate(id_i(max_numCntcts),    Stat=iErr6)
    allocate(Head_Cp(this%mBucket),  Stat=iErr7)
    allocate(Next_cp(max_numCntcts), Stat=ierr8)
    iErrSum=abs(iErr1)+abs(iErr2)+abs(iErr3)+abs(iErr4)+abs(iErr5)+abs(iErr6)+abs(iErr7)+abs(iErr8)
    if(iErrSum/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"CL_InitContactList","Allocation failed ")
        
    Bucket = 0
    do i = 1, max_numCntcts
      Next(i)=-i-1
      CntctStatus(i)=-1
    end do
    Next(max_numCntcts) = 0
    this%nextInsert = 1

    Head_Cp = -1
    Next_cp = -1
   
  end subroutine CL_InitContactList

  !**********************************************************************
  ! reallocate contact list
  !**********************************************************************
  subroutine CL_reallocateCL(this,nCL_new)
    implicit none
    class(ContactList):: this
    integer,intent(in):: nCL_new

    ! locals
    integer::sizep,sizen
    integer,dimension(:),allocatable:: IntVec

    sizep= this%mBucket
    sizen= int(1.2_RK*real(sizep,kind=RK))
    sizen= max(sizen, nCL_new+1)
    sizen= min(sizen,DEM_Opt%numPrtcl)
    this%mBucket = sizen

    ! ======= integer vector part =======
    call move_alloc(Bucket, IntVec)
    allocate(Bucket(sizen))
    Bucket(1:sizep)=IntVec
    Bucket(sizeP+1:sizen)=0    ! Added at 11:24, 2020-09-11, Gong Zheng

    call move_alloc(Head_Cp,IntVec)
    allocate(Head_Cp(sizen))
    Head_Cp(1:sizep)=IntVec
    Head_Cp(sizep+1:sizen)=-1  ! Added at 11:24, 2020-09-11, Gong Zheng
    deallocate(IntVec)

   ! here nothing is done for id_j, next, CntctStatus, tang_del
   ! considering that Max_numCntcts is big enough.
   ! If not, a fatal error will occur in CL_AddContactPP/CL_AddContactPW/CL_AddContactPPG/CL_AddContactPPFix

  end subroutine CL_reallocateCL

  !**********************************************************************
  ! Adding a contact pair to the contact list (particle & particle) 
  !**********************************************************************
  subroutine CL_AddContactPP(this,id1,id2,ovrlp)
    implicit none
    class(ContactList):: this
    integer,intent(in):: id1,id2
    real(RK),intent(in)::ovrlp
    integer::item1,item2,gid2,nextI

    integer::i,j,k

    gid2=GPrtcl_id(id2)
    call this%Find_Insert(id1,gid2,item1,item2)
    
    ! item1 is the status of the item insertion. 1:old, 2:new 
    ! item2 is the container index
    if(item1>0) then
      CntctStatus(item2)=item1
      this%numCntcts(1) = this%numCntcts(1) + 1
      call ContactForce_PP(id1,id2,item2,ovrlp)

      id_i(item2)= GPrtcl_id(id1)
      Next_Cp(item2)=Head_Cp(id2)
      Head_Cp(id2)=item2

    elseif(item1 == -1 ) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPP","The inserted item's id is greater than allowed value:"//num2str(id1))
    else
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPP","The container is full and there is no space for new item" )      
    endif
  end subroutine CL_AddContactPP

  !**********************************************************************
  ! Adding a contact pair to the contact list (particle & ghost particle) 
  !**********************************************************************
  subroutine CL_AddContactPPG(this,id1,id2,ovrlp)
    implicit none
    class(ContactList):: this
    integer,intent(in):: id1,id2
    real(RK),intent(in)::ovrlp
    integer::item1,item2,gid2

    integer::i,j,k

    gid2=GhostP_id(id2)
    call this%Find_Insert(id1,gid2,item1,item2)
    
    ! item1 is the status of the item insertion. 1:old, 2:new 
    ! item2 is the container index
    if(item1>0) then
      CntctStatus(item2)=item1
      this%numCntcts(2) = this%numCntcts(2) + 1
      call ContactForce_PPG(id1,id2,item2,ovrlp)

    elseif(item1 == -1 ) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPPG","The inserted item's id is greater than allowed value:"//num2str(id1))
    else
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPPG","The container is full and there is no space for new item" )      
    endif
  end subroutine CL_AddContactPPG

  !**********************************************************************
  ! Adding a contact pair to the contact list (particle & fixed particle) 
  !**********************************************************************
  subroutine CL_AddContactPPFix(this,id1,id2,ovrlp)
    implicit none
    class(ContactList):: this
    integer,intent(in):: id1,id2
    real(RK),intent(in)::ovrlp
    integer::item1,item2,gid2

    integer::i,j,k

    gid2=GPFix_id(id2)
    call this%Find_Insert(id1,gid2,item1,item2)
    
    ! item1 is the status of the item insertion. 1:old, 2:new 
    ! item2 is the container index
    if(item1>0) then
      CntctStatus(item2)=item1
      this%numCntcts(3) = this%numCntcts(3) + 1
      call ContactForce_PPFix(id1,id2,item2,ovrlp)

    elseif(item1 == -1 ) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPPFix","The inserted item's id is greater than allowed value:"//num2str(id1))
    else
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPPFix","The container is full and there is no space for new item" )      
    endif
  end subroutine CL_AddContactPPFix

  !**********************************************************************
  ! Adding a contact pair to the contact list 
  !**********************************************************************
  subroutine CL_AddContactPW(this,id1,id2,ovrlp,nv)
    implicit none
    class(ContactList):: this
    integer,intent(in):: id1,id2
    integer::wid,item1,item2
    real(RK)::ovrlp
    type(real3)::nv

    integer::i,j,k

    wid= DEMGeometry%pWall(id2)%wall_id
    call this%Find_Insert(id1,wid,item1,item2)
    
    ! item1 is the status of the item insertion. 1:old, 2:new 
    ! item2 is the container index
    if(item1>0) then
      CntctStatus(item2)=item1
      this%numCntcts(4) = this%numCntcts(4) + 1
      call ContactForce_PW(id1,id2,item2,ovrlp,nv)
    elseif(item1 == -1 ) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPW","The inserted item's id is greater than allowed value:"//num2str(id1))
    else
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPW","The container is full and there is no space for new item")      
    endif
  end subroutine CL_AddContactPW
    
  !***********************************************************************************************
  !* searching for an item, if it exists, it would return (/1, index of container that contain the item/),
  !     if it is new, pushing the item into the list and return (/2, index of container that
  !     contains the new item /), otherwise (-2,-2). 
  !********************************************************************************************
  subroutine Find_insert(this,bkt_id,value,item1,item2)
    implicit none
    class(ContactList) this
    integer,intent(in) :: bkt_id,value
    integer,intent(out):: item1, item2

    ! locals
    integer::n,nextI
        
    ! The bucket id is not in the range and must return (-1,-1), this is an error for the linked list
    if(bkt_id>GPrtcl_list%mlocal) then
      item1 =-1; item2 = -1
      return
    end if
        
     n = Bucket(bkt_id) 
     ! the item already exists in the list; returning the proper code and the index of container in which the item exists
    do while(n>0)
      if(id_j(n)==value) then
        item1 = 1; item2 = n
        return
      end if
      n = Next(n)
    end do
        
    !the container is full and there is no more space for this item
    if(this%nextInsert==0) then
      item1 =-2; item2 = -2
      return
    endif
        
    ! the item is new and it should be pushed into the list, inserting new item in the list 
    nextI = this%nextInsert     
    this%nextInsert = -Next(nextI)
    id_j(nextI)=value
    Next(nextI)=Bucket(bkt_id)
    Bucket(bkt_id) = nextI
    Tang_del(nextI)= zero_r3
    item1 =2; item2 = nextI
  end subroutine Find_insert
 
  !**********************************************************************
  ! Removing all released contacts (those which are not in contact in this time step) 
  ! from contact list.
  !**********************************************************************  
  subroutine CL_RemvReleased(this)
    class(ContactList)::this

    ! locals
    integer:: i,n

    do i=1, this%Max_numCntcts
      n=Next(i); if(n<=0) cycle
      if(CntctStatus(n) == -2) then
        CntctStatus(n) = -1
        Next(i)=Next(n)
                
        Next(n) = -this%nextInsert
        this%nextInsert = n
      endif
    enddo
        
    do i=1, this%mBucket
      n=Bucket(i); if(n<=0) cycle
      if(CntctStatus(n)==-2) then
        CntctStatus(n) =-1
        Bucket(i) = Next(n)
                
        Next(n) = - this%nextInsert
        this%nextInsert = n
      endif
    enddo

  end subroutine CL_RemvReleased

  !**********************************************************************
  ! CL_PreIteration
  !**********************************************************************
  subroutine CL_PreIteration(this)
    class(ContactList) this

    ! locals
    integer:: i

    this%numCntcts = 0
    do i=1,this%Max_numCntcts
       if(CntctStatus(i)>0) CntctStatus(i) = -2  ! flag previous contact
    enddo
    Head_Cp = -1
    Next_Cp = -1

  end subroutine CL_PreIteration

  !**********************************************************************
  ! CL_copy
  !**********************************************************************
  subroutine CL_copy(this,id1,id2)
    implicit none
    class(ContactList)::this
    integer,intent(in)::id1,id2

    ! locals
    integer::n,NextI
   
    n=bucket(id1)
    do while(n>0)
      CntctStatus(n)=-1
      NextI = Next(n)
      Next(n) = -this%nextInsert
      this%nextInsert = n      
      n = NextI
    enddo

    bucket(id1) = bucket(id2)
    Head_Cp(id1)= Head_Cp(id2)
    bucket(id2) = 0
    Head_Cp(id2)=-1

  end subroutine CL_copy

  !**********************************************************************
  ! CL_getPrtcl_nlink
  !**********************************************************************
  function CL_getPrtcl_nlink(this,pid) result(res)
    implicit none
    class(ContactList)::this
    integer,intent(in)::pid

    ! locals
    integer::res,n

    res = 0
    n = Bucket(pid)
    do while(n>0)
      res = res + 1
      n = Next(n)
    enddo

    n = Head_cp(pid)
    do while(n .ne. -1)
      res = res + 1
      n = Next_Cp(n)
    enddo
  end function CL_getPrtcl_nlink

  !**********************************************************************
  ! CL_IsCntct
  !**********************************************************************
  function CL_IsCntct(this,pid) result(res)
    implicit none
    class(ContactList)::this
    integer,intent(in)::pid

    ! locals
    integer::res

    if(Bucket(pid)>0 .or. Head_cp(pid)>0) then
      res= 1
    else
      res= 0
    endif
  end function CL_IsCntct

  !**********************************************************************
  ! CL_Prepare_Restart
  !**********************************************************************
  subroutine CL_Prepare_Restart(this,nlink_ind)
    implicit none
    class(ContactList)::this
    integer,intent(in)::nlink_ind

    ! locals
    integer::i,CLid

    CLid = nlink_ind
    do i=1,this%Max_numCntcts
      if(CntctStatus(i)>0) then
        CLid = CLid + 1
        CntctStatus(i) = CLid 
      endif
    enddo
  end subroutine CL_Prepare_Restart

  !**********************************************************************
  ! CL_Get_numCntcts
  !**********************************************************************
  function CL_Get_numCntcts(this) result(res)
    implicit none
    class(ContactList)::this
    integer:: res

    ! locals
    integer:: i,n

    res=0
    DO i=1,GPrtcl_list%nlocal
      n=Bucket(i)
      do while(n>0)
        res=res+1; n=Next(n)
      enddo
     
      n=Head_cp(i)
      do while(n>0)
        res=res+1; n=Next_cp(n)
      enddo
    ENDDO
  end function CL_Get_numCntcts

  !**********************************************************************
  ! CL_Get_numTanDel
  !**********************************************************************
  function CL_Get_numTanDel(this) result(res)
    implicit none
    class(ContactList)::this
    integer:: res

    ! locals
    integer::i
    
    res=0
    do i=1,this%Max_numCntcts
      if(CntctStatus(i)>0) then
        res = res + 1
      endif
    enddo
  end function CL_Get_numTanDel

  !**********************************************************************
  ! CL_GetNext_TanDel
  !**********************************************************************
  subroutine CL_GetNext_TanDel(this,TanDel,prev,now)
    implicit none
    class(ContactList)::this 
    type(real3),intent(out)::TanDel
    integer,intent(in)::prev
    integer,intent(out)::now

    ! locals
    integer::i
    
    do i=prev,this%Max_numCntcts
      if(CntctStatus(i)>0) then
        TanDel = Tang_del(i)
        now=i;exit
      endif
    enddo

  end subroutine CL_GetNext_TanDel

  !**********************************************************************
  ! CL_Gather_Cntctlink
  !**********************************************************************
  subroutine CL_Gather_Cntctlink(this,pid,buf_send,m)
    implicit none
    class(ContactList)::this
    integer,intent(in)::pid
    real(RK),dimension(*),intent(out)::buf_send
    integer,intent(inout)::m

    ! locals
    integer::n

    n = Bucket(pid)
    do while(n>0)
      buf_send(m)=real(id_j(n));  m=m+1
      buf_send(m)=tang_del(n)%x;  m=m+1
      buf_send(m)=tang_del(n)%y;  m=m+1
      buf_send(m)=tang_del(n)%z;  m=m+1
      n = Next(n)
    enddo        

    n = Head_cp(pid)
    do while(n .ne. -1)
      buf_send(m)=real(id_i(n));  m=m+1
      buf_send(m)=tang_del(n)%x;  m=m+1
      buf_send(m)=tang_del(n)%y;  m=m+1
      buf_send(m)=tang_del(n)%z;  m=m+1     
      n = Next_Cp(n)
    enddo

  end subroutine CL_Gather_Cntctlink

  !**********************************************************************
  ! CL_Add_RestartCntctlink
  !**********************************************************************
  subroutine CL_Add_RestartCntctlink(this,pid1,ncv,CntctVec,TanDelVel)
    implicit none
    class(ContactList)::this
    integer,intent(in)::pid1,ncv
    integer,dimension(20),intent(in)::CntctVec
    type(real3),dimension(20),intent(in)::TanDelVel

    ! locals
    integer::i,j,nlocal,pid2,gid1,gid2,nextI
  
    nlocal=GPrtcl_list%nlocal
    gid1= GPrtcl_id(pid1)
    DO j=1,ncv
      pid2=0;  gid2=CntctVec(j)
      DO i=1,nlocal
        if(GPrtcl_id(i)==gid2) then
          pid2=i; exit
        endif
      ENDDO

      nextI = this%nextInsert     
      this%nextInsert = -Next(nextI)
      IF(this%nextInsert==0) THEN
        call DEMLogInfo%CheckForError(ErrT_Abort,"CL_Add_RestartCntctlink","The container is full and there is no space for new item") 
      ENDIF

      ! here gid2 can be a particle_within_this_processor/ghost_particle/wall
      id_j(nextI) = gid2
      Next(nextI) = Bucket(pid1)
      Bucket(pid1) = nextI
      CntctStatus(nextI)=2
      Tang_del(nextI) = TanDelVel(j)
      IF(pid2>0 .and. gid1<gid2) THEN  ! gid2 is also within this process
        id_i(nextI)= gid1
        Next_Cp(nextI)=Head_Cp(pid2)
        Head_Cp(pid2)=nextI
      ENDIF
    ENDDO
  end subroutine CL_Add_RestartCntctlink

  !**********************************************************************
  ! CL_Gather_Cntctlink_Restart
  !**********************************************************************
  subroutine CL_Gather_Cntctlink_Restart(this,pid,CntctVec,ncv)
    implicit none
    class(ContactList)::this
    integer,intent(in)::pid
    integer,dimension(40),intent(out)::CntctVec
    integer,intent(out)::ncv

    ! locals
    integer::n
    
    ncv =0
    n = Bucket(pid)
    do while(n>0)
      ncv = ncv + 1
      ! for monosize particles, a particle can contact with NO MORE THAN 12 neighbor particles.
      if(ncv>20) call DEMLogInfo%CheckForError(ErrT_Abort,"CL_Gather_Cntctlink_Restart","so big ncv")
      CntctVec(2*ncv-1) =  id_j(n)
      CntctVec(2*ncv)   =  CntctStatus(n)
      n = Next(n)
    enddo        

    n = Head_cp(pid)
    do while(n .ne. -1)
      ncv = ncv + 1
      if(ncv>20) call DEMLogInfo%CheckForError(ErrT_Abort,"CL_Gather_Cntctlink_Restart","so big ncv")
      CntctVec(2*ncv-1) =  id_i(n)
      CntctVec(2*ncv)   =  CntctStatus(n)  
      n = Next_Cp(n)
    enddo

  end subroutine CL_Gather_Cntctlink_Restart
  !**********************************************************************
  ! CL_printCL
  !********************************************************************** 
  subroutine CL_printCL(this)
    implicit none
    class(ContactList)::this
   
    ! locals
    integer::i,n

    open(nrank+50,file='global_id'//trim(num2str(nrank))//'.txt',status='replace')
    open(nrank+60,file='cntctlist'//trim(num2str(nrank))//'.txt',status='replace')
    DO i=1,GPrtcl_list%nlocal
      write(nrank+50,*)i,GPrtcl_id(i)
    ENDDO

    DO i=1,GPrtcl_list%nlocal
      n = Bucket(i)
      do while(n>0)
        write(nrank+60,*)i,GPrtcl_id(i),id_j(n) 
        n = Next(n)
      enddo
     
      n = Head_cp(i)
      do while(n>0)
        write(nrank+60,*)i,GPrtcl_id(i),id_i(n) 
        n = Next_cp(n)
      enddo
    ENDDO
    close(nrank+50)
    close(nrank+60)
  end subroutine CL_printCL
  
  !**********************************************************************
  ! CL_Add_Cntctlink
  !**********************************************************************
  subroutine CL_Add_Cntctlink(this,pid,buf_recv,m)
    implicit none
    class(ContactList)::this
    integer,intent(in)::pid
    real(RK),dimension(*),intent(in)::buf_recv
    integer,intent(inout)::m   

    ! locals
    integer::value,nextI
    real(RK)::realt

    do
      realt = buf_recv(m); m=m+1
      if(abs(realt-END_OF_PRTCL)<1.00E-10_RK) return

      nextI = this%nextInsert     
      this%nextInsert = -Next(nextI)
      if(this%nextInsert==0) then
        call DEMLogInfo%CheckForError(ErrT_Abort,"CL_Add_Cntctlink","The container is full and there is no space for new item") 
      endif
      id_j(nextI) = int(realt+0.2)
      Next(nextI) = Bucket(pid)
      Bucket(pid) = nextI
      CntctStatus(nextI)=2

      Tang_del(nextI)%x = buf_recv(m); m=m+1
      Tang_del(nextI)%y = buf_recv(m); m=m+1
      Tang_del(nextI)%z = buf_recv(m); m=m+1
    enddo
   
  end subroutine CL_Add_Cntctlink

  !**********************************************************************
  ! ContactForce_PP
  !**********************************************************************
  subroutine ContactForce_PP(pid,pjd,ind,ovrlp)
    implicit none
    integer,intent(in)::pid,pjd,ind
    real(RK),intent(in)::ovrlp
    
    integer:: jid,pt_i,pt_j
    real(RK)::Ri,Rj,fn,ft,vrn,ft_fric,kt
    type(real3)::rveli,rvelj,norm_v,ovlp_t,fnij,ftij,Vrij,Vij_n,Vij_t,Mij,Mji,Mri,Mrj,veli,velj
    type(BinaryProperty)::Prop_ij
    type(real4)::posi,posj    
      
    ! type of particle to get physical properties
    pt_i = GPrtcl_pType(pid)
    pt_j = GPrtcl_pType(pjd)
    prop_ij = DEMProperty%Prtcl_BnryProp(pt_i, pt_j)
       
    veli = GPrtcl_linVel(1,pid)
    velj = GPrtcl_linVel(1,pjd)
    rveli = GPrtcl_rotVel(1,pid)
    rvelj = GPrtcl_rotVel(1,pjd)
    posi = GPrtcl_PosR(pid)
    posj = GPrtcl_PosR(pjd)
    Ri = posi%w  
    Rj = posj%w  
   ! ovrlp = posj .ovlp. posi ! normal overlap
    norm_v = posj .nv. posi  ! normal vector, posj-posi
    
    ! relative velocity at contact point
    !(vi - vj ) + cross( (Ri*wi + Rj*wj), nij) ; 
    Vrij = veli-velj + ((Ri*rveli+ Rj*rvelj).cross.norm_v)
    
    ! normal and tangential velocities vectors 
    Vij_n = (Vrij .dot. norm_v)*norm_v
    Vij_t = Vrij - Vij_n
    ! normal relative velocity 
    vrn = Vrij .dot. norm_v
    
    ! tangential overlap vector
    ovlp_t= (Vij_t*DEM_opt%dt)+Tang_del(ind)
   
    ! computing the normal and tangential contact forces 
    if(DEM_opt%CF_Type ==CFT_LSD_nl .or. DEM_opt%CF_Type ==CFT_LSD_l) then
      fn =  -prop_ij%LinearSpringStiffness_n*ovrlp -prop_ij%dmp_n* vrn         ! 2.34, p29
      fnij = fn * norm_v
      ftij = ((-prop_ij%LinearSpringStiffness_t)*ovlp_t)-(prop_ij%dmp_t*Vij_t) ! 2.48, p33
    else
      fn = -(four/three*prop_ij%YoungsModEff*sqrt(prop_ij%RadEff)*ovrlp**1.5_RK)-(prop_ij%dmp_n*ovrlp**0.25_RK * vrn) ! 2.62, p39
      fnij = fn * norm_v
      ftij = (-sixteen/three*prop_ij%ShearModEff*sqrt(prop_ij%RadEff*ovrlp))*ovlp_t  ! 2.72, p44
    endif
    
    ! Coulomb's friction law
    ft = norm(ftij)
    ft_fric = prop_ij%FrictionCoe * abs(fn)
    if(abs(ft).gt.ft_fric) then  
       if(norm(ovlp_t)>zero ) then
         ftij =ftij*ft_fric/ft   ! 2.49, p33
         if(DEM_opt%CF_Type==CFT_LSD_l) then
           ovlp_t =(-one)*(ftij/prop_ij%LinearSpringStiffness_t )   ! 2.50, p33
         elseif(DEM_opt%CF_Type == CFT_nLin_l) then
           kt = sixteen/three*prop_ij%ShearModEff*sqrt(prop_ij%RadEff*ovrlp)
           ovlp_t =(-one)*(ftij/kt)
         endif
       else
       ftij = zero_r3   
     endif
    endif
    
    ! computing torque acting on spheres it also includes rolling torque 
    call clc_Torque(fn,ftij,prop_ij%FrictionCoe_Roll, norm_v, rveli, rvelj, Ri, Rj, Mij, Mji, Mri, Mrj)
    
    ! updating the contact force and torques of particles i and j
    GPrtcl_cntctForce(pid) = GPrtcl_cntctForce(pid) + (fnij+ftij)
    GPrtcl_cntctForce(pjd) = GPrtcl_cntctForce(pjd) - (fnij+ftij) 
    GPrtcl_torque(pid) = GPrtcl_torque(pid) + Mij+Mri
    GPrtcl_torque(pjd) = GPrtcl_torque(pjd) + Mji+Mrj     
    
    ! setting the updated contact info pair into the contact list 
    Tang_del(ind) = ovlp_t
  end subroutine ContactForce_PP

  !**********************************************************************
  ! ContactForce_PPG
  !**********************************************************************
  subroutine ContactForce_PPG(pid,pjd,ind,ovrlp)
    implicit none
    integer,intent(in)::pid,pjd,ind
    real(RK),intent(in)::ovrlp
    
    integer:: jid,pt_i,pt_j
    real(RK)::Ri,Rj,fn,ft,vrn,ft_fric,kt
    type(real3)::rveli,rvelj,norm_v,ovlp_t,fnij,ftij,Vrij,Vij_n,Vij_t,Mij,Mji,Mri,Mrj,veli,velj
    type(BinaryProperty)::Prop_ij
    type(real4)::posi,posj    
      
    ! type of particle to get physical properties
    pt_i = GPrtcl_pType(pid)
    pt_j = GhostP_pType(pjd)
    prop_ij = DEMProperty%Prtcl_BnryProp(pt_i, pt_j)
       
    veli = GPrtcl_linVel(1,pid)
    velj = GhostP_linVel(pjd)
    rveli = GPrtcl_rotVel(1,pid)
    rvelj = GhostP_rotVel(pjd)
    posi = GPrtcl_PosR(pid)
    posj = GhostP_PosR(pjd)
    Ri = posi%w  
    Rj = posj%w  
   ! ovrlp = posj .ovlp. posi ! normal overlap
    norm_v = posj .nv. posi  ! normal vector, posj-posi
    
    ! relative velocity at contact point
    !(vi - vj ) + cross( (Ri*wi + Rj*wj), nij) ; 
    Vrij = veli-velj + ((Ri*rveli+ Rj*rvelj).cross.norm_v)
    
    ! normal and tangential velocities vectors 
    vrn = Vrij .dot. norm_v
    Vij_n = vrn*norm_v
    Vij_t = Vrij - Vij_n
    
    ! tangential overlap vector
    ovlp_t= (Vij_t*DEM_opt%dt)+Tang_del(ind)
   
    ! computing the normal and tangential contact forces 
    if(DEM_opt%CF_Type ==CFT_LSD_nl .or. DEM_opt%CF_Type ==CFT_LSD_l) then
      fn =  -prop_ij%LinearSpringStiffness_n*ovrlp -prop_ij%dmp_n* vrn         ! 2.34, p29
      fnij = fn * norm_v
      ftij = ((-prop_ij%LinearSpringStiffness_t)*ovlp_t)-(prop_ij%dmp_t*Vij_t) ! 2.48, p33
    else
      fn = -(four/three*prop_ij%YoungsModEff*sqrt(prop_ij%RadEff)*ovrlp**1.5_RK)-(prop_ij%dmp_n*ovrlp**0.25_RK * vrn) ! 2.62, p39
      fnij = fn * norm_v
      ftij = (-sixteen/three*prop_ij%ShearModEff*sqrt(prop_ij%RadEff*ovrlp))*ovlp_t  ! 2.72, p44
    endif
    
    ! Coulomb's friction law
    ft = norm(ftij)
    ft_fric = prop_ij%FrictionCoe * abs(fn)
    if(abs(ft).gt.ft_fric) then  
       if(norm(ovlp_t)>zero ) then
         ftij =ftij*ft_fric/ft   ! 2.49, p33
         if(DEM_opt%CF_Type==CFT_LSD_l) then
           ovlp_t =(-one)*(ftij/prop_ij%LinearSpringStiffness_t )   ! 2.50, p33
         elseif(DEM_opt%CF_Type == CFT_nLin_l) then
           kt = sixteen/three*prop_ij%ShearModEff*sqrt(prop_ij%RadEff*ovrlp)

           ovlp_t =(-one)*(ftij/kt)
         endif
       else
       ftij = zero_r3   
     endif
    endif
    
    ! computing torque acting on spheres it also includes rolling torque 
    call clc_Torque(fn,ftij,prop_ij%FrictionCoe_Roll, norm_v, rveli, rvelj, Ri, Rj, Mij, Mji, Mri, Mrj)
    
    ! updating the contact force and torques of particles i and j
    GPrtcl_cntctForce(pid) = GPrtcl_cntctForce(pid) + (fnij+ftij) 
    GPrtcl_torque(pid) = GPrtcl_torque(pid) + Mij+Mri     
    
    ! setting the updated contact info pair into the contact list 
    Tang_del(ind) = ovlp_t
  end subroutine ContactForce_PPG

  !**********************************************************************
  ! ContactForce_PPFix
  !**********************************************************************
  subroutine ContactForce_PPFix(pid,pjd,ind,ovrlp)
    implicit none
    integer,intent(in)::pid,pjd,ind
    real(RK),intent(in)::ovrlp
    
    integer:: jid,pt_i,pt_j
    real(RK)::Ri,Rj,fn,ft,vrn,ft_fric,kt
    type(real3)::rveli,norm_v,ovlp_t,fnij,ftij,Vrij,Vij_n,Vij_t,Mij,Mji,Mri,Mrj,veli
    type(BinaryProperty)::Prop_ij
    type(real4)::posi,posj    
      
    ! type of particle to get physical properties
    pt_i = GPrtcl_pType(pid)
    pt_j = GPFix_pType(pjd)
    prop_ij = DEMProperty%Prtcl_BnryProp(pt_i, pt_j)
       
    veli = GPrtcl_linVel(1,pid)
    rveli = GPrtcl_rotVel(1,pid)
    posi = GPrtcl_PosR(pid)
    posj = GPFix_PosR(pjd)
    Ri = posi%w  
    Rj = posj%w  
   ! ovrlp = posj .ovlp. posi ! normal overlap
    norm_v = posj .nv. posi  ! normal vector, posj-posi
    
    ! relative velocity at contact point
    !(vi - vj ) + cross( (Ri*wi + Rj*wj), nij) ; 
    Vrij = veli + ((Ri*rveli).cross.norm_v)
    
    ! normal and tangential velocities vectors 
    vrn = Vrij .dot. norm_v
    Vij_n = vrn*norm_v
    Vij_t = Vrij - Vij_n
    
    ! tangential overlap vector
    ovlp_t= (Vij_t*DEM_opt%dt)+Tang_del(ind)
   
    ! computing the normal and tangential contact forces 
    if(DEM_opt%CF_Type ==CFT_LSD_nl .or. DEM_opt%CF_Type ==CFT_LSD_l) then
      fn =  -prop_ij%LinearSpringStiffness_n*ovrlp -prop_ij%dmp_n* vrn         ! 2.34, p29
      fnij = fn * norm_v
      ftij = ((-prop_ij%LinearSpringStiffness_t)*ovlp_t)-(prop_ij%dmp_t*Vij_t) ! 2.48, p33
    else
      fn = -(four/three*prop_ij%YoungsModEff*sqrt(prop_ij%RadEff)*ovrlp**1.5_RK)-(prop_ij%dmp_n*ovrlp**0.25_RK * vrn) ! 2.62, p39
      fnij = fn * norm_v
      ftij = (-sixteen/three*prop_ij%ShearModEff*sqrt(prop_ij%RadEff*ovrlp))*ovlp_t  ! 2.72, p44
    endif
    
    ! Coulomb's friction law
    ft = norm(ftij)
    ft_fric = prop_ij%FrictionCoe * abs(fn)
    if(abs(ft).gt.ft_fric) then  
       if(norm(ovlp_t)>zero ) then
         ftij =ftij*ft_fric/ft   ! 2.49, p33
         if(DEM_opt%CF_Type==CFT_LSD_l) then
           ovlp_t =(-one)*(ftij/prop_ij%LinearSpringStiffness_t )   ! 2.50, p33
         elseif(DEM_opt%CF_Type == CFT_nLin_l) then
           kt = sixteen/three*prop_ij%ShearModEff*sqrt(prop_ij%RadEff*ovrlp)

           ovlp_t =(-one)*(ftij/kt)
         endif
       else
       ftij = zero_r3   
     endif
    endif
    
    ! computing torque acting on spheres it also includes rolling torque 
    call clc_Torque(fn,ftij,prop_ij%FrictionCoe_Roll, norm_v, rveli, zero_r3, Ri, Rj, Mij, Mji, Mri, Mrj)
    
    ! updating the contact force and torques of particles i and j
    GPrtcl_cntctForce(pid) = GPrtcl_cntctForce(pid) + (fnij+ftij) 
    GPrtcl_torque(pid) = GPrtcl_torque(pid) + Mij+Mri     
    
    ! setting the updated contact info pair into the contact list 
    Tang_del(ind) = ovlp_t
  end subroutine ContactForce_PPFix

  !**********************************************************************
  ! ContactForce_PW
  !**********************************************************************
  subroutine ContactForce_PW(pid,mwi,ind,ovrlp,norm_v )
    implicit none
    integer,intent(in)::pid,mwi,ind
    real(RK),intent(in):: ovrlp
    type(real3),intent(inout):: norm_v
    
    integer:: jid,pt_i, w_pt
    real(RK)::Ri,Rj,fn,ft,ft_fric,vrn,kt
    type(real3)::rveli,rvelj,Vrij,Vij_n,Vij_t,ovlp_t,fnij,ftij,Mij,Mji,Mri,Mrj,veli,velj
    type(BinaryProperty):: Prop_ij
    type(real4):: posi
        
    ! particle property type
    pt_i = GPrtcl_pType(pid)
    
    veli = GPrtcl_linVel(1,pid)
    rveli= GPrtcl_rotVel(1,pid)
    posi = GPrtcl_PosR(pid)
    
    ! getting normal vector, normal distance, velocity at contact point and wall property type
    velj = DEMGeometry%pWall(mwi)%trans_vel
    w_pt = DEMGeometry%pWall(mwi)%wall_Type
    
    ! since the normal vector points from particle i to j we must negate the normal vector
    norm_v = (-one) * norm_v
    Ri = posi%w
    
    ! relative velocity at contact point 
    !(vi - vj ) + cross( (Ri*wi + Rj*wj), nij) ; 
    Vrij = veli-velj + ((Ri*rveli) .cross. norm_v )
    
    ! normal and tangential velocities vectors
    vrn = Vrij .dot. norm_v   
    Vij_n = vrn*norm_v
    Vij_t = Vrij - Vij_n
     
    ! getting particle-wall property
    prop_ij = DEMProperty%PrtclWall_BnryProp(pt_i,w_pt)

    ! tangential overlap vector
    ovlp_t = (Vij_t*DEM_opt%dt)+Tang_del(ind)
          
    ! computing the normal and tangential contact forces
    if(DEM_opt%CF_Type ==CFT_LSD_nl .or.  DEM_opt%CF_Type ==CFT_LSD_l) then
      fn = -prop_ij%LinearSpringStiffness_n * ovrlp - prop_ij%dmp_n* vrn
      fnij =fn * norm_v
      ftij =((-prop_ij%LinearSpringStiffness_t)*ovlp_t)-(prop_ij%dmp_t * Vij_t)  
    else
      fn = -(four/three*prop_ij%YoungsModEff*sqrt(prop_ij%RadEff)*ovrlp**1.5_RK)-(prop_ij%dmp_n*ovrlp**0.25_RK * vrn)
      fnij = fn * norm_v
      ftij = (- sixteen/three * prop_ij%ShearModEff * sqrt(prop_ij%RadEff*ovrlp) ) * ovlp_t        
    endif
    
    ! Coulomb's friction law
    ft = norm(ftij)
    ft_fric = prop_ij%FrictionCoe * abs(fn)
    if( ft > ft_fric ) then
      if(norm(ovlp_t) > zero ) then
        ftij = ftij*ft_fric/ft
        if(DEM_opt%CF_Type == CFT_LSD_l ) then
          ovlp_t = (-one) * (ftij/prop_ij%LinearSpringStiffness_t )
        elseif(DEM_opt%CF_Type == CFT_nLin_l) then
          kt = sixteen/three * prop_ij%ShearModEff * sqrt(prop_ij%RadEff*ovrlp)
          ovlp_t =  (-one) * (ftij/kt)
        endif
      else
        ftij = zero_r3    
      endif
    endif
    
    ! computing torque acting on spheres it also includes rolling torque 
    Rj=1.00e20_RK*Ri
    call clc_Torque(fn,ftij, prop_ij%FrictionCoe_Roll,norm_v,rveli,zero_r3,Ri,Rj,Mij,Mji,Mri,Mrj)
    
    ! updating the contact force and torques of particle i   
    GPrtcl_cntctForce(pid) = GPrtcl_cntctForce(pid) + fnij+ftij    
    GPrtcl_torque(pid) = GPrtcl_torque(pid) + Mij+Mri
    
    ! setting the updated contact pair into the contact list 
    Tang_del(ind) = ovlp_t

  end subroutine ContactForce_PW

  !**********************************************************************
  ! clc_Torque
  !**********************************************************************
  subroutine clc_Torque(fn,ftij,roll_fric,nij,wi,wj,Ri,Rj,Mij,Mji,Mri,Mrj)
    implicit none
    real(RK),intent(in):: fn,roll_fric
    type(real3),intent(in):: ftij,nij,wi,wj 
    real(RK),intent(in):: Ri,Rj
    type(real3),intent(out):: Mij,Mji,Mri,Mrj
    real(RK)::RadEff,w_hat_mag
    type(real3)::M,w_hat
    
    ! tangential torque
    M = nij .cross. ftij
    Mij = Ri*M
    Mji = Rj*M
    
    ! rolling torque
    w_hat = wi-wj
    w_hat_mag = norm(w_hat) ! 2.123, p56
    if(w_hat_mag .gt. 0.000001_RK) then
      w_hat = w_hat/w_hat_mag
    else
      w_hat = zero_r3
    end if
    
    RadEff=(Ri*Rj)/(Ri+Rj)
    if(DEM_opt%CT_Model == CTM_ConstantTorque) then    
      Mri=(-roll_fric*abs(fn)*RadEff)*w_hat  ! 2.122, p56
    else
      Mri=(-roll_fric*abs(fn)*RadEff*norm((Ri*wi+Rj*wj).cross.nij))*w_hat  ! 2.124, p57
    end if
    Mrj = (-one)*Mri
  end subroutine clc_Torque

end module Prtcl_CL_and_CF
