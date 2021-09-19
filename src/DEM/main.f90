program main_DEM
  use MPI
  use Prtcl_decomp_2d
  use Prtcl_DEMSystem
  use Prtcl_Parameters
  use Prtcl_IOAndVisu
  use Prtcl_LogInfo
  use Prtcl_Variables
#ifdef DEM_DEBUG
  use Prtcl_CL_and_CF
#endif
  implicit none
  character(len=64)::chDEMPrm
  integer :: ierror,i
  
  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)

  ! read DEM options
  if(command_argument_count()/=1 .and. nrank==0) write(*,*)'command argument wrong!'
  call get_command_argument(1,chDEMPrm)
  call DEM_opt%ReadDEMOption( chDEMPrm)
  call DEM_decomp%Init_DECOMP(chDEMPrm)
  call DEM%Initialize(chDEMPrm) ! Topest level initialing for DEM body
#ifdef DEM_DEBUG
  call GPPW_CntctList%printCL(DEM_opt%ifirst-1)
#endif

  call DEM_IO%dump_visu(DEM_opt%ifirst-1)
  print*, nrank,GPrtcl_list%nlocal,GPrtcl_list%mlocalFix
  do i= DEM_opt%ifirst, DEM_opt%ilast
    call DEM%iterate(i)
  enddo
  call DEM_IO%Final_visu()

  if(nrank==0)call DEMLogInfo%OutInfo("Good job! DEM finished successfully!",1)
  call MPI_FINALIZE(ierror)
end program main_DEM
