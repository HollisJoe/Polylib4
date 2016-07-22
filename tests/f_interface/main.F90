
!////////////////////////////////////////////////////////////////////////////
!///
!/// Fortranインターフェーステスト
!///
!////////////////////////////////////////////////////////////////////////////


!//----------------------------------------------------
!//  メインルーチン
!//----------------------------------------------------
program main

#ifdef MPI_PL
  use mpi
#endif
  implicit none

  ! Fortran用型、構造体、精度の定義
  !    gfortranではインクルードファイル内で#ifdefが使えないため
  !    FPolylib_User.incから呼ばれるFPolylib_precision.inc内で
  !    精度を切り替えられない
  !    インクルードファイルを個別に呼び出す形にする

#if 1
#include "FPolylib_precision_def.inc"
#else
#ifdef _REAL_IS_DOUBLE_
  integer,  parameter :: PL_REAL_PN = 8
#else
  integer,  parameter :: PL_REAL_PN = 4
#endif
#endif

  include 'FPolylib_define.inc'
  include 'FPolylib_struct.inc'
  include 'FPolylib_prototype.inc'

  integer :: ret,ierr
  integer :: num_rank,myrank
  integer :: i
  character(len=PL_FILE_PATH_LEN) :: config_file_name
  character(len=PL_FILE_PATH_LEN) :: config_file_name_out
  character(len=PL_FORMAT_LEN)   :: fmt_out
  character(len=PL_STR_LEN)      :: extend
  integer num_pg
  integer(PL_TAG_PN), allocatable :: pg_tags(:)
  character(PL_GRP_PATH_LEN) :: pg_name
  integer movable
  integer num_tri
  integer(PL_TAG_PN), allocatable :: tags_tri(:)
  integer pl_type
  real(PL_REAL_PN) vertex(9)
  real(PL_REAL_PN),save :: scale = 1.0

#ifdef MPI_PL
!
! モデル空間
!   想定空間 x: -50〜150
!            y: -50〜150
!            z:   0〜150
!
!  (初期値）
!   plane    x:   0〜100
!            y:   0〜100
!            z:   0〜  0
!   spehre   x: -50〜 50
!            y: -50〜 50
!            z:  50〜150
!  (移動）
!     Fortranは移動関数が登録出来ないため、移動なし
!

#if 0
! 並列数１：分割なし
  type (FParallelBboxStruct) :: bboxes(1)
  bboxes(1)%bpos  (1) = -50.0
  bboxes(1)%bpos  (2) = -50.0
  bboxes(1)%bpos  (3) = -50.0
  bboxes(1)%bbsize(1) =  20
  bboxes(1)%bbsize(2) =  20
  bboxes(1)%bbsize(3) =  20
  bboxes(1)%gcsize(1) =  2
  bboxes(1)%gcsize(2) =  2
  bboxes(1)%gcsize(3) =  2
  bboxes(1)%dx(1)     =  10.0
  bboxes(1)%dx(2)     =  10.0
  bboxes(1)%dx(3)     =  10.0
#endif
#if 0
! 並列数２：Y方向を +50 の位置で２分割
  type (FParallelBboxStruct) bboxes(2)
  bboxes(1)%bpos  (1) = -50.0
  bboxes(1)%bpos  (2) = -50.0
  bboxes(1)%bpos  (3) = -50.0
  bboxes(1)%bbsize(1) =  20
  bboxes(1)%bbsize(2) =  10
  bboxes(1)%bbsize(3) =  20
  bboxes(1)%gcsize(1) =  2
  bboxes(1)%gcsize(2) =  2
  bboxes(1)%gcsize(3) =  2
  bboxes(1)%dx(1)     =  10.0
  bboxes(1)%dx(2)     =  10.0
  bboxes(1)%dx(3)     =  10.0

  bboxes(2)%bpos  (1) = -50.0
  bboxes(2)%bpos  (2) =  50.0
  bboxes(2)%bpos  (3) = -50.0
  bboxes(2)%bbsize(1) =  20
  bboxes(2)%bbsize(2) =  10
  bboxes(2)%bbsize(3) =  20
  bboxes(2)%gcsize(1) =  2
  bboxes(2)%gcsize(2) =  2
  bboxes(2)%gcsize(3) =  2
  bboxes(2)%dx(1)     =  10.0
  bboxes(2)%dx(2)     =  10.0
  bboxes(2)%dx(3)     =  10.0
#endif
#if 1
! 並列数４：X,Y方向を +50 の位置で２分割
  type (FParallelBboxStruct) bboxes(4)
  bboxes(1)%bpos  (1) = -50.0
  bboxes(1)%bpos  (2) = -50.0
  bboxes(1)%bpos  (3) = -50.0
  bboxes(1)%bbsize(1) =  10
  bboxes(1)%bbsize(2) =  10
  bboxes(1)%bbsize(3) =  20
  bboxes(1)%gcsize(1) =  2
  bboxes(1)%gcsize(2) =  2
  bboxes(1)%gcsize(3) =  2
  bboxes(1)%dx(1)     =  10.0
  bboxes(1)%dx(2)     =  10.0
  bboxes(1)%dx(3)     =  10.0

  bboxes(2)%bpos  (1) =  50.0
  bboxes(2)%bpos  (2) = -50.0
  bboxes(2)%bpos  (3) = -50.0
  bboxes(2)%bbsize(1) =  10
  bboxes(2)%bbsize(2) =  10
  bboxes(2)%bbsize(3) =  20
  bboxes(2)%gcsize(1) =  2
  bboxes(2)%gcsize(2) =  2
  bboxes(2)%gcsize(3) =  2
  bboxes(2)%dx(1)     =  10.0
  bboxes(2)%dx(2)     =  10.0
  bboxes(2)%dx(3)     =  10.0

  bboxes(3)%bpos  (1) = -50.0
  bboxes(3)%bpos  (2) =  50.0
  bboxes(3)%bpos  (3) = -50.0
  bboxes(3)%bbsize(1) =  10
  bboxes(3)%bbsize(2) =  10
  bboxes(3)%bbsize(3) =  20
  bboxes(3)%gcsize(1) =  2
  bboxes(3)%gcsize(2) =  2
  bboxes(3)%gcsize(3) =  2
  bboxes(3)%dx(1)     =  10.0
  bboxes(3)%dx(2)     =  10.0
  bboxes(3)%dx(3)     =  10.0

  bboxes(4)%bpos  (1) =  50.0
  bboxes(4)%bpos  (2) =  50.0
  bboxes(4)%bpos  (3) = -50.0
  bboxes(4)%bbsize(1) =  10
  bboxes(4)%bbsize(2) =  10
  bboxes(4)%bbsize(3) =  20
  bboxes(4)%gcsize(1) =  2
  bboxes(4)%gcsize(2) =  2
  bboxes(4)%gcsize(3) =  2
  bboxes(4)%dx(1)     =  10.0
  bboxes(4)%dx(2)     =  10.0
  bboxes(4)%dx(3)     =  10.0
#endif

#endif

  !-------------------------------------------
  config_file_name = 'polylib_config.tp'

  !-------------------------------------------
  !  初期化
  !-------------------------------------------

  ! MPI初期化
#ifdef MPI_PL
  call mpi_init( ierr )
  call mpi_comm_size( mpi_comm_world, num_rank, ierr )
  call mpi_comm_rank( mpi_comm_world, myrank , ierr )
  write(*,*) 'num_rank=',num_rank,' myrank=',myrank
#else
  num_rank = 1
  myrank   = 0
#endif

  ! Polylib(Fortran)初期化
  call fpolylib_instance( ret )

  ! 並列計算関連情報の設定と初期化
#ifdef MPI_PL
  call fpolylib_init_parallel_info (  &
                    bboxes(myrank+1)%bpos,    &
                    bboxes(myrank+1)%bbsize,  &
                    bboxes(myrank+1)%gcsize,  &
                    bboxes(myrank+1)%dx,      &
                    ret &
             )
  if( ret .ne. 0 ) then
    write(0,*) '(rk:',myrank,') # ERROR: fpolylib_init_parallel_info()'
    stop 1
  endif
#endif

  !-------------------------------------------
  !  ロード
  !-------------------------------------------

  ! 初期化ファイルを指定してデータロード
  ! write(*,*) '(rk:',myrank,') load() start'
  call fpolylib_load( config_file_name, scale, ret );
  ! write(*,*) '(rk:',myrank,') load() end  ret=',ret
  if( ret .ne. 0 ) then
    write(0,*) '(rk:',myrank,') # ERROR: fpolylib_load()'
    stop 1
  endif

  call fpolylib_get_root_groups_tags_num ( num_pg, ret );
  write(*,*) '(rk:',myrank,')  fpolylib_get_root_groups_tags_num() num_pg=',num_pg

  allocate  ( pg_tags(num_pg) )

  call fpolylib_get_root_groups_tags ( num_pg, pg_tags, ret );
  write(*,*) '(rk:',myrank,')  fpolylib_get_root_groups_tags() num_pg=',num_pg
  if( ret .ne. 0 ) then
    write(0,*) '(rk:',myrank,') # ERROR: fpolylib_get_root_groups_tags()'
    stop 1
  endif

  do i=1, num_pg

    call fpolylib_group_get_name( pg_tags(i), pg_name )
    call fpolylib_group_get_movable( pg_tags(i), movable )
    call fpolylib_group_get_num_triangles( pg_tags(i), num_tri, ret )
    write(*,*) '(rk:',myrank,')  i=',i,' pg_name(1:16)=[',pg_name(1:16),']  movable=',movable,' num_tri=',num_tri
    allocate  ( tags_tri(num_tri) )

    call fpolylib_group_get_triangles( pg_tags(i), num_tri, tags_tri, ret );
    write(*,*) '(rk:',myrank,')  num_tri=',num_tri
    if( ret .ne. 0 ) then
      write(0,*) '(rk:',myrank,') # ERROR: fpolylib_group_get_triangles()'
      stop 1
    endif

    pl_type = fpolylib_triangle_get_pl_type( tags_tri(1) );
    write(*,*) '(rk:',myrank,')  pl_type=',pl_type

    !call fpolylib_triangle_get_vertexes( tags_tri(1), vertex );
    !write(*,*) '(rk:',myrank,')  p1=',vertex(1),' ',vertex(2),' ',vertex(3)
    !write(*,*) '(rk:',myrank,')  p2=',vertex(4),' ',vertex(5),' ',vertex(6)
    !write(*,*) '(rk:',myrank,')  p3=',vertex(7),' ',vertex(8),' ',vertex(9)

    deallocate( tags_tri );

  enddo


    !-------------------------------------------
    !  セーブ
    !-------------------------------------------

  fmt_out = FILE_FMT_STL_A
  write(*,*) '(rk:',myrank,')  fmt_out=',fmt_out
  extend  = ''
  write(*,*) '(rk:',myrank,')  save()  start'
  call fpolylib_save( config_file_name_out, fmt_out, extend, ret );
  write(*,*) '(rk:',myrank,')  save()  end  config_file_name_out=',config_file_name_out(1:32)
  if( ret .ne. 0 ) then
    write(0,*) '(rk:',myrank,') # ERROR: fpolylib_save()()'
    stop 1
  endif

  !-------------------------------------------
  !  終了化
  !-------------------------------------------
  deallocate  ( pg_tags )

  ! MPI終了化
#ifdef MPI_PL
  write(*,*) '(rk:',myrank,')  mpi_finalize()  start'
  call mpi_finalize( ierr )
  write(*,*) '(rk:',myrank,')  mpi_finalize()  end'
#endif

  if( myrank .eq. 0 )  then
    write(0,*) '------------------------------------------'
    write(0,*) '   PASS :  f_interface (Normal End)'
    write(0,*) '------------------------------------------'
  endif


end program main
