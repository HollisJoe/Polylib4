!/*
!* Polylib - Polygon Management Library
!*
!* Copyright (c) 2010-2011 VCAD System Research Program, RIKEN.
!* All rights reserved.
!*
!* Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
!* All rights reserved.
!*
!*/

!*********************************************************
!
! Fortranインクルードファイル
!    文字列長の定義
!    Fortranユーザアプリ用のインクルードファイル
!    (Polylib内部はC++で実装されており当インクルードファイルは使用していない)
!
!*********************************************************

!----------------------------------------------------
! タグの精度　実体はアドレスなので8byte
!   integer(kind=PL_TAG_PN) :: tag   // fort90
!   integer(PL_TAG_PN)      :: tag   // fort90
!   integer*8 tag // fort77

  integer, parameter :: PL_TAG_PN = 8


!----------------------------------------------------
! Fortran文字列長の定義
!     FPolylib.h と合致していること
!
!    (Fortranアプリ使用方法)
!
!      character(len=PL_FILE_PATH_LEN) :: str  // fort90
!      character(PL_FILE_PATH_LEN)     :: str  // fort90
!      character*NNN str   // fort77   NNNは文字列の長さ
!      //character str(NNN)  // fort77   NNNは文字列の長さ  文字の配列となる
!

! ファイルパス、ファイル名のfortran文字列長
  integer, parameter :: PL_FILE_PATH_LEN = 256
  integer, parameter :: PL_FILE_NAME_LEN =  64

! ポリゴングループ名のfortran文字列長
  integer, parameter :: PL_GRP_PATH_LEN = 256
  integer, parameter :: PL_GRP_NAME_LEN =  64

! 属性用文字列のfortran文字列長
!     属性キー, 属性値 etc.
   integer, parameter :: PL_GRP_ATR_LEN = 32

! ファイルフォーマットのfortran文字列長
   integer, parameter :: PL_FORMAT_LEN = 8

! 文字列のfortran文字列長
!    その他
   integer, parameter :: PL_STR_LEN = 32

!----------------------------------------------------
! ファイルフォーマット
!
   character(len=PL_FORMAT_LEN), parameter :: FILE_FMT_STL_A  = 'stl_a'
   character(len=PL_FORMAT_LEN), parameter :: FILE_FMT_STL_AA = 'stl_aa'
   character(len=PL_FORMAT_LEN), parameter :: FILE_FMT_STL_B  = 'stl_b'
   character(len=PL_FORMAT_LEN), parameter :: FILE_FMT_STL_BB = 'stl_bb'
   character(len=PL_FORMAT_LEN), parameter :: FILE_FMT_NPT_A  = 'npt_a'
   character(len=PL_FORMAT_LEN), parameter :: FILE_FMT_NPT_B  = 'npt_b'

!----------------------------------------------------
