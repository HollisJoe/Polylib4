c/*
c* Polylib - Polygon Management Library
c*
c* Copyright (c) 2010-2011 VCAD System Research Program, RIKEN.
c* All rights reserved.
c*
c* Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
c* All rights reserved.
c*
c*/

c*********************************************************
c
c Fortranインクルードファイル
c    文字列長の定義
c    Fortran77ユーザアプリ用のインクルードファイル
c    (Polylib内部はC++で実装されており当インクルードファイルは使用していない)
c
c*********************************************************

c----------------------------------------------------
c タグの精度　実体はアドレスなので8byte
c   integer*8 tag // fort77

      integer PL_TAG_PN
      parameter( PL_TAG_PN = 8 )


c----------------------------------------------------
c Fortran文字列長の定義
c     FPolylib.h と合致していること
c
c    (Fortranアプリ使用方法)
c
c      character*NNN str   // fort77   NNNは文字列の長さ
c

c ファイルパス、ファイル名のfortran文字列長
      integer PL_FILE_PATH_LEN, PL_FILE_NAME_LEN
      parameter( PL_FILE_PATH_LEN = 256, PL_FILE_NAME_LEN =  64 )

c ポリゴングループ名のfortran文字列長
      integer PL_GRP_PATH_LEN, PL_GRP_NAME_LEN
      parameter( PL_GRP_PATH_LEN = 256, PL_GRP_NAME_LEN =  64 )

c 属性用文字列のfortran文字列長
c     属性キー, 属性値 etc.
      integer PL_GRP_ATR_LEN
      parameter( PL_GRP_ATR_LEN = 32 )

c ファイルフォーマットのfortran文字列長
      integer PL_FORMAT_LEN
      parameter( PL_FORMAT_LEN = 8 )

c 文字列のfortran文字列長
c    その他
      integer PL_STR_LEN
      parameter( PL_STR_LEN = 32 )

!----------------------------------------------------
! ファイルフォーマット
!
      character*8 FILE_FMT_STL_A
      parameter ( FILE_FMT_STL_A   = 'stl_a' )
      character*8 FILE_FMT_STL_AA
      parameter ( FILE_FMT_STL_AA  = 'stl_aa' )
      character*8 FILE_FMT_STL_B
      parameter ( FILE_FMT_STL_B   = 'stl_b' )
      character*8 FILE_FMT_STL_BB
      parameter ( FILE_FMT_STL_BB  = 'stl_bb' )
      character*8 FILE_FMT_NPT_A
      parameter ( FILE_FMT_NPT_A   = 'npt_a' )
      character*8 FILE_FMT_NPT_B
      parameter ( FILE_FMT_NPT_B   = 'npt_b' )

c----------------------------------------------------
