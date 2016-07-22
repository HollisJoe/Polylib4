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
!    Fortran90ユーザアプリ用のインクルードファイル
!    (Polylib内部はC++で実装されており当インクルードファイルは使用していない)
!
!*********************************************************

!----------------------------------------------------
! Fortran実数型の指定
!    以下のようにユーザアプリ側で記載する
!----------------------------------------------------
!  include 'FPolylib_precision.inc'
!    or
!#include "FPolylib_precision_def.inc"

!----------------------------------------------------
! Fortran文字列長
!----------------------------------------------------
  include 'FPolylib_define.inc'

!----------------------------------------------------
! Fortran構造体
!----------------------------------------------------
  include 'FPolylib_struct.inc'

!----------------------------------------------------
! Fortranプロトタイプ
!----------------------------------------------------
  include 'FPolylib_prototype.inc'
