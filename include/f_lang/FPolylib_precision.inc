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
!    実数の型
!    Fortranユーザアプリ用のインクルードファイル
!    (Polylib内部はC++で実装されており当インクルードファイルは使用していない)
!
!    Fortranの include 'FPolylib_precision.inc' では
!    内部で#ifdefで切り替えられない。
!      -> #ifdefで切り替える用として、FPolylib_precision_def.inc を別途提供する
!    　   #include "FPolylib_precision_def.inc" とソースに記述する 　
!     
!
!*********************************************************

!----------------------------------------------------
! Fortran実数型の指定
!----------------------------------------------------

!  実数の精度は以下のパラメータ文を切り替えることにより行うものとする
!  コンパイルオプションでは制御しない

  integer,  parameter :: DP      = 8    ! 倍精度
  integer,  parameter :: SP      = 4    ! 単精度

  integer,  parameter :: PL_REAL_PN   = SP       ! ソースコード　精度
  ! integer,  parameter :: PL_REAL_PN    = DP         ! ソースコード　精度

!
!    (Fortranアプリ使用方法)
!
!      real(kind=PL_REAL_PN) :: x
!        or
!      real(PL_REAL_PN) :: x
!        or
!      real*4  x       ! fortran77等で明に記載する場合
!      real*8  x
!
