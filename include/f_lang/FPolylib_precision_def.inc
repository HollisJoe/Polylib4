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
!    プリプロセッサで実数の型を切り替える時に使用する
!    (Polylib内部はC++で実装されており当インクルードファイルは使用していない)
!
!    使用方法：
!         #include "FPolylib_precision_def.inc"
!
!         ※　include 'FPolylib_precision_def.inc' ではコンパイルエラーとなるので注意
!         ※　プリプロセッサで切り替える必要がなければ
!             include 'FPolylib_precision.inc' を使用する
!     
!
!*********************************************************

!----------------------------------------------------
! Fortran実数型の指定
!     PL_REAL_PN は小文字に変換されていることを前提としている
!----------------------------------------------------

#ifdef _REAL_IS_DOUBLE_
  integer,  parameter :: pl_real_pn = 8
#else
  integer,  parameter :: pl_real_pn = 4
#endif

! [備考]
!   当インクルードファイルを使わずソースに直接、記述する場合は
!   以下のように大文字を使うこと
!
!   #ifdef _REAL_IS_DOUBLE_
!     integer,  parameter :: PL_REAL_PN = 8
!   #else
!     integer,  parameter :: PL_REAL_PN = 4
!   #endif


