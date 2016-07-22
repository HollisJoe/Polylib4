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
!    構造体の定義
!    Fortranユーザアプリ用のインクルードファイル
!    Fortran77用はFpolylib_struct_f77.incを使用してください。
!    (Polylib内部はC++で実装されており当インクルードファイルは使用していない)
!
!*********************************************************

!
! 領域情報構造体(Fortran用)
!

type FParallelBboxStruct

    !cc sequence
    real(PL_REAL_PN)   bpos(3)
    integer            bbsize(3)
    integer            gcsize(3)
    real(PL_REAL_PN)   dx(3)

end type FParallelBboxStruct

!
! 長田パッチパラメータ（曲面補間用制御点情報）
!
type NpatchParamStruct

    real(PL_REAL_PN)   cp_side1_1(3)     !  p1p2辺の３次ベジェ制御点1
    real(PL_REAL_PN)   cp_side1_2(3)     !  p1p2辺の３次ベジェ制御点2
    real(PL_REAL_PN)   cp_side2_1(3)     !  p2p3辺の３次ベジェ制御点1
    real(PL_REAL_PN)   cp_side2_2(3)     !  p2p3辺の３次ベジェ制御点2
    real(PL_REAL_PN)   cp_side3_1(3)     !  p3p1辺の３次ベジェ制御点1
    real(PL_REAL_PN)   cp_side3_2(3)     !  p3p1辺の３次ベジェ制御点2
    real(PL_REAL_PN)   cp_center (3)     !  三角形中央の３次ベジェ制御点

end type NpatchParamStruct


!
! 構造体:PolylibMoveParamsStruct
! polylib_move()の引数として利用するパラメタセットの構造体です
!

type PolylibMoveParamsStruct

    ! 現在の計算ステップ番号
    integer m_current_step  

    ! 移動後の計算ステップ番号
    integer m_next_step

    ! １計算ステップあたりの時間変異
    real(PL_REAL_PN) m_delta_t

    ! ユーザ定義パラメータ
    real(PL_REAL_PN) m_params(10)

end type PolylibMoveParamsStruct

