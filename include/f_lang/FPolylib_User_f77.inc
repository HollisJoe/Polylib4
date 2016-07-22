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
c Fortran77インクルードファイル
c    Fortran77ユーザアプリ用のインクルードファイル
c    (Polylib内部はC++で実装されており当インクルードファイルは使用していない)
c    構造体は使用出来ない
c
c*********************************************************

c----------------------------------------------------
c Fortran実数型の指定
c     定義しても使う方法がないため設定しない
c----------------------------------------------------
cc      include 'FPolylib_precision_f77.inc'

c----------------------------------------------------
c Fortran文字列長
c----------------------------------------------------
      include 'FPolylib_define_f77.inc'

c----------------------------------------------------
c Fortran構造体
c     構造体は使用出来ない
c----------------------------------------------------
cc      include 'FPolylib_struct_f77_float.inc'

c----------------------------------------------------
c Fortranプロトタイプ
c     実数の型により切り替えること
c----------------------------------------------------
      include 'FPolylib_prototype_f77_float.inc'
c      include 'FPolylib_prototype_f77_double.inc'

      include 'FPolylib_prototype_f77_integer.inc'
