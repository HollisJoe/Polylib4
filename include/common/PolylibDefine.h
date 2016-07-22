/*
* Polylib - Polygon Management Library
*
* Copyright (c) 2010-2011 VCAD System Research Program, RIKEN.
* All rights reserved.
*
* Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
* All rights reserved.
*
*/

#ifndef polylib_define_h
#define polylib_define_h

#include <limits.h>
#include <float.h>
#ifdef MPI
#include "mpi.h"
#endif

/** 実数型の指定
 * - デフォルトでは、PL_REAL=float
 * - コンパイル時オプション-D_REAL_IS_DOUBLE_を付与することで
 *   PL_REAL=doubleになる
 */
#ifdef _REAL_IS_DOUBLE_
#define PL_REAL double
#else
#define PL_REAL float
#endif

#ifdef MPI_PL
    #ifdef _REAL_IS_DOUBLE_
    #define PL_MPI_REAL    MPI_DOUBLE
    #else
    #define PL_MPI_REAL    MPI_FLOAT
    #endif
#endif

//
// 最大値／最小値
#define PL_INT_MAX   INT_MAX
#define PL_INT_MIN   INT_MIN
#define PL_REAL_MAX  FLT_MAX
                  // FLT_MINは正の最小値
                  // 負の最大値は -FLT_MAXである
#define PL_REAL_MIN  FLT_MIN

// PolygonGroup用のタグ
//      タグとしてm_internal_id を使う場合
//#define PL_GRP_TAG    int
//
//      タグとしてPolygonGroup用のポインタを使う場合
//          Fortranで使うことも考慮し、符号なしは使わない
//          PolygonGroup* ptr;
//          PL_GRP_TAG foo = reinterpret_cast<PL_GRP_TAG>(ptr); 
//      注意：vectorで管理されているため、
//              PolygonGroupの追加/削除があると無効になってしまう
#define PL_GRP_TAG   long long int 

// Triangle用のタグ
//      タグとしてTriangle用のポインタを使う
//      注意：vectorで管理されているため、
//              Triangleの追加/削除/ソートがあると無効になってしまう
//              並列実行時のmigrate処理で無効になることに注意
//              有効期間が短いことに注意
#define PL_ELM_TAG   long long int 

// 無効タグ
#define PL_NULL_TAG  0


// 形状（ポリゴン）のタイプ
#define PL_TYPE_UNKNOWN     0
        // 3角形(STL)
#define PL_TYPE_TRIANGLE    1
        // 長田パッチ
#define PL_TYPE_NPT         2


// 集合演算のタイプ
//    MPI_ReduceのMPI_SUM,MPI_MAX,MPI_MINに相当する
typedef enum {
    PL_OP_SUM = 1,    ///< 総和
    PL_OP_MAX = 2,    ///< MAX
    PL_OP_MIN = 3,    ///< MIN
} PL_OP_TYPE;


// 形状（ポリゴン）ファイルの形式
//      C言語用 （PolygonIOクラスに準拠）

                                //  STL アスキーファイル
#define FILE_FMT_STL_A  "stl_a"
                                //  STL アスキーファイル
#define FILE_FMT_STL_AA "stl_aa"
                                //  STL バイナリファイル
#define FILE_FMT_STL_B  "stl_b"
                                //  STL バイナリファイル
#define FILE_FMT_STL_BB "stl_bb"

                                //  長田パッチ アスキーファイル
#define FILE_FMT_NPT_A  "npt_a"
                                //  長田パッチ バイナリファイル
#define FILE_FMT_NPT_B  "npt_b"

                                //  デフォルトファイル形式
#define FILE_FMT_DEFAULT    FILE_FMT_STL_B

#endif // polylib_define_h
