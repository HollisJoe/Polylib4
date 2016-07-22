#ifndef _CALC_GEO_MATRIX_H_
#define _CALC_GEO_MATRIX_H_

/*
 * Npatch - Nagata Patch Library
 *
 * Copyright (c) 2015-2016 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

////////////////////////////////////////////////////////////////////////////
///
/// マトリック演算
/// （行ベクトル系）
///
////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {  // for C++
#else
#endif


///
/// ４×１マトリクスと４×４マトリクスの積を求める
///
/// @param [in]   matA   ４×１マトリクス Ａ
/// @param [in]   matB   ４×４マトリクス Ｂ
/// @param [out]  matC   １×４マトリクス Ｃ
/// @return 終了コード   なし
///

INLINE void
Calc_3dMat4Multi14(
            GEO_REAL matA[4],     // [in]  ４×１マトリクス Ａ
            GEO_REAL matB[4][4],  // [in]  ４×４マトリクス Ｂ
            GEO_REAL matC[4]      // [out] １×４マトリクス Ｃ （演算結果）
        )
{
    int l_cnt_b;        /* 演算カウンタ  matB   */

    for ( l_cnt_b = 0; l_cnt_b < 4; l_cnt_b++ ) {
        matC[l_cnt_b] = ( matA[0] * matB[0][l_cnt_b] +
                          matA[1] * matB[1][l_cnt_b] +
                          matA[2] * matB[2][l_cnt_b] +
                          matA[3] * matB[3][l_cnt_b] );
    }
}


///
/// ４×４マトリクスと４×１マトリクスの積を求める
///
/// @param [in]   matA   ４×４マトリクス Ａ
/// @param [in]   matB   ４×１マトリクス Ｂ
/// @param [out]  matC   ４×１マトリクス Ｃ
/// @return 終了コード   なし
///
INLINE void
Calc_3dMat4Multi41( GEO_REAL matA[4][4], GEO_REAL matB[4], GEO_REAL matC[4] )
{
    int l_cnt_a;        /* 演算カウンタ　matA	*/

    for ( l_cnt_a = 0; l_cnt_a < 4; l_cnt_a++ ) {
        matC[l_cnt_a] =	( matA[l_cnt_a][0] * matB[0] +
                          matA[l_cnt_a][1] * matB[1] +
                          matA[l_cnt_a][2] * matB[2] +
                          matA[l_cnt_a][3] * matB[3] );
    }
}


//
/// ４×４マトリクスと４×４マトリクスの積を求める
///
/// @param [in]   matA   ４×４マトリクス Ａ
/// @param [in]   matB   ４×４マトリクス Ｂ
/// @param [out]  matC   ４×４マトリクス Ｃ (演算結果）
/// @return 終了コード   なし
///
INLINE void
Calc_3dMat4Multi44(
        GEO_REAL matA[4][4],
        GEO_REAL matB[4][4],
        GEO_REAL matC[4][4]
     )
{
    int l_cnt_a;        /* 演算カウンタ  matA   */
    int l_cnt_b;        /* 演算カウンタ  matB   */


    for ( l_cnt_a = 0; l_cnt_a < 4; l_cnt_a++ ) {
        for ( l_cnt_b = 0; l_cnt_b < 4; l_cnt_b++ ) {
            matC[l_cnt_a][l_cnt_b] = ( matA[l_cnt_a][0] * matB[0][l_cnt_b] +
                                       matA[l_cnt_a][1] * matB[1][l_cnt_b] +
                                       matA[l_cnt_a][2] * matB[2][l_cnt_b] +
                                       matA[l_cnt_a][3] * matB[3][l_cnt_b] );
        }
    }
}


///
/// 平行移動用の４×４変換マトリクスを求める
///
/// @param [in]   mov_x   Ｘ軸平行移動量
/// @param [in]   mov_y   Ｙ軸平行移動量
/// @param [in]   mov_z   Ｚ軸平行移動量
/// @param [out]  mat4      ４×４変換マトリクス（行ベクトル系）
/// @return 終了コード   なし
///

INLINE void
Calc_3dMat4Mov(
        GEO_REAL mov_x,
        GEO_REAL mov_y,
        GEO_REAL mov_z,
        GEO_REAL matA[4][4]
     )
{
    int l_cnt_c;

    for ( l_cnt_c = 0; l_cnt_c < 3; l_cnt_c++ ) {
        matA[l_cnt_c][0] = 0.0;
        matA[l_cnt_c][1] = 0.0;
        matA[l_cnt_c][2] = 0.0;
        matA[l_cnt_c][3] = 0.0;
    }

    matA[0][0] = 1.0;
    matA[1][1] = 1.0;
    matA[2][2] = 1.0;
    matA[3][3] = 1.0;

    /* 平行移動量セット	*/
    matA[3][0] = mov_x;
    matA[3][1] = mov_y;
    matA[3][2] = mov_z;
}


///
/// 座標軸回転用４×４変換マトリクスを求める
///
/// @param [in]   x_axis    新座標系のX軸となるベクトル
/// @param [in]   y_axis    新座標系のY軸となるベクトル
/// @param [in]   z_axis    新座標系のZ軸となるベクトル
/// @param [out]  matA      ４×４変換マトリクス（行ベクトル系）
/// @return 終了コード   なし
///

INLINE void
Calc_3dMat4RotAxis(
        GEO_REAL x_axis[3],
        GEO_REAL y_axis[3],
        GEO_REAL z_axis[3],
        GEO_REAL matA[4][4]
	)
{
	matA[0][0] = x_axis[0];
	matA[1][0] = x_axis[1];
	matA[2][0] = x_axis[2];
	matA[3][0] = 0.0;

	matA[0][1] = y_axis[0];
	matA[1][1] = y_axis[1];
	matA[2][1] = y_axis[2];
	matA[3][1] = 0.0;

	matA[0][2] = z_axis[0];
	matA[1][2] = z_axis[1];
	matA[2][2] = z_axis[2];
	matA[3][2] = 0.0;

	matA[0][3] = 0.0;
	matA[1][3] = 0.0;
	matA[2][3] = 0.0;
	matA[3][3] = 1.0;

}


///
/// 任意座標軸設定用４×４変換マトリクスを求める
///        （ワールド --> ローカル座標変換)
///
/// @param [in]   orig      新座標軸の原点となる点
/// @param [in]   x_axis    新座標系のX軸となるベクトル
/// @param [in]   y_axis    新座標系のY軸となるベクトル
/// @param [in]   z_axis    新座標系のZ軸となるベクトル
/// @param [out]  mat4      ４×４変換マトリクス（行ベクトル系）
/// @return 終了コード   なし
///

INLINE void
Calc_3dMat4TranAxis(
        GEO_REAL orig[3],
        GEO_REAL x_axis[3],
        GEO_REAL y_axis[3],
        GEO_REAL z_axis[3],
        GEO_REAL matA[4][4]
     )
{
	matA[0][0] = x_axis[0];
	matA[1][0] = x_axis[1];
	matA[2][0] = x_axis[2];
	matA[3][0] = -( x_axis[0]*orig[0] + x_axis[1]*orig[1] + x_axis[2]*orig[2] );

	matA[0][1] = y_axis[0];
	matA[1][1] = y_axis[1];
	matA[2][1] = y_axis[2];
	matA[3][1] = -( y_axis[0]*orig[0] + y_axis[1]*orig[1] + y_axis[2]*orig[2] );

	matA[0][2] = z_axis[0];
	matA[1][2] = z_axis[1];
	matA[2][2] = z_axis[2];
	matA[3][2] = -( z_axis[0]*orig[0] + z_axis[1]*orig[1] + z_axis[2]*orig[2] );

	matA[0][3] = 0.0;
	matA[1][3] = 0.0;
	matA[2][3] = 0.0;
	matA[3][3] = 1.0;
}


///
/// 任意座標軸設定用４×４変換マトリクスを求める
///        （ローカル --> ワールド座標変換)
///
/// @param [in]   orig      任意座標軸の原点となる点（ワールド座標系で指定）
/// @param [in]   x_axis    任意座標系のX軸となるベクトル（ワールド座標系で指定）
/// @param [in]   y_axis    任意座標系のY軸となるベクトル（ワールド座標系で指定）
/// @param [in]   z_axis    任意座標系のZ軸となるベクトル（ワールド座標系で指定）
/// @param [out]  mat4      ４×４変換マトリクス（行ベクトル系）
/// @return 終了コード   なし
///

INLINE void
Calc_3dMat4TranAxisInv(
        GEO_REAL orig[3],
        GEO_REAL x_axis[3],
        GEO_REAL y_axis[3],
        GEO_REAL z_axis[3],
        GEO_REAL matA[4][4]
     )
{
    GEO_REAL  axis_mat[4][4];     /* 回転用座標変換マトリクス */
    GEO_REAL  mov_mat[4][4];      /* 移動用座標変換マトリクス */
    GEO_REAL  ex_work;            /* 転置用ワーク             */

    /* 回転マトリクス */
    Calc_3dMat4RotAxis( x_axis, y_axis, z_axis, axis_mat );

    /* 回転マトリクスの回転部分（３×３）を転置させる */
    ex_work = axis_mat[0][1];   /* a(1,2) <--> a(2,1) */
    axis_mat[0][1] = axis_mat[1][0];
    axis_mat[1][0] = ex_work;

    ex_work = axis_mat[0][2];   /* a(1,3) <--> a(3,1) */
    axis_mat[0][2] = axis_mat[2][0];
    axis_mat[2][0] = ex_work;

    ex_work = axis_mat[1][2];   /* a(2,3) <--> a(3,2) */
    axis_mat[1][2] = axis_mat[2][1];
    axis_mat[2][1] = ex_work;

    /* 移動マトリクス */
    Calc_3dMat4Mov( orig[0], orig[1], orig[2], mov_mat );

    Calc_3dMat4Multi44( axis_mat, mov_mat, matA );
}


///
/// 回転用４×４変換マトリクスを求める（任意軸）
///
/// @param [in]   pos       回転軸上の１点の座標
/// @param [in]   vec       回転軸のベクトル
/// @param [in]   th_rad    回転角 (rad)
/// @param [out]  matA      回転用４×４変換マトリクス（行ベクトル系）
/// @return 終了コード   なし
///

INLINE void
Calc_3dMat4Rot2(
        GEO_REAL pos[3],    // [in]  回転軸上の１点の座標
        GEO_REAL vec[3],    // [in]  回転軸のベクトル
        GEO_REAL th_rad,    // [in]  回転角(rad)
        GEO_REAL matA[4][4] // [out] 回転用４×４変換マトリクス（行ベクトル系）
     )
{
    GEO_REAL  value_sin;          // sinθ 値
    GEO_REAL  value_cos;          // cosθ 値

    value_sin = sin( th_rad );
    value_cos = cos( th_rad );

    matA[0][0] = vec[0]*vec[0] + (1.0 - (vec[0]*vec[0]))*value_cos;
    matA[0][1] = vec[0]*vec[1]*(1.0 - value_cos) + vec[2]*value_sin;
    matA[0][2] = vec[0]*vec[2]*(1.0 - value_cos) - vec[1]*value_sin;
    matA[0][3] = 0.0;

    matA[1][0] = vec[0]*vec[1]*(1.0 - value_cos) - vec[2]*value_sin;
    matA[1][1] = vec[1]*vec[1] + (1.0 - vec[1]*vec[1] )*value_cos;
    matA[1][2] = vec[1]*vec[2]*(1.0 - value_cos) + vec[0]*value_sin;
    matA[1][3] = 0.0;

    matA[2][0] = vec[0]*vec[2]*(1.0 - value_cos) + vec[1]*value_sin;
    matA[2][1] = vec[1]*vec[2]*(1.0 - value_cos) - vec[0]*value_sin;
    matA[2][2] = vec[2]*vec[2] + (1.0 - vec[2]*vec[2])*value_cos;
    matA[2][3] = 0.0;

    matA[3][0] = pos[0]-pos[0]*matA[0][0]-pos[1]*matA[1][0]-pos[2]*matA[2][0];
    matA[3][1] = pos[1]-pos[0]*matA[0][1]-pos[1]*matA[1][1]-pos[2]*matA[2][1];
    matA[3][2] = pos[2]-pos[0]*matA[0][2]-pos[1]*matA[1][2]-pos[2]*matA[2][2];
    matA[3][3] = 1.0;
}

#ifdef __cplusplus
} // extern "C" or extern
#else
#endif


#endif // _CALC_GEO_MATRIX_H_
