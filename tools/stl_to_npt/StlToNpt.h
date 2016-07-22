/* -- Mode: c++ --*/
/*
* Polylib - Polygon Management Library
*
* Copyright (c) 2015-2016 Advanced Institute for Computational Science, RIKEN.
* All rights reserved.
*
*/

#ifndef _STL_TO_NPT_H_
#define _STL_TO_NPT_H_

#include "Polylib.h"

using namespace std;
using namespace PolylibNS;

////////////////////////////////////////////////////////////////////////////
///
/// プログラム:StlToNpt
/// STL→長田パッチファイル変換プログラム
///
////////////////////////////////////////////////////////////////////////////

//----------------------------------------------------
//  プロトタイプ宣言
//----------------------------------------------------

///
/// 各ランクにあるポリゴン情報を集約する
///     ポリゴンの重複を削除する
///     長田パッチを生成する
///
/// @param [in]    num_rank         ランク数
/// @param [in]    myrank           自身のランクID
/// @param [in]    num_tri_rank     各ランク内の３角形の数
/// @param [in]    tri_list_rank    各ランク内の３角形情報
///                                    IDの重複あり
/// @param [in]    p1_vertex_norm   頂点1の法線ベクトル
/// @param [in]    p1_setted_flag   頂点1の法線ベクトルの設定フラグ
/// @param [in]    p2_vertex_norm   頂点2の法線ベクトル
/// @param [in]    p2_setted_flag   頂点2の法線ベクトルの設定フラグ
/// @param [in]    p3_vertex_norm   頂点3の法線ベクトル
/// @param [in]    p3_setted_flag   頂点3の法線ベクトルの設定フラグ
/// @param [in]    num_npt_alloc    長田パッチポリゴン領域獲得数
///                                   （出力長田パッチ数と同一のはず）
/// @param [out]   npt_list         長田パッチポリゴン情報
/// @return リターンコード   =0 正常  !=0 異常
/// @attention
///     出力長田パッチ数とnum_npt_allocが合致しない場合はエラーとする
///     NptTriangleオブジェクトは内部で生成される（使用後deleteが必要）

int
gather_polygons_to_npt(
            int num_rank,
            int myrank,
            int num_tri_rank,
            std::vector<Triangle*> *tri_list_rank,
            PL_REAL** p1_vertex_norm,
            int*      p1_setted_flag,
            PL_REAL** p2_vertex_norm,
            int*      p2_setted_flag,
            PL_REAL** p3_vertex_norm,     
            int*      p3_setted_flag,
            int num_npt_alloc,
            std::vector<NptTriangle*>&  npt_list
        );


///
/// 頂点ベクトル決定
///
/// @param [in]    pg           ポリゴングループ
/// @param [in]    min          担当領域Min座標
/// @param [in]    max          担当領域Max座標
/// @param [in]    id           自身の３角形のinternal_id
/// @param [in]    norm         自身の法線ベクトル
/// @param [in]    pos          頂点座標
/// @param [in]    edge_degree  エッジ判定角度
/// @param [out]   vertex_norm  頂点法線ベクトル（単位ベクトル）
/// @param [out]   isetted_flg  頂点設定フラグ
///                                =0 設定なし
///                                =1 設定あり
//cc///                                =0 自身の３角形以外で設定なし
//cc///                                =1 自身の３角形以外で設定あり
/// @return リターンコード   =0 正常  !=0 異常
/// @attention
///     
int
get_vertex_normal (
        PolygonGroup& pg,
        PL_REAL  min[3],
        PL_REAL  max[3],
        long long int id,
        PL_REAL  norm[3],
        PL_REAL  pos[3],
        PL_REAL  edge_degree,
        PL_REAL  vertex_norm[3],
        int& isetted_flg
    );

#endif // _STL_TO_NPT_H_
