/* -- Mode: c++ --*/
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

#ifndef polylib_func_h
#define polylib_func_h
#include <string>
#include <vector>
#include <iostream>
#include "polygons/Triangle.h"
#include "polygons/NptTriangle.h"
#include "groups/PolygonGroup.h"
#include "common/PolylibStat.h"
#include "common/PolylibCommon.h"
#include "common/BBox.h"
#include "common/Vec3.h"

#define INLINE inline

using namespace Vec3class;

namespace PolylibNS {

///
/// ポリゴンのデシリアリズ
///   シリアライズされたデータよりオブジェクトの生成を行う
///
/// @param[in]  pl_type  ポリゴンタイプ(PL_TYPE_TRIANGLE/PL_TYPE_NPT)
/// @param[in]  pbuff    バッファ格納位置先頭ポインタ
///                        シリアリズされたデータ
/// @return 生成したオブジェクト(Triangle/NptTriangle)
///
INLINE Triangle* deserialize_polygon(
           int         pl_type,
           const char* pbuff
         )
{
    if( pl_type == PL_TYPE_NPT ) {
        NptTriangle* pTri = new NptTriangle( pbuff );  // deserialize
        return pTri;
    } else {   // PL_TYPE_TRIANGLE
        Triangle* pTri = new Triangle( pbuff );  // deserialize
        return pTri;
    }
}

///
/// ポリゴンの複製
///   ポリゴンの種別(Triangle/NptTrinangle)を判別し、適正なポリゴンを生成する
///
///  @param[in]  tri       複製元ポリゴン
///  @param[out] copy_tri  複製ポリゴン
///  @attention  関数内でアロケーションするので、copy_tri使用後deleteしてください
///
INLINE void  copy_polygon(
           Triangle* tri,
           Triangle* &copy_tri
         )
{
    int pl_type = tri->get_pl_type();
    if( pl_type == PL_TYPE_NPT ) {
       // NptTriangle
        NptTriangle* pNpt = dynamic_cast<NptTriangle*>(tri);
        copy_tri = new NptTriangle(*pNpt);
    } else {
        // Triangle
        copy_tri = new Triangle(*tri);
    }
}

///
/// ポリゴン（複数）の複製・追加
///   ポリゴンの種別(Triangle/NptTrinangle)を判別し、適正なポリゴンを生成する
///
///  @param[in]  tri_list          複製元ポリゴン（複数）
///  @param[in/out] copy_tri_list  複製ポリゴン  （複数）
///                                   複製したものが追加される   
///  @attention  関数内でアロケーションするので、
///       使用後 copy_trias内のポリゴンはdeleteしてください
///       ディープコピーしています。
///
INLINE void  copy_polygons(
           const std::vector<Triangle*>& tri_list,
           std::vector<Triangle*>& copy_tri_list
         )
{
    if( tri_list.size()>0 ) {
        copy_tri_list.reserve( copy_tri_list.size() + tri_list.size() );
        int pl_type = tri_list[0]->get_pl_type();
        if( pl_type == PL_TYPE_NPT ) {
            // NptTriangle
            for(int i=0; i<tri_list.size(); i++ ) {
                NptTriangle* pNpt = dynamic_cast<NptTriangle*>(tri_list[i]);
                copy_tri_list.push_back( new NptTriangle(*pNpt) );
            }
        } else {
            // Triangle
            for(int i=0; i<tri_list.size(); i++ ) {
                copy_tri_list.push_back( new Triangle( *(tri_list[i]) ) );
            }
        }
    }
}

///
/// ポリゴンvectorの型変換
///   vector<Triangle*> から vector<NptTriangle*> に変換する
///
///  @param[in]  tri_list   Triangleのリスト
///  @param[out] npt_list   NptTriangleのリスト
///  @return POLYLIB_STATで定義される値が返る
///
INLINE POLYLIB_STAT  convert_polygons_to_npt (
           std::vector<Triangle*>&    tri_list,
           std::vector<NptTriangle*>& npt_list
         )
{
    if( tri_list.size() == 0 ) {
        return PLSTAT_OK;
    }

    int pl_type = tri_list[0]->get_pl_type();
    if( pl_type != PL_TYPE_NPT ) {
        return PLSTAT_NOT_NPT;
    } 

    // Triangle* -> NptTriangle* 変換
    for( int i=0; i<tri_list.size(); i++ ) {
        NptTriangle* pNpt = dynamic_cast<NptTriangle*>(tri_list[i]);
        if( pNpt == 0 ) {
            return PLSTAT_NOT_NPT;
        }
        npt_list.push_back( pNpt );
    }

    return PLSTAT_OK;
}

///
/// ポリゴンvectorの型変換
///   vector<NptTriangle*> から vector<Triangle*> に変換する
///
///  @param[in]   npt_list   NptTriangleのリスト
///  @param[out]  tri_list   Triangleのリスト
///      
///
INLINE void convert_polygons_to_tri (
           std::vector<NptTriangle*>& npt_list,
           std::vector<Triangle*>&    tri_list
         )
{
    // NptTriangle* -> Triangle* 変換
    //    Triangleはベースクラスなのでそのまま変換される
    for( int i=0; i<npt_list.size(); i++ ) {
        tri_list.push_back( npt_list[i] );
    }
}

///
/// ２次元配列（整数）の領域確保
///   領域解放は free_array_2d() を使用すること
///
///  @param[in]   n1         1次元目サイズ
///  @param[in]   n2         2次元目サイズ
///  @return ２次元配列のポインタ
///  @attention  x[i][j]でアクセスする
///
INLINE int** alloc_array_2d_int ( int n1, int n2 )
{
    int*  data = (int* )malloc(n1*n2*sizeof(int));
    int** x    = (int**)malloc(n1*sizeof(int*));
    for( int i=0; i<n1; i++ ) {
        x[i] = &(data[i*n2]);
    }
    return x;
}

///
/// ２次元配列（実数）の領域確保
///   領域解放は free_array_2d() を使用すること
///
///  @param[in]   n1         1次元目サイズ
///  @param[in]   n2         2次元目サイズ
///  @return ２次元配列のポインタ
///  @attention  x[i][j]でアクセスする
///
INLINE PL_REAL** alloc_array_2d_real ( int n1, int n2 )
{
    PL_REAL*  data = (PL_REAL* )malloc(n1*n2*sizeof(PL_REAL));
    PL_REAL** x    = (PL_REAL**)malloc(n1*sizeof(PL_REAL*));
    for( int i=0; i<n1; i++ ) {
        x[i] = &(data[i*n2]);
    }
    return x;
}

///
/// ２次元配列の領域解放
///   alloc_array_2d_*() でアロケーションした領域の解放
///
///  @param[in]   x     ２次元配列のポインタ
///  @return 戻り値なし
///
INLINE void free_array_2d ( void** x )
{
    if( x == NULL ) return;
    free( x[0] ); 
    free( x ); 
}


} //namespace PolylibNS

#endif // polylib_func_h

