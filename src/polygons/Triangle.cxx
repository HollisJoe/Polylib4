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

//  <string> for windows, <string.h> for linux
#include <string>
#include <string.h>
#include <stdio.h>

#include "polygons/Triangle.h"


using namespace std;
using namespace PolylibNS;


/************************************************************************
 *
 * Triangleクラス
 *
 ***********************************************************************/

/// コンストラクタ (deserialize)
Triangle::Triangle(
            const char* pbuff
        )
{
    m_AtrI = NULL;
    m_AtrR = NULL;

    char* p = const_cast<char*>(pbuff);

    int* pAtrI;
    PL_REAL* pAtrR;

    // ポリゴンタイプ取得 確認用
    int pl_type;
    memcpy( &pl_type, p, sizeof(int) );
    p += sizeof(int);
    //PL_DBGOSH << "   pl_type="<<pl_type << endl;

    // ポリゴンID取得
    memcpy( &m_id, p, sizeof(long long int) );
    p += sizeof(long long int);
    //PL_DBGOSH << "   id="<<id << endl;

    // ３頂点座標
    PL_REAL vtx[9];
    memcpy( vtx, p, 9*sizeof(PL_REAL) );
    REAL9_TO_VEC3_3( vtx, m_vertex );
    p += 9*sizeof(PL_REAL);
    //PL_DBGOSH << "   vertex[0]= "<<vertex[0].x<<" "<<vertex[0].y<<" "<<vertex[0].z << endl;
    //PL_DBGOSH << "   vertex[1]= "<<vertex[1].x<<" "<<vertex[1].y<<" "<<vertex[1].z << endl;
    //PL_DBGOSH << "   vertex[2]= "<<vertex[2].x<<" "<<vertex[2].y<<" "<<vertex[2].z << endl;

    // ユーザ定義ID
    memcpy( &m_exid, p, sizeof(short int) );
    p += sizeof(short int);
    //PL_DBGOSH << "   exid= "<<exid << endl;

    // ユーザ定義属性
    memcpy( &m_numAtrI, p, sizeof(unsigned char) );     // 整数属性数
    p += sizeof(unsigned char);
    memcpy( &m_numAtrR, p, sizeof(unsigned char) );     // 実数属性数
    p += sizeof(unsigned char);
    //PL_DBGOSH << "   numAtrI= "<<numAtrI "   numAtrR= "<<numAtrR << endl;

    if( m_numAtrI > 0 ) {
        m_AtrI = (int*)malloc( m_numAtrI*sizeof(int) );
        memcpy( m_AtrI, p, m_numAtrI*sizeof(int) );     // 整数属性数
        p += m_numAtrI*sizeof(int);
    }

    if( m_numAtrR > 0 ) {
        m_AtrR = (PL_REAL*)malloc( m_numAtrR*sizeof(PL_REAL) );
        memcpy( m_AtrR, p, m_numAtrR*sizeof(PL_REAL) );     // 実数属性数
        p += m_numAtrR*sizeof(PL_REAL);
    }

    // 法線ベクトルと面積の更新
    calc_normal();
    calc_area();
}


// ポリゴンのシリアリズ
char* Triangle::serialize( const char* pbuff )
{
    char* p = const_cast<char*>(pbuff);

    // ポリゴンタイプを設定
    int pl_type = get_pl_type();
    memcpy( p, &pl_type, sizeof(int) );
    p += sizeof(int);

    // ポリゴンID
    memcpy( p, &m_id, sizeof(long long int) );
    p += sizeof(long long int);

    // ３頂点座標
    PL_REAL vtx[9];
    VEC3_3_TO_REAL9( m_vertex, vtx );
    memcpy( p, vtx, 9*sizeof(PL_REAL) );
    p += 9*sizeof(PL_REAL);

    // ユーザ定義ID
    memcpy( p, &m_exid, sizeof(short int) );
    p += sizeof(short int);

    // ユーザ定義属性
    memcpy( p, &m_numAtrI, sizeof(unsigned char) );     // 整数属性数
    p += sizeof(unsigned char);

    memcpy( p, &m_numAtrR, sizeof(unsigned char) );     // 実数属性数
    p += sizeof(unsigned char);

    if( m_numAtrI > 0 ) {
        memcpy( p, m_AtrI, m_numAtrI*sizeof(int) ); // 整数属性設定
        p += m_numAtrI*sizeof(int);
    }

    if( m_numAtrR > 0 ) {
        memcpy( p, m_AtrR, m_numAtrR*sizeof(PL_REAL) ); // 実数属性設定
        p += m_numAtrR*sizeof(PL_REAL);
    }

    return p;
}


///
/// Bounding box of this triangle
///
/// @param[in] detail       曲面補間可能要素の時に曲面補間したBounding boxを返すか否か
///                             false: ３角形の頂点にて決定
///                             長田パッチ用であるので、３角形では無視する
/// @return Bounding box
///
//   VTree関連で使用する　特に長田パッチのbbox取得用
BBox Triangle::get_bbox( bool detail )
{

    BBox bbox;
    bbox.init();
    bbox.add( m_vertex[0] );
    bbox.add( m_vertex[1] );
    bbox.add( m_vertex[2] );
    return bbox;
}


///
/// 頂点を設定
///
/// @param[in] vertex           三角形の3頂点
/// @param[in] update_param     法線ベクトルを再計算するか？
/// @param[in] calc_area        面積を再計算するか？
/// @attention  長田パッチの場合、update_param=trueで
///             長田パッチのパラメータ更新を更新する
///
///  仮想関数のため、ヘッダでのinline展開は止める
///
void Triangle::set_vertexes(
        const Vec3<PL_REAL> vertex[3], 
        bool    update_param, 
        bool    calc_area
     )
{
    m_vertex[0] = vertex[0];
    m_vertex[1] = vertex[1];
    m_vertex[2] = vertex[2];
    if(update_param) this->calc_normal();
    if(calc_area) this->calc_area();
}



// protected //////////////////////////////////////////////////////////////////
// システムで一意のポリゴンIDを作成する
long long int Triangle::create_unique_id() 
{
    static long long int unique_id = 0;
    return ++unique_id;
}

