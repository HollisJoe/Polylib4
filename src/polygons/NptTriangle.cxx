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

#include "polygons/NptTriangle.h"
#ifdef USE_NPATCH_LIB
#include "Npt.h"
#endif


using namespace std;
using namespace PolylibNS;

//  NpatchParam p  -> PL_REAL r[21]
#define NPATCH_PARAM_TO_REAL21(p,r) \
            {   \
                (r)[ 0]=(p).cp_side1_1.x;  (r)[ 1]=(p).cp_side1_1.y;  (r)[ 2]=(p).cp_side1_1.z; \
                (r)[ 3]=(p).cp_side1_2.x;  (r)[ 4]=(p).cp_side1_2.y;  (r)[ 5]=(p).cp_side1_2.z; \
                (r)[ 6]=(p).cp_side2_1.x;  (r)[ 7]=(p).cp_side2_1.y;  (r)[ 8]=(p).cp_side2_1.z; \
                (r)[ 9]=(p).cp_side2_2.x;  (r)[10]=(p).cp_side2_2.y;  (r)[11]=(p).cp_side2_2.z; \
                (r)[12]=(p).cp_side3_1.x;  (r)[13]=(p).cp_side3_1.y;  (r)[14]=(p).cp_side3_1.z; \
                (r)[15]=(p).cp_side3_2.x;  (r)[16]=(p).cp_side3_2.y;  (r)[17]=(p).cp_side3_2.z; \
                (r)[18]=(p).cp_center.x;   (r)[19]=(p).cp_center.y;   (r)[20]=(p).cp_center.z;  \
            }


//  PL_REAL r[21]  ->  NpatchParam p
#define REAL21_TO_NPATCH_PARAM(r,p) \
            {   \
                (p).cp_side1_1.x=(r)[ 0];  (p).cp_side1_1.y=(r)[ 1];  (p).cp_side1_1.z=(r)[ 2]; \
                (p).cp_side1_2.x=(r)[ 3];  (p).cp_side1_2.y=(r)[ 4];  (p).cp_side1_2.z=(r)[ 5]; \
                (p).cp_side2_1.x=(r)[ 6];  (p).cp_side2_1.y=(r)[ 7];  (p).cp_side2_1.z=(r)[ 8]; \
                (p).cp_side2_2.x=(r)[ 9];  (p).cp_side2_2.y=(r)[10];  (p).cp_side2_2.z=(r)[11]; \
                (p).cp_side3_1.x=(r)[12];  (p).cp_side3_1.y=(r)[13];  (p).cp_side3_1.z=(r)[14]; \
                (p).cp_side3_2.x=(r)[15];  (p).cp_side3_2.y=(r)[16];  (p).cp_side3_2.z=(r)[17]; \
                (p).cp_center.x =(r)[18];  (p).cp_center.y =(r)[19];  (p).cp_center.z =(r)[20]; \
            }


/************************************************************************
 *
 * NptTriangleクラス
 *
 ***********************************************************************/

/// コンストラクタ (deserialize)
NptTriangle::NptTriangle(
            const char* pbuff
      ) : Triangle ( pbuff )
{
    char* p = const_cast<char*>(pbuff);

    size_t size = Triangle::serialized_size();

    p += size;

    // 長田パッチパラメータ
    PL_REAL   cp[21];
    memcpy( cp, p, 21*sizeof(PL_REAL) );
    REAL21_TO_NPATCH_PARAM( cp, m_npatchParam );
    p += 21*sizeof(PL_REAL);
}


// ポリゴンのシリアリズ
char* NptTriangle::serialize( const char* pbuff )
{
    char* p = Triangle::serialize( pbuff );

    // 長田パッチ
    PL_REAL   cp[21]; 
    NPATCH_PARAM_TO_REAL21( m_npatchParam, cp );
    memcpy( p, cp, 21*sizeof(PL_REAL) );
    p += 21*sizeof(PL_REAL);

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
BBox NptTriangle::get_bbox( bool detail )
{
    if( detail ) {
        // 曲面補正
        BBox bbox;
        bbox.init();
        // 3頂点のMinMax
        bbox.add( m_vertex[0] );
        bbox.add( m_vertex[1] );
        bbox.add( m_vertex[2] );
        // 補正点のMinMax
        Vec3<PL_REAL> pos_c;    // 補正点座標
            // 辺1の中点 eta  = 0.5; xi = 0.0; の補正点
        correct( 0.5, 0.0, pos_c );  bbox.add( pos_c );
            // 辺2の中点 eta  = 1.0; xi = 0.5; の補正点
        correct( 1.0, 0.5, pos_c );  bbox.add( pos_c );
            // 辺3の中点 eta  = 0.5; xi = 0.5; の補正点
        correct( 0.5, 0.5, pos_c );  bbox.add( pos_c );
            // 3角形の重心   eta  = 2.0/3.0, xi = 0.5*2.0/3.0; の補正点
        correct( (2.0/3.0), (0.5*2.0/3.0), pos_c );  bbox.add( pos_c );
           
        return bbox;
    
    } else {
        return Triangle::get_bbox( detail );
    }
}


#ifdef USE_NPATCH_LIB
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
void NptTriangle::set_vertexes(
        const Vec3<PL_REAL> vertex[3], 
        bool    update_param, 
        bool    calc_area
     )
{
    if(update_param) {
        // 長田パッチパラメータ（制御点）の更新
            // 旧座標
        PL_REAL p1[3], p2[3], p3[3];
            // 旧長田パッチパラメータ
        PL_REAL side1_1[3], side1_2[3], side2_1[3], side2_2[3], side3_1[3], side3_2[3], center[3];
            // 新座標
        PL_REAL p1_n[3], p2_n[3], p3_n[3];
            // 新長田パッチパラメータ
        PL_REAL side1_1_n[3], side1_2_n[3], side2_1_n[3], side2_2_n[3], side3_1_n[3], side3_2_n[3], center_n[3];
        
        p1[0] = m_vertex[0].x;  p1[1] = m_vertex[0].y;  p1[2] = m_vertex[0].z;
        p2[0] = m_vertex[1].x;  p2[1] = m_vertex[1].y;  p2[2] = m_vertex[1].z;
        p3[0] = m_vertex[2].x;  p3[1] = m_vertex[2].y;  p3[2] = m_vertex[2].z;

        p1_n[0] = vertex[0].x;  p1_n[1] = vertex[0].y;  p1_n[2] = vertex[0].z;
        p2_n[0] = vertex[1].x;  p2_n[1] = vertex[1].y;  p2_n[2] = vertex[1].z;
        p3_n[0] = vertex[2].x;  p3_n[1] = vertex[2].y;  p3_n[2] = vertex[2].z;
        
        Vec3<PL_REAL> *pside1_1 = &m_npatchParam.cp_side1_1;
        Vec3<PL_REAL> *pside1_2 = &m_npatchParam.cp_side1_2;
        Vec3<PL_REAL> *pside2_1 = &m_npatchParam.cp_side2_1;
        Vec3<PL_REAL> *pside2_2 = &m_npatchParam.cp_side2_2;
        Vec3<PL_REAL> *pside3_1 = &m_npatchParam.cp_side3_1;
        Vec3<PL_REAL> *pside3_2 = &m_npatchParam.cp_side3_2;
        Vec3<PL_REAL> *pcenter  = &m_npatchParam.cp_center;
        side1_1[0] = pside1_1->x;  side1_1[1] = pside1_1->y;  side1_1[2] = pside1_1->z;
        side1_2[0] = pside1_2->x;  side1_2[1] = pside1_2->y;  side1_2[2] = pside2_1->z;
        side2_1[0] = pside2_1->x;  side2_1[1] = pside2_1->y;  side2_1[2] = pside2_1->z;
        side2_2[0] = pside2_2->x;  side2_2[1] = pside2_2->y;  side2_2[2] = pside2_2->z;
        side3_1[0] = pside3_1->x;  side3_1[1] = pside3_1->y;  side3_1[2] = pside3_1->z;
        side3_2[0] = pside3_2->x;  side3_2[1] = pside3_2->y;  side3_2[2] = pside3_2->z;
        center [0] = pcenter->x;   center[1]  = pcenter->y;   center[2]  = pcenter->z;
        
        // 頂点移動後の長田パッチパラメータ取得
        npt_move_vertex (
                p1,   p2,   p3,   side1_1,   side1_2,   side2_1,   side2_2,   side3_1,   side3_2,   center,
                p1_n, p2_n, p3_n, side1_1_n, side1_2_n, side2_1_n, side2_2_n, side3_1_n, side3_2_n, center_n
            );
        
        pside1_1->x = side1_1_n[0];  pside1_1->y = side1_1_n[1];  pside1_1->z = side1_1_n[2];
        pside1_2->x = side1_2_n[0];  pside1_2->y = side1_2_n[1];  pside1_2->z = side1_2_n[2];
        pside2_1->x = side2_1_n[0];  pside2_1->y = side2_1_n[1];  pside2_1->z = side2_1_n[2];
        pside2_2->x = side2_2_n[0];  pside2_2->y = side2_2_n[1];  pside2_2->z = side2_2_n[2];
        pside3_1->x = side3_1_n[0];  pside3_1->y = side3_1_n[1];  pside3_1->z = side3_1_n[2];
        pside3_2->x = side3_2_n[0];  pside3_2->y = side3_2_n[1];  pside3_2->z = side3_2_n[2];
        pcenter->x  = center_n[0];   pcenter->y  = center_n[1];   pcenter->z  = center_n[2];
    }
    
    m_vertex[0] = vertex[0];
    m_vertex[1] = vertex[1];
    m_vertex[2] = vertex[2];
    if(update_param) this->calc_normal();
    if(calc_area) this->calc_area();
}
#endif

///
/// 頂点・長田パッチパラメータを設定
///
/// @param[in] vertex           三角形の3頂点
/// @param[in] param            長田パッチパラメータ
/// @param[in] update_param     法線ベクトルを再計算するか？
/// @param[in] calc_area        面積を再計算するか？
/// @attention  長田パッチの場合、update_param=trueで
///             長田パッチのパラメータ更新を更新する
///
///  仮想関数のため、ヘッダでのinline展開は止める
///
void NptTriangle::set_vertexes(
        const Vec3<PL_REAL> vertex[3], 
        const NpatchParam&  param,
        bool    update_param, 
        bool    calc_area
     )
{
    m_vertex[0] = vertex[0];
    m_vertex[1] = vertex[1];
    m_vertex[2] = vertex[2];
    m_npatchParam = param;
    if(update_param) this->calc_normal();
    if(calc_area) this->calc_area();
}


#ifdef USE_NPATCH_LIB
///
/// 点の近似曲面補正
///
///  @param[in]     pos     3角形平面内の点
///  @param[out]    pos_o   曲面上に補正された座標
///  @return なし
///
void  NptTriangle::correct( 
           const Vec3<PL_REAL>&   pos,
           Vec3<PL_REAL>&         pos_o
         )
{
    PL_REAL  pos_wk[3];
    PL_REAL  p1_wk[3];
    PL_REAL  p2_wk[3];
    PL_REAL  p3_wk[3];
    PL_REAL  eta;
    PL_REAL  xi;
    
    pos_wk[0] = pos[0];  pos_wk[1] = pos[1];  pos_wk[2] = pos[2];
    p1_wk[0] = m_vertex[0].x;  p1_wk[1] = m_vertex[0].y;  p1_wk[2] = m_vertex[0].z;
    p2_wk[0] = m_vertex[1].x;  p2_wk[1] = m_vertex[1].y;  p2_wk[2] = m_vertex[1].z;
    p3_wk[0] = m_vertex[2].x;  p3_wk[1] = m_vertex[2].y;  p3_wk[2] = m_vertex[2].z;

    // 入力座標 -> η、ξパラメータ変換
    npt_cvt_pos_to_eta_xi( 
                pos_wk, p1_wk, p2_wk, p3_wk,
                &eta, &xi
            );
    
    // 点の近似曲面補正
    correct( eta, xi, pos_o );

}
#endif

///
/// 点の近似曲面補正
///
///  @param[in]     eta     ηパラメータ
///  @param[in]     xi      ξパラメータ
///  @param[out]    pos_o   曲面上に補正された座標
///  @return なし
///  @attention
///     ηとξのパラメータで３角形上の座標が決まる
///     (参考)
///         頂点1         eta  = 0.0; xi = 0.0;
///         頂点2         eta  = 1.0; xi = 0.0;
///         頂点3         eta  = 1.0; xi = 1.0;
///         辺1の中点     eta  = 0.5; xi = 0.0;
///         辺2の中点     eta  = 1.0; xi = 0.5;
///         辺3の中点     eta  = 0.5; xi = 0.5;
///         3角形の重心   eta  = 2.0/3.0, xi = 0.5*2.0/3.0;
///
void  NptTriangle::correct( 
           PL_REAL          eta,
           PL_REAL          xi,
           Vec3<PL_REAL>&   pos_o
         )
{

    Vec3<PL_REAL> *side1_1 = &m_npatchParam.cp_side1_1;
    Vec3<PL_REAL> *side1_2 = &m_npatchParam.cp_side1_2;
    Vec3<PL_REAL> *side2_1 = &m_npatchParam.cp_side2_1;
    Vec3<PL_REAL> *side2_2 = &m_npatchParam.cp_side2_2;
    Vec3<PL_REAL> *side3_1 = &m_npatchParam.cp_side3_1;
    Vec3<PL_REAL> *side3_2 = &m_npatchParam.cp_side3_2;
    Vec3<PL_REAL> *center  = &m_npatchParam.cp_center;


    // x(u,v,w) =    p1*w*w*w + cp_side1_1*3*u*w*w + cp_side1_2*3*u*u*w 
    //            +  p2*u*u*u + cp_side2_1*3*u*u*v + cp_side2_2*3*u*v*v
    //            +  p3*v*v*v + cp_side3_1*3*v*v*w + cp_dide3_2*3*v*w*w
    //            +  cp_center*6*u*v*w
    //
    //     u = eta - xi
    //     v = xi
    //     w = 1 - eta
    //     u + v + w = 1
    //
    PL_REAL u,v,w;

    u = eta - xi;
    v = xi;
    w = 1.0 - eta;

    pos_o.x =    m_vertex[0].x*w*w*w + side1_1->x*3.0*u*w*w + side1_2->x*3.0*u*u*w
               + m_vertex[1].x*u*u*u + side2_1->x*3.0*u*u*v + side2_2->x*3.0*u*v*v
               + m_vertex[2].x*v*v*v + side3_1->x*3.0*v*v*w + side3_2->x*3.0*v*w*w
               + center->x*6.0*u*v*w;

    pos_o.y =    m_vertex[0].y*w*w*w + side1_1->y*3.0*u*w*w + side1_2->y*3.0*u*u*w
               + m_vertex[1].y*u*u*u + side2_1->y*3.0*u*u*v + side2_2->y*3.0*u*v*v
               + m_vertex[2].y*v*v*v + side3_1->y*3.0*v*v*w + side3_2->y*3.0*v*w*w
               + center->y*6.0*u*v*w;

    pos_o.z =    m_vertex[0].z*w*w*w + side1_1->z*3.0*u*w*w + side1_2->z*3.0*u*u*w
               + m_vertex[1].z*u*u*u + side2_1->z*3.0*u*u*v + side2_2->z*3.0*u*v*v
               + m_vertex[2].z*v*v*v + side3_1->z*3.0*v*v*w + side3_2->z*3.0*v*w*w
               + center->z*6.0*u*v*w;
}

