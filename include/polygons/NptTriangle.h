// -*- Mode: c++ -*-
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

#ifndef polylib_npatchtriangle_h
#define polylib_npatchtriangle_h

#include "common/Vec3.h"
#include "polygons/Triangle.h"

using namespace Vec3class;

namespace PolylibNS{


////////////////////////////////////////////////////////////////////////////
///
/// クラス:NpatchParam
/// 長田パッチパラメータ（制御点情報）
///
////////////////////////////////////////////////////////////////////////////

struct NpatchParam {

    /// 長田パッチ曲面補間用制御点
    //      三角形の頂点座標と頂点の法線ベクトルより計算可能

    ///   p1->p2辺の３次ベジェ制御点1
    Vec3<PL_REAL>   cp_side1_1;

    ///   p1->p2辺の３次ベジェ制御点2
    Vec3<PL_REAL>   cp_side1_2;

    ///   p2->p3辺の３次ベジェ制御点1
    Vec3<PL_REAL>   cp_side2_1;
    
    ///  p2->p3辺の３次ベジェ制御点2
    Vec3<PL_REAL>   cp_side2_2;
    
    ///  p3->p1辺の３次ベジェ制御点1
    Vec3<PL_REAL>   cp_side3_1;
    
    ///  p3->p1辺の３次ベジェ制御点2
    Vec3<PL_REAL>   cp_side3_2;  

    ///  三角形中央の３次ベジェ制御点
    Vec3<PL_REAL>   cp_center;
};

////////////////////////////////////////////////////////////////////////////
///
/// クラス:NptTriangleクラス
///
////////////////////////////////////////////////////////////////////////////


class NptTriangle : public Triangle {
public:

    ///
    /// コンストラクタ
    ///
    NptTriangle()
    {
    }

    ///
    /// コピーコンストラクタ
    ///
    NptTriangle( const NptTriangle  &tria )
        : Triangle ( tria )
    {
        m_npatchParam = tria.m_npatchParam;
    }


    ///
    /// コンストラクタ
    ///
    /// @param[in] vertex   ポリゴンの頂点
    /// @param[in] param     長田パッチパラメータ
    /// @param[in] id        三角形ポリゴンID（ユニークな識別子）
    ///                          システム全体でユニークな識別子であること
    /// @param[in] num_atrI    ユーザ定義属性数（整数型）
    /// @param[in] num_atrR    ユーザ定義属性数（実数型）
    /// @param[in] atrI        ユーザ定義属性（整数型）
    /// @param[in] atrR        ユーザ定義属性（実数型）
    /// @attention
    ///     id=0の場合、IDは内部で採番される
    ///
    NptTriangle(
            const Vec3<PL_REAL> vertex[3],
            const NpatchParam&  param,
            long long int       id = 0,
            int                 num_atrI = 0,
            int                 num_atrR = 0,
            const int           *atrI    = NULL,
            const PL_REAL       *atrR    = NULL
        ) : Triangle ( vertex, id, num_atrI, num_atrR, atrI, atrR  )
    {
        m_npatchParam = param;
    }

    ///
    /// コンストラクタ (deserialize)
    ///    シリアライズされた通信バッファよりオブジェクトの生成を行う
    ///
    /// @param[in] pbuff    バッファ格納先頭位置ポインタ
    ///                        シリアリズされたデータ
    /// @attention シリアライズはserialize()で行う
    ///
    NptTriangle(
            const char* pbuff
        );

    ///
    /// デストラクタ
    ///
    virtual ~NptTriangle() {}

    ///
    /// ポリゴンのリスケール
    ///
    /// @param[in] scale        スケール
    /// @return 戻り値なし
    ///
    virtual void rescale( PL_REAL scale )
    {
        m_vertex[0] *= scale;
        m_vertex[1] *= scale;
        m_vertex[2] *= scale;

        m_npatchParam.cp_side1_1 *= scale;
        m_npatchParam.cp_side1_2 *= scale;
        m_npatchParam.cp_side2_1 *= scale;
        m_npatchParam.cp_side2_2 *= scale;
        m_npatchParam.cp_side3_1 *= scale;
        m_npatchParam.cp_side3_2 *= scale;
        m_npatchParam.cp_center  *= scale;
    } 

    ///
    /// ポリゴンの使用メモリサイズ取得
    ///    Polylib::used_memory_size()で使用する
    ///
    /// @return 使用メモリサイズ
    ///
    virtual size_t used_memory_size( void )
    {
        size_t size = 0;
        //size += sizeof(this);               // ユーザ定義属性を除くサイズ
        size += sizeof(NptTriangle);               // ユーザ定義属性を除くサイズ
        size += m_numAtrI*sizeof(int);      // ユーザ定義属性（整数）
        size += m_numAtrR*sizeof(PL_REAL);  // ユーザ定義属性（実数）

        return size;
    }


    ///
    /// ポリゴンのシリアリズした時のサイズ取得
    ///
    /// @return シリアリズ後のサイズ
    ///
    virtual size_t serialized_size( void )
    {
        size_t size = Triangle::serialized_size();
        size += 21*sizeof(PL_REAL);     // 長田パッチパラメータ

        return size;
    }

   ///
    /// ポリゴンのシリアリズ
    ///   １ポリゴンをバッファに格納する
    ///
    /// @param[out] pbuff  バッファ格納位置先頭ポインタ
    ///                        シリアリズされたデータ
    /// @return バッファ格納位置Nextポインタ
    /// @attention  デシリアリズはPolylibNS::deserialize_polygon()で行う
    ///     最終的にはコンストラクタ NptTriangle( pbuff ) が呼び出される
    ///
    virtual char* serialize( const char* pbuff );

    ///
    /// Bounding box of this triangle
    ///
    /// @param[in] detail       曲面補間したBounding boxを返すか否か
    ///                             false: ３角形の頂点にて決定
    ///                             true:  曲面補間して決定
    /// @return Bounding box
    ///
    // VTree関連で使用する　特に長田パッチのbbox取得用
    virtual BBox get_bbox( bool detail = false );

    //=======================================================================
    // Setter/Getter
    //=======================================================================

    ///
    /// ポリゴンタイプ取得
    ///
    /// @return ポリゴンタイプ PL_TYPE_NPT
    ///
    virtual int get_pl_type( void ) 
    {
        return PL_TYPE_NPT;
    }

#ifdef USE_NPATCH_LIB
    ///
    /// 頂点を設定
    ///
    /// @param[in] vertex           三角形の3頂点
    /// @param[in] update_param     面の法線ベクトル、長田パッチのパラメータ更新するか？
    /// @param[in] calc_area        面積を再計算するか？
    /// @attention 長田パッチのパラメータ更新は遅いため
    ///    長田パッチとわかっている場合は、長田パッチを同時に
    ///    更新するタイプは良い
    ///
    virtual void set_vertexes(
        const Vec3<PL_REAL> vertex[3],
        bool    update_param, 
        bool    calc_area
        ) ;
#endif

    ///
    /// 頂点・長田パッチパラメータを設定
    ///
    /// @param[in] vertex           三角形の3頂点
    /// @param[in] param            長田パッチパラメータ
    /// @param[in] update_param     面の法線ベクトルを更新するか？
    /// @param[in] calc_area        面積を再計算するか？
    /// @attention 長田パッチのパラメータは絶対座標の制御点であるため
    ///        移動の時は、頂点と同様に移動させる必要がある
    ///
    void set_vertexes(
        const Vec3<PL_REAL> vertex[3],
        const NpatchParam&  param,
        bool  update_param, 
        bool  calc_area
        ) ;

    ///
    /// 長田パッチパラメータの設定
    ///
    /// @param[in] param        長田パッチパラメータ
    ///
    void set_npatch_param(
        const NpatchParam&      param
        )
    {
        m_npatchParam = param;
    }


    ///
    /// 長田パッチパラメータの取得
    ///
    /// @param[in] param        長田パッチパラメータ
    ///
    NpatchParam* get_npatch_param( void ) const
    {
        // const_cast:const宣言を外して返す
        return const_cast< NpatchParam* >(&m_npatchParam);
    }

#ifdef USE_NPATCH_LIB
    ///
    /// 点の近似曲面補正
    ///
    ///  @param[in]     pos     3角形平面内の点
    ///  @param[out]    pos_o   曲面上に補正された座標
    ///  @return なし
    ///
    void  correct( 
           const Vec3<PL_REAL>&   pos,
           Vec3<PL_REAL>&         pos_o
         );
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
    void  correct( 
           PL_REAL          eta,
           PL_REAL          xi,
           Vec3<PL_REAL>&   pos_o
         );


protected:
    //=======================================================================
    // クラス変数
    //=======================================================================

    /// 長田パッチパラメータ
    NpatchParam     m_npatchParam;
    
};


} //namespace PolylibNS

#endif  // polylib_npatchtriangle_h

