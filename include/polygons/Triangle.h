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

#ifndef polylib_triangle_h
#define polylib_triangle_h

#include "common/Vec3.h"
#include "common/BBox.h"

using namespace Vec3class;

namespace PolylibNS{

////////////////////////////////////////////////////////////////////////////
///
/// クラス:Triangle
///
////////////////////////////////////////////////////////////////////////////

class Triangle {
public:
    ///
    /// コンストラクタ。
    ///
    Triangle()
    {
        m_numAtrI = 0;
        m_numAtrR = 0;
        m_AtrI = NULL;  // デストラクタのため設定
        m_AtrR = NULL;
    }

    ///
    /// コンストラクタ（コピーコンストラクタ）
    ///
    /// @param[in] tria コピー元。
    ///
    Triangle(
        const Triangle  &tria
    )
    {
        m_vertex[0] = tria.m_vertex[0];
        m_vertex[1] = tria.m_vertex[1];
        m_vertex[2] = tria.m_vertex[2];
        m_normal  = tria.m_normal;
        m_area    = tria.m_area;
        m_id      = tria.m_id;
        m_exid    = tria.m_exid;
        m_numAtrI = tria.m_numAtrI;
        m_numAtrR = tria.m_numAtrR;
        if( m_numAtrI > 0 ) {
            m_AtrI = (int*)malloc( m_numAtrI*sizeof(int) );
            for( int i=0; i<m_numAtrI; i++ ) {
                m_AtrI[i] = tria.m_AtrI[i];
            }
        } else {
            m_AtrI = NULL;
        }
        if( m_numAtrR > 0 ) {
            m_AtrR = (PL_REAL*)malloc( m_numAtrR*sizeof(PL_REAL) );
            for( int i=0; i<m_numAtrR; i++ ) {
                m_AtrR[i] = tria.m_AtrR[i];
            }
        } else {
            m_AtrR = NULL;
        }
    }


    ///
    /// コンストラクタ。
    ///
    /// @param[in] vertex ポリゴンの頂点
    /// @param[in] id        三角形ポリゴンID（ユニークな識別子）
    ///                          システム全体でユニークな識別子であること
    /// @param[in] num_atrI    ユーザ定義属性数（整数型）
    /// @param[in] num_atrR    ユーザ定義属性数（実数型）
    /// @param[in] atrI        ユーザ定義属性（整数型）
    /// @param[in] atrR        ユーザ定義属性（実数型）
    /// @attention 面積と法線はvertexを元に自動計算される
    ///    atrI,atrR=NULLの場合、属性領域は初期値が設定される
    ///     id=0の場合、IDは内部で採番される
    ///
    Triangle(
        const Vec3<PL_REAL> vertex[3],
        long long int       id = 0,
        int                 num_atrI = 0,
        int                 num_atrR = 0,
        const int           *atrI    = NULL,
        const PL_REAL       *atrR    = NULL
    )
    {
        m_vertex[0] = vertex[0];
        m_vertex[1] = vertex[1];
        m_vertex[2] = vertex[2];
        calc_normal();
        calc_area();
        if ( id != 0 ) {
            m_id = id;
        } else {
            m_id = create_unique_id();
        }
        m_exid    = 0;
        m_numAtrI = (unsigned char)num_atrI;
        m_numAtrR = (unsigned char)num_atrR;
        if( m_numAtrI > 0 ) {
            m_AtrI = (int*)malloc( m_numAtrI*sizeof(int) );
            if( atrI != NULL ) {
                for( int i=0; i<m_numAtrI; i++ ) {
                    m_AtrI[i] = atrI[i];
                }
            } else {
                for( int i=0; i<m_numAtrI; i++ ) {
                    m_AtrI[i] = 0;
                }
            }
        } else {
            m_AtrI = NULL;
        }
        if( m_numAtrR > 0 ) {
            m_AtrR = (PL_REAL*)malloc( m_numAtrR*sizeof(PL_REAL) );
            if( atrR != NULL ) {
                for( int i=0; i<m_numAtrR; i++ ) {
                    m_AtrR[i] = atrR[i];
                }
            } else {
                for( int i=0; i<m_numAtrR; i++ ) {
                    m_AtrR[i] = 0;
                }
            }
        } else {
            m_AtrR = NULL;
        }
    }

    ///
    /// コンストラクタ。
    ///
    /// @param[in] vertex   ポリゴンの頂点
    /// @param[in] normal   法線
    /// @param[in] id        三角形ポリゴンID（ユニークな識別子）
    ///                          システム全体でユニークな識別子であること
    /// @param[in] num_atrI    ユーザ定義属性数（整数型）
    /// @param[in] num_atrR    ユーザ定義属性数（実数型）
    /// @param[in] atrI        ユーザ定義属性（整数型）
    /// @param[in] atrR        ユーザ定義属性（実数型）
    /// @attention 面積はvertexを元に自動計算される
    ///     id=0の場合、IDは内部で採番される
    ///
    Triangle(
        const Vec3<PL_REAL> vertex[3], 
        const Vec3<PL_REAL>&    normal,
        long long int       id = 0,
        int                 num_atrI = 0,
        int                 num_atrR = 0,
        const int           *atrI    = NULL,
        const PL_REAL       *atrR    = NULL
        )
    {
        m_vertex[0] = vertex[0];
        m_vertex[1] = vertex[1];
        m_vertex[2] = vertex[2];
        m_normal    = normal;
        calc_area();
        if ( id != 0 ) {
            m_id = id;
        } else {
            m_id = create_unique_id();
        }
        m_exid    = 0;
        m_numAtrI = (unsigned char)num_atrI;
        m_numAtrR = (unsigned char)num_atrR;
        if( m_numAtrI > 0 ) {
            m_AtrI = (int*)malloc( m_numAtrI*sizeof(int) );
            if( atrI != NULL ) {
                for( int i=0; i<m_numAtrI; i++ ) {
                    m_AtrI[i] = atrI[i];
                }
            } else {
                for( int i=0; i<m_numAtrI; i++ ) {
                    m_AtrI[i] = 0;
                }
            }
        } else {
            m_AtrI = NULL;
        }
        if( m_numAtrR > 0 ) {
            m_AtrR = (PL_REAL*)malloc( m_numAtrR*sizeof(PL_REAL) );
            if( atrR != NULL ) {
                for( int i=0; i<m_numAtrR; i++ ) {
                    m_AtrR[i] = atrR[i];
                }
            } else {
                for( int i=0; i<m_numAtrR; i++ ) {
                    m_AtrR[i] = 0;
                }
            }
        } else {
            m_AtrR = NULL;
        }
    }

    ///
    /// コンストラクタ (deserialize)
    ///    シリアライズされた通信バッファよりオブジェクトの生成を行う
    ///
    /// @param[in] pbuff    バッファ格納先頭位置ポインタ
    ///                        シリアリズされたデータ
    /// @attention シリアライズはserialize()で行う
    ///
    Triangle(
            const char* pbuff
        );

    ///
    /// デストラクタ
    ///
    virtual ~Triangle()
    {
        if( m_AtrI != NULL )  free(m_AtrI);
        if( m_AtrR != NULL )  free(m_AtrR);
    }

    ///
    /// ユーザ定義属性数の設定
    ///     （ユーザ限定関数）
    ///     PolygonGroup内で同一属性数である必要があるため
    ///     一般的にはPolygonGroup::set_num_atr()にて設定する
    ///
    /// @param[in] num_atrI    ユーザ定義属性数（整数型）
    /// @param[in] num_atrR    ユーザ定義属性数（実数型）
    /// @return    なし
    /// @attention 属性格納領域をallocationする
    ///              新規アロケーション時は初期値を設定する
    ///              既にアロケーションされていた場合、reallocする
    ///
    void set_num_atr(
        int num_atrI, 
        int num_atrR
        )
    {
        if( num_atrI == m_numAtrI ) {
            // 何もしない
        } else if( num_atrI == 0 ) {
            // 解放する
            m_numAtrI = 0;
            if( m_AtrI != NULL )  {
                free(m_AtrI);
                m_AtrI = NULL;
            }
        } else {    // num_atrI > 0
            if( m_numAtrI == 0 ) {
                // 新規allocation
                m_numAtrI = num_atrI;
                m_AtrI = (int*)malloc( m_numAtrI*sizeof(int) );
            } else if( num_atrI < m_numAtrI ) {
                // realloc サイズ縮小
                int *p = m_AtrI;
                m_AtrI = (int*)realloc( p, num_atrI*sizeof(int) );
                m_numAtrI = num_atrI;
            } else {    // num_atrI > m_numAtrI
                // realloc サイズ拡大  拡大したところの初期化
                int *p = m_AtrI;
                m_AtrI = (int*)realloc( p, num_atrI*sizeof(int) );
                for( int i=m_numAtrI; i<num_atrI; i++ ) {
                    m_AtrI[i] = 0;
                }
                m_numAtrI = num_atrI;
            }
        }
        
        if( num_atrR == m_numAtrR ) {
            // 何もしない
        } else if( num_atrR == 0 ) {
            // 解放する
            m_numAtrR = 0;
            if( m_AtrR != NULL )  {
                free(m_AtrR);
                m_AtrR = NULL;
            }
        } else {    // num_atrR > 0
            if( m_numAtrR == 0 ) {
                // 新規allocation
                m_numAtrR = num_atrR;
                m_AtrR = (PL_REAL*)malloc( m_numAtrR*sizeof(PL_REAL) );
            } else if( num_atrR < m_numAtrR ) {
                // realloc サイズ縮小
                PL_REAL *p = m_AtrR;
                m_AtrR = (PL_REAL*)realloc( p, num_atrR*sizeof(PL_REAL) );
                m_numAtrR = num_atrR;
            } else {    // num_atrI > m_numAtrI
                // realloc サイズ拡大  拡大したところの初期化
                PL_REAL *p = m_AtrR;
                m_AtrR = (PL_REAL*)realloc( p, num_atrR*sizeof(PL_REAL) );
                for( int i=m_numAtrR; i<num_atrR; i++ ) {
                    m_AtrR[i] = 0.0;
                }
                m_numAtrR = num_atrR;
            }
        }
    }

    ///
    /// ユーザ定義属性数（整数型）の取得
    ///
    /// @return    ユーザ定義属性数（整数型）
    ///
    int get_num_atrI( void )
    {
        return (int)m_numAtrI;
    }


    ///
    /// ユーザ定義属性数（実数型）の取得
    ///
    /// @return    ユーザ定義属性数（実数型）
    ///
    int get_num_atrR( void )
    {
        return (int)m_numAtrR;
    }

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
        size += sizeof(Triangle);           // ユーザ定義属性を除くサイズ
        size += m_numAtrI*sizeof(int);      // ユーザ定義属性（整数）
        size += m_numAtrR*sizeof(PL_REAL);  // ユーザ定義属性（実数）

        return size;
    }

    ///
    /// ポリゴンのシリアリズした時のサイズ取得
    ///     通信バッファにシリアリズした時のサイズ
    ///
    /// @return シリアリズ後のサイズ
    ///
    virtual size_t serialized_size( void )
    {
        size_t size = 0;
        size += sizeof(int);                // ポリゴンタイプを設定
        size += sizeof(long long int);      // ポリゴンID
        size += 9*sizeof(PL_REAL);          // ３角形頂点座標
        size += sizeof(short int);          // ユーザ定義ID
        size += sizeof(unsigned char);      // ユーザ定義属性数（整数）
        size += sizeof(unsigned char);      // ユーザ定義属性数（実数）
        size += m_numAtrI*sizeof(int);      // ユーザ定義属性（整数）
        size += m_numAtrR*sizeof(PL_REAL);  // ユーザ定義属性（実数）

        return size;
    }

    ///
    /// ポリゴンのシリアリズ
    ///   １ポリゴンを通信バッファに格納する
    ///
    /// @param[out]  pbuff  バッファ格納位置先頭ポインタ
    ///                        シリアリズされたデータ
    /// @return バッファ格納位置Nextポインタ
    /// @attention  デシリアリズはPolylibNS::deserialize_polygon()で行う
    ///     最終的にはコンストラクタ Triangle( pbuff ) が呼び出される
    ///
    virtual char* serialize( const char* pbuff );

    ///
    /// Bounding box of this triangle
    ///
    /// @param[in] detail       曲面補間可能要素の時に曲面補間したBounding boxを返すか否か
    ///                             false: ３角形の頂点にて決定
    /// @return Bounding box
    ///
    //   VTree関連で使用する　特に長田パッチのbbox取得用
    virtual BBox get_bbox( bool detail = false );

    //=======================================================================
    // Setter/Getter
    //=======================================================================

    ///
    /// ポリゴンタイプ取得
    ///
    /// @return ポリゴンタイプ PL_TYPE_TRIANGLE
    ///
    virtual int get_pl_type( void ) 
    {
        return PL_TYPE_TRIANGLE;
    }


    ///
    /// 頂点を設定
    ///
    /// @param[in] vertex           三角形の3頂点
    /// @param[in] update_param     面の法線ベクトルを再計算するか？
    /// @param[in] calc_area        面積を再計算するか？
    /// @attention  長田パッチの場合、update_param=trueで
    ///    長田パッチのパラメータ更新も更新する
    ///    但し、この長田パッチのパラメータ更新は遅いため
    ///    長田パッチとわかっている場合は、長田パッチを同時に
    ///    更新するタイプは良い         
    ///
    virtual void set_vertexes(
        const Vec3<PL_REAL> vertex[3], 
        bool    update_param, 
        bool    calc_area
        ) ;

    ///
    /// 法線ベクトル・面積の更新
    ///
    /// @param[in] update_normal      法線ベクトルを再計算するか？
    /// @param[in] update_area        面積を再計算するか？
    /// @attention  get_vertexes()で座標のポインタを取得して
    ///    直接座標を書き換えた時などに、法線ベクトル・面積を更新するのに使用する
    ///
    virtual void update(
            bool    update_normal, 
            bool    update_area
        )
    {
        if( update_normal ) calc_normal();
        if( update_area   ) calc_area();
    }

    ///
    /// 頂点の取得
    ///
    /// @return 三角形の3頂点へのポインタ
    ///
    Vec3<PL_REAL>* get_vertexes() const 
    {
        return const_cast< Vec3<PL_REAL>* >(m_vertex);
    }

    ///
    /// 法線ベクトルを取得。
    ///
    /// @return 法線ベクトル
    /// @attention  return値はポインタではありません
    ///             メンバー変数でなく計算により求める方法に変更するかもしれないため
    ///
    Vec3<PL_REAL> get_normal() const
    {
        return m_normal;
    }

    ///
    /// 面積を取得。
    ///
    /// @return 面積。
    ///
    PL_REAL get_area() const
    {
        return m_area;
    }

    ///
    /// ポリゴンIDを設定。
    ///
    ///  @param[in] id  ポリゴンID。
    ///  @attention 一般ユーザ使用不可
    ///
    void set_id(long long int id)
    {
        m_id = id;
    }

    ///
    /// ポリゴンIDを返す。
    ///
    ///  @return ポリゴンID。
    ///
    long long int get_id() const
    {
        return m_id;
    }


    ///
    /// ユーザ定義IDを設定
    ///    FFV-C 専用
    ///
    void set_exid( int exid ) {
        m_exid = exid;
    }

    ///
    /// ユーザ定義IDを取得
    ///    FFV-C 専用
    ///
    /// @return ユーザ定義ID。
    ///
    int get_exid() const {
        return m_exid;
    }

    ///
    /// ユーザ定義属性数（整数型）のポインタ取得
    ///
    ///  @return    ユーザ定義属性数（整数型）のポインタ
    ///                  =NULL の場合、属性なし
    ///  @attention  ユーザ定義属性の設定は以下で行う
    ///                int* p = get_polygons_pAtrI();
    ///                p[2] = 1  属性種類の3番目に１を設定
    ///
    int* get_pAtrI() const 
    {
        return const_cast<int*>(m_AtrI); // const_cast:const宣言を外して返す
    }

    ///
    /// ユーザ定義属性数（実数型）のポインタ取得
    ///
    ///  @return    ユーザ定義属性数（実数型）のポインタ
    ///                  =NULL の場合、属性なし
    ///  @attention  ユーザ定義属性の設定は以下で行う
    ///                PL_REAL* p = get_polygons_pAtrR();
    ///                p[2] = 1.0  属性種類の3番目に1.0を設定
    ///
    PL_REAL* get_pAtrR() const
    {
        return const_cast<PL_REAL*>(m_AtrR); // const_cast:const宣言を外して返す
    }


protected:
    ///
    /// 法線ベクトル算出。
    ///
    void calc_normal()
    {
        Vec3<PL_REAL> vd[3];
        vd[0].assign( m_vertex[0].x, m_vertex[0].y, m_vertex[0].z );
        vd[1].assign( m_vertex[1].x, m_vertex[1].y, m_vertex[1].z );
        vd[2].assign( m_vertex[2].x, m_vertex[2].y, m_vertex[2].z );
        Vec3<PL_REAL> ad = vd[1] - vd[0];
        Vec3<PL_REAL> bd = vd[2] - vd[0];

        Vec3<PL_REAL> normald = (cross(ad,bd)).normalize();
        m_normal[0] = normald[0];
        m_normal[1] = normald[1];
        m_normal[2] = normald[2];
    }

    ///
    /// 面積算出。
    ///
    void calc_area()
    {
        Vec3<PL_REAL> a = m_vertex[1] - m_vertex[0];
        Vec3<PL_REAL> b = m_vertex[2] - m_vertex[0];
        PL_REAL al = a.length();
        PL_REAL bl = b.length();
        PL_REAL ab = dot(a,b);
        PL_REAL f = al*al*bl*bl - ab*ab;
        if(f<0.0) f=0.0;
        m_area = 0.5*sqrtf(f);
    }


    ///
    /// システムで一意のポリゴンIDを作成する
    ///     MPI環境の場合、ランク０以外で当ルーチンが呼ばれてはならない
    ///     ポリゴン生成時にユーザ指定でポリゴンIDを指定することも可能であるが
    ///      全体で整合性のあるものに制御すること
    ///
    ///  @return    ポリゴンID
    ///
    long long int create_unique_id();      //  ポリゴンID m_id採番用



    //=======================================================================
    // クラス変数
    //=======================================================================
protected:

    // ---- 8byte -------------------------------

    /// 一意となるポリゴンID（内部識別用）
    //              ファイルには出力しない
    //    Polylib環境下でユニークとするため4byteでは不足するので8byteとしている
    //    検索等で返されるのはTriangleクラスなので、このクラスにポリゴンIDを
    //    保持させておく必要がある
    long long int       m_id;      // 並列実行時の重複を解消する、
                                   // save時のソートなどに使用

    /// ユーザ定義属性（整数型）
    //   属性数はm_numAtrIに設定
    //   何番目の属性かはユーザにて管理する
    int*     m_AtrI;

    /// ユーザ定義属性（実数型）
    //   属性数はm_numAtrRに設定
    //   何番目の属性かはユーザにて管理する
    PL_REAL*     m_AtrR;


    // ---- PL_REAL系 4byte or 8byte ------------------------

    /// 三角形の頂点座標（反時計回りで並んでいる）
    //                  STLファイル, NPTファイルに出力
    Vec3<PL_REAL>   m_vertex[3];

    /// 三角形の法線ベクトル
    //                  STLファイルに出力
    Vec3<PL_REAL>   m_normal;

    /// 三角形の面積
    //     曲面補間可能なものでも補間しないものとする
    //              ファイルには出力しない
    PL_REAL m_area;


    // ---- 2byte -------------------------------------------

    /// 三角形のユーザ定義ID（FFV-Cでのみ使用すること）
    //    ユーザ定義属性側に移動させることも考えたが、
    //    STLバイナリファイルの時のみユーザ定義属性の先頭に属性が追加されているのは良くなため
    //    当属性は残すものとする
    short int m_exid;   // 現状、STLバイナリファイルの未定義域 2byteに格納されるのみ
                        // FFV-C 専用


    // ---- 1byte -------------------------------------------

    /// ユーザ定義属性数（整数型）
    unsigned char  m_numAtrI;

    /// ユーザ定義属性数数（実数型）
    unsigned char  m_numAtrR;

};


} //namespace PolylibNS

#endif  // polylib_triangle_h

