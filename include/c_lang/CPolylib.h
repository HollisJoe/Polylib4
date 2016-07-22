/* -*- Mode: c++ -*- */
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

#ifndef c_polylib_h
#define c_polylib_h

#ifdef MPI_PL
#include "mpi.h"
#endif

#ifdef __cplusplus
extern "C" {  // for C++
#else
#endif

//#ifdef MPI_PL
//#include "mpi.h"   // compile errorとなる
//#endif

#include <stdbool.h>

#include "common/PolylibStat.h"
#include "common/PolylibDefine.h"



// #ifdef HAVE_CONFIG_H
// #include "config.h"
// #endif 

///
/// C言語用Polylib
///
/// 注意：
///       C言語APIではタグを操作（ハンドル）用として使用します。
///       タグには有効期間があります。セッション中で永続的に有効というわけではありません。
///         PolygonGroupタグ
///             PolygonGroupの追加、削除、挿入があると該当する箇所以外も含めて全て無効となります。
///         Triangleタグ
///             Triangleの追加、削除、挿入、ソートがあると該当する箇所以外も含めて全て無効となります。
///             並列環境でTriangleを移動した場合、migrate処理を実行しますが
///             この際に削除、挿入、ソートが行われるのでタグが全て無効となります。

///
/// 領域情報構造体
///
typedef struct {
    PL_REAL bpos[3];
    unsigned int bbsize[3];
    unsigned int gcsize[3];
    PL_REAL dx[3];
} ParallelBboxStruct;

///
/// 三角形ポリゴン情報構造体
///
typedef struct {
    PL_REAL vertex[9];  ///< ３頂点座標
    PL_REAL normal[3];  ///< 面法線ベクトル
    //PL_REAL area;     ///< 面積
} TriangleStruct;

///
/// 長田パッチパラメータ（曲面補間用制御点情報）
///
typedef struct {
    PL_REAL   cp_side1_1[3];    ///<  p1p2辺の３次ベジェ制御点1
    PL_REAL   cp_side1_2[3];    ///<  p1p2辺の３次ベジェ制御点2
    PL_REAL   cp_side2_1[3];    ///<  p2p3辺の３次ベジェ制御点1
    PL_REAL   cp_side2_2[3];    ///<  p2p3辺の３次ベジェ制御点2
    PL_REAL   cp_side3_1[3];    ///<  p3p1辺の３次ベジェ制御点1
    PL_REAL   cp_side3_2[3];    ///<  p3p1辺の３次ベジェ制御点2
    PL_REAL   cp_center[3];     ///<  三角形中央の３次ベジェ制御点
} NpatchParamStruct;

///
/// 長田パッチポリゴン情報構造体
///
typedef struct {
    PL_REAL vertex[9];          ///< ３頂点座標
    NpatchParamStruct   param;  ///< 長田パッチパラメータ
} NptTriangleStruct;

////////////////////////////////////////////////////////////////////////////
///
/// 構造体:PolylibMoveParamsStruct
/// polylib_move()の引数として利用するパラメタセットの構造体です。
/// 本構造体メンバ変数ではパラメタが不足する場合は、C++側のPolylibMoveParamsも含めて
/// ユーザ定義する必要がある
///      PolylibMoveParamsの継承クラスを作成する、構造体もそれに合わせて新規に作成するなど
///
////////////////////////////////////////////////////////////////////////////
typedef struct {
    /// 現在の計算ステップ番号
    int m_current_step;

    /// 移動後の計算ステップ番号
    int m_next_step;

    /// １計算ステップあたりの時間変異
    PL_REAL m_delta_t;

    /// ユーザ定義パラメータ
    //    m_current_step, m_next_step, m_delta_t以外のパラメータが必要な時に任意に使用する
    PL_REAL m_params[10];
} PolylibMoveParamsStruct;

//**********************************************************
//  Polylibクラス用ラッパー関数
//**********************************************************

///
/// C言語用Polylib環境の構築
///     Polylibインスタンス生成
/// Polylib::get_instanceメソッドの代替関数。
///
///  @return POLYLIB_STATで定義される値が返る。
///  @attention 最初に呼び出すこと
///
POLYLIB_STAT  polylib_instance( void );


#ifdef MPI_PL
///
/// 並列計算関連情報の設定と初期化
///     (各ランクが１領域を担当している場合）
/// Polylib::init_parallel_infoメソッドのラッパー関数。
///
///  @param[in] comm    MPIコミュニケーター
///  @param[in] bpos    自PE担当領域の基点座標
///  @param[in] bbsize  同、計算領域のボクセル数
///  @param[in] gcsize  同、ガイドセルのボクセル数
///  @param[in] dx      同、ボクセル１辺の長さ
///  @return POLYLIB_STATで定義される値が返る。
///  @attention 全rankで各々設定を行い、その領域情報を全rankへ配信する
///
POLYLIB_STAT
polylib_init_parallel_info(
                              MPI_Comm comm,
                              PL_REAL bpos[3],
                              unsigned int bbsize[3],
                              unsigned int gcsize[3],
                              PL_REAL dx[3]
);

///
/// 並列計算関連情報の設定と初期化
/// (各ランクが複数領域を担当している場合）
/// 全rankで各々設定を行い、その領域情報を全rankへ配信する。
///
///  @param[in] comm    MPIコミュニケーター
///  @param[in] num     自PEが担当する領域
///  @param[in] bbox    担当するboundary box情報（複数）
///  @return POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT
polylib_init_parallel_info2(
    MPI_Comm comm,
    int      num,
    ParallelBboxStruct *bboxes
);
#endif

///
/// Polylib::loadメソッドのラッパー関数。
/// 引数で指定された設定ファイルを読み込み、グループツリーを作成する。
/// 続いて設定ファイルで指定されたSTLファイルを読み込み、KD木を作成する。
///  @param[in] config_name  設定ファイル名
///                             NULLの場合、デフォルト値 polylib_config.tpとする
///  @param[in] scale        縮尺率
///  @return    POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT polylib_load(
                char*   config_name,
                PL_REAL scale
            );

///
/// PolygoGroupツリー、三角形ポリゴン情報の保存。
///   Polylib::saveメソッドのラッパー関数。
/// グループツリーの情報を設定ファイルへ出力。三角形ポリゴン情報をSTL/NPT
/// ファイルへ出力
///
///  @param[out]    p_fname 設定ファイル名返却用
///                             内部領域のアドレスを返す
///                             free不可
///  @param[in]     format  形状ファイルのフォーマット
///                             PolylibDefine.hで定義されているFILE_FMT_* 参照
///  @param[in]     extend  ファイル名に付加する文字列。NULLを指定した
///                         場合は、付加文字列として本メソッド呼び出し時の
///                         年月日時分秒(YYYYMMDD24hhmmss)を用いる
///  @return    POLYLIB_STATで定義される値が返る。
///  @attention ファイル名命名規約は次の通り。
///         設定ファイル : polylib_config_付加文字.tpp
///         STL/NPTファイル  : ポリゴングループ名_付加文字.拡張子
///
POLYLIB_STAT polylib_save(
                          char  **p_fname,
                          char  *format,
                          char  *extend
                          );

///
/// 三角形ポリゴン座標の移動
///   Polylib::moveメソッドのラッパー関数
/// 本クラスインスタンス配下の全PolygonGroupのmoveメソッドが呼び出される。
/// move関数は、polylib_group_set_move_func_c関数で登録する
///
///  @param[in]     param   移動計算パラメータセット
///  @return    POLYLIB_STATで定義される値が返る
///
POLYLIB_STAT polylib_move(
                          PolylibMoveParamsStruct* param
                          );


#ifdef MPI_PL
///
/// ポリゴンデータのPE間移動
///   Polylib::mograteメソッドのラッパー関数
/// オブジェクトのインスタンス毎に登録が必要
///  @return POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT polylib_migrate( void );
#endif

///
/// ポリゴングループタグの取得
///      Polylib::get_groupメソッドのラッパー関数
///  @param[in]     group_name  ポリゴングループ名
///  @param[out]    tag         ポリゴングループ
///  @return    POLYLIB_STATで定義される値が返る
///
POLYLIB_STAT polylib_get_group_tag(
        char*       group_name,
        PL_GRP_TAG* grp_tag
    );
    
///
/// PolygonGroupツリーの最上位ノードの取得（Cインターフェース用）
///      Polylib::get_root_groupsメソッドのラッパー関数
///
///  @param[out]        n       ポリゴングループ数
///  @param[out]        tags    ポリゴングループ　タグ
///  @return    POLYLIB_STATで定義される値が返る
///  @attention   tagsは使用後freeして下さい
///
POLYLIB_STAT polylib_get_root_groups_tags( int* n, PL_GRP_TAG** tags );

///
/// リーフPolygonGroupリストの取得（Cインターフェース用）
///      Polylib::get_leaf_groupsメソッドのラッパー関数
///
///  @param[out]        n       PolygonGroup数
///  @param[out]        tags    PolygonGroupタグ
///  @return    POLYLIB_STATで定義される値が返る
///  @attention   tagsは使用後freeして下さい
///
POLYLIB_STAT  polylib_get_leaf_groups_tags( int* n, PL_GRP_TAG** tags );

///
/// Polylib::search_polygonsメソッドのラッパー関数。
/// 位置ベクトルmin_posとmax_posにより特定される矩形領域に含まれる、
/// 特定のグループとその子孫グループに属する三角形ポリゴンをKD探索に
/// より抽出する。Polylib内でメモリ領域が確保される
///  @param[out]    num         抽出された三角形ポリゴン数
///  @param[out]    tags        三角形ポリゴンのタグ（ハンドル)
///  @param[in]     group_name  抽出グループ名。
///  @param[in]     min_pos     抽出する矩形領域の最小値。(x,y,z順の配列)
///  @param[in]     max_pos     抽出する矩形領域の最大値。(x,y,z順の配列)
///  @param[in]     every       抽出オプション。
///   1：3頂点が全て検索領域に含まれるポリゴンを抽出する。
///   0：三角形のBBoxが一部でも検索領域と交差するものを抽出する
///  @return    POLYLIB_STATで定義される値が返る
///  @attention tags,trisはfreeしてください。
///             MPI並列計算時は,min_pos, max_posは各ランクの矩形領域を
///             超えないようにして下さい (各ランク内の担当領域内のみ検索するため）
///
POLYLIB_STAT polylib_search_polygons(
    int         *num,
    PL_ELM_TAG  **tags,
    char*       group_name,
    PL_REAL     min_pos[3],
    PL_REAL     max_pos[3],
    int         every 
    );

///
/// Polylib::search_polygonsメソッドのラッパー関数。
/// 位置ベクトルmin_posとmax_posにより特定される矩形領域に含まれる、
/// 特定のグループとその子孫グループに属する三角形ポリゴンをKD探索に
/// より抽出する。Polylib内でメモリ領域が確保される。
///  @param[out]    num         抽出された三角形ポリゴン数
///  @param[out]    tags        三角形ポリゴンのタグ（ハンドル)
///  @param[out]    tris        三角形ポリゴン
///  @param[in]     group_name  抽出グループ名。
///  @param[in]     min_pos     抽出する矩形領域の最小値。(x,y,z順の配列)
///  @param[in]     max_pos     抽出する矩形領域の最大値。(x,y,z順の配列)
///  @param[in]     every       抽出オプション。
///   1：3頂点が全て検索領域に含まれるポリゴンを抽出する。
///   0：三角形のBBoxが一部でも検索領域と交差するものを抽出する。
///  @return    POLYLIB_STATで定義される値が返る
///  @attention tags,trisはfreeしてください。
///             MPI並列計算時は,min_pos, max_posは各ランクの矩形領域を
///             超えないようにして下さい (各ランク内の担当領域内のみ検索するため）
///             ポリゴンに長田パッチが含まれている場合エラーになります。
///
POLYLIB_STAT polylib_search_polygons_triangle(
    int         *num,
    PL_ELM_TAG  **tags,
    TriangleStruct **tris,
    char*       group_name,
    PL_REAL     min_pos[3],
    PL_REAL     max_pos[3],
    int         every 
    );

///
/// Polylib::search_polygonsメソッドのラッパー関数。
/// 位置ベクトルmin_posとmax_posにより特定される矩形領域に含まれる、
/// 特定のグループとその子孫グループに属する三角形ポリゴンをKD探索に
/// より抽出する。Polylib内でメモリ領域が確保される
///  @param[out]    num         抽出された三角形(長田パッチ）ポリゴン数
///  @param[out]    tags        三角形(長田パッチ）ポリゴンのタグ（ハンドル)
///  @param[out]    tris        三角形(長田パッチ）ポリゴン
///  @param[in]     group_name  抽出グループ名。
///  @param[in]     min_pos     抽出する矩形領域の最小値。(x,y,z順の配列)
///  @param[in]     max_pos     抽出する矩形領域の最大値。(x,y,z順の配列)
///   1：3頂点が全て検索領域に含まれるポリゴンを抽出する。
///   0：三角形のBBoxが一部でも検索領域と交差するものを抽出する。
///  @param[in]     every       抽出オプション。
///  @return    POLYLIB_STATで定義される値が返る
///  @attention tags,trisはfreeしてください。
///             MPI並列計算時は,min_pos, max_posは各ランクの矩形領域を
///             超えないようにして下さい (各ランク内の担当領域内のみ検索するため）
///             ポリゴンに長田パッチ以外が含まれている場合エラーになります。
///
POLYLIB_STAT polylib_search_polygons_npt(
    int         *num,
    PL_ELM_TAG  **tags,
    NptTriangleStruct** tris,
    char*       group_name,
    PL_REAL     min_pos[3],
    PL_REAL     max_pos[3],
    int         every 
    );


///
///  指定した点に最も近いポリゴンの検索
///     Polylib::search_nearest_polygonメソッドのラッパー関数。
///  @param[out]    tag         ポリゴンのタグ（ハンドル)
///                               PL_NULL_TAG 検索なし 
///  @param[in]     group_name  グループ名。
///  @param[in]     pos         指定点
///  @return    POLYLIB_STATで定義される値が返る
///             MPI並列計算時は,posは各ランクの矩形領域を
///             超えないようにして下さい (各ランク内の担当領域内のみ検索するため）
///
POLYLIB_STAT polylib_search_nearest_polygon(
        PL_ELM_TAG  *tag,
        char*       group_name,
        PL_REAL     pos[3]
    );

///
/// Polylib::show_group_hierarchyメソッドのラッパー関数。
/// グループ階層構造リストを標準出力に出力する。
///
void polylib_show_group_hierarchy();

///
/// グループの情報を出力する。(親グループ名、自身の名前、ファイル名、
///      Polylib::show_group_infoメソッドのラッパー関数。
///   登録三角形数、3頂点ベクトルの座標、法線ベクトルの座標、面積)
///  @param[in] group_name グループ名
///  @return POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT polylib_show_group_info(char* group_name);


///
///  Polylibが利用中の概算メモリ量を返す
///
/// @return 利用中のメモリ量(byte)
///
size_t polylib_used_memory_size();


///
/// Polylibが利用中の概算メモリ量(MB)を返す
///
/// @return 利用中のメモリ量(Mbyte)
///
size_t polylib_used_memory_size_mb();


//**********************************************************
//  PolygonGroupクラス用ラッパー関数
//**********************************************************

///
/// PolygonGroupのポリゴンを求める
///    PolygonGroup::get_trianglesのラッパー関数
///  @param[in]   tag_pg       PolygonGroupを操作するためのタグ
///  @param[out]  num_tri      三角形ポリゴン数
///  @param[out]  tags_tri     三角形ポリゴンのタグ（ハンドル)
///  @return    POLYLIB_STATで定義される値が返る
///  @attention tags_ trisはfreeしてください。
///
POLYLIB_STAT polylib_group_get_triangles( 
            PL_GRP_TAG tag_pg, 
            int        *num_tri,
            PL_ELM_TAG  **tags_tri
         );

///
/// PolygonGroupのポリゴン要素数を求める
///    PolygonGroup::get_group_num_triaのラッパー関数
///  @param[in]   tag_pg       PolygonGroupを操作するためのタグ
///  @return    三角形ポリゴン数
///  @attention 全プロセス通した要素数が必要な場合は、
///     polylib_group_get_num_global_tria()を使用する
///
int polylib_group_get_num_triangles( 
            PL_GRP_TAG tag_pg 
         );

#ifdef MPI_PL
///
/// PolygonGroupのポリゴン要素数(global)を求める
///    PolygonGroup::get_group_num_global_triaのラッパー関数
///  @param[in]   tag_pg       PolygonGroupを操作するためのタグ
///  @return    三角形ポリゴン数(global)
///  @attention 並列環境用
///     ポリゴンの重複を削除するための通信あり
///     polylib_group_get_num_tria()よりも大幅に処理時間がかかることに注意
///
int polylib_group_get_num_global_triangles( 
            PL_GRP_TAG tag_pg 
         );
#endif

///
/// PolygonGroupのポリゴンの面積を求める
///    PolygonGroup::get_group_num_triaのラッパー関数
///  @param[in]   tag_pg       PolygonGroupを操作するためのタグ
///  @return    ポリゴンの面積
///  @attention 全プロセス通したポリゴンの面積が必要な場合は、
///     polylib_group_get_global_area()を使用する
///
PL_REAL polylib_group_get_area( 
            PL_GRP_TAG tag_pg 
         );

#ifdef MPI_PL
///
/// PolygonGroupのポリゴンの面積(global)を求める
///    PolygonGroup::get_group_num_global_triaのラッパー関数
///  @param[in]   tag_pg       PolygonGroupを操作するためのタグ
///  @return    ポリゴンの面積(global)
///  @attention 並列環境用
///     ポリゴンの重複を削除するための通信あり
///     polylib_group_get_area()よりも大幅に処理時間がかかることに注意
///
PL_REAL polylib_group_get_global_area( 
            PL_GRP_TAG tag_pg 
         );
#endif

/// グループ内のポリゴン属性（整数）の集合演算値を返す
///     並列化されている場合は全プロセスを通した値
///     (PL_OP_SUM：重複ポリゴン分は無視される）
///     全プロセスに同じ値が返る
///  @param[in]   tag_pg       PolygonGroupを操作するためのタグ
///  @param[in]  op     演算種類　PL_OP_SUM/PL_OP_MAX/PL_OP_MIN
///  @param[in]  atr_no ポリゴン整数属性の何番目か　0〜
///  @param[out] val    属性値
///  @return    POLYLIB_STATで定義される値が返る。
///                ポリゴンが存在しない
///                ポリゴン属性が存在しないなど
///
POLYLIB_STAT polylib_group_get_polygons_reduce_atrI(
            PL_GRP_TAG tag_pg, 
            PL_OP_TYPE op,
            int        atr_no,
            int*       val
    );

/// グループ内のポリゴン属性（実数）の集合演算値を返す
///     並列化されている場合は全プロセスを通した値
///     (PL_OP_SUM：重複ポリゴン分は無視される）
///     全プロセスに同じ値が返る
///  @param[in]   tag_pg       PolygonGroupを操作するためのタグ
///  @param[in]  op     演算種類　PL_OP_SUM/PL_OP_MAX/PL_OP_MIN
///  @param[in]  atr_no ポリゴン実数属性の何番目か　0〜
///  @param[out] val    属性値
///  @return    POLYLIB_STATで定義される値が返る。
///                ポリゴンが存在しない
///                ポリゴン属性が存在しないなど
///
POLYLIB_STAT polylib_group_get_polygons_reduce_atrR(
            PL_GRP_TAG tag_pg, 
            PL_OP_TYPE op,
            int        atr_no,
            PL_REAL*   val
    );


///
/// PolygonGroupに移動関数を登録するためのラッパー関数。
/// オブジェクトのインスタンス毎に登録が必要
///  @param[in] tag       PolygonGroupを操作するためのタグ
///  @param[in] func      移動関数へのポインタ
///  @return POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT polylib_group_set_move_func_c( 
            PL_GRP_TAG tag, 
            void (*func)(PL_GRP_TAG,PolylibMoveParamsStruct*)
         );

///
/// KD木の再構築フラグの設定
///     ユーザ定義の移動関数内の最後で呼び出す
///
///  @param[in] tag       PolygonGroupを操作するためのタグ
///  @return    戻り値なし
///
void polylib_group_set_need_rebuild( PL_GRP_TAG tag );


//----------------------------------------------
// setter / getter
//----------------------------------------------

/// PolygonGroup グループ名設定
///  @param[in] tag       PolygonGroupを操作するためのタグ
///  @param[in] name      グループ名
///  @return   戻り値なし
void polylib_group_set_name( 
            PL_GRP_TAG tag, 
            char* name
         );

/// PolygonGroup グループ名取得
///  @param[in]  tag       PolygonGroupを操作するためのタグ
///  @param[out] name      グループ名
///  @return   戻り値なし
void polylib_group_get_name( 
            PL_GRP_TAG tag, 
            char* name
         );

/// PolygonGroup 移動対象フラグ設定
///  @param[in]  tag       PolygonGroupを操作するためのタグ
///  @param[in]  moval     移動対象フラグ 1:true  0:false
///  @return   戻り値なし
void polylib_group_set_movable( 
            PL_GRP_TAG tag, 
            int        movable
         );

/// PolygonGroup 移動対象フラグ取得
///  @param[in]  tag       PolygonGroupを操作するためのタグ
///  @param[out] moval     移動対象フラグ 1:true  0:false
void polylib_group_get_movable( 
            PL_GRP_TAG tag, 
            int*       movable
         );

/// PolygonGroup ユーザ定義属性の設定
///  @param[in]     tag     PolygonGroupを操作するためのタグ
///  @param[in]     key     キー
///  @param[in]     val     属性値
///  @return なし
///  @attention 既に登録されていた場合、上書きする
void polylib_group_set_atr( 
            PL_GRP_TAG tag,
            char* key,
            char* val
         );

/// PolygonGroup ユーザ定義属性の取得
///  @param[in]     tag     PolygonGroupを操作するためのタグ
///  @param[in]     key     キー
///  @param[out]    val     属性値
///  @return OK/NG  NG:キーと属性のペアが登録されていない
POLYLIB_STAT polylib_group_get_atr( 
            PL_GRP_TAG tag,
            char* key,
            char* val
         );

/// ポリゴングループ内のポリゴンのユーザ定義属性数の設定
///  @param[in]     tag      PolygonGroupを操作するためのタグ integer*8
///                            integer*8
///  @param[in]     num_atrI 整数属性数
///  @param[in]     num_atrR 実数属性数
///  @return 戻り値なし
void polylib_group_set_num_polygon_atr (
                PL_GRP_TAG tag,
                int num_atrI, 
                int num_atrR
            );

///
/// PolygonGroup 親グループを取得
///
///  @param[in]     tag     PolygonGroupを操作するためのタグ
///  @return 親グループのタグ
///
PL_GRP_TAG polylib_group_get_parent( PL_GRP_TAG tag );

///
/// PolygonGroup 子グループを取得
///  @param[in]     tag     PolygonGroupを操作するためのタグ
///  @param[out]    num     ポリゴングループ数
///  @param[out]    child_tags ポリゴングループ　タグ
///  @return 戻り値なし
///  @attention   tagsは使用後freeして下さい
void polylib_group_get_children(
                  PL_GRP_TAG tag,
                  int* num,
                  PL_GRP_TAG** child_tags
              );


//**********************************************************
//  Triangleクラス用ラッパー関数
//      三角形オブジェクトと長田パッチオブジェクトでは動作が違う場合あり
//**********************************************************

/// ポリゴンタイプ取得
///  @param[in]     tag     Triangle/NptTriangleを操作するためのタグ
///  @return PL_TYPE_TRIANGLE / PL_TYPE_NPT（長田パッチ）
int polylib_triangle_get_pl_type(
            PL_ELM_TAG tag
         );

//----------------------------------------------
// setter / getter
//----------------------------------------------

/// 頂点座標設定
///   基本はTriangle用
///  @param[in]     tag     Triangle/NptTriangleを操作するためのタグ
///  @param[in]     vertex  ３頂点の座標
///  @return 戻り値なし
///  @attention 頂点設定時、法線ベクトル,面積も内部で設定する
///             長田パッチのパラメータは更新されないので注意
void polylib_triangle_set_vertexes( 
            PL_ELM_TAG tag, 
            PL_REAL vertex[9]
         );

/// 頂点座標・長田パッチパラメータ設定
///  @param[in]     tag     Triangle/NptTriangleを操作するためのタグ
///  @param[in]     vertex  ３頂点の座標
///  @param[in]     param   長田パッチパラメータ
///  @return 戻り値なし
///  @attention 頂点設定時、法線ベクトル,面積も内部で設定する
void polylib_triangle_set_vertexes_npatch( 
            PL_ELM_TAG tag, 
            PL_REAL vertex[9],
            NpatchParamStruct* param
         );

/// 頂点座標取得
///  @param[in]     tag     Triangle/NptTriangleを操作するためのタグ
///  @param[out]    vertex  ３頂点の座標
///  @return 戻り値なし
void polylib_triangle_get_vertexes( 
            PL_ELM_TAG tag, 
            PL_REAL vertex[9]
         );

/// 法線ベクトル取得
///  @param[in]     tag     Triangle/NptTriangleを操作するためのタグ
///  @param[out]    norm    法線ベクトル
///  @return 戻り値なし
void polylib_triangle_get_normal( 
            PL_ELM_TAG tag, 
            PL_REAL norm[3]
         );

/// 長田パッチパラメータ設定
///  @param[in]     tag     NptTriangleを操作するためのタグ
///  @param[in]     param   長田パッチパラメータ
///  @return POLYLIB_STATで定義される値が返る。
///      オブジェクトが長田パッチでない時はエラーを返す
POLYLIB_STAT polylib_triangle_set_npatchParam( 
            PL_ELM_TAG tag, 
            NpatchParamStruct *param
         );

/// 長田パッチパラメータ取得
///  @param[in]     tag     NptTriangleを操作するためのタグ
///  @param[out]    param   長田パッチパラメータ
///  @return POLYLIB_STATで定義される値が返る。
///      オブジェクトが長田パッチでない時はエラーを返す
POLYLIB_STAT polylib_triangle_get_npatchParam( 
            PL_ELM_TAG tag, 
            NpatchParamStruct *param
         );

/// ユーザ定義属性数（整数型）の取得
///  @param[in]     tag      Triangle/NptTriangleを操作するためのタグ
///  @return 整数属性数
int polylib_triangle_get_num_atrI(
                PL_ELM_TAG tag
            );

/// ユーザ定義属性数（実数型）の取得
///  @param[in]     tag      Triangle/NptTriangleを操作するためのタグ
///  @return 実数属性数
int polylib_triangle_get_num_atrR(
                PL_ELM_TAG tag
            );

/// ユーザ定義属性（整数型）のポインタ取得
///  @param[in]     tag      Triangle/NptTriangleを操作するためのタグ
///  @return ユーザ定義属性（整数型）へのポインタ
///          freeしてはならない
int* polylib_triangle_get_pAtrI(
                PL_ELM_TAG tag
            );


/// ユーザ定義属性（実数型）のポインタ取得
///  @param[in]     tag      Triangle/NptTriangleを操作するためのタグ
///  @return ユーザ定義属性（実数型）へのポインタ
///          freeしてはならない
PL_REAL* polylib_triangle_get_pAtrR(
                PL_ELM_TAG tag
            );

#ifdef __cplusplus
} // extern "C" or extern
#else
#endif

#endif //c_polylib_h
