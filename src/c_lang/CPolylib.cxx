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

#ifdef MPI_PL
#include "mpi.h"
#endif
#include <string>
#include "c_lang/CPolylib.h"
#include "Polylib.h"

using namespace std;
using namespace PolylibNS;

/************************************************************************
 *
 * C言語用Polylib
 *
 ***********************************************************************/

///
///   Polylibインスタンス(静的変数）
//static Polylib* p_polylib_instance = NULL;
Polylib* p_polylib_instance = NULL;  // Fortranで使用するためstaticを外す


/// C言語用Polylib環境の構築
///     Polylibインスタンス生成
///
///  @return POLYLIB_STATで定義される値が返る。
///  @attention 最初に呼び出すこと
///
POLYLIB_STAT  polylib_instance( void )
{
    p_polylib_instance = Polylib::get_instance();
    if( p_polylib_instance == NULL ) {
        return PLSTAT_NG;
    }
    return PLSTAT_OK;
}


#ifdef MPI_PL
///
/// MPIPolylib::init_parallel_infoメソッドのラッパー関数。
/// 並列計算関連情報の設定と初期化を行う。
    /// (各ランクが１領域を担当している場合）
/// 全rankで各々設定を行い、その領域情報を全rankへ配信する。
///
///  @param[in] comm    MPIコミュニケーター
///  @param[in] bpos    自PE担当領域の基点座標
///  @param[in] bbsize  同、計算領域のボクセル数
///  @param[in] gcsize  同、ガイドセルのボクセル数
///  @param[in] dx      同、ボクセル１辺の長さ
///  @return POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT
polylib_init_parallel_info(
                              MPI_Comm comm,
                              PL_REAL bpos[3],
                              unsigned int bbsize[3],
                              unsigned int gcsize[3],
                              PL_REAL dx[3]
    )
{
    return p_polylib_instance->init_parallel_info( comm, bpos, bbsize, gcsize, dx );
}

///
/// 並列計算関連情報の設定と初期化を行う。
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
                 ParallelBboxStruct *bboxes_s
    )
{
    std::vector<ParallelBbox>  bboxes;
    for(int i=0; i<num; i++ ) {
        ParallelBbox bbox;
        bbox.bpos  [0]=bboxes_s[i].bpos  [0]; 
        bbox.bpos  [1]=bboxes_s[i].bpos  [1];
        bbox.bpos  [2]=bboxes_s[i].bpos  [2];
        bbox.bbsize[0]=bboxes_s[i].bbsize[0];
        bbox.bbsize[1]=bboxes_s[i].bbsize[1];
        bbox.bbsize[2]=bboxes_s[i].bbsize[2];
        bbox.gcsize[0]=bboxes_s[i].gcsize[0];
        bbox.gcsize[1]=bboxes_s[i].gcsize[1];
        bbox.gcsize[2]=bboxes_s[i].gcsize[2];
        bbox.dx    [0]=bboxes_s[i].dx    [0];
        bbox.dx    [1]=bboxes_s[i].dx    [1];
        bbox.dx    [2]=bboxes_s[i].dx    [2];

        bboxes.push_back( bbox );
    }

    return p_polylib_instance->init_parallel_info( comm, bboxes );
}
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
POLYLIB_STAT
polylib_load(
                char*   config_name,
                PL_REAL scale
            )
{
    std::string config_file_name;
    if( config_name == NULL ) {
        config_file_name = "polylib_config.tp";
    } else {
        config_file_name = config_name;
    }

    return p_polylib_instance->load( config_file_name, scale );
}


///
/// PolygoGroupツリー、三角形ポリゴン情報の保存。
///    Polylib::saveメソッドのラッパー関数。
/// グループツリーの情報を設定ファイルへ出力。三角形ポリゴン情報をSTL/NPT
/// ファイルへ出力
///
///  @param[out]    p_fname 設定ファイル名
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
POLYLIB_STAT
polylib_save(
               char  **p_fname,
               char  *format,
               char  *extend
             )
{
    POLYLIB_STAT ret;
    std::string  config_name_out;
    std::string  format_wk = format;
    std::string  extend_wk = extend;
    static char file_name[128]; 

    file_name[0] = '\0';

    ret = p_polylib_instance->save( config_name_out, format_wk, extend_wk );

    if( ret == PLSTAT_OK ) {
        strcpy( file_name, config_name_out.c_str() );
        *p_fname = file_name;
    }

    return ret;
}


///
/// 三角形ポリゴン座標の移動。
///   Polylib::moveメソッドのラッパー関数
/// 本クラスインスタンス配下の全PolygonGroupのmoveメソッドが呼び出される。
/// moveメソッドは、polygongroup_set_move_func_c関数で登録する
///
///  @param[in]     param   移動計算パラメータセット
///  @return    POLYLIB_STATで定義される値が返る
///
POLYLIB_STAT 
polylib_move(
               PolylibMoveParamsStruct* param
             )
{
    PolylibMoveParams params;
    params.m_current_step = param->m_current_step;
    params.m_next_step    = param->m_next_step;
    params.m_delta_t      = param->m_delta_t;
    memcpy( params.m_params, param->m_params, 10*sizeof(PL_REAL) );

    return p_polylib_instance->move( params );
}

#ifdef MPI_PL
///
/// Polylib::mograteメソッドのラッパー関数
/// オブジェクトのインスタンス毎に登録が必要
///  @return POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT polylib_migrate( void )
{
    return p_polylib_instance->migrate();
}
#endif

///
/// ポリゴングループタグの取得
///  @param[in]     group_name  ポリゴングループ名
///  @param[out]    tag         ポリゴングループ
///  @return    POLYLIB_STATで定義される値が返る
///
POLYLIB_STAT polylib_get_group_tag(
        char*       group_name,
        PL_GRP_TAG* grp_tag
    )
{
    std::string  name = group_name;

    PolygonGroup* pg=p_polylib_instance->get_group( name );
    if( pg == NULL )  {
       *grp_tag = PL_NULL_TAG;
       return PLSTAT_NG;
    }

    *grp_tag = reinterpret_cast<PL_GRP_TAG>(pg);

    return PLSTAT_OK;
}

///
/// PolygonGroupツリーの最上位ノードの取得（Cインターフェース用）
///      Polylib::get_root_groupsメソッドのラッパー関数
///
///  @param[out]        n       ポリゴングループ数
///  @param[out]        tags    ポリゴングループ　タグ
///  @return    POLYLIB_STATで定義される値が返る
///  @attention   tagsは使用後freeして下さい
///
POLYLIB_STAT polylib_get_root_groups_tags( int* n, PL_GRP_TAG** tags )
{
    std::vector<PolygonGroup *> *grp_list = p_polylib_instance->get_root_groups();
    if( grp_list == NULL )  {
        *n = 0;
        *tags = NULL;
        return PLSTAT_OK;
    }

    *n = grp_list->size();
    *tags = (PL_GRP_TAG*) malloc( (*n)*sizeof(PL_GRP_TAG) );
    
    for(int i=0; i<(*n); i++ ) {
        (*tags)[i] = reinterpret_cast<PL_GRP_TAG>( (*grp_list)[i] );
    }

    delete grp_list;

    return PLSTAT_OK;
}

///
/// リーフPolygonGroupリストの取得（Cインターフェース用）
///      Polylib::get_leaf_groupsメソッドのラッパー関数
///
///  @param[out]        n       PolygonGroup数
///  @param[out]        tags    PolygonGroupタグ
///  @return    POLYLIB_STATで定義される値が返る
///  @attention   tagsは使用後freeして下さい
///
POLYLIB_STAT  polylib_get_leaf_groups_tags( int* n, PL_GRP_TAG** tags )
{
    std::vector<PolygonGroup *> *grp_list = p_polylib_instance->get_leaf_groups();
    if( grp_list == NULL )  {
        *n = 0;
        *tags = NULL;
        return PLSTAT_OK;
    }

    *n = grp_list->size();
    if( *n > 0 ) {
        *tags = (PL_GRP_TAG*) malloc( (*n)*sizeof(PL_GRP_TAG) );
    } else {
        *tags = NULL;
    }
    
    for(int i=0; i<(*n); i++ ) {
        (*tags)[i] = reinterpret_cast<PL_GRP_TAG>( (*grp_list)[i] );
    }

    delete grp_list;

    return PLSTAT_OK;
}

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
    )
{
    POLYLIB_STAT ret;
    std::vector<Triangle*> tri_list_tmp;
    std::string      group_name_tmp = group_name;
    Vec3<PL_REAL>    min_pos_tmp( min_pos );
    Vec3<PL_REAL>    max_pos_tmp( max_pos );
    bool every_tmp;
    if( every == 1 ) {
        every_tmp = true;
    } else {
        every_tmp = false;
    }

    ret = p_polylib_instance->search_polygons (
                      tri_list_tmp, 
                      group_name_tmp, min_pos_tmp, max_pos_tmp, every_tmp
                );
    if( ret != PLSTAT_OK ) {
        *num = 0;
        *tags = NULL;
        return ret;
    }

    *num = tri_list_tmp.size();
    if( *num > 0 ) {
        *tags = (PL_ELM_TAG*) malloc( (*num)*sizeof(PL_ELM_TAG) );
    } else {
        *tags = NULL;
    }
    
    for(int i=0; i<(*num); i++ ) {
        (*tags)[i] = reinterpret_cast<PL_ELM_TAG>( tri_list_tmp[i] );
    }

    return PLSTAT_OK;
}


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
    )
{
    POLYLIB_STAT ret;
    std::vector<Triangle*> tri_list_tmp;
    std::string      group_name_tmp = group_name;
    Vec3<PL_REAL>    min_pos_tmp( min_pos );
    Vec3<PL_REAL>    max_pos_tmp( max_pos );
    bool every_tmp;
    if( every == 1 ) {
        every_tmp = true;
    } else {
        every_tmp = false;
    }

    ret = p_polylib_instance->search_polygons (
                      tri_list_tmp, 
                      group_name_tmp, min_pos_tmp, max_pos_tmp, every_tmp
                );
    if( ret != PLSTAT_OK ) {
        *num = 0;
        *tags = NULL;
        *tris = NULL;
        return ret;
    }

    *num = tri_list_tmp.size();
    if( *num > 0 ) {
        *tags = (PL_ELM_TAG*) malloc( (*num)*sizeof(PL_ELM_TAG) );
        *tris = (TriangleStruct*) malloc( (*num)*sizeof(TriangleStruct) );
    } else {
        *tags = NULL;
        *tris = NULL;
    }
    
    for(int i=0; i<(*num); i++ ) {
        (*tags)[i] = reinterpret_cast<PL_ELM_TAG>( tri_list_tmp[i] );
        Vec3<PL_REAL>* vertexes = tri_list_tmp[i]->get_vertexes();
        Vec3<PL_REAL>  normal   = tri_list_tmp[i]->get_normal();
        // クラス->構造体変換
        VEC3_3_TO_REAL9( vertexes, (*tris)[i].vertex );
        VEC3_TO_REAL   ( normal,   (*tris)[i].normal );
    }

    return PLSTAT_OK;
}

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
    )
{
    POLYLIB_STAT ret;
    std::vector<NptTriangle*> tri_list_tmp;
    std::string      group_name_tmp = group_name;
    Vec3<PL_REAL>    min_pos_tmp( min_pos );
    Vec3<PL_REAL>    max_pos_tmp( max_pos );
    bool every_tmp;
    if( every == 1 ) {
        every_tmp = true;
    } else {
        every_tmp = false;
    }

    ret = p_polylib_instance->search_polygons (
                      tri_list_tmp, 
                      group_name_tmp, min_pos_tmp, max_pos_tmp, every_tmp
                );
    if( ret != PLSTAT_OK ) {
        *num = 0;
        *tags = NULL;
        *tris = NULL;
        return ret;
    }

    *num = tri_list_tmp.size();
    if( *num > 0 ) {
        *tags = (PL_ELM_TAG*) malloc( (*num)*sizeof(PL_ELM_TAG) );
        *tris = (NptTriangleStruct*) malloc( (*num)*sizeof(NptTriangleStruct) );
    } else {
        *tags = NULL;
        *tris = NULL;
    }
    
    for(int i=0; i<(*num); i++ ) {
        (*tags)[i] = reinterpret_cast<PL_ELM_TAG>( tri_list_tmp[i] );
        Vec3<PL_REAL>* vertexes = tri_list_tmp[i]->get_vertexes();
        // クラス->構造体変換
        VEC3_3_TO_REAL9( vertexes, (*tris)[i].vertex );
        polylib_triangle_get_npatchParam(  (*tags)[i],  &((*tris)[i].param) ); 
    }

    return PLSTAT_OK;
}

///
///  指定した点に最も近いポリゴンの検索
///     Polylib::search_nearest_polygonメソッドのラッパー関数。
///  @param[out]    tag         ポリゴンのタグ（ハンドル)
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
    )
{
    POLYLIB_STAT ret;
    Triangle*        pTri;
    std::string      group_name_tmp = group_name;
    Vec3<PL_REAL>    pos_tmp( pos );

    ret = p_polylib_instance->search_nearest_polygon (
                      pTri, 
                      group_name_tmp, pos_tmp
                );
    if( ret != PLSTAT_OK ) {
        *tag = PL_NULL_TAG;
        return ret;
    }
    if( pTri == 0 ) {
        *tag = PL_NULL_TAG;
        return PLSTAT_OK;
    }

    *tag = reinterpret_cast<PL_ELM_TAG>( pTri );

    return PLSTAT_OK;
}

///
/// Polylib::show_group_hierarchyメソッドのラッパー関数。
/// グループ階層構造リストを標準出力に出力する。
///
void polylib_show_group_hierarchy()
{
    p_polylib_instance->show_group_hierarchy();
}


///
/// グループの情報を出力する。(親グループ名、自身の名前、ファイル名、
///      Polylib::show_group_infoメソッドのラッパー関数。
///   登録三角形数、3頂点ベクトルの座標、法線ベクトルの座標、面積)
///  @param[in] group_name グループ名
///  @return POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT polylib_show_group_info(char* group_name)
{
    std::string group_name_tmp = group_name;

    return p_polylib_instance->show_group_info( group_name_tmp );
}


///
///  Polylibが利用中の概算メモリ量を返す
///
/// @return 利用中のメモリ量(byte)
///
size_t polylib_used_memory_size()
{
    return p_polylib_instance->used_memory_size();
}


///
/// Polylibが利用中の概算メモリ量(MB)を返す
///
/// @return 利用中のメモリ量(Mbyte)
///
size_t polylib_used_memory_size_mb()
{
    return ( polylib_used_memory_size()/(1024*1024) );
}


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
         )
{
    *num_tri = 0;
    *tags_tri = NULL;

    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(tag_pg);

    std::vector<Triangle*>* tri_list = pg->get_triangles();
    if( tri_list == NULL ) {
        return PLSTAT_OK;
    }

    *num_tri = tri_list->size();
    *tags_tri = (PL_ELM_TAG*) malloc( (*num_tri)*sizeof(PL_ELM_TAG) );

    for(int i=0; i<(*num_tri); i++ ) {
        (*tags_tri)[i] = reinterpret_cast<PL_ELM_TAG>( (*tri_list)[i] );
    }

    return PLSTAT_OK;
}

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
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(tag_pg);

    return pg->get_group_num_tria();
}

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
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(tag_pg);

    return pg->get_group_num_global_tria();
}
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
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(tag_pg);

    return pg->get_group_area();
}

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
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(tag_pg);

    return pg->get_group_global_area();
}
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
    )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(tag_pg);

    return pg->get_polygons_reduce_atrI( op, atr_no, *val );
}

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
    )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(tag_pg);

    return pg->get_polygons_reduce_atrR( op, atr_no, *val );
}


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
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(tag);
    return pg->set_move_func_c( func );
}

///
/// KD木の再構築フラグの設定
///     ユーザ定義の移動関数内の最後で呼び出す
/// 
///  @param[in] tag       PolygonGroupを操作するためのタグ
///  @return    戻り値なし
/// 
void polylib_group_set_need_rebuild( PL_GRP_TAG tag )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(tag);
    return pg->set_need_rebuild();
}

//----------------------------------------------
// setter / getter
//----------------------------------------------

/// グループ名設定
void polylib_group_set_name(
            PL_GRP_TAG tag,
            char* name
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(tag);
    std::string name_tmp = name;
    pg->set_name( name );
}

/// グループ名取得
void polylib_group_get_name(
            PL_GRP_TAG tag,
            char* name
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(tag);
    std::string name_tmp = pg->get_name();
    strcpy( name, name_tmp.c_str() );
}

/// 移動対象フラグ設定
///  @param[in]  tag       PolygonGroupを操作するためのタグ
///  @param[in]  moval     移動対象フラグ 1:true  0:false
///  @return   戻り値なし
void polylib_group_set_movable(
            PL_GRP_TAG tag,
            int        movable
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(tag);
    bool moval_tmp;
    if( movable == 1 ) {
        moval_tmp = true;
    } else {
        moval_tmp = false;
    }

    pg->set_movable( moval_tmp );
}

/// 移動対象フラグ取得
///  @param[in]  tag       PolygonGroupを操作するためのタグ
///  @param[out] moval     移動対象フラグ 1:true  0:false
void polylib_group_get_movable(
            PL_GRP_TAG tag,
            int*       movable
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(tag);
    bool moval_tmp;

    moval_tmp = pg->get_movable();
    if( moval_tmp == true ) {
        *movable = 1;
    } else {
        *movable = 0;
    }
}

/// PolygonGroup ユーザ定義属性の設定
///  @param[in]     key     キー
///  @param[in]     val     属性値
///  @return なし
///  @attention 既に登録されていた場合、上書きする
void polylib_group_set_atr(
            PL_GRP_TAG tag,
            char* key,
            char* val
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(tag);
    std::string key_tmp = key;
    std::string val_tmp = val;

    pg->set_atr( key_tmp, val_tmp );
}

/// PolygonGroup ユーザ定義属性の取得
///  @param[in]     key     キー
///  @param[out]    val     属性値
///  @return OK/NG  NG:キーと属性のペアが登録されていない
POLYLIB_STAT polylib_group_get_atr(
            PL_GRP_TAG tag,
            char* key,
            char* val
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(tag);
    std::string key_tmp;
    std::string val_tmp;

    POLYLIB_STAT ret = pg->get_atr( key_tmp, val_tmp );
    
    if( ret == PLSTAT_OK ) {
        strcpy( val, val_tmp.c_str() );
    }
    return ret;
}

/// ポリゴングループ内のポリゴンのユーザ定義属性数の設定
void polylib_group_set_num_polygon_atr(
            PL_GRP_TAG tag,
            int num_atrI,
            int num_atrR
            )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(tag);
    pg->set_num_polygon_atr( num_atrI,num_atrR );
}

///
/// PolygonGroup 親グループを取得
/// 
///  @param[in]     tag     PolygonGroupを操作するためのタグ
///  @return 親グループのタグ
/// 
PL_GRP_TAG polylib_group_get_parent( PL_GRP_TAG tag )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(tag);
    PolygonGroup* parent = pg->get_parent();
    return  reinterpret_cast<PL_GRP_TAG>(parent);
}

///
/// 子グループを取得
///  @param[in]     tag     PolygonGroupを操作するためのタグ
///  @param[out]    num     ポリゴングループ数
///  @param[out]    child_tags ポリゴングループ　タグ
///  @return 戻り値なし
///  @attention   tagsは使用後freeして下さい
///
void polylib_group_get_children(
                  PL_GRP_TAG tag,
                  int* num,
                  PL_GRP_TAG** child_tags
              ) 

{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(tag);
    std::vector<PolygonGroup*>& pg_child = pg->get_children();

    *num = pg_child.size();
    if( (*num) == 0 ) {
        *child_tags = NULL;
        return;
    }

    *child_tags = (PL_GRP_TAG*) malloc( (*num)*sizeof(PL_GRP_TAG) );
    
    for(int i=0; i<(*num); i++ ) {
        (*child_tags)[i] = reinterpret_cast<PL_GRP_TAG>( pg_child[i] );
    }
}

//**********************************************************
//  Triangleクラス用ラッパー関数
//      三角形オブジェクトと長田パッチオブジェクトでは動作が違う場合あり
//**********************************************************

/// ポリゴンタイプ取得
///  @param[in]     tag     Triangle/NptTriangleを操作するためのタグ
///  @return PL_TYPE_TRIANGLE / PL_TYPE_NPT（長田パッチ）
int polylib_triangle_get_pl_type(
            PL_ELM_TAG tag
         )
{
    Triangle* tri = reinterpret_cast<Triangle*>(tag);

    return tri->get_pl_type();
}

/// 頂点座標設定
///   基本はTriangle用
            // 頂点設定時、法線ベクトル,面積も内部で設定する
            // 長田パッチのパラメータは更新されないので注意
void polylib_triangle_set_vertexes(
            PL_ELM_TAG tag,
            PL_REAL vertex[9]
         )
{
    Triangle* tri = reinterpret_cast<Triangle*>(tag);
    Vec3<PL_REAL> vertex_tmp[3];
    REAL9_TO_VEC3_3(vertex, vertex_tmp);

    tri->set_vertexes( vertex_tmp, true, true );
}

/// 頂点座標・長田パッチパラメータ設定
            // 法線ベクトル,面積も内部で設定する
void polylib_triangle_set_vertexes_npatch(
            PL_ELM_TAG tag,
            PL_REAL vertex[9],
            NpatchParamStruct* param
         )
{
    NptTriangle* tri = reinterpret_cast<NptTriangle*>(tag);
    Vec3<PL_REAL> vertex_tmp[3];
    NpatchParam param_tmp;

    REAL9_TO_VEC3_3(vertex, vertex_tmp);
    REAL_TO_VEC3(param->cp_side1_1, param_tmp.cp_side1_1 );
    REAL_TO_VEC3(param->cp_side1_2, param_tmp.cp_side1_2 );
    REAL_TO_VEC3(param->cp_side2_1, param_tmp.cp_side2_1 );
    REAL_TO_VEC3(param->cp_side2_2, param_tmp.cp_side2_2 );
    REAL_TO_VEC3(param->cp_side3_1, param_tmp.cp_side3_1 );
    REAL_TO_VEC3(param->cp_side3_2, param_tmp.cp_side3_2 );
    REAL_TO_VEC3(param->cp_center , param_tmp.cp_center  );
    
    tri->set_vertexes( vertex_tmp, param_tmp, true, true );
}

/// 頂点座標取得
void polylib_triangle_get_vertexes(
            PL_ELM_TAG tag,
            PL_REAL vertex[9]
         )
{
    Triangle* tri = reinterpret_cast<Triangle*>(tag);
    Vec3<PL_REAL>* vertex_tmp = tri->get_vertexes();

    VEC3_3_TO_REAL9( vertex_tmp, vertex );
}

/// 法線ベクトル取得
void polylib_triangle_get_normal(
            PL_ELM_TAG tag,
            PL_REAL norm[3]
         )
{
    Triangle* tri = reinterpret_cast<Triangle*>(tag);
    Vec3<PL_REAL> normal_tmp = tri->get_normal();

    VEC3_TO_REAL( normal_tmp, norm );
}

/// 長田パッチパラメータ設定
/// オブジェクトが長田パッチでない時はエラーを返す
POLYLIB_STAT polylib_triangle_set_npatchParam(
            PL_ELM_TAG tag,
            NpatchParamStruct *param
         )
{
    NptTriangle* tri = reinterpret_cast<NptTriangle*>(tag);
    int pl_type = tri->get_pl_type();
    if( pl_type != PL_TYPE_NPT ) {
        return PLSTAT_NG;
    }

    NpatchParam param_tmp;

    REAL_TO_VEC3(param->cp_side1_1, param_tmp.cp_side1_1 );
    REAL_TO_VEC3(param->cp_side1_2, param_tmp.cp_side1_2 );
    REAL_TO_VEC3(param->cp_side2_1, param_tmp.cp_side2_1 );
    REAL_TO_VEC3(param->cp_side2_2, param_tmp.cp_side2_2 );
    REAL_TO_VEC3(param->cp_side3_1, param_tmp.cp_side3_1 );
    REAL_TO_VEC3(param->cp_side3_2, param_tmp.cp_side3_2 );
    REAL_TO_VEC3(param->cp_center , param_tmp.cp_center  );
    
    tri->set_npatch_param( param_tmp );
    return PLSTAT_OK;
}

/// 長田パッチパラメータ取得
            // オブジェクトが長田パッチでない時はエラーを返す
POLYLIB_STAT polylib_triangle_get_npatchParam(
            PL_ELM_TAG tag,
            NpatchParamStruct *param
         )
{
    NptTriangle* tri = reinterpret_cast<NptTriangle*>(tag);
    int pl_type = tri->get_pl_type();
    if( pl_type != PL_TYPE_NPT ) {
        return PLSTAT_NG;
    }

    NpatchParam* param_tmp = tri->get_npatch_param();

    VEC3_TO_REAL( param_tmp->cp_side1_1, param->cp_side1_1 );
    VEC3_TO_REAL( param_tmp->cp_side1_2, param->cp_side1_2 );
    VEC3_TO_REAL( param_tmp->cp_side2_1, param->cp_side2_1 );
    VEC3_TO_REAL( param_tmp->cp_side2_2, param->cp_side2_2 );
    VEC3_TO_REAL( param_tmp->cp_side3_1, param->cp_side3_1 );
    VEC3_TO_REAL( param_tmp->cp_side3_2, param->cp_side3_2 );
    VEC3_TO_REAL( param_tmp->cp_center , param->cp_center  );

    return PLSTAT_OK;
}


/// ユーザ定義属性数（整数型）の取得
int polylib_triangle_get_num_atrI(
                PL_ELM_TAG tag
            )
{
    Triangle* tri = reinterpret_cast<Triangle*>(tag);
    return tri->get_num_atrI();
}

/// ユーザ定義属性数（実数型）の取得
int polylib_triangle_get_num_atrR(
                PL_ELM_TAG tag
            )
{
    Triangle* tri = reinterpret_cast<Triangle*>(tag);
    return tri->get_num_atrR();
}

/// ユーザ定義属性（整数型）のポインタ取得
int* polylib_triangle_get_pAtrI(
            PL_ELM_TAG tag
            )
{
    Triangle* tri = reinterpret_cast<Triangle*>(tag);
    return tri->get_pAtrI();
}


/// ユーザ定義属性（実数型）のポインタ取得
PL_REAL* polylib_triangle_get_pAtrR(
            PL_ELM_TAG tag
            )
{
    Triangle* tri = reinterpret_cast<Triangle*>(tag);
    return tri->get_pAtrR();
}

