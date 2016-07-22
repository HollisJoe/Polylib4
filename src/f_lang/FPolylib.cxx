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

#include "f_lang/FPolylib.h"
#include "c_lang/CPolylib.h"
#include "Polylib.h"

using namespace std;
using namespace PolylibNS;

/************************************************************************
 *
 * Fortran用Polylib
 *
 ***********************************************************************/

///
///   Polylibインスタンス(静的変数）
///     CPolylib.cxx で実態を定義
extern Polylib* p_polylib_instance;

//--------------------------------------------------
//  内部関数
//--------------------------------------------------

///
/// Fortranの文字列をC++の文字列に変換する
///   日本語非対応
///
///  @param[in]  fstring      Fortran文字列
///                              charcter*len_f
///                              文字列の有効長以降はspaceが詰められてことが前提
///  @param[in]  len_f        Fortran文字列の長さ
///  @param[out] cstring      C++文字列
///  @param[out] len_c        C++文字列の長さ
///                                終端の'\0' 含まず
///
void fpolylib_fstring_to_cstring (
               char*        fstring,
               int          len_f,
               std::string& cstring,
               int&         len_c
       )
{
    int ip;
    char* cstring_tmp = new char[len_f+1];

    cstring_tmp[0] = '\0';

    for( ip=0; ip<len_f; ip++ ) { 
        if( fstring[ip] == ' ' ) {
            break;
        } else if ( fstring[ip] == '\0' ) {
            break;
        }
        cstring_tmp[ip] = fstring[ip];
    }
    cstring_tmp[ip] = '\0';

    cstring = cstring_tmp;
    len_c   = ip;

    delete[] cstring_tmp;
}

///
/// C++の文字列をFortranの文字列に変換する
///   日本語非対応
///
///  @param[in]  cstring      C++文字列
///  @param[out] fstring      Fortran文字列
///                              charcter*len_f
///                              文字列の有効長以降はspaceを詰める
///  @param[in]  len_f        Fortran文字列の長さ
///
void fpolylib_cstring_to_fstring (
               std::string& cstring,
               char*        fstring,
               int          len_f 
       )
{
    const char* p     = cstring.c_str();
    int   len_c = strlen( p );

    memcpy( fstring, p, len_c );
    memset( &fstring[len_c], ' ', (len_f-len_c) );
}



//--------------------------------------------------
//  外部公開関数
//--------------------------------------------------

/// Fortran用Polylib環境の構築
///     Polylibインスタンス生成
///
///  @param[out] ret   POLYLIB_STATで定義される値が返る。
///  @attention 最初に呼び出すこと
///
void  fpolylib_instance_ ( POLYLIB_STAT* ret )
{
    *ret = polylib_instance();  // C用呼び出し
}


#ifdef MPI_PL
///
/// 並列計算関連情報の設定と初期化を行う。
/// Polylib::init_parallel_infoメソッドのFortranラッパー関数。
    /// (各ランクが１領域を担当している場合）
/// 全rankで各々設定を行い、その領域情報を全rankへ配信する。
///
///  @param[in]  bpos    自PE担当領域の基点座標
///  @param[in]  bbsize  同、計算領域のボクセル数
///  @param[in]  gcsize  同、ガイドセルのボクセル数
///  @param[in]  dx      同、ボクセル１辺の長さ
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void
fpolylib_init_parallel_info_ (
                              PL_REAL bpos[3],
                              int bbsize[3],
                              int gcsize[3],
                              PL_REAL dx[3],
                              POLYLIB_STAT* ret
    )
{
    unsigned int bbsize_uint[3];
    unsigned int gcsize_uint[3];
    bbsize_uint[0]=bbsize[0];  bbsize_uint[1]=bbsize[1];  bbsize_uint[2]=bbsize[2];
    gcsize_uint[0]=gcsize[0];  gcsize_uint[1]=gcsize[1];  gcsize_uint[2]=gcsize[2];

    *ret = p_polylib_instance->init_parallel_info( MPI_COMM_WORLD, bpos, bbsize_uint, gcsize_uint, dx );
}

///
/// 並列計算関連情報の設定と初期化を行う。
/// (各ランクが複数領域を担当している場合）
/// 全rankで各々設定を行い、その領域情報を全rankへ配信する。
///
///  @param[in] num     自PEが担当する領域
///  @param[in] bbox    担当するboundary box情報（複数）
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void
fpolylib_init_parallel_info2_ (
                 int*      num,
                 FParallelBboxStruct *bboxes_s,
                 POLYLIB_STAT* ret
    )
{
    //
    std::vector<ParallelBbox>  bboxes;
    for(int i=0; i<*num; i++ ) {
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

    *ret = p_polylib_instance->init_parallel_info( MPI_COMM_WORLD, bboxes );
}
#endif

///
/// 引数で指定された設定ファイルを読み込み、グループツリーを作成する。
/// Polylib::loadメソッドのFortran用ラッパー関数。
/// 続いて設定ファイルで指定されたSTLファイルを読み込み、KD木を作成する。
///  @param[in] config_name  設定ファイル名
///                     Fortran型の文字列  (\0で終了しない）
///                     character*256 長さPL_FILE_PATH_LEN 
///                     Fortran側から呼ぶときは、文字列の最初の空白の１個前までを
///                     有効なファイル名と見なします。
///                     すべて空白の時は"polylib_config.tp"とみなす
///  @param[in] scale        縮尺率
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void
fpolylib_load_ (
                char*   config_name,
                PL_REAL* scale,
                POLYLIB_STAT* ret
            )
{
    std::string config_file_name;
    int len_c;
    char* config_name_tmp = new char[PL_FILE_PATH_LEN+1];

    fpolylib_fstring_to_cstring( 
             config_name, PL_FILE_PATH_LEN,
             config_file_name, len_c
          );
     
    if( len_c  == 0 ) {
        config_file_name = "polylib_config.tp";
    }

    *ret = p_polylib_instance->load( config_file_name, *scale );

    delete[] config_name_tmp;
}


///
/// PolygoGroupツリー、三角形ポリゴン情報の保存。
///    Polylib::saveメソッドのラッパー関数。
/// グループツリーの情報を設定ファイルへ出力。三角形ポリゴン情報をSTL/NPT
/// ファイルへ出力
///
///  @param[out]    o_fname 設定ファイル名
///                             character*256
///  @param[in]     format  形状ファイルのフォーマット
///                             PolylibDefine.hで定義されているFILE_FMT_* 参照
///                             character*8
///  @param[in]     extend  ファイル名に付加する文字列。NULLを指定した
///                         場合は、付加文字列として本メソッド呼び出し時の
///                         年月日時分秒(YYYYMMDD24hhmmss)を用いる
///                             character*32
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///  @attention ファイル名命名規約は次の通り。
///         設定ファイル : polylib_config_付加文字.tpp
///         STL/NPTファイル  : ポリゴングループ名_付加文字.拡張子
///
void
fpolylib_save_ (
               char  *o_fname,
               char  *format,
               char  *extend,
               POLYLIB_STAT* ret
             )
{
    std::string  config_name_out;
    std::string  format_wk;
    std::string  extend_wk;
    int len_c;

    fpolylib_fstring_to_cstring ( format, PL_FORMAT_LEN, format_wk, len_c );
    fpolylib_fstring_to_cstring ( extend, PL_STR_LEN, extend_wk, len_c );
    
    *ret = p_polylib_instance->save( config_name_out, format_wk, extend_wk );

    if( *ret == PLSTAT_OK ) {
        fpolylib_cstring_to_fstring ( config_name_out,
                                      o_fname, PL_FILE_PATH_LEN );
    }
}


///
/// 三角形ポリゴン座標の移動。
///   Polylib::moveメソッドのラッパー関数
/// 本クラスインスタンス配下の全PolygonGroupのmoveメソッドが呼び出される。
///
///  @param[in]     param   移動計算パラメータセット
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void
fpolylib_move_ (
               PolylibMoveParamsStruct* param,
                          POLYLIB_STAT* ret
             )
{
    PolylibMoveParams params;
    params.m_current_step = param->m_current_step;
    params.m_next_step    = param->m_next_step;
    params.m_delta_t      = param->m_delta_t;
    memcpy( params.m_params, param->m_params, 10*sizeof(PL_REAL) );

    *ret = p_polylib_instance->move( params );
}

#ifdef MPI_PL
///
/// Polylib::mograteメソッドのラッパー関数
/// オブジェクトのインスタンス毎に登録が必要
///  @return POLYLIB_STATで定義される値が返る。
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_migrate_ ( POLYLIB_STAT* ret )
{
    *ret = p_polylib_instance->migrate();
}
#endif

///
/// ポリゴングループタグの取得
///  @param[in]     group_name  ポリゴングループ名
///  @param[out]    tag         ポリゴングループ
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_get_group_tag_ (
        char*       group_name,
        PL_GRP_TAG* grp_tag,
        POLYLIB_STAT* ret
    )
{
    std::string  name;
    int          len_c;

    fpolylib_fstring_to_cstring (
               group_name, PL_GRP_PATH_LEN,
               name,       len_c
         );

    PolygonGroup* pg=p_polylib_instance->get_group( name );
    if( pg == NULL )  {
       *grp_tag = PL_NULL_TAG;
       *ret = PLSTAT_NG;
       return;
    }

    *grp_tag = reinterpret_cast<PL_GRP_TAG>(pg);

    *ret = PLSTAT_OK;
}

///
/// PolygonGroupツリーの最上位ノード数の取得
///      Polylib::get_root_groupsメソッドのラッパー関数
///
///  @param[out]        n       ポリゴングループ数
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_get_root_groups_tags_num_ ( int* n, POLYLIB_STAT* ret )
{
    std::vector<PolygonGroup *> *grp_list = p_polylib_instance->get_root_groups();
    if( grp_list == NULL )  {
        *n = 0;
        *ret = PLSTAT_OK;
        return;
    }

    *n = grp_list->size();

    delete grp_list;

    *ret = PLSTAT_OK;
}

///
/// PolygonGroupツリーの最上位ノードの取得
///      Polylib::get_root_groupsメソッドのラッパー関数
///
///  @param[out]        n       ポリゴングループ数
///  @param[out]        tags    ポリゴングループ　タグ
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_get_root_groups_tags_ ( int* n, PL_GRP_TAG* tags, POLYLIB_STAT* ret )
{
    std::vector<PolygonGroup *> *grp_list = p_polylib_instance->get_root_groups();
    if( grp_list == NULL )  {
        *n = 0;
        *ret = PLSTAT_OK;
        return;
    }

    *n = grp_list->size();
    
    for(int i=0; i<(*n); i++ ) {
        tags[i] = reinterpret_cast<PL_GRP_TAG>( (*grp_list)[i] );
    }

    delete grp_list;

    *ret = PLSTAT_OK;
}

///
/// リーフPolygonGroupリストの個数取得
///
///  @param[out]        n       PolygonGroup数
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_get_leaf_groups_tags_num_ ( int* n, POLYLIB_STAT* ret )
{
    std::vector<PolygonGroup *> *grp_list = p_polylib_instance->get_leaf_groups();
    if( grp_list == NULL )  {
        *n = 0;
        *ret = PLSTAT_OK;
        return;
    }

    *n = grp_list->size();
    delete grp_list;

    *ret = PLSTAT_OK;
}

///
/// リーフPolygonGroupリストの取得（Cインターフェース用）
///      Polylib::get_leaf_groupsメソッドのラッパー関数
///
///  @param[out]        n       PolygonGroup数
///  @param[out]        tags    PolygonGroupタグ
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_get_leaf_groups_tags_ ( int* n, PL_GRP_TAG* tags, POLYLIB_STAT* ret )
{
    std::vector<PolygonGroup *> *grp_list = p_polylib_instance->get_leaf_groups();
    if( grp_list == NULL )  {
        *n = 0;
        *ret = PLSTAT_OK;
        return;
    }

    *n = grp_list->size();
    
    for(int i=0; i<(*n); i++ ) {
        tags[i] = reinterpret_cast<PL_GRP_TAG>( (*grp_list)[i] );
    }

    delete grp_list;

    *ret = PLSTAT_OK;
}

///
/// 位置ベクトルmin_posとmax_posにより特定される矩形領域に含まれる、
/// ポリゴン数を取得する
///  @param[out]    num         抽出された三角形ポリゴン数
///  @param[in]     group_name  抽出グループ名。
///  @param[in]     min_pos     抽出する矩形領域の最小値。(x,y,z順の配列)
///  @param[in]     max_pos     抽出する矩形領域の最大値。(x,y,z順の配列)
///  @param[in]     every       抽出オプション。
///   1：3頂点が全て検索領域に含まれるポリゴンを抽出する。
///   0：三角形のBBoxが一部でも検索領域と交差するものを抽出する
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_search_polygons_num_ (
    int         *num,
    char*       group_name,
    PL_REAL     min_pos[3],
    PL_REAL     max_pos[3],
    int         every,
    POLYLIB_STAT* ret
    )
{
    std::vector<Triangle*> tri_list_tmp;
    std::string      group_name_tmp; int len_c;
    fpolylib_fstring_to_cstring ( group_name, PL_GRP_PATH_LEN, 
                                  group_name_tmp, len_c );

    Vec3<PL_REAL>    min_pos_tmp( min_pos );
    Vec3<PL_REAL>    max_pos_tmp( max_pos );
    bool every_tmp;
    if( every == 1 ) {
        every_tmp = true;
    } else {
        every_tmp = false;
    }

    *ret = p_polylib_instance->search_polygons (
                      tri_list_tmp, 
                      group_name_tmp, min_pos_tmp, max_pos_tmp, every_tmp
                );
    if( *ret != PLSTAT_OK ) {
        *num = 0;
        return;
    }

    *num = tri_list_tmp.size();
    
    *ret = PLSTAT_OK;
}


///
/// Polylib::search_polygonsメソッドのラッパー関数。
/// 位置ベクトルmin_posとmax_posにより特定される矩形領域に含まれる、
/// 特定のグループとその子孫グループに属する三角形ポリゴンをKD探索に
/// より抽出する。
///  @param[out]    num         抽出された三角形ポリゴン数
///  @param[out]    tags        三角形ポリゴンのタグ（ハンドル)
///  @param[in]     group_name  抽出グループ名。
///  @param[in]     min_pos     抽出する矩形領域の最小値。(x,y,z順の配列)
///  @param[in]     max_pos     抽出する矩形領域の最大値。(x,y,z順の配列)
///  @param[in]     every       抽出オプション。
///   1：3頂点が全て検索領域に含まれるポリゴンを抽出する。
///   0：三角形のBBoxが一部でも検索領域と交差するものを抽出する
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///  @attention MPI並列計算時は,min_pos, max_posは各ランクの矩形領域を
///             超えないようにして下さい (各ランク内の担当領域内のみ検索するため）
///
void fpolylib_search_polygons_ (
    int         *num,
    PL_ELM_TAG  *tags,
    char*       group_name,
    PL_REAL     min_pos[3],
    PL_REAL     max_pos[3],
    int         every,
    POLYLIB_STAT* ret
    )
{
    std::vector<Triangle*> tri_list_tmp;
    std::string      group_name_tmp; int len_c;
    fpolylib_fstring_to_cstring ( group_name, PL_GRP_PATH_LEN, 
                                  group_name_tmp, len_c );
    Vec3<PL_REAL>    min_pos_tmp( min_pos );
    Vec3<PL_REAL>    max_pos_tmp( max_pos );
    bool every_tmp;
    if( every == 1 ) {
        every_tmp = true;
    } else {
        every_tmp = false;
    }

    *ret = p_polylib_instance->search_polygons (
                      tri_list_tmp, 
                      group_name_tmp, min_pos_tmp, max_pos_tmp, every_tmp
                );
    if( *ret != PLSTAT_OK ) {
        *num = 0;
        return;
    }

    *num = tri_list_tmp.size();
    
    for(int i=0; i<(*num); i++ ) {
        tags[i] = reinterpret_cast<PL_ELM_TAG>( tri_list_tmp[i] );
    }

    *ret = PLSTAT_OK;
}

///
///  指定した点に最も近いポリゴンの検索
///     Polylib::search_nearest_polygonメソッドのラッパー関数。
///  @param[out]    tag         ポリゴンのタグ（ハンドル)
///  @param[in]     group_name  グループ名。
///  @param[in]     pos         指定点
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_search_nearest_polygon_ (
        PL_ELM_TAG  *tag,
        char*       group_name,
        PL_REAL     pos[3],
        POLYLIB_STAT* ret
    )
{
    Triangle*        pTri;
    std::string      group_name_tmp; int len_c;
    fpolylib_fstring_to_cstring ( group_name, PL_GRP_PATH_LEN, 
                                  group_name_tmp, len_c );
    Vec3<PL_REAL>    pos_tmp( pos );

    *ret = p_polylib_instance->search_nearest_polygon (
                      pTri, 
                      group_name_tmp, pos_tmp
                );
    if( *ret != PLSTAT_OK ) {
        *tag = PL_NULL_TAG;
        return;
    }
    if( pTri == 0 ) {
        *tag = PL_NULL_TAG;
        *ret = PLSTAT_OK;
        return;
    }

    *tag = reinterpret_cast<PL_ELM_TAG>( pTri );

    *ret = PLSTAT_OK;
}

///
/// Polylib::show_group_hierarchyメソッドのラッパー関数。
/// グループ階層構造リストを標準出力に出力する。
///
void fpolylib_show_group_hierarchy_ ()
{
    p_polylib_instance->show_group_hierarchy();
}


///
/// グループの情報を出力する。(親グループ名、自身の名前、ファイル名、
///      Polylib::show_group_infoメソッドのラッパー関数。
///   登録三角形数、3頂点ベクトルの座標、法線ベクトルの座標、面積)
///  @param[in] group_name グループ名
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_show_group_info_ (char* group_name, POLYLIB_STAT* ret )
{
    std::string group_name_tmp; int len_c;
    fpolylib_fstring_to_cstring ( group_name, PL_GRP_PATH_LEN, 
                                  group_name_tmp, len_c );

    *ret = p_polylib_instance->show_group_info( group_name_tmp );
}


///
///  Polylibが利用中の概算メモリ量を返す
///
/// @return 利用中のメモリ量(byte)
///
int fpolylib_used_memory_size_ ()
{
    return p_polylib_instance->used_memory_size();
}


///
/// Polylibが利用中の概算メモリ量(MB)を返す
///
/// @return 利用中のメモリ量(Mbyte)
///
int fpolylib_used_memory_size_mb_ ()
{
    unsigned int sz = p_polylib_instance->used_memory_size();
    return ( sz/(1024*1024) );
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
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_group_get_triangles_ (
            PL_GRP_TAG* tag_pg,
            int        *num_tri,
            PL_ELM_TAG  *tags_tri,
            POLYLIB_STAT* ret
         )
{
    *num_tri = 0;

    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(*tag_pg);

    std::vector<Triangle*>* tri_list = pg->get_triangles();
    if( tri_list == NULL ) {
        *ret = PLSTAT_OK;
        return;
    }

    *num_tri = tri_list->size();

    for(int i=0; i<(*num_tri); i++ ) {
        tags_tri[i] = reinterpret_cast<PL_ELM_TAG>( (*tri_list)[i] );
    }

    *ret = PLSTAT_OK;
}

///
/// PolygonGroupのポリゴン数を求める
///    PolygonGroup::get_trianglesのラッパー関数
///  @param[in]   tag_pg       PolygonGroupを操作するためのタグ
///  @param[out]  num_tri      三角形ポリゴン数
///  @param[out]  tags_tri     三角形ポリゴンのタグ（ハンドル)
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_group_get_num_triangles_ (
            PL_GRP_TAG* tag_pg,
            int        *num_tri,
            POLYLIB_STAT* ret
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(*tag_pg);
    *num_tri = pg->get_group_num_tria();
    *ret = PLSTAT_OK;
}

#ifdef MPI_PL
/// PolygonGroupのポリゴン数(global)を求める
///    PolygonGroup::get_trianglesのラッパー関数
///  @param[in]   tag_pg       PolygonGroupを操作するためのタグ
///                                integer*8
///  @param[out]  num_tri      三角形ポリゴン数(global)
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_group_get_num_global_triangles_ (
            PL_GRP_TAG *tag_pg,
            int        *num_tri,
            POLYLIB_STAT* ret
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(*tag_pg);
    *num_tri = pg->get_group_num_global_tria();
    *ret = PLSTAT_OK;
}
#endif

/// PolygonGroupのポリゴンの面積を求める
///    PolygonGroup::get_trianglesのラッパー関数
///  @param[in]   tag_pg       PolygonGroupを操作するためのタグ
///                                integer*8
///  @param[out]  area         ポリゴンの面積
///  @param[out]  ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_group_get_area_ (
            PL_GRP_TAG *tag_pg,
            PL_REAL    *area,
            POLYLIB_STAT* ret
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(*tag_pg);
    *area = pg->get_group_area();
    *ret = PLSTAT_OK;
}

#ifdef MPI_PL
/// PolygonGroupのポリゴンの面積(global)を求める
///    PolygonGroup::get_trianglesのラッパー関数
///  @param[in]   tag_pg       PolygonGroupを操作するためのタグ
///                                integer*8
///  @param[out]  area         ポリゴンの面積(global)
///  @param[out]  ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_group_get_global_area_ (
            PL_GRP_TAG *tag_pg,
            PL_REAL    *area,
            POLYLIB_STAT* ret
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(*tag_pg);
    *area = pg->get_group_global_area();
    *ret = PLSTAT_OK;
}
#endif

/// グループ内のポリゴン属性（整数）の集合演算値を返す
///     並列化されている場合は全プロセスを通した値
///     (PL_OP_SUM：重複ポリゴン分は無視される）
///     全プロセスに同じ値が返る
///  @param[in]  tag_pg       PolygonGroupを操作するためのタグ
///  @param[in]  op     演算種類　PL_OP_SUM/PL_OP_MAX/PL_OP_MIN
///  @param[in]  atr_no ポリゴン整数属性の何番目か　0〜
///  @param[out] val    属性値
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_group_get_polygons_reduce_atrI_ (
            PL_GRP_TAG*   tag_pg,
            PL_OP_TYPE*   op,
            int*          atr_no,
            int*          val,
            POLYLIB_STAT* ret
    )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(*tag_pg);
    *ret = pg->get_polygons_reduce_atrI( *op, *atr_no, *val );
}

/// グループ内のポリゴン属性（実数）の集合演算値を返す
///     並列化されている場合は全プロセスを通した値
///     (PL_OP_SUM：重複ポリゴン分は無視される）
///     全プロセスに同じ値が返る
///  @param[in]  tag_pg       PolygonGroupを操作するためのタグ
///  @param[in]  op     演算種類　PL_OP_SUM/PL_OP_MAX/PL_OP_MIN
///  @param[in]  atr_no ポリゴン実数属性の何番目か　0〜
///  @param[out] val    属性値
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_group_get_polygons_reduce_atrR_ (
            PL_GRP_TAG*   tag_pg,
            PL_OP_TYPE*   op,
            int*          atr_no,
            PL_REAL*      val,
            POLYLIB_STAT* ret
    )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(*tag_pg);
    *ret = pg->get_polygons_reduce_atrR( *op, *atr_no, *val );
}

///
/// KD木の再構築フラグの設定
///     ユーザ定義の移動関数内の最後で呼び出す
/// 
///  @param[in] tag       PolygonGroupを操作するためのタグ
///  @return    戻り値なし
/// 
void fpolylib_group_set_need_rebuild_ ( PL_GRP_TAG* tag_pg )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(*tag_pg);
    pg->set_need_rebuild();
}


//----------------------------------------------
// setter / getter
//----------------------------------------------

/// グループ名設定
void fpolylib_group_set_name_ (
            PL_GRP_TAG* tag,
            char* name
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(*tag);
    std::string group_name_tmp; int len_c;
    fpolylib_fstring_to_cstring ( name, PL_GRP_PATH_LEN, 
                                  group_name_tmp, len_c );

    pg->set_name( group_name_tmp );
}

/// グループ名取得
void fpolylib_group_get_name_ (
            PL_GRP_TAG* tag,
            char* name
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(*tag);
    std::string name_tmp = pg->get_name();

    fpolylib_cstring_to_fstring ( name_tmp, name, PL_GRP_PATH_LEN );
}

/// 移動対象フラグ設定
///  @param[in]  tag       PolygonGroupを操作するためのタグ
///  @param[in]  moval     移動対象フラグ 1:true  0:false
///  @return   戻り値なし
void fpolylib_group_set_movable_ (
            PL_GRP_TAG* tag,
            int*        movable
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(*tag);
    bool moval_tmp;
    if( *movable == 1 ) {
        moval_tmp = true;
    } else {
        moval_tmp = false;
    }

    pg->set_movable( moval_tmp );
}

/// 移動対象フラグ取得
///  @param[in]  tag       PolygonGroupを操作するためのタグ
///  @param[out] moval     移動対象フラグ 1:true  0:false
void fpolylib_group_get_movable_ (
            PL_GRP_TAG* tag,
            int*        movable
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(*tag);
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
void fpolylib_group_set_atr_ (
            PL_GRP_TAG* tag,
            char* key,
            char* val
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(*tag);
    std::string key_tmp;
    std::string val_tmp;
    int len_c;
    fpolylib_fstring_to_cstring ( key, PL_GRP_ATR_LEN, 
                                  key_tmp, len_c );
    fpolylib_fstring_to_cstring ( val, PL_GRP_ATR_LEN, 
                                  val_tmp, len_c );

    pg->set_atr( key_tmp, val_tmp );
}

/// PolygonGroup ユーザ定義属性の取得
///  @param[in]     key     キー
///  @param[out]    val     属性値
///  @param[out]    ret     POLYLIB_STAT(integer)で定義される値が返る
///                          NG:キーと属性のペアが登録されていない
void fpolylib_group_get_atr_ (
            PL_GRP_TAG* tag,
            char* key,
            char* val,
            POLYLIB_STAT* ret
         )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(*tag);
    std::string key_tmp;
    std::string val_tmp;
    int len_c;
    fpolylib_fstring_to_cstring ( key, PL_GRP_ATR_LEN, 
                                  key_tmp, len_c );

    *ret = pg->get_atr( key_tmp, val_tmp );
    
    if( *ret == PLSTAT_OK ) {
        fpolylib_cstring_to_fstring ( val_tmp, val, PL_GRP_ATR_LEN );
    }
}

/// ポリゴングループ内のポリゴンのユーザ定義属性数の設定
void fpolylib_group_set_num_polygon_atr_ (
            PL_GRP_TAG* tag,
            int* num_atrI,
            int* num_atrR
            )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(*tag);
    pg->set_num_polygon_atr( *num_atrI,*num_atrR );
}


///
/// PolygonGroup 親グループを取得
/// 
///  @param[in]  tag         PolygonGroupを操作するためのタグ
///  @param[out] parent_tag  親グループのタグ  integer*8
///  @return 親グループのタグ
/// 
void fpolylib_group_get_parent_ ( PL_GRP_TAG* tag, PL_GRP_TAG* parent_tag )
{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(*tag);
    PolygonGroup* parent = pg->get_parent();
    *parent_tag = reinterpret_cast<PL_GRP_TAG>(parent);
}

///
/// 子グループ数を取得
///  @param[in]     tag     PolygonGroupを操作するためのタグ
///  @param[out]    num     ポリゴングループ数
///  @return 戻り値なし
///
void fpolylib_group_get_children_num_ (
                  PL_GRP_TAG* tag,
                  int* num
              ) 

{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(*tag);
    std::vector<PolygonGroup*>& pg_child = pg->get_children();

    *num = pg_child.size();
}


///
/// 子グループを取得
///  @param[in]     tag     PolygonGroupを操作するためのタグ
///  @param[out]    num     ポリゴングループ数
///  @param[out]    child_tags ポリゴングループ　タグ
///  @return 戻り値なし
///
void fpolylib_group_get_children_ (
                  PL_GRP_TAG* tag,
                  int* num,
                  PL_GRP_TAG* child_tags
              ) 

{
    PolygonGroup* pg = reinterpret_cast<PolygonGroup*>(*tag);
    std::vector<PolygonGroup*>& pg_child = pg->get_children();

    *num = pg_child.size();
    if( (*num) == 0 ) {
        return;
    }
    
    for(int i=0; i<(*num); i++ ) {
        child_tags[i] = reinterpret_cast<PL_GRP_TAG>( pg_child[i] );
    }
}

//**********************************************************
//  Triangleクラス用ラッパー関数
//      三角形オブジェクトと長田パッチオブジェクトでは動作が違う場合あり
//**********************************************************

/// ポリゴンタイプ取得
///  @param[in]     tag     Triangle/NptTriangleを操作するためのタグ
///  @return PL_TYPE_TRIANGLE / PL_TYPE_NPT（長田パッチ）
int fpolylib_triangle_get_pl_type_ (
            PL_ELM_TAG* tag
         )
{
    Triangle* tri = reinterpret_cast<Triangle*>(*tag);

    return tri->get_pl_type();
}

/// 頂点座標設定
void fpolylib_triangle_set_vertexes_ (
            PL_ELM_TAG* tag,
            PL_REAL vertex[9]
         )
{
    Triangle* tri = reinterpret_cast<Triangle*>(*tag);
    Vec3<PL_REAL> vertex_tmp[3];
    REAL9_TO_VEC3_3(vertex, vertex_tmp);

    tri->set_vertexes( vertex_tmp, true, true );
}

/// 頂点座標・長田パッチパラメータ設定
            // 法線ベクトル,面積も内部で設定する
void fpolylib_triangle_set_vertexes_npatch_ (
            PL_ELM_TAG* tag,
            PL_REAL vertex[9],
            NpatchParamStruct* param
         )
{
    NptTriangle* tri = reinterpret_cast<NptTriangle*>(*tag);
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


/// 頂点座標・長田パッチパラメータ設定
///   引数に構造体を使用しない版
            // 法線ベクトル,面積も内部で設定する
void fpolylib_triangle_set_vertexes_npatch2_ (
            PL_ELM_TAG* tag,
            PL_REAL vertex[9],
            PL_REAL cp_side1_1[3],
            PL_REAL cp_side1_2[3],
            PL_REAL cp_side2_1[3],
            PL_REAL cp_side2_2[3],
            PL_REAL cp_side3_1[3],
            PL_REAL cp_side3_2[3],
            PL_REAL cp_center[3]
         )
{
    NptTriangle* tri = reinterpret_cast<NptTriangle*>(*tag);
    Vec3<PL_REAL> vertex_tmp[3];
    NpatchParam param_tmp;

    REAL9_TO_VEC3_3(vertex, vertex_tmp);
    REAL_TO_VEC3(cp_side1_1, param_tmp.cp_side1_1 );
    REAL_TO_VEC3(cp_side1_2, param_tmp.cp_side1_2 );
    REAL_TO_VEC3(cp_side2_1, param_tmp.cp_side2_1 );
    REAL_TO_VEC3(cp_side2_2, param_tmp.cp_side2_2 );
    REAL_TO_VEC3(cp_side3_1, param_tmp.cp_side3_1 );
    REAL_TO_VEC3(cp_side3_2, param_tmp.cp_side3_2 );
    REAL_TO_VEC3(cp_center , param_tmp.cp_center  );
    
    tri->set_vertexes( vertex_tmp, param_tmp, true, true );
}


/// 頂点座標取得
void fpolylib_triangle_get_vertexes_ (
            PL_ELM_TAG* tag,
            PL_REAL vertex[9]
         )
{
    Triangle* tri = reinterpret_cast<Triangle*>(*tag);
    Vec3<PL_REAL>* vertex_tmp = tri->get_vertexes();

    VEC3_3_TO_REAL9( vertex_tmp, vertex );
}

/// 法線ベクトル取得
void fpolylib_triangle_get_normal_ (
            PL_ELM_TAG* tag,
            PL_REAL norm[3]
         )
{
    Triangle* tri = reinterpret_cast<Triangle*>(*tag);
    Vec3<PL_REAL> normal_tmp = tri->get_normal();

    VEC3_TO_REAL( normal_tmp, norm );
}

/// 長田パッチパラメータ設定
/// オブジェクトが長田パッチでない時はエラーを返す
void fpolylib_triangle_set_npatchParam_ (
            PL_ELM_TAG* tag,
            NpatchParamStruct *param,
            POLYLIB_STAT* ret
         )
{
    NptTriangle* tri = reinterpret_cast<NptTriangle*>(*tag);
    int pl_type = tri->get_pl_type();
    if( pl_type != PL_TYPE_NPT ) {
        *ret = PLSTAT_NG;
        return;
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
    *ret = PLSTAT_OK;
}

/// 長田パッチパラメータ設定
///   引数に構造体を使用しない版
/// オブジェクトが長田パッチでない時はエラーを返す
void fpolylib_triangle_set_npatchParam2_ (
            PL_ELM_TAG* tag,
            PL_REAL cp_side1_1[3],
            PL_REAL cp_side1_2[3],
            PL_REAL cp_side2_1[3],
            PL_REAL cp_side2_2[3],
            PL_REAL cp_side3_1[3],
            PL_REAL cp_side3_2[3],
            PL_REAL cp_center[3],
            POLYLIB_STAT* ret
         )
{
    NptTriangle* tri = reinterpret_cast<NptTriangle*>(*tag);
    int pl_type = tri->get_pl_type();
    if( pl_type != PL_TYPE_NPT ) {
        *ret = PLSTAT_NG;
        return;
    }

    NpatchParam param_tmp;

    REAL_TO_VEC3(cp_side1_1, param_tmp.cp_side1_1 );
    REAL_TO_VEC3(cp_side1_2, param_tmp.cp_side1_2 );
    REAL_TO_VEC3(cp_side2_1, param_tmp.cp_side2_1 );
    REAL_TO_VEC3(cp_side2_2, param_tmp.cp_side2_2 );
    REAL_TO_VEC3(cp_side3_1, param_tmp.cp_side3_1 );
    REAL_TO_VEC3(cp_side3_2, param_tmp.cp_side3_2 );
    REAL_TO_VEC3(cp_center , param_tmp.cp_center  );
    
    tri->set_npatch_param( param_tmp );
    *ret = PLSTAT_OK;
}

/// 長田パッチパラメータ取得
/// オブジェクトが長田パッチでない時はエラーを返す
void fpolylib_triangle_get_npatchParam_ (
            PL_ELM_TAG* tag,
            NpatchParamStruct *param,
            POLYLIB_STAT* ret
         )
{
    NptTriangle* tri = reinterpret_cast<NptTriangle*>(*tag);
    int pl_type = tri->get_pl_type();
    if( pl_type != PL_TYPE_NPT ) {
        *ret = PLSTAT_NG;
        return;
    }

    NpatchParam* param_tmp = tri->get_npatch_param();

    VEC3_TO_REAL( param_tmp->cp_side1_1, param->cp_side1_1 );
    VEC3_TO_REAL( param_tmp->cp_side1_2, param->cp_side1_2 );
    VEC3_TO_REAL( param_tmp->cp_side2_1, param->cp_side2_1 );
    VEC3_TO_REAL( param_tmp->cp_side2_2, param->cp_side2_2 );
    VEC3_TO_REAL( param_tmp->cp_side3_1, param->cp_side3_1 );
    VEC3_TO_REAL( param_tmp->cp_side3_2, param->cp_side3_2 );
    VEC3_TO_REAL( param_tmp->cp_center , param->cp_center  );

    *ret = PLSTAT_OK;
}

/// 長田パッチパラメータ取得
///   引数に構造体を使用しない版
/// オブジェクトが長田パッチでない時はエラーを返す
void fpolylib_triangle_get_npatchParam2_ (
            PL_ELM_TAG* tag,
            PL_REAL cp_side1_1[3],
            PL_REAL cp_side1_2[3],
            PL_REAL cp_side2_1[3],
            PL_REAL cp_side2_2[3],
            PL_REAL cp_side3_1[3],
            PL_REAL cp_side3_2[3],
            PL_REAL cp_center[3],
            POLYLIB_STAT* ret
         )
{
    NptTriangle* tri = reinterpret_cast<NptTriangle*>(*tag);
    int pl_type = tri->get_pl_type();
    if( pl_type != PL_TYPE_NPT ) {
        *ret = PLSTAT_NG;
        return;
    }

    NpatchParam* param_tmp = tri->get_npatch_param();

    VEC3_TO_REAL( param_tmp->cp_side1_1, cp_side1_1 );
    VEC3_TO_REAL( param_tmp->cp_side1_2, cp_side1_2 );
    VEC3_TO_REAL( param_tmp->cp_side2_1, cp_side2_1 );
    VEC3_TO_REAL( param_tmp->cp_side2_2, cp_side2_2 );
    VEC3_TO_REAL( param_tmp->cp_side3_1, cp_side3_1 );
    VEC3_TO_REAL( param_tmp->cp_side3_2, cp_side3_2 );
    VEC3_TO_REAL( param_tmp->cp_center , cp_center  );

    *ret = PLSTAT_OK;
}


/// ポリゴン（形状）のユーザ定義属性数（整数型）の取得
int fpolylib_triangle_get_num_atrI_ (
                PL_ELM_TAG* tag
            )
{
    Triangle* tri = reinterpret_cast<Triangle*>(*tag);
    return tri->get_num_atrI();
}

/// ポリゴン（形状）のユーザ定義属性数（整数型）の取得
int fpolylib_triangle_get_num_atrR_ (
                PL_ELM_TAG* tag
            )
{
    Triangle* tri = reinterpret_cast<Triangle*>(*tag);
    return tri->get_num_atrR();
}

/// ポリゴン（形状）のユーザ定義属性（整数型）の取得
int fpolylib_triangle_get_atrI_ (
                PL_ELM_TAG* tag,
                int*        atr_no  // 1〜
            )
{
    Triangle* tri = reinterpret_cast<Triangle*>(*tag);
    int* p = tri->get_pAtrI();
    return p[*atr_no-1]; 
}

/// ポリゴン（形状）のユーザ定義属性（整数型）の設定
void fpolylib_triangle_set_atrI_ (
                PL_ELM_TAG* tag,
                int*        atr_no,  // 1〜
                int*        val
            )
{
    Triangle* tri = reinterpret_cast<Triangle*>(*tag);
    int* p = tri->get_pAtrI();
    p[*atr_no-1] = *val; 
}


/// ポリゴン（形状）のユーザ定義属性（整数型）の取得
PL_REAL fpolylib_triangle_get_atrR_ (
                PL_ELM_TAG* tag,
                int*        atr_no  // 1〜
            )
{
    Triangle* tri = reinterpret_cast<Triangle*>(*tag);
    PL_REAL* p = tri->get_pAtrR();
    return p[*atr_no-1]; 
}

/// ポリゴン（形状）のユーザ定義属性（整数型）の設定
void fpolylib_triangle_set_atrR_ (
                PL_ELM_TAG* tag,
                int*        atr_no,  // 1〜
                PL_REAL*        val
            )
{
    Triangle* tri = reinterpret_cast<Triangle*>(*tag);
    PL_REAL* p = tri->get_pAtrR();
    p[*atr_no-1] = *val; 
}



