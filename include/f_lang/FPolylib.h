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

///
/// Fortranインターフェース実装用のヘッダ
///    Fortranアプリ側から使用するヘッダではありません
///

#ifndef f_polylib_h
#define f_polylib_h

#ifdef MPI_PL
#include "mpi.h"
#endif

#define POLYLIB_FALSE 0
#define POLYLIB_TRUE  1

/// ポリゴングループ名のfortran文字列長
#define PL_GRP_PATH_LEN 256
#define PL_GRP_NAME_LEN 64

/// 属性用文字列のfortran文字列長
///     属性キー, 属性値 etc.
#define PL_GRP_ATR_LEN  32


/// ファイルパス、ファイル名のfortran文字列長
#define PL_FILE_PATH_LEN    256
#define PL_FILE_NAME_LEN    64

/// ファイルフォーマットのfortran文字列長
#define PL_FORMAT_LEN       8

/// 文字列のfortran文字列長
///    その他
#define PL_STR_LEN          32


#ifdef __cplusplus
extern "C" {  // for C++
#else
#endif

#include <stdbool.h>

#include "c_lang/CPolylib.h"

///
/// 領域情報構造体(Fortran用)
///
typedef struct {
    PL_REAL bpos[3];
    int bbsize[3];   // fortranにはunsigned intがないためintとする
    int gcsize[3];   // fortranにはunsigned intがないためintとする
    PL_REAL dx[3];
} FParallelBboxStruct;
//} fparallelbboxstruct;


///
/// Fortran用Polylib
///
/// 注意：Fortran用インターフェースのためのC言語インターフェースを提供する
///       Fortran言語APIではタグを操作（ハンドル）用として使用します。
///       タグには有効期間があります。セッション中で永続的に有効というわけではありません。
///         PolygonGroupタグ
///             PolygonGroupの追加、削除、挿入があると該当する箇所以外も含めて全て無効となります。
///         Triangleタグ
///             Triangleの追加、削除、挿入、ソートがあると該当する箇所以外も含めて全て無効となります。
///             並列環境でTriangleを移動した場合、migrate処理を実行しますが
///             この際に削除、挿入、ソートが行われるのでタグが全て無効となります。
///
///         Fortran用インターフェース基本方針
///             サブルーチン名・関数名：C言語APIの名前に先頭に'f' 末尾に'_'を付加する
///             引数は全て、アドレス渡しとする
///             戻り値は基本的に引数に変更。
///             unsigned int は　すべて int に変更
///             C言語の char*型(NULL termination) はFortranではNULL文字以降を
///             すべてspaceとする
///             当ラッパー関数内ではCの構造体を使用しているが
///             Fortran側では、Fortranの構造体のインクルードファイルを使用すること
///
///         Fortran用インターフェース制限事項
///             (1) FortranからPolylib環境を構築する場合は
///                 MPIのコミュニケータはMPI_COMM_WORLDとなります。
///                 他のコミュニケータが必要な場合は、Polylib環境をC++/Cで構築する
///                 必要があります。
///             (2) Fortranからは移動関数の登録は出来ません。
///                 C++/Cより移動関数を登録してください。


///
/// F言語用Polylib環境の構築
///     Polylibインスタンス生成
/// Polylib::get_instanceメソッドの代替関数。
///
///  @param[out]  ret   POLYLIB_STATで定義される値が返る。
///  @attention 最初に呼び出すこと
///
void  fpolylib_instance_ (
                              POLYLIB_STAT* ret 
                    );



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
///  @attention  C++/Cと違いコミュニケーターの指定は出来ません。
///      固定でMPI_COMM_WORLDとなります。
///      その他のコミュニケータを使う場合は、C++/Cのメソッドを使用してください。
///
void
fpolylib_init_parallel_info_ (
                              PL_REAL bpos[3],
                              int bbsize[3],
                              int gcsize[3],
                              PL_REAL dx[3],
                              POLYLIB_STAT* ret
);

///
/// 並列計算関連情報の設定と初期化を行う。
/// (各ランクが複数領域を担当している場合）
/// 全rankで各々設定を行い、その領域情報を全rankへ配信する。
///
///  @param[in]  bpos    自PE担当領域の基点座標
///  @param[in]  num     自PEが担当する領域
///  @param[in]  bbox    担当するboundary box情報（複数）
///                         Fortranの構造体で受けること
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///  @attention  C++/Cと違いコミュニケーターの指定は出来ません。
///      固定でMPI_COMM_WORLDとなります。
///      その他のコミュニケータを使う場合は、C++/Cのメソッドを使用してください。
///      構造体に対応していないFortranでは使用出来ません。C++/Cのメソッドを使用してください。
///
void
fpolylib_init_parallel_info2_ (
        int      *num,
        FParallelBboxStruct *bbox,
        POLYLIB_STAT* ret
    );
#endif

///
///  Polylib::loadメソッドのラッパー関数。
///  引数で指定された設定ファイルを読み込み、グループツリーを作成する。
///  続いて設定ファイルで指定されたSTLファイルを読み込み、KD木を作成する。
///  @param[in]  fname 設定ファイル名
///                     Fortran型の文字列  (\0で終了しない）
///                     character(len=PL_FILE_PATH_LEN), character*256
///                     Fortran側から呼ぶときは、文字列の最初の空白の１個前までを
///                     有効なファイル名と見なします。
///                     すべて空白の時は"polylib_config.tp"とみなす
///  @param[in]  scale 縮尺率
///                      PL_REALが倍精度の時に 単精度定数1.0等を指定するのはNG
///                      1.0d0等で倍精度定数を指定すること
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_load_ (
                   char* config_name,
                   PL_REAL* scale,
                   POLYLIB_STAT* ret
                );

///
/// Polylib::saveメソッドのラッパー関数。
/// PolygoGroupツリー、三角形ポリゴン情報の保存。
/// グループツリーの情報を設定ファイルへ出力。三角形ポリゴン情報をSTL
/// ファイル or NPTファイルへ出力
///
///  @param[out]    o_fname 設定ファイル名
///                             Fortran型の文字列(長さPL_FILE_PATH_LEN)  (\0で終了しない）
///                             character(len=PL_FILE_PATH_LEN), character*256
///  @param[in]     format  形状ファイルのフォーマット
///                             Fortran型の文字列(長さPL_FORMAT_LEN)  (\0で終了しない）
///                             character(len=PL_FORMAT_LEN), character*8
///                             FPolylib_define.incで定義されているFILE_FMT_* 参照
///                               'stl_a   '
///                               'stl_aa  '
///                               'stl_b   '
///                               'stl_bb  '
///                               'npt_a   '
///                               'npt_b   '
///  @param[in]     extend  ファイル名に付加する文字列
///                             Fortran型の文字列(長さPL_STR_LEN)
///                             character(len=PL_STR_LEN), character*32
//                              空白を指定した場合は、付加文字列として本メソッド呼び出し時の
///                             年月日時分秒(YYYYMMDD24hhmmss)を用いる。
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///  @attention ファイル名命名規約は次の通り。
///         設定ファイル : polylib_config_付加文字.tpp
///         STL or NPTファイル  : ポリゴングループ名_付加文字.拡張子
///
void fpolylib_save_ (
                          char* o_fname,
                          char* format,
                          char* extend,
                          POLYLIB_STAT* ret
                    );

///
/// 三角形ポリゴン座標の移動。
/// Polylib::moveメソッドのラッパー関数
/// 本クラスインスタンス配下の全PolygonGroupのmoveメソッドが呼び出される。
/// move関数は、Fortran関数では登録できない。CまたｈC++から行う。
///    C言語：polygongroup_set_move_func_c関数
///    C++  ：PolygonGroup::pg->set_move_func関数
///
///  @param[in]     param   移動計算パラメータセット
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_move_ (
                          PolylibMoveParamsStruct* param,
                          POLYLIB_STAT* ret
                    );


#ifdef MPI_PL
///
/// Polylib::mograteメソッドのラッパー関数
/// オブジェクトのインスタンス毎に登録が必要
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_migrate_ ( POLYLIB_STAT* ret );
#endif

///
/// ポリゴングループタグの取得
///  @param[in]     group_name  ポリゴングループ名
///                                character*256
///  @param[out]    tag         ポリゴングループ
///                                integer*8
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_get_group_tag_ (
        char*       group_name,
        PL_GRP_TAG* grp_tag,
        POLYLIB_STAT* ret
    );

///
/// ルートのポリゴングループの個数取得
///
///  @param[out] n       ポリゴングループ数
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_get_root_groups_tags_num_ (
                    int* n,
                    POLYLIB_STAT* ret
                  );

///
/// PolygonGroupツリーの最上位ノードの取得（Cインターフェース用）
/// Polylib::get_root_groups_tagsメソッドのラッパー関数
///
///  @param[out] n       ポリゴングループ数
///  @param[out] tags    ポリゴングループ　タグ
///                            integer*8
///       get_root_groups_tags_num_()で得られる個数以上をallocate済みであること
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_get_root_groups_tags_ (
                    int* n, 
                    PL_GRP_TAG* tags,
                    POLYLIB_STAT* ret
                  );

///
/// リーフのポリゴングループの個数取得
///
///  @param[out] n       ポリゴングループ数
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_get_leaf_groups_tags_num_ (
                    int* n,
                    POLYLIB_STAT* ret
                  );

///
/// リーフPolygonGroupリストの取得
///
///  @param[out] n       ポリゴングループ数
///  @param[out] tags    ポリゴングループ　タグ
///                          integer*8
///         get_leaf_groups_tags_num_()で得られる個数以上をallocate済みであること
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_get_leaf_groups_tags_ (
                    int* n, 
                    PL_GRP_TAG* tags,
                    POLYLIB_STAT* ret
                  );

///
/// 位置ベクトルmin_posとmax_posにより特定される矩形領域に含まれる、
/// ポリゴン数取得
///  @param[out]    num         抽出されたポリゴン数
///  @param[in]     group_name  抽出グループ名
///                             Fortran型の文字列(長さPL_GRP_PATH_LEN)
///                                character(len=PL_GRP_PATH_LEN), character*256
///  @param[in]     min_pos     抽出する矩形領域の最小値。(x,y,z順の配列)
///  @param[in]     max_pos     抽出する矩形領域の最大値。(x,y,z順の配列)
///  @param[in]     every       抽出オプション。
///   1：3頂点が全て検索領域に含まれるポリゴンを抽出する。
///   0：三角形のBBoxが一部でも検索領域と交差するものを抽出する。
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_search_polygons_num_ (
        int         *num,
        char*       group_name,
        PL_REAL     min_pos[3],
        PL_REAL     max_pos[3],
        int         *every,
        POLYLIB_STAT* ret
    );

///
/// Polylib::search_polygonsメソッドのラッパー関数。
/// 位置ベクトルmin_posとmax_posにより特定される矩形領域に含まれる、
/// 特定のグループとその子孫グループに属する三角形ポリゴンをKD探索に
/// より抽出する。
///  @param[out]    num         抽出された三角形ポリゴン数
///  @param[out]    tags        ポリゴンのタグ（ハンドル)
///                                integer*8
///       polylib_search_polygons_num_()で得られる個数以上をallocate済みであること
///  @param[in]     group_name  抽出グループ名
///                             Fortran型の文字列(長さPL_GRP_PATH_LEN)
///                                 character(len=PL_GRP_PATH_LEN)
///  @param[in]     min_pos     抽出する矩形領域の最小値。(x,y,z順の配列)
///  @param[in]     max_pos     抽出する矩形領域の最大値。(x,y,z順の配列)
///  @param[in]     every       抽出オプション。
///   1：3頂点が全て検索領域に含まれるポリゴンを抽出する。
///   0：三角形のBBoxが一部でも検索領域と交差するものを抽出する。
///  @param[in]     max         三角形ポリゴンMAX数（領域確保数）

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
    int         *every,
    POLYLIB_STAT* ret
    );

///
/// 指定した点に最も近いポリゴンの検索
///  @param[out] tag         ポリゴンのタグ（ハンドル)
///                                integer*8
///                                 ポリゴンがない場合、
///                                 PL_NULL_TAGi(=0)が返される
///  @param[in]  group_name  グループ名
///                             Fortran型の文字列(長さPL_GRP_PATH_LEN)
///                             character(len=PL_GRP_PATH_LEN), character*256
///  @param[in]  pos         指定点
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///  @attention MPI並列計算時は,posは各ランクの矩形領域を
///             超えないようにして下さい (各ランク内の担当領域内のみ検索するため）
///
void fpolylib_search_nearest_polygon_ (
        PL_ELM_TAG  *tag,
        char*       group_name,
        PL_REAL     pos[3],
        POLYLIB_STAT* ret
    );

///
/// Polylib::show_group_hierarchyメソッドのラッパー関数。
/// グループ階層構造リストを標準出力に出力する。
///
void fpolylib_show_group_hierarchy_ ();

///
/// Polylib::show_group_infoメソッドのラッパー関数。
/// グループの情報を出力する。(親グループ名、自信の名前、ファイル名、
///   登録三角形数、3頂点ベクトルの座標、法線ベクトルの座標、面積)
///  @param[in] group_name グループ名
///                             Fortran型の文字列(長さPL_GRP_PATH_LEN)
///                             character(len=PL_GRP_PATH_LEN)
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_show_group_info_ (char* group_name, POLYLIB_STAT* ret );

///
///  Polylibが利用中の概算メモリ量を返す
///
/// @return 利用中のメモリ量(byte)
///            2GB以上返せないので注意
///            2GB以上が予想される場合はfpolylib_used_memory_size_mb()を使用すること
///
int fpolylib_used_memory_size_ ();

///
///  Polylibが利用中の概算メモリ量(MB)を返す
///
/// @return 利用中のメモリ量(Mbyte)
///
int fpolylib_used_memory_size_mb_ ();


//**********************************************************
//  PolygonGroupクラス用ラッパー関数
//**********************************************************


/// PolygonGroupのポリゴンを求める
///    PolygonGroup::get_trianglesのラッパー関数
///  @param[in]   tag_pg       PolygonGroupを操作するためのタグ
///                                integer*8
///  @param[out]  num_tri      三角形ポリゴン数
///  @param[out]  tags_tri     三角形ポリゴンのタグ（ハンドル)
///                                integer*8
///                                tags_triはpolylib_group_get_triangles_num()で
///                                求めた個数分allocateされていること
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_group_get_triangles_ (
            PL_GRP_TAG *tag_pg,
            int        *num_tri,
            PL_ELM_TAG  *tags_tri,
            POLYLIB_STAT* ret
         );

/// PolygonGroupのポリゴン数を求める
///    PolygonGroup::get_trianglesのラッパー関数
///  @param[in]   tag_pg       PolygonGroupを操作するためのタグ
///                                integer*8
///  @param[out]  num_tri      三角形ポリゴン数
///  @param[out] ret     POLYLIB_STAT(integer)で定義される値が返る
///
void fpolylib_group_get_num_triangles_ (
            PL_GRP_TAG *tag_pg,
            int        *num_tri,
            POLYLIB_STAT* ret
         );

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
         );
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
         );

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
         );
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
    );

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
    );

///
/// KD木の再構築フラグの設定
///     ユーザ定義の移動関数内の最後で呼び出す
/// 
///  @param[in] tag_pg       PolygonGroupを操作するためのタグ
///  @return    戻り値なし
/// 
void fpolylib_group_set_need_rebuild_ ( PL_GRP_TAG* tag_pg );


//----------------------------------------------
// setter / getter
//----------------------------------------------
void fpolylib_group_set_name_ ( 
            PL_GRP_TAG *tag,    // PolygonGroupを操作するためのタグ integer*8 
            char* name          // Fortran型文字列(長さPL_GRP_PATH_LEN)
         );

void fpolylib_group_get_name_ ( 
            PL_GRP_TAG *tag,    // PolygonGroupを操作するためのタグ integer*8
            char* name          // Fortran型文字列(長さPL_GRP_PATH_LEN)
         );


/// 移動対象フラグ設定
void fpolylib_group_set_movable_ ( 
            PL_GRP_TAG* tag,     // PolygonGroupを操作するためのタグ integer*8
            int*        movable
         );

/// 移動対象フラグ取得
void fpolylib_group_get_movable_ ( 
            PL_GRP_TAG* tag,     // PolygonGroupを操作するためのタグ integer*8
            int*        movable
         );

/// ユーザ定義属性の設定
void fpolylib_group_set_atr_ ( 
            PL_GRP_TAG* tag,    // PolygonGroupを操作するためのタグ integer*8
            char* key,          // Fortran型文字列(長さPL_GRP_ATR_LEN)
                                //    character(len=PL_GRP_ATR_LEN), character*32
            char* str           // Fortran型文字列(長さPL_GRP_ATR_LEN)
                                //    character(len=PL_GRP_ATR_LEN), character*32
         );

/// ユーザ定義属性の取得
void fpolylib_group_get_atr_ ( 
            PL_GRP_TAG* tag,    // PolygonGroupを操作するためのタグ integer*8
            char* key,          // Fortran型文字列(長さPL_GRP_ATR_LEN)
                                //    character(len=PL_GRP_ATR_LEN), character*32
            char* str,          // Fortran型文字列(長さPL_GRP_ATR_LEN)
                                //    character(len=PL_GRP_ATR_LEN), character*32
            POLYLIB_STAT* ret   // POLYLIB_STAT(integer)で定義される値が返る
         );

/// ポリゴングループ内のポリゴンのユーザ定義属性数の設定
///  @param[in]     tag      PolygonGroupを操作するためのタグ integer*8
///                            integer*8
///  @param[in]     num_atrI 整数属性数
///  @param[in]     num_atrR 実数属性数
///  @return 戻り値なし
void fpolylib_group_set_num_polygon_atr_ (
            PL_GRP_TAG *tag,
            int *num_atrI, 
            int *num_atrR
            );

///
/// PolygonGroup 親グループを取得
/// 
///  @param[in]  tag         PolygonGroupを操作するためのタグ
///                            integer*8
///  @param[out] parent_tag  親グループのタグ  integer*8
/// 
void fpolylib_group_get_parent_ (
              PL_GRP_TAG* tag,
              PL_GRP_TAG* parent_tag
          );

///
/// 子グループ数を取得
///  @param[in]     tag     PolygonGroupを操作するためのタグ
///                            integer*8
///  @param[out]    num     ポリゴングループ数
///  @return 戻り値なし
///
void fpolylib_group_get_children_num_ (
                  PL_GRP_TAG* tag,
                  int* num
              );

///
/// 子グループを取得
///  @param[in]     tag     PolygonGroupを操作するためのタグ
///                            integer*8
///  @param[out]    num     ポリゴングループ数
///  @param[out]    child_tags ポリゴングループ　タグ
///                            integer*8
///  @return 戻り値なし
///
void fpolylib_group_get_children_ (
                  PL_GRP_TAG* tag,
                  int* num,
                  PL_GRP_TAG* child_tags
              );



//**********************************************************
//  Triangleクラス用ラッパー関数
//      三角形オブジェクトと長田パッチオブジェクトでは動作が違う場合あり
//**********************************************************

/// ポリゴンタイプ取得
///  @param[in]     tag     Triangle/NptTriangleを操作するためのタグ
///                            integer*8
///  @return PL_TYPE_TRIANGLE / PL_TYPE_NPT（長田パッチ）
int fpolylib_triangle_get_pl_type_ (
            PL_ELM_TAG* tag
         );

//----------------------------------------------
// setter / getter
//----------------------------------------------

/// 頂点座標設定
/// 頂点設定時、法線ベクトル,面積も内部で設定する
/// 長田パッチのパラメータは更新されないので注意
///
///  @param[in]     tag     Triangle/NptTriangleを操作するためのタグ
///                            integer*8
///  @param[in]     vertex  ３頂点の座標
void fpolylib_triangle_set_vertexes_ ( 
            PL_ELM_TAG *tag, 
            PL_REAL vertex[9]
         );

/// 頂点座標・長田パッチパラメータ設定
///  @param[in]     tag     Triangle/NptTriangleを操作するためのタグ
///                            integer*8
///  @param[in]     vertex  ３頂点の座標
///  @param[in]     param   長田パッチパラメータ
///  @return 戻り値なし
///  @attention 法線ベクトル,面積も内部で設定する
void fpolylib_triangle_set_vertexes_npatch_ (
            PL_ELM_TAG* tag,
            PL_REAL vertex[9],
            NpatchParamStruct* param
         );

/// 頂点座標・長田パッチパラメータ設定
///   引数に構造体を使用しない版
///  @param[in]     tag         Triangle/NptTriangleを操作するためのタグ
///                               integer*8
///  @param[in]     vertex      ３頂点の座標
///  @param[in]     cp_side1_1  長田パッチパラメータ p1p2辺の３次ベジェ制御点1 
///  @param[in]     cp_side1_2  長田パッチパラメータ p1p2辺の３次ベジェ制御点2 
///  @param[in]     cp_side2_1  長田パッチパラメータ p2p3辺の３次ベジェ制御点1 
///  @param[in]     cp_side2_2  長田パッチパラメータ p2p3辺の３次ベジェ制御点2 
///  @param[in]     cp_side3_1  長田パッチパラメータ p3p1辺の３次ベジェ制御点1 
///  @param[in]     cp_side3_2  長田パッチパラメータ p3p1辺の３次ベジェ制御点2 
///  @param[in]     cp_center   長田パッチパラメータ 三角形中央の３次ベジェ制御点
///  @return 戻り値なし
///  @attention 法線ベクトル,面積も内部で設定する
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
         );


/// 頂点座標取得
///  @param[in]     tag     Triangle/NptTriangleを操作するためのタグ
///                            integer*8
///  @param[out]    vertex  ３頂点の座標
///  @return 戻り値なし
void fpolylib_triangle_get_vertexes_ ( 
            PL_ELM_TAG *tag, 
            PL_REAL vertex[9]
         );

/// 法線ベクトル取得
///  @param[in]     tag     Triangle/NptTriangleを操作するためのタグ
///                            integer*8
///  @param[out]    norm    法線ベクトル
///  @return 戻り値なし
void fpolylib_triangle_get_normal_ ( 
            PL_ELM_TAG *tag, 
            PL_REAL norm[3]
         );

/// 長田パッチパラメータ設定
///  @param[in]     tag     NptTriangleを操作するためのタグ
///                            integer*8
///  @param[in]     param   長田パッチパラメータ
///  @param[out]    ret     POLYLIB_STAT(integer)で定義される値が返る
///      オブジェクトが長田パッチでない時はエラーを返す
void fpolylib_triangle_set_npatchParam_ ( 
            PL_ELM_TAG *tag, 
            NpatchParamStruct *param,
            POLYLIB_STAT* ret
         );

/// 長田パッチパラメータ設定
///   引数に構造体を使用しない版
///  @param[in]     tag         NptTriangleを操作するためのタグ
///                               integer*8
///  @param[in]     cp_side1_1  長田パッチパラメータ p1p2辺の３次ベジェ制御点1 
///  @param[in]     cp_side1_2  長田パッチパラメータ p1p2辺の３次ベジェ制御点2 
///  @param[in]     cp_side2_1  長田パッチパラメータ p2p3辺の３次ベジェ制御点1 
///  @param[in]     cp_side2_2  長田パッチパラメータ p2p3辺の３次ベジェ制御点2 
///  @param[in]     cp_side3_1  長田パッチパラメータ p3p1辺の３次ベジェ制御点1 
///  @param[in]     cp_side3_2  長田パッチパラメータ p3p1辺の３次ベジェ制御点2 
///  @param[in]     cp_center   長田パッチパラメータ 三角形中央の３次ベジェ制御点
///  @param[out]    ret     POLYLIB_STAT(integer)で定義される値が返る
///      オブジェクトが長田パッチでない時はエラーを返す
void fpolylib_triangle_set_npatchParam2_ ( 
            PL_ELM_TAG *tag, 
            PL_REAL cp_side1_1[3],    
            PL_REAL cp_side1_2[3],    
            PL_REAL cp_side2_1[3],    
            PL_REAL cp_side2_2[3],    
            PL_REAL cp_side3_1[3],    
            PL_REAL cp_side3_2[3],    
            PL_REAL cp_center[3],
            POLYLIB_STAT* ret
         );

/// 長田パッチパラメータ取得
///  @param[in]     tag     NptTriangleを操作するためのタグ
///                            integer*8
///  @param[out]    param   長田パッチパラメータ
///  @param[out]    ret     POLYLIB_STAT(integer)で定義される値が返る
///      オブジェクトが長田パッチでない時はエラーを返す
void fpolylib_triangle_get_npatchParam_ ( 
            PL_ELM_TAG *tag, 
            NpatchParamStruct *param,
            POLYLIB_STAT* ret
         );

/// 長田パッチパラメータ取得
///   引数に構造体を使用しない版
///  @param[in]     tag         NptTriangleを操作するためのタグ
///                               integer*8
///  @param[out]    cp_side1_1  長田パッチパラメータ p1p2辺の３次ベジェ制御点1 
///  @param[out]    cp_side1_2  長田パッチパラメータ p1p2辺の３次ベジェ制御点2 
///  @param[out]    cp_side2_1  長田パッチパラメータ p2p3辺の３次ベジェ制御点1 
///  @param[out]    cp_side2_2  長田パッチパラメータ p2p3辺の３次ベジェ制御点2 
///  @param[out]    cp_side3_1  長田パッチパラメータ p3p1辺の３次ベジェ制御点1 
///  @param[out]    cp_side3_2  長田パッチパラメータ p3p1辺の３次ベジェ制御点2 
///  @param[out]    cp_center   長田パッチパラメータ 三角形中央の３次ベジェ制御点
///  @param[out]    ret     POLYLIB_STAT(integer)で定義される値が返る
///      オブジェクトが長田パッチでない時はエラーを返す
void fpolylib_triangle_get_npatchParam2_ ( 
            PL_ELM_TAG *tag, 
            PL_REAL cp_side1_1[3],    
            PL_REAL cp_side1_2[3],    
            PL_REAL cp_side2_1[3],    
            PL_REAL cp_side2_2[3],    
            PL_REAL cp_side3_1[3],    
            PL_REAL cp_side3_2[3],    
            PL_REAL cp_center[3],
            POLYLIB_STAT* ret
         );

/// ユーザ定義属性数性（整数型）の取得
///  @param[in]     tag      Triangle/NptTriangleを操作するためのタグ
///                            integer*8
///  @return 整数属性数
int fpolylib_triangle_get_num_atrI_ (
                PL_ELM_TAG *tag
            );

/// ユーザ定義属性数性（実数型）の取得
///  @param[in]     tag      Triangle/NptTriangleを操作するためのタグ
///                            integer*8
///  @return 実数属性数
int fpolylib_triangle_get_num_atrR_ (
                PL_ELM_TAG *tag
            );

/// ポリゴン（形状）のユーザ定義属性（整数型）の取得
///  @param[in]     tag      Triangle/NptTriangleを操作するためのタグ
///                            integer*8
///  @param[in]     atr_no   整数属性の何番目か　開始:1
///  @return ユーザ定義属性（整数型）値
int fpolylib_triangle_get_atrI_ (
            PL_ELM_TAG* tag,
            int*        atr_no 
            );

/// ポリゴン（形状）のユーザ定義属性（整数型）の設定
///  @param[in]     tag      Triangle/NptTriangleを操作するためのタグ
///                            integer*8
///  @param[in]     atr_no   整数属性の何番目か　開始:1
///  @param[in]     val      整数属性値
///  @return 戻り値なし
void fpolylib_triangle_set_atrI_ (
            PL_ELM_TAG* tag,
            int*        atr_no,
            int*        val
            );

/// ポリゴン（形状）のユーザ定義属性（実数型）の取得
///  @param[in]     tag      Triangle/NptTriangleを操作するためのタグ
///                            integer*8
///  @param[in]     atr_no   実数属性の何番目か　開始:1
///  @return ユーザ定義属性（実数型）値
PL_REAL fpolylib_triangle_get_atrR_ (
            PL_ELM_TAG* tag,
            int*        atr_no
            );

/// ポリゴン（形状）のユーザ定義属性（実数型）の設定
///  @param[in]     tag      Triangle/NptTriangleを操作するためのタグ
///                            integer*8
///  @param[in]     atr_no   実数属性の何番目か　開始:1
///  @param[in]     val      実数属性値
///  @return 戻り値なし
void fpolylib_triangle_set_atrR_ (
            PL_ELM_TAG* tag,
            int*        atr_no,
            PL_REAL*    val
            );

#ifdef __cplusplus
} // extern "C" or extern
#else
#endif

#endif //f_polylib_h
