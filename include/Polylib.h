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

#ifndef polylib_h
#define polylib_h
#include <string>
#include <vector>
#include <iostream>
#include "polygons/Triangle.h"
#include "polygons/NptTriangle.h"
#include "groups/PolygonGroup.h"
#include "file_io/PolygonIO.h"
#include "file_io/FileIO_func.h"
#include "common/PolylibStat.h"
#include "common/PolylibCommon.h"
#include "common/BBox.h"
#include "common/Vec3.h"
#include "Polylib_func.h"

#include "TextParser.h"
#include "polyVersion.h"

//#define TIME_MEASURE
#ifdef TIME_MEASURE
#include <stdio.h>
#include <time.h>
#endif // TIME_MEASURE

#ifdef WIN32
#include <stdio.h>
#include <time.h>
#endif

#ifdef MPI_PL
#include "mpi.h"
#endif

using namespace Vec3class;

namespace PolylibNS {

#ifdef MPI_PL

////////////////////////////////////////////////////////////////////////////
///
/// クラス:ParallelAreaInfo
/// 並列プロセス領域情報。
///
////////////////////////////////////////////////////////////////////////////

//
// Fortranではunsignedは付かないが、Fortran側で対応する
//   Fortran用に別途作成するという手もある
struct ParallelBbox {
    PL_REAL bpos[3];         // 基点座標
    unsigned int bbsize[3];  // 計算領域のボクセル数
    unsigned int gcsize[3];  // ガイドセルのボクセル数
    PL_REAL dx[3];           // ボクセル１辺の長さ
};

// 内部処理用(各ランク毎）

/// 計算領域情報。
struct CalcAreaInfo {
    /// 基点座標
    Vec3<PL_REAL> m_bpos;

    /// 計算領域のボクセル数
    Vec3<PL_REAL> m_bbsize;

    /// ガイドセルのボクセル数
    Vec3<PL_REAL> m_gcsize;

    /// ボクセル１辺の長さ
    Vec3<PL_REAL> m_dx;

    /// ガイドセルを含めた担当領域の最小位置
    Vec3<PL_REAL> m_gcell_min;

    /// ガイドセルを含めた担当領域の最大位置
    Vec3<PL_REAL> m_gcell_max;

    /// ガイドセルを含めたBounding Box
    ///     m_gcell_min,m_gcell_maxをBBox化したのみ
    BBox m_gcell_bbox;

};

/// 並列プロセス領域情報。
//struct ParallelInfo {
struct ParallelAreaInfo {

    /// ランクNo
    int m_rank;

    /// 各ランク内の計算領域情報
    std::vector<CalcAreaInfo> m_areas;

};



#endif

////////////////////////////////////////////////////////////////////////////
///
/// クラス:PolylibMoveParams
/// Polylib::move()の引数として利用するパラメタセットクラスです。
///
////////////////////////////////////////////////////////////////////////////
class PolylibMoveParams {
public:
    /// 現在の計算ステップ番号
    int m_current_step;

    /// 移動後の計算ステップ番号
    int m_next_step;

    /// １計算ステップあたりの時間変異
    PL_REAL m_delta_t;
    
    /// ユーザ定義パラメータ
    //    m_current_step, m_next_step, m_delta_t以外のパラメータが必要な時に任意に使用する
    PL_REAL m_params[10];
};

////////////////////////////////////////////////////////////////////////////
///
/// クラス:Polylib
/// ポリゴンを管理する為のクラスライブラリです。
///
////////////////////////////////////////////////////////////////////////////

class Polylib
{
public:
    ///
    /// singletonのPolylibインスタンス取得。
    /// デフォルトのFactoryクラスであるPolygonGroupFactoryを使用してインスタンス
    /// を生成する。
    ///
    ///  @return    Polylibクラスのインスタンス。
    ///  @attention 呼び出し側でdeleteはできません。
    ///
    static Polylib* get_instance();

    ///
    /// 検索モードの指定
    /// ポリゴン検索時にポリゴンを曲面補正して検索させるかどうかを指定する
    ///  @param[in] detail   詳細検索指定
    ///                          true  : 曲面補正が可能な場合、曲面補正して検索
    ///                          false : 3角形平面で検索
    static void set_srch_mode( bool detail );

    ///
    /// 検索モードの取得
    /// ポリゴン検索時にポリゴンを曲面補正して検索させるかどうかを指定する
    ///
    ///  @return   検索モード
    ///                          true  : 曲面補正が可能な場合、曲面補正して検索
    ///                          false : 3角形平面で検索
    static bool get_srch_mode();

#ifdef MPI_PL
    ///
    /// MPIのコミュニケータを返す
    ///
    ///  @return MPIのコミュニケータ
    ///      PolygonGroup内のポリゴンの集合演算を行うときに使用する
    MPI_Comm get_MPI_Comm( void )
    {
        return m_comm;
    }

    ///
    /// MPIのランクＮｏを返す
    ///
    ///  @return MPIのランクＮｏ
    ///      PolygonGroup内のscatter/gatherで使用する
    int get_MPI_myrank( void )
    {
        return m_myrank;
    }

    ///
    /// MPIの並列プロセス数を返す
    ///
    ///  @return MPIの並列プロセス数
    ///      PolygonGroup内のscatter/gatherで使用する
    int get_MPI_numproc( void )
    {
        return m_numproc;
    }

    ///
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
        init_parallel_info(
        MPI_Comm comm,
        PL_REAL bpos[3],
        unsigned int bbsize[3],
        unsigned int gcsize[3],
        PL_REAL dx[3]
    );

    ///
    /// 並列計算関連情報の設定と初期化を行う。
    /// (各ランクが複数領域を担当している場合）
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
        init_parallel_info(
        MPI_Comm comm,
        const std::vector<ParallelBbox>&  bboxes
    );

    ///
    /// 自PE担当領域情報取得
    /// @return 自PE領域情報
    ///
    ParallelAreaInfo* get_myproc_area()
    {
        return &m_myproc_area;
    }

    ///
    /// 自PEを除く全PE担当領域情報リスト取得
    /// @return 自PEを除く全PE担当領域情報
    ///
    std::vector<ParallelAreaInfo>* get_other_procs_area()
    {
        return &m_other_procs_area;
    }

#endif

    ///
    /// PolygoGroup、三角形ポリゴン情報の読み込み。
    /// 引数で指定された設定ファイル (TextParser 形式) を読み込み、グループツリーを作成する。
    /// 続いて設定ファイルで指定されたSTL(or NPT)ファイルを読み込み、KD木を作成する。
    ///
    ///  @param[in] config_name 設定ファイル名
    ///  @param[in] scale       縮尺率
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///  @attention set_max_memory_size_mb関数にてメモリ使用サイズMAXを
    ///             指定すると、必要に応じてポリゴンファイルを分割して
    //              ロードします。
    ///
    POLYLIB_STAT load(
        const std::string   config_name = "polylib_config.tp",
        PL_REAL             scale = 1.0
        );

    ///
    /// PolygoGroupツリー、三角形ポリゴン情報の保存。
    /// グループツリーの情報を設定ファイルへ出力。三角形ポリゴン情報をSTL/NPT
    /// ファイルへ出力。
    ///
    ///  @param[out] config_name_out  保存した設定ファイル名の返却用。
    ///  @param[in]  file_format      PolygonIO.hクラスで定義されているSTL/NPTファイルの
    ///                               フォーマット
    ///  @param[in]  extend           ファイル名に付加する文字列。省略可。省略した
    ///                               場合は、付加文字列として本メソッド呼び出し時
    ///                               の年月日時分秒(YYYYMMDD24hhmmss)を用いる。
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///  @attention ファイル名命名規約は次の通り。
    ///         設定ファイル : polylib_config_ランク番号_付加文字.tpp。
    ///         STL/NPTファイル  : ポリゴングループ名_ランク番号_付加文字.拡張子。
    ///         set_max_memory_size_mb関数にてメモリ使用サイズMAXを
    ///         指定しても、メモリ削減版は動作しません
    ///
    POLYLIB_STAT save(
        std::string&            config_name_out,
        const std::string&      file_format,
        std::string             extend = ""
        );

#ifdef MPI_PL
    ///
    /// 全rank並列でのデータ保存。
    /// 各rankの本クラスインスタンスが保持するグループ階層構造を設定ファイルに各rank毎に書き出す。
    /// 同時にポリゴンデータも指定されたフォーマットのSTL/NPTファイルに各rank毎に書き出す。
    /// 設定ファイル命名規則は以下の通り
    ///   polylib_config_ランク番号_付加文字列.tpp
    /// STL/NPTファイル命名規則は以下の通り
    ///   ポリゴングループ名称_ランク番号_付加文字列.拡張子
    ///
    /// @param[out] config_name_out   設定ファイル名返却用
    /// @param[in] file_format  STL/NPTファイルフォーマット。 "stl_a":アスキー形式　"stl_b":バイナリ形式 "obj_a":アスキー形式　"obj_b","obj_bb":バイナリ形式,"obj_bb"は、頂点法線付き。
    /// @param[in]  extend              ファイル名に付加する文字列。省略可。省略
    ///                                 した場合は、付加文字列として本メソッド呼
    ///                                 び出し時の年月日時分秒(YYYYMMDD24hhmmss)
    ///                                 を用いる。
    /// @return POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT
        save_parallel(
        std::string&  config_name_out,
        const std::string& file_format,
        std::string extend = ""
        );

#endif

    ///
    /// 三角形ポリゴン座標の移動
    /// 本クラスインスタンス配下の全PolygonGroupのmoveメソッドが呼び出される。
    /// moveメソッドは、PolygonGroupクラスを拡張したクラスに利用者が記述する。
    ///
    ///  @param[in] params  Polylib.hで宣言された移動計算パラメータセット。
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT move(
        PolylibMoveParams   &params
        );

#ifdef MPI_PL
    ///
    /// ポリゴンデータのPE間移動
    /// 本クラスインスタンス配下の全PolygonGroupのポリゴンデータについて、
    /// moveメソッドにより移動した三角形ポリゴン情報を隣接PE間でやり取りする。
    ///
    /// @return POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT
        migrate();
#endif

    ///
    /// PolygonGroupツリーの最上位ノードの取得。
    ///
    ///  @return    最上位ノードのvector。
    ///  @attention 返却したPolygonGroupは、削除不可。vectorは要削除。
    ///
    std::vector<PolygonGroup *> *get_root_groups() const;


    /// リーフPolygonGroupリストの取得。
    /// PolygonGroupツリーの末端ノード（リーフ）をリスト化する。
    ///
    /// @return    リーフPolygonGroupのvector.
    /// @attension  返却したPolygonGroupは削除不可。vectorは要削除。
    ///
    std::vector<PolygonGroup *> *get_leaf_groups() const;


    ///
    /// ポリゴンの検索
    /// 位置ベクトルmin_posとmax_posにより特定される矩形領域に含まれる、
    /// ポリゴンをgroup_nameで指定されたグループの下から探索する
    ///
    ///  @param[in,out] tri_list    検索されたポリゴンリスト(Triangle/NptTriangle)
    ///                                 ポリゴンが追加されて返される
    ///  @param[in]     group_name  抽出グループ名
    ///  @param[in]     min_pos     抽出する矩形領域の最小値
    ///  @param[in]     max_pos     抽出する矩形領域の最大値
    ///  @param[in]     every       true:3頂点が全て検索領域に含まれるものを抽出
    ///                             false:3頂点の一部でも検索領域と重なるものを抽出
    ///  @return    POLYLIB_STATで定義される値
    ///  @attention tri_list内のポリゴン要素は、削除不可
    ///             tri_list内のポリゴン要素はサブクラスのNptTriangleである可能性あり
    ///             MPI並列計算時は,min_pos, max_posは各ランクの矩形領域を
    ///             超えないようにして下さい (各ランク内の担当領域内のみ検索するため）
    ///
    POLYLIB_STAT search_polygons(
        std::vector<Triangle*>& tri_list,
        const std::string&      group_name, 
        const Vec3<PL_REAL>&    min_pos, 
        const Vec3<PL_REAL>&    max_pos, 
        const bool              every
        ) const;

    ///
    /// 長田パッチポリゴンの検索
    /// 位置ベクトルmin_posとmax_posにより特定される矩形領域に含まれる、
    /// ポリゴンをgroup_nameで指定されたグループの下から探索する
    ///
    ///  @param[in,out] tri_list    検索された長田パッチポリゴンリスト(NptTriangle)
    ///                                 ポリゴンが追加されて返される
    ///  @param[in] group_name  抽出グループ名
    ///  @param[in] min_pos     抽出する矩形領域の最小値
    ///  @param[in] max_pos     抽出する矩形領域の最大値
    ///  @param[in] every       true:3頂点が全て検索領域に含まれるものを抽出
    ///                         false:3頂点の一部でも検索領域と重なるものを抽出
    ///  @return    POLYLIB_STATで定義される値
    ///  @attention tri_list内のポリゴン要素は、削除不可
    ///             MPI並列計算時は,min_pos, max_posは各ランクの矩形領域を
    ///             超えないようにして下さい (各ランク内の担当領域内のみ検索するため）
    ///             
    ///
    POLYLIB_STAT search_polygons(
        std::vector<NptTriangle*>&  tri_list,
        const std::string&          group_name, 
        const Vec3<PL_REAL>&        min_pos,
        const Vec3<PL_REAL>&        max_pos, 
        const bool                  every
        ) const;

    ///
    /// 指定した点に最も近い三角形ポリゴンの検索
    ///
    ///  @param[out] tri        検索されたポリゴン(Triangle/NptTriangle)
    ///                             != 0 検索されたポリゴン
    ///  @param[in]  group_name 抽出グループ名
    ///  @param[in]  pos        指定した点
    ///  @return    POLYLIB_STATで定義される値
    ///  @attention triは削除不可
    ///             MPI並列計算時は,posは各ランクの矩形領域を
    ///             超えないようにして下さい (各ランク内の担当領域内のみ検索するため）
    ///
    POLYLIB_STAT search_nearest_polygon(
        Triangle*&              tri,
        const std::string&      group_name,
        const Vec3<PL_REAL>&    pos
        ) const;


    ///
    /// 指定した点に最も近い長田パッチポリゴンの検索 
    ///
    ///  @param[out] tri        検索されたポリゴン(NptTriangle)
    ///                             != 0 検索されたポリゴン
    ///  @param[in]  group_name 抽出グループ名
    ///  @param[in]  pos            指定した点
    ///  @return    POLYLIB_STATで定義される値
    ///  @attention triは削除不可
    ///             MPI並列計算時は,posは各ランクの矩形領域を
    ///             超えないようにして下さい (各ランク内の担当領域内のみ検索するため）
    ///
    POLYLIB_STAT  search_nearest_polygon(
        NptTriangle*&           tri,
        const std::string&      group_name,
        const Vec3<PL_REAL>&    pos
        ) const;

    ///
    /// 引数のグループ名が既存グループと重複しないかチェック。
    ///
    ///  @param[in] pg_name     グループ名
    ///  @param[in] parent_path 親グループまでのフルパス
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///  @attention Polylib内部で使用する関数であり、通常は利用者が用いるもの
    ///             ではない。
    ///
    POLYLIB_STAT check_group_name(
        const std::string   &pg_name, 
        const std::string   &parent_path
        );

    ///
    /// PolygonGroupの追加。
    /// 本クラスが管理しているPolygonGroupのリストにPolygonGroupを追加する
    ///
    ///  @param[in] pg      PolygonGroup
    ///
    void add_pg_list(
        PolygonGroup    *pg
        );

    ///
    /// グループ階層構造を標準出力に出力。
    /// 2010.10.20 引数FILE *追加。
    ///  @param[in] fp  出力先ファイル。指定されていれれば、標準出力へ出力する。
    ///
    ///  @attention テスト用の関数であり、通常は利用者が用いるものではない。
    ///
    void show_group_hierarchy(
        FILE    *fp = NULL
        );

    ///
    /// グループの情報と配下の三角形ポリゴン情報を標準出力に出力。
    ///     親グループ名、自身の名前、STLファイル名、属性、面積、座標値
    ///
    ///  @param[in] group_name グループ名
    ///  @param[in] detail     ポリゴンの座標値・法線ベクトルを出力するか否か
    ///  @return    POLYLIB_STATで定義される値が返る
    ///  @attention テスト用の関数であり、通常は利用者が用いるものではない
    ///
    POLYLIB_STAT show_group_info(
            const std::string&      group_name,
            bool detail = false
        );

    ///
    /// 全てのグループの情報と配下の三角形ポリゴン情報を標準出力に出力。
    ///     親グループ名、自身の名前、STLファイル名、属性、面積、座標値
    ///
    ///  @param[in] detail     ポリゴンの座標値・法線ベクトルを出力するか否か
    ///  @return    戻り値なし
    ///  @attention テスト用の関数であり、通常は利用者が用いるものではない
    ///
    void show_all_group_info( bool detail=false );

    ///
    /// Polylibが利用中の概算メモリ量を返す
    ///
    /// @return 利用中のメモリ量(byte)
    ///

    //unsigned int used_memory_size();  // 4GBまでしか使えないので変更
    size_t used_memory_size();


    ///
    /// Polylibが利用中の概算メモリ量(MB)を返す
    ///
    /// @return 利用中のメモリ量(Mbyte)
    ///
    size_t used_memory_size_mb()
    {
        return ( used_memory_size()/(1024*1024) );
    }

#ifdef MPI_PL
    ///
    /// Polylibが利用する最大メモリサイズ(MB)を設定する
    ///
    /// @param[in] max_size_mb 最大メモリサイズ(Mbyte)
    ///               >= 0  0を設定するとメモリ制限なしとなる
    /// @return 戻り値なし
    /// @attention 現状、MPI環境でPolygonGroupのポリゴンのload処理に関係するのみ
    ///            何個に分割してloadするかがこの値で決まる
    ///            設定がない場合は、メモリ制限なしとみなす
    ///            MPI内部で使用する分は考慮しないので、余裕を持って設定すること
    ///
    void set_max_memory_size_mb( int max_size_mb )
    {
        m_max_memory_size_mb = max_size_mb;
    }
#endif

    ///
    /// グループの取得。
    /// nameで与えられた名前のPolygonGroupを返す。
    ///  @param[in] name グループ名
    ///  @return    ポリゴングループクラスのポインタ。エラー時はNULLが返る。
    ///  @attention オーバーロードメソッドあり。
    ///
    PolygonGroup* get_group(
        const std::string&      name
        ) const;

    /**
    * @brief バージョン番号の文字列を返す
    */
    std::string getVersionInfo()
    {
      std::string str(PL_VERSION_NO);
      return str;
    }


protected:

    ///
    /// コンストラクタ
    ///
    /// @attention
    ///   singletonのため、子クラス以外からの呼び出し不可とする
    ///
    Polylib();

    ///
    /// デストラクタ
    ///
    ~Polylib();

#ifdef MPI_PL
    ///
    /// 全PEの担当領域の収集
    ///
    ///  @param[in]  myproc_area    自PEの担当領域
    ///  @param[out] all_procs_area 全PEの担当領域
    ///  @return POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT
        allgather_ParallelAreaInfo(
                 ParallelAreaInfo& myproc_area,
                 std::vector<ParallelAreaInfo>& all_procs_area
             );
#endif

    ///
    /// グループツリー作成。
    /// TextParser クラスを使い、
    /// PolygonGroupを作成し、グループツリーに登録する。
    ///
    ///  @param[in] TextParser のインスタンス
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///  //@attention オーバーロードメソッドあり。
    ///
    ///

    POLYLIB_STAT make_group_tree(
        TextParser *  tp_ptr
        );

    ///
    /// STL/NPTファイルの読み込み。
    /// グループツリーの全リーフについて、設定されているSTL/NPTファイルから
    /// ポリゴン情報を読み込む。読み込んだ後、KD木の生成、法線の計算、面積の
    /// 計算を行う。
    ///
    ///  @param[in] scale       縮尺率
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT load_polygons(
        PL_REAL     scale = 1.0
    );


    ///
    /// 設定ファイルの保存。
    /// メモリに展開しているグループツリー情報から設定ファイルを生成する。
    ///
    ///  @param[in] rank_no ランク番号
    ///  @param[in] extend  ファイル名に付加する文字列
    ///  @param[in] format  PolygonIOクラスで定義されているSTL/NPTファイルのフォー
    ///                     マット。
    ///  @return    作成した設定ファイルの名称。エラー時はNULLが返る。
    ///
    char* save_config_file(
        const std::string& rank_no,
        const std::string& extend,
        const std::string& format
        );

    /// TextParser 内部データから　"filepath" "filepath[*]" というリーフを
    /// すべて削除する.
    ///
    /// recursiveの動作の為、引数にtp_ptrが必要
    ///
    /// @param[in] tp_ptr TextParser　へのポインタ.
    ///  @return    POLYLIB_STATで定義される値が返る。


    POLYLIB_STAT clearfilepath(TextParser* tp_ptr);

    /// TextParser 内部データに　saveしたstl ファイルの　"filepath"を書き込む。
    ///
    ///　saveしたSTLファイルとPolygonGroupの階層は、save_stl_file に
    ///　map を渡し保持してもらう。その　map の内容に基づき、TextParser内部のデータを
    ///　変更する.
    ///
    /// @param[in] polygons_fname_map saveしたSTL/NPTファイルとその階層のmap型データ
    /// @return POLYLIB_STATで定義される値が返る。

    POLYLIB_STAT setfilepath( std::map<std::string,std::string>& polygons_fname_map);

#ifdef MPI_PL
    /// ランク毎にPolygoGroupツリー、ポリゴン情報の保存
    /// グループツリー情報を設定ファイルへ出力。ポリゴン情報をSTL/NPTファイル
    /// へ出力。
    ///
    ///  @param[out] config_name_out  保存した設定ファイル名の返却用。
    ///  @param[in]  myrank         自ランク番号。
    ///  @param[in]  maxrank        最大ランク番号。
    ///  @param[in]  extend         ファイ名に付加される文字列。
    ///  @param[in]  file_format    ファイルフォーマット指定。
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///  @attention ファイル名命名規約は次の通り。
    ///         定義ファイル : polylib_config_ランク番号_付加文字.xml。
    ///         STLファイル  : ポリゴングループ名_ランク番号_付加文字.拡張子。
    ///  @attention save_parallel()の内部実装用メソッド
    POLYLIB_STAT save_at_rank(
        //std::string     *p_config_name,
        std::string&     config_name_out,
        int             myrank,
        int             maxrank,
        const std::string&     extend,
        const std::string&     file_format
        );
#endif

    ///
    /// グループ名の表示。
    /// 指定されたグループ以下の階層構造をツリー形式で標準出力に出力する。
    /// 2010.10.20 引数FILE *追加。
    ///
    ///  @param[in] p   検索の基点となるPolygonGroupのポインタ
    ///  @param[in] tab 階層の深さを示すスペース
    ///  @param[in] fp  出力先ファイル。指定されて行ければ、標準出力へ出力する。
    ///
    void show_group_name(
        PolygonGroup    *p, 
        std::string     tab,
        FILE            *fp
        );

    ///
    /// グループの取得。
    /// internal_idで与えられたm_internal_idを持つPolygonGroupを返す。
    ///  @param[in] internal_id ポリゴングループID
    ///  @return    ポリゴングループクラスのポインタ。エラー時はNULLが返る。
    ///
    PolygonGroup* get_group(
        int internal_id
        ) const;

#ifdef MPI_PL

    ///
    /// プロセス担当領域クラスのポインタを返す
    ///  @param[in] rank ランクNo（自ランク以外）
    ///  @return プロセス担当領域クラスのポインタ
    ///
    ParallelAreaInfo* get_proc_area(int rank);

#endif


private:

    ///
    /// グループの検索。
    /// 基点となるポリゴングループに連なる子孫ポリゴングループを全て抽出する。
    ///  @param[in]  p  探索の基点となるポリゴングループへのポインタ
    ///  @param[out] pg 抽出した子孫ポリゴングループのリスト
    ///
    void search_group(
        PolygonGroup                *p, 
        std::vector<PolygonGroup*>  *pg
        ) const;

    ///
    /// 設定ファイルの保存。 PolylibConfig 内部にあったものをここへ。
    //  暫定措置
    //  file name を作ってsave
    ///
    ///  @param[in] rank_no ランク番号
    ///  @param[in] extend  ファイル名に付加する文字列
    ///  @param[in] format  TriMeshIOクラスで定義されているSTLファイルのフォー
    ///                     マット。
    ///  @return    作成した設定ファイルの名称。エラー時はNULLが返る。
    ///
    char * polylib_config_save_file(
        const std::string&  rank_no,
        const std::string&  extend
    );


#ifdef MPI_PL

    ///
    /// 移動除外三角形IDリストの作成
    ///
    /// @return POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT
        //select_excluded_trias( PolygonGroup *p_pg );
        select_excluded_trias( void );

#endif


protected:
    //=======================================================================
    // クラス変数
    //=======================================================================

    /// 全てのポリゴングループリスト
    std::vector<PolygonGroup*>  m_pg_list;

    /// TextParser へのポインタ
    TextParser* tp;

#ifdef MPI_PL

    /// 自PE担当領域情報
    //ParallelInfo               m_myproc;
    ParallelAreaInfo               m_myproc_area;

    /// 自PEを除く全PE担当領域情報リスト
    //std::vector<ParallelInfo*> m_other_procs;
    std::vector<ParallelAreaInfo> m_other_procs_area; // ポインタを止める

    /// 隣接PE担当領域情報リスト
    //std::vector<ParallelInfo*> m_neibour_procs;
    std::vector<ParallelAreaInfo> m_neibour_procs_area; // ポインタを止める
    
    /// migrate除外三角形IDマップ(k:グループID, v:三角形IDリスト)
    ///    一時情報  m_neibour_procs_area数分あり(順番も同一）
    //     ここに必要かどうか実装時に再検討。削除可能であれば削除する
    std::vector< std::map< int, std::vector<long long int> > > m_exclusion_map_procs;
           //              grp_id           polygon_id

    /// 自プロセスのランクNo
    int m_myrank;

    /// 全プロセス数
    int m_numproc;

    /// MPIコミュニケーター
    MPI_Comm m_comm;

#endif

    /// MAXメモリーサイズ(MB)
    int m_max_memory_size_mb;


};



} //namespace PolylibNS

#endif // polylib_h

