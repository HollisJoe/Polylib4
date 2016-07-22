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

#ifndef polylib_polygongroup_h
#define polylib_polygongroup_h

#include "common/PolylibCommon.h"
#include "common/PolylibStat.h"
#include "common/Vec3.h"
#include "groups/VTree.h"
#include "TextParser.h"
#include "c_lang/CPolylib.h"
#include <vector>
#include <map>

using namespace Vec3class;

namespace PolylibNS {

class Polylib;
//class Polygons;
class PolylibMoveParams;
class Triangle;
class NptTriangle;
class BBox;

// ユーザ定義属性（内部でのみ使用）
struct UsrAtr {
    /// キー
    std::string  key;

    /// 値
    std::string  value;
};

////////////////////////////////////////////////////////////////////////////
///
/// クラス:PolygonGroup
/// ポリゴングループを管理するクラスです。
///
////////////////////////////////////////////////////////////////////////////

class PolygonGroup {
public:
    ///
    /// コンストラクタ
    ///
    PolygonGroup();

    ///
    /// デストラクタ
    ///
    virtual ~PolygonGroup();

    ///
    /// 引数で与えられるポリゴンリストを複製し、KD木の生成を行う。
    ///
    ///  @param[in] tri_list    設定するポリゴンリスト(Triangle)
    ///  @param[in] clear       true:ポリゴン複製、面積計算、KD木生成を行う。
    ///                         false:面積計算、KD木生成だけを行う。
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///  @attention  ポリゴンがNptTriangleに限定できるときは
    ///                引数をstd::vector<NptTriangle*>とした方が高速です。
    ///
    POLYLIB_STAT init(
        const std::vector<Triangle*>    *tri_list, 
        bool                                clear = true
        );


    ///
    /// 引数で与えられるポリゴンリストを複製し、KD木の生成を行う。
    ///
    ///  @param[in] tri_list    設定する長田パッチポリゴンリスト(NptTriangle)
    ///  @param[in] clear       true:ポリゴン複製、面積計算、KD木生成を行う。
    ///                         false:面積計算、KD木生成だけを行う。
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT init(
        const std::vector<NptTriangle*>     *tri_list,
        bool                                clear = true
        );

    ///
    /// ポリゴンのポインタの登録
    ///
    ///  @param[in] tri_list    設定するポリゴンリスト(Triangle)
    ///  @param[in] build_tree  KD木の作成フラグ
    ///                            build_polygon_tree()実行
    ///  @return    なし
    ///  @attention  vectorのポインタがコピーされる。
    //               vector内要素のTriangleの管理はPolygonGroupに移る
    ///              以降は、利用者側でvectorおよびvector内要素(Triangle)をdeleteしないこと
    ///              auto変数等で自動で解放される変数に注意。
    ///              使用法が難しいので一般ユーザは使用不可
    //               
    ///
    void set_triangles_ptr(
             std::vector<Triangle*>  *tri_list,
             bool build_tree = true
        );

    //
    /// PolygonGroupツリーの作成。
    /// 設定ファイルの内容を再帰的に呼び出し、PolygonGroupツリーを作成する。
    ///
    ///  @param[in] polylib     Polygonクラスのインスタンス
    ///  @param[in] parent      親グループ
    ///  @param[in] tp          TextParser のインスタンス
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT build_group_tree(
        Polylib                 *polylib,
        PolygonGroup            *parent,
        TextParser* tp
        );

    ///
    /// ポリゴンの法線ベクトルの計算、面積の計算、KD木の生成を行う。
    ///
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT build_polygon_tree();

    ///
    /// STL/NPTファイルからポリゴン情報を読み込む
    ///    (非メモリ削減版）
    ///
    ///  @return POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT load_polygons_file(PL_REAL scale=1.0);

#ifdef MPI_PL
    ///
    /// STL/NPTファイルの読み込み(MPI メモリ削減版）
    /// ポリゴングループに設定されているSTL/NPTファイルから
    /// ポリゴン情報を読み込む。読み込んだ後、KD木の生成、法線の計算、面積の
    /// 計算を行う。
    ///
    ///  @param[in] scale       縮尺率
    ///  @param[in] size_mb     使用可能メモリ(Mbyte）
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT load_polygons_mem_reduced(
        PL_REAL        scale,
        int            size_mb
    );
#endif

    ///
    /// ポリゴン情報をSTL/NPTファイルに出力する。
    /// TextParser 対応版
    ///  @param[in] rank_no ファイル名に付加するランク番号。
    ///  @param[in] extend  ファイル名に付加する自由文字列。
    ///  @param[in] format  ファイルフォーマット。
    ///  @param[in,out] polygons_fname_map STL/NPT ファイル名とポリゴングループのパス
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///  @attention  ポリゴングループが持つポリゴンリストを使用する
    ///              MPI環境の場合、各ランク内のポリゴンリストを使う
    ///
    POLYLIB_STAT save_polygons_file(
        const std::string&     rank_no,
        const std::string&     extend,
        const std::string&     format,
        std::map<std::string,std::string>& polygons_fname_map
        );

    ///
    /// ポリゴン情報をSTL/NPTファイルに出力する
    /// TextParser 対応版
    ///  @param[in] rank_no ファイル名に付加するランク番号。
    ///  @param[in] extend  ファイル名に付加する自由文字列。
    ///  @param[in] format  ファイルフォーマット。
    ///  @param[in] tri_list 三角形ポリゴンリストの領域
    ///  @param[in,out] polygons_fname_map STL/NPT ファイル名とポリゴングループのパス
    ///  @return    POLYLIB_STATで定義される値が返る
    ///  @attention  ポリゴングループが持つポリゴンリストではなく、
    ///              外部から与えたポリゴンリストを使うことに注意
    ///              MPI環境においてRank0に集約したポリゴンリストを使う想定
    //               ポリゴングループのパスは使用する
    ///
    POLYLIB_STAT save_polygons_file(
        const std::string&     rank_no,
        const std::string&     extend,
        const std::string&     format,
        std::vector<Triangle*> *tri_list,
        std::map<std::string,std::string>& polygons_fname_map
        );

#ifdef MPI_PL
    ///
    /// ポリゴン情報分散
    /// ポリゴングループ毎にRank0に集約されているポリゴン情報を他ランクに分散する
    /// （非メモリ削減版）
    ///
    /// @return POLYLIB_STATで定義される値が返る。
    /// @attention 
    ///   (IN)  ランク0    : ポリゴングループに全ポリゴン設定
    ///   (IN)  ランク0以外: ポリゴングループにポリゴン設定なし
    ///   (OUT) 全ランク   : ポリゴングループに担当領域内のポリゴン設定
    ///
    POLYLIB_STAT
        scatter_polygons( void );

    ///
    /// ポリゴン情報分散
    /// ポリゴングループ毎にRank0に集約されているポリゴン情報を他ランクに分散する
    /// （メモリ削減版）
    ///     メモリ制限により分割されたポリゴン毎の処理となる
    ///
    /// @param[in,out] tri_list_div_all   ポリゴングループの分割された全ポリゴン
    ///                                  rank0 : (in)  分割された全ポリゴン
    ///                                          (out) 空リスト
    ///                                  rank0以外 : (in/out) 空リスト
    /// @param[out] tri_list_div_local ポリゴングループの分割されたポリゴン
    ///                                  (担当領域内のポリゴン）
    ///                                     (in) 空リスト
    /// @return POLYLIB_STATで定義される値が返る。
    /// @attention 
    ///   (IN)  全ランク  : ポリゴングループにポリゴン設定なし
    ///   (OUT) 全ランク  : ポリゴングループにポリゴン設定なし
    ///   ランク０において担当領域内検索を行うため一時的に
    ///   ポリゴングループにポリゴンtri_list_div_allが設定される
    ///
    POLYLIB_STAT
        scatter_polygons(
            std::vector<Triangle*>  &tri_list_div_all,
            std::vector<Triangle*>  &tri_list_div_local
        );

    ///
    /// ポリゴン情報集約
    /// ポリゴングループ(複数）の各ランクに分散されているポリゴン情報をRank0に集約する
    ///
    ///  @param[in] pg     ポリゴングループ
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT
        gather_polygons(
            std::vector<Triangle*>  &tri_list
        );


#endif

    ///
    /// 三角形ポリゴン移動メソッド
    ///     カスタマイズのためにはset_move_func()で移動関数を登録しておくこと
    ///
    ///  @param[in] params  Polylib.hで宣言しているパラメタセットクラス。
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT move(
        PolylibMoveParams   &params
        );

    ///
    /// 三角形ポリゴン移動関数登録
    ///     move()と違い、オブジェクトのインスタンス毎に登録が必要
    ///
    ///  @param[in] func    移動関数へのポインタ
    ///                         引数のPolygonGroup*はthisが設定されて呼ばれる
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT set_move_func(
            void (*func)(PolygonGroup*,PolylibMoveParams*) 
        )
        {
            m_move_func = func;
            return PLSTAT_OK;
        };

    ///
    /// 三角形ポリゴン移動関数登録
    ///     move()と違い、オブジェクトのインスタンス毎に登録が必要
    ///
    ///  @param[in] func_c  移動関数(Cの関数）へのポインタ
    ///                         引数のPolygonGroup*はthisが設定されて呼ばれる
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT set_move_func_c(
            //extern "C" void (*func)( PL_GRP_TAG, PolylibMoveParamsStruct* ) 
                        //  Cの関数ポインタの登録はうまく動作しないかも？
            void (*func)( PL_GRP_TAG, ::PolylibMoveParamsStruct* ) 
        )
        {
            m_move_func_c = func;
            return PLSTAT_OK;
        };

    ///
    /// KD木の再構築フラグの設定
    ///     ユーザ定義の移動関数内の最後で呼び出す
    ///
    ///  @return    戻り値なし
    ///
    void set_need_rebuild( void )
        {
            m_need_rebuild = true;
        };


    ///
    /// 指定矩形領域に含まれるポリゴンを抽出する。
    ///
    ///  @param[in,out] tri_list    検索されたポリゴンリスト(Triangle/NptTriangle)
    ///  @param[in]     bbox        矩形領域。
    ///  @param[in]     every       true:3頂点が全て検索領域に含まれるものを抽出。
    ///                             false:1頂点でも検索領域に含まれるものを抽出。
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///  @attention tri_list内のポリゴン要素は、削除不可
    ///             tri_list内のポリゴン要素はサブクラスのNptTriangleである可能性あり
    ///             内部でKD木探索実施
    POLYLIB_STAT search(
        std::vector<Triangle*>&     tri_list,
        const BBox&                bbox, 
        bool                        every 
        ) const;



    ///
    /// 指定矩形領域（複数）に含まれるポリゴンを抽出する
    ///
    ///  @param[in,out] tri_list    検索されたポリゴンリスト(Triangle/NptTriangle)
    ///  @param[in]     bbox        矩形領域。
    ///  @param[in]     every       true:3頂点が全て検索領域に含まれるものを抽出。
    ///                             false:1頂点でも検索領域に含まれるものを抽出。
    ///  @param[in]     duplicate  ポリゴンID 重複指定 false:重複なし true:重複の可能性あり
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///  @attention tri_list内のポリゴン要素は、削除不可
    ///             tri_list内のポリゴン要素はサブクラスのNptTriangleである可能性あり
    ///             tri_listはクリアされず、検索されたポリゴンが追加される。
    ///             内部でKD木探索実施
    ///             各ランクの担当領域が複数対応となったため追加
    ///             複数領域のため、単純に加算するとIDが重複する可能性あり
    POLYLIB_STAT search(
        std::vector<Triangle*>&     tri_list,
        const std::vector<BBox>&    bboxes, 
        bool                        every, 
        bool                        duplicate = false
        ) const;

    ///
    /// 指定矩形領域に含まれるポリゴンを抽出する。
    ///
    ///  @param[in,out] tri_list    検索されたポリゴンリスト(NptTriangle)
    ///  @param[in]     bbox        矩形領域。
    ///  @param[in]     every       true:3頂点が全て検索領域に含まれるものを抽出。
    ///                             false:1頂点でも検索領域に含まれるものを抽出。
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///  @attention tri_list内のポリゴン要素は、削除不可
    ///             tri_listはクリアされず、検索されたポリゴンが追加される
    ///             内部でKD木探索実施
    POLYLIB_STAT search(
        std::vector<NptTriangle*>&  tri_list,
        const BBox &                bbox, 
        bool                        every 
        ) const;


    /// 
    /// 指定位置に最も近いポリゴンを検索する
    ///
    ///  @param[out] tri      検索されたポリゴン(Triangle/NptTriangle)
    ///                             != 0 検索されたポリゴン
    ///  @param[in]   pos     指定位置
    ///  @return    POLYLIB_STATで定義される値
    ///  @attention triは削除不可
    ///             内部でKD木探索実施
    ///
    POLYLIB_STAT  search_nearest(
        Triangle*&              tri,
        const Vec3<PL_REAL>&    pos
        ) const;

    ///
    /// PolygonGroupのフルパス名を取得する。
    ///
    ///  @return フルパス名。
    ///
    std::string acq_fullpath();

    ///
    /// カンマ区切りでSTLファイル名リストを取得。
    ///
    ///  @return ファイル名リスト。
    ///
    std::string acq_file_name();

#ifdef MPI_PL

    ///
    /// 自領域内ポリゴンのみ抽出してポリゴンを返す
    ///     ポリゴン自体の複製は行っていない
    /// 
    /// @param[out] tri_list  自身の担当領域内のポリゴン
    /// @return POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT
        get_inbounded_polygons(
                std::vector<Triangle*>  &tri_list
            );

    ///
    /// 自領域内ポリゴンのみ抽出してポリゴン情報を再構築
    /// 
    /// @return POLYLIB_STATで定義される値が返る。
    /// @attention  ポリゴンのload処理内およびmigrate処理後に実行する
    ///
    POLYLIB_STAT
        erase_outbounded_polygons( void );

    // 以下はMPI環境のPolylib実装用

    ///
    /// PE領域間移動する三角形ポリゴンリストの取得。
    ///
    ///  @param[in,out] tri_list        検索されたポリゴンリスト(Triangle)
    ///  @param[in] neibour_bboxes      隣接PE領域バウンディングボックス
    ///  @param[in] exclude_tria_ids    領域移動対象外三角形IDリスト
    ///  @return    検索結果三角形リスト
    ///  @attention 戻り値は使用後領域を解放すること
    ///
    //const std::vector<Triangle*>* search_outbounded(
        //const BBox &        neibour_bbox,
    POLYLIB_STAT search_outbounded(
        std::vector<Triangle*>&     tri_list,
        std::vector<BBox>&          neibour_bboxes,
        std::vector<long long int>& exclude_tria_ids
        );

    ///
    /// 三角形リストの追加（ID重複は追加されない）
    ///
    ///  @param[in] tri_list    三角形ポリゴンリストのポインタ。
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///  @attention 内部IDが重複した三角形は追加しない。KD木の再構築はしない。
    ///
    POLYLIB_STAT add_triangles(
        const std::vector<Triangle*>&  tri_list
        );
#endif

    ///
    /// ポリゴン情報を再構築する。（KD木の再構築をおこなう）
    ///
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT rebuild_polygons();

    ///
    /// グループ情報（ランク番号、親グループ名、自分のグループ名、ファイル名、
    /// 頂点数、各頂点のXYZ座標値、法線ベクトルのXYZ座標値、面積）を出力する。
    ///
    ///  @param[in] irank ランクＮｏ(この値を出力しているのみ）
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT show_group_info(
            int irank   = -1,
            bool detail = false
        );

    // add keno 20120331
    /// ポリゴングループの要素数を返す
    ///  @return    ポリゴングループの要素数
    ///  @attention 並列環境の場合、各ランク担当分の要素数
    ///     全プロセス通した要素数が必要な場合は、
    ///     get_group_num_global_tria()を使用する
    int get_group_num_tria( void );

#ifdef MPI_PL
    /// ポリゴングループの要素数(global)を返す
    ///  @return    ポリゴングループの要素数(global)
    ///  @attention 並列環境用
    ///     ポリゴンの重複を削除するための通信あり
    ///     get_group_num_tria()よりも大幅に処理時間がかかることに注意
    int get_group_num_global_tria( void );
#endif

    /// グループ内のポリゴンの面積を積算して返す
    ///     並列化されている場合は該当ランク内のみの面積
    /// @attention
    ///     全プロセス通した面積が必要な場合は、
    ///     get_group_global_area()を使用する
    ///     ポリゴンがない場合は0.0を返す
    PL_REAL get_group_area( void );

#ifdef MPI_PL
    /// グループ内のポリゴンの面積(global)を積算して返す
    ///     全プロセスを通算した面積(重複ポリゴン分は無視される）
    ///     全プロセスに同じ値が返る
    ///     ポリゴンがない場合は0.0を返す
    PL_REAL get_group_global_area( void );
#endif

    /// グループ内のポリゴン属性（整数）の集合演算値を返す
    ///     並列化されている場合は全プロセスを通した値
    ///     (PL_OP_SUM：重複ポリゴン分は無視される）
    ///     全プロセスに同じ値が返る
    ///  @param[in]  op     演算種類　PL_OP_SUM/PL_OP_MAX/PL_OP_MIN
    ///  @param[in]  atr_no ポリゴン整数属性の何番目か　0〜
    ///  @param[out] val    属性値
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///                ポリゴンが存在しない
    ///                ポリゴン属性が存在しないなど
    ///
    POLYLIB_STAT get_polygons_reduce_atrI(
            PL_OP_TYPE op,
            int        atr_no,
            int&       val
        );

    /// グループ内のポリゴン属性（実数）の集合演算値を返す
    ///     並列化されている場合は全プロセスを通した値
    ///     (PL_OP_SUM：重複ポリゴン分は無視される）
    ///     全プロセスに同じ値が返る
    ///  @param[in]  op     演算種類　PL_OP_SUM/PL_OP_MAX/PL_OP_MIN
    ///  @param[in]  atr_no ポリゴン実数属性の何番目か　0〜
    ///  @param[out] val    属性値
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///                ポリゴンが存在しない
    ///                ポリゴン属性が存在しないなど
    ///
    POLYLIB_STAT get_polygons_reduce_atrR(
            PL_OP_TYPE op,
            int        atr_no,
            PL_REAL&   val
        );


    /// ポリゴンの縮尺変換＆KD木再構築
    POLYLIB_STAT rescale_polygons( PL_REAL scale );

/*  当メソッドは削除。各アプリにて実装するものとする
    このメソッドはFFV-Cで使われているもの
    ///
    /// グループ配下の全Triangleオブジェクトのm_exidを更新する。
    ///
    /// @param[in] id 更新するID番号。
    /// @return    ステータスコード。
    ///
    POLYLIB_STAT set_all_exid_of_trias(
        int id
        );
*/

    ///
    /// move()メソッド実行により、頂点が隣接セルよりも遠くへ移動した三角形情報
    /// を報告（前処理）
    ///
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///  @attention 本メソッドはデバッグ用です。
    ///         set_move_func()で登録した移動関数内で座標移動処理前に呼ぶこと。
    ///
    POLYLIB_STAT init_check_leaped();

    ///
    /// move()メソッド実行により、頂点が隣接セルよりも遠くへ移動した三角形情報
    /// を報告（後処理）。該当する三角形について、以下の情報をcerrへ出力する。
    ///     ・ポリゴングループID
    ///     ・三角形ID
    ///     ・移動前/後の頂点座標
    ///
    ///  @param[in] origin      計算領域起点座標
    ///  @param[in] cell_size   ボクセルサイズ
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///  @attention 本メソッドはデバッグ用です。
    ///         set_move_func()で登録した移動関数内で座標移動処理後に呼ぶこと。
    ///
    POLYLIB_STAT check_leaped(
            std::vector< Vec3<PL_REAL> >& origin,
            std::vector< Vec3<PL_REAL> >& cell_size
        );


    //=======================================================================
    // Setter/Getter
    //=======================================================================

    ///
    /// ポリゴングループIDを取得
    ///   システム内でユニークなID
    /// メンバー名修正( m_id -> m_internal_id) 2010.10.20
    ///
    ///  @return ポリゴングループID。
    ///
    int get_internal_id()
    {
        return m_internal_id;
    };

    ///
    /// グループ名を取得。
    ///
    /// @return グループ名。
    ///
    std::string get_name(void)
    {
        return m_name;
    };

    ///
    /// グループ名を設定。
    ///
    /// @param[in] name グループ名。
    ///
    void set_name(const std::string& name)
    {
        m_name = name;
    };

    ///
    /// 移動対象フラグを取得。
    ///
    ///  @return 移動対象フラグ。
    ///
    bool get_movable()
    {
        return m_movable;
    };


    ///
    /// 移動対象フラグを取得。
    ///
    ///  @param[in] 移動対象フラグ。
    ///
    void set_movable( bool movable )
    {
        m_movable = movable;
    };

    ///
    /// ポリゴングループのユーザ定義属性取得。
    ///
    ///  @param[in]     key     キー
    ///  @param[out]    val     属性値
    ///  @return OK/NG  NG:キーと属性のペアが登録されていない
    ///
    POLYLIB_STAT get_atr( std::string& key, std::string& val ) const;

    ///
    /// ポリゴングループのユーザ定義属性設定
    ///
    ///  @param[in]     key     キー
    ///  @param[in]     val     属性値
    ///  @return なし
    ///  @attention 既に登録されていた場合、上書きする
    ///
    void set_atr( std::string& key, std::string& val );

    ///
    /// ポリゴン(Triangle/NptTriangle)のユーザ定義属性数（整数型）の取得
    ///
    ///  @return    ユーザ定義属性数（整数型）
    ///  @attention ポリゴングループにポリゴンが存在しない場合は０が返る
    ///      並列環境でポリゴンが存在しないランクがあるかもしれないので注意
    ///
    int get_num_polygon_atrI( void ) {
        if( m_tri_list != NULL && m_tri_list->size() > 0 ) {
            return (*m_tri_list)[0]->get_num_atrI();
        } else {
            return 0;
        }
    }

    ///
    /// ポリゴン(Triangle/NptTriangle)のユーザ定義属性数（実数型）の取得
    ///
    ///  @return    ユーザ定義属性数（実数型）
    ///  @attention ポリゴングループにポリゴンが存在しない場合は０が返る
    ///      並列環境でポリゴンが存在しないランクがあるかもしれないので注意
    ///
    int get_num_polygon_atrR( void )
    {
        if( m_tri_list != NULL && m_tri_list->size() > 0 ) {
            return (*m_tri_list)[0]->get_num_atrR();
        } else {
            return 0;
        }
    }

    ///
    /// ポリゴン(Triangle/NptTriangle)のユーザ定義属性数の設定
    ///    PolygonGroup内の全ポリゴンに属性数を設定する
    ///
    /// @param[in] num_atrI    ユーザ定義属性数（整数型）
    /// @param[in] num_atrR    ユーザ定義属性数（実数型）
    /// @return    なし
    /// @attention ポリゴングループにポリゴンが存在しない場合は何もしない
    ///      並列環境でポリゴンが存在しないランクがあるかもしれないので注意
    ///      エラーとはしていない
    ///
    void set_num_polygon_atr(
            int num_atrI,
            int num_atrR
        )
    {
        if( m_tri_list == NULL ) return;

        for(int i=0; i<m_tri_list->size(); i++ ) {
            (*m_tri_list)[i]->set_num_atr( num_atrI,num_atrR );
        }
    }

    ///
    /// 親グループのフルパス名を設定。
    ///
    /// @param[in] ppath 親グループのフルパス名。
    ///
    void set_parent_path(std::string ppath)
    {
        m_parent_path = ppath;
    };
    
    ///
    /// 親グループのフルパス名を取得。
    ///
    /// @return 親グループのフルパス名。
    ///
    std::string get_parent_path(void)
    {
        return m_parent_path;
    };

    ///
    /// 親グループを設定。
    ///
    /// @param[in] p 親グループのポインタ。
    ///
    void set_parent(PolygonGroup* p) 
    {
        m_parent = p;
    };

    ///
    /// 親グループを取得
    ///
    /// @return 親グループのポインタ。
    ///
    PolygonGroup* get_parent(void)
    {
        return m_parent;
    };

    ///
    /// 子グループを設定
    ///
    /// @param[in] p    子グループのリスト。
    ///
    void set_children(std::vector<PolygonGroup*>& p)
    {
        m_children = p;
    };

    ///
    /// 子グループを取得
    ///
    /// @return 子グループのリスト（参照型）
    ///
    std::vector<PolygonGroup*>& get_children(void)
    {
        return m_children;
    };

    ///
    /// 子グループを追加
    ///
    /// @param[in] p    子グループ。
    ///
    void add_children(PolygonGroup* p) 
    {
        m_children.push_back(p);
    };

    ///
    /// 子グループを削除
    ///
    /// @param[in] p    子グループ。
    ///
    void remove_child(PolygonGroup* p)
    {
        std::vector<PolygonGroup*>::iterator itr;
        for( itr=m_children.begin(); itr!=m_children.end(); itr++ ) {
            if( *itr == p ) {
                itr = m_children.erase(itr);
                break;
            }
        }
    }

    ///
    /// STL/NPTファイル名とファイルフォーマットを設定。
    ///
    ///  @param[in] fname STLファイル名とファイルフォーマットの対応マップ。
    ///
    void set_file_name(std::map<std::string, std::string> fname) 
    {
        m_polygon_files = fname;
    };

    ///
    /// STL/NPTファイル名とファイルフォーマットの対応マップ取得。
    ///
    ///  @return STLファイル名とファイルフォーマットの対応マップ。
    ///
    std::map<std::string, std::string> get_file_name() const 
    {
        return m_polygon_files;
    };

    ///
    /// ポリゴンリストを取得。
    ///
    /// @return 三角形ポリゴンリスト
    ///             Triangle/NptTriangle
    ///
    std::vector<Triangle*>* get_triangles()
    {
        //return m_polygons->get_tri_list();
        return m_tri_list;
    };

    ///
    /// KD木オブジェクトを取得。
    ///
    /// @return KD木ポリゴンリスト
    /// @attention  ユーザは使用不可。Polylib::used_memory_size()のみ使用
    ///
    VTree *get_vtree()
    {
        //return m_polygons->get_vtree();
        return m_vtree;
    };


    ///
    /// move()による移動前三角形一時保存リストの個数を取得。
    ///
    ///  @return 一時保存リストサイズ。
    ///
    size_t get_num_of_trias_before_move()
    {
        if (m_trias_before_move == NULL)    return 0;
        else                                return m_trias_before_move->size();
    };

    ///
    /// configファイルに記述するParamタグのクラス名(value="...")。
    ///
    //  static const char *ATT_NAME_config;
    static const char *ATT_NAME_CLASS;

protected:

    ///
    /// 設定ファイルから取得したPolygonGroupの属性情報をインスタンスにセットする。
    ///       ユーザ定義属性以外を設定する
    /// 
    /// "filepath" に関して、先にfilepathが複数　(filepath[0])が存在するかどうか
    ///  をチェックして、複数ならばその処理を行い、filepath の処理は終了する。
    ///  複数でないことが分かったら、filepath が単体で存在するかをチェックして、
    ///  存在するならば、処理を行う。
    ///
    ///  @param[in] polylib     Polygonクラスのインスタンス。
    ///  @param[in] parent      親グループ。
    ///  @param[in] tp              TextParserクラスのインスタンス
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT setup_attribute (
        Polylib                 *polylib,
        PolygonGroup            *parent,
        TextParser *tp
        );

    ///
    /// 設定ファイルから取得したPolygonGroupのユーザ定義属性情報をインスタンスにセットする。
    ///
    ///  @param[in] node_label      ユーザ定義属性ノードの名前
    ///  @param[in] tp              TextParserクラスのインスタンス
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT setup_user_attribute (
        const std::string&                 node_label,
        TextParser *tp
        );

    ///
    /// 2点が隣接ボクセルよりも離れているか？
    ///
    ///  @param[in] origin      計算領域起点座標。
    ///  @param[in] cell_size   ボクセルサイズ。
    ///  @param[in] pos1            点(1)。
    ///  @param[in] pos2            点(2)。
    ///  @return    true:2点が隣接ボクセルよりも離れている。
    ///
    bool is_far(
        Vec3<PL_REAL>& origin,
        Vec3<PL_REAL>& cell_size,
        Vec3<PL_REAL>& pos1,
        Vec3<PL_REAL>& pos2
        );


private:
    ///
    /// 三角形ポリゴンリストの初期化
    ///
    void init_tri_list();

    ///
    /// 三角形ポリゴンリストの削除
    ///         リスト毎削除
    ///
    void delete_tri_list();

    ///
    /// STLファイル名を作成。ファイル名は、以下の通り。
    /// グループ名のフルパス_ランク番号_自由文字列.フォーマット文字列。
    ///
    ///  @param[in] rank_no ファイル名に付加するランク番号。
    ///  @param[in] extend  ファイル名に付加する自由文字列。
    ///  @param[in] format  STL/NPTファイルフォーマット。
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    char *mk_polygons_fname(
        const std::string&      rank_no,
        const std::string&      extend,
        const std::string&      format
        );


    ///
    /// ポリゴンファイル名を作成。ファイル名は、以下の通り。
    /// グループ名のフルパス_ランク番号_自由文字列.フォーマット文字列。
    /// TextParser 対応版
    ///
    ///  @param[in] rank_no ファイル名に付加するランク番号。
    ///  @param[in] extend  ファイル名に付加する自由文字列。
    ///  @param[in] format  STL/NPTファイルフォーマット。
    ///  @param[in,out] polygons_fname_map ポリゴンファイル名とポリゴングループのパス
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    char *mk_polygons_fname(
        const std::string&      rank_no,
        const std::string&      extend,
        const std::string&      format,
        std::map<std::string,std::string>& polygons_fname_map
        );


    ///
    /// 全PolygonGroupに一意のグループIDを作成する。
    ///
    ///  @return    グループID。
    ///
    //int create_global_id();       //fj><  m_internal_id採番用
    int create_global_unique_id();      //  メソッド名変更, m_internal_id採番用

#ifdef MPI_PL
    /// グループ内のポリゴン属性（整数）の総和を返す
    ///     private & 並列化専用関数
    ///     並列化環境で重複ポリゴンを省く処理を実装している
    ///  @param[in]  ids    ポリゴンID
    ///  @param[in]  atrs   ポリゴン属性（整数）
    ///  @param[out] sum    属性総和値（整数）
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT get_polygons_reduce_sum_atrI(
            std::vector<long long int>& ids,
            std::vector<int>&           atrs,
            int&   sum
        );

    /// グループ内のポリゴン属性（実数）の総和を返す
    ///     private & 並列化専用関数
    ///     並列化環境で重複ポリゴンを省く処理を実装している
    ///  @param[in]  ids    ポリゴンID
    ///  @param[in]  atrs   ポリゴン属性（実数）
    ///  @param[out] sum    属性総和値（実数）
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT get_polygons_reduce_sum_atrR(
            std::vector<long long int>& ids,
            std::vector<PL_REAL>&       atrs,
            PL_REAL&   sum
        );

    ///
    /// ポリゴン情報分散
    /// ランク０の処理： ランク０以外に該当ポリゴンを送信する
    ///
    /// @param[in]  pg     ポリゴングループ
    /// @return POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT
        scatter_send_polygons( void );

    ///
    /// ポリゴン情報分散
    /// ランク０以外の処理： ランク０より該当ポリゴンを受信する
    ///
    /// @param[out] tri_list  担当領域内のポリゴン
    /// @return POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT
        scatter_receive_polygons(
            std::vector<Triangle*>  &tri_list
        );

    ///
    /// ポリゴン情報集約
    /// ランク０の処理： 分割データをランク０に集約する
    ///
    /// @param[out] tri_list  ポリゴングループ内の全ポリゴン
    /// @return POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT
        gather_recv_polygons(
            std::vector<Triangle*>  &tri_list
        );

    ///
    /// ポリゴン情報集約
    /// ランク０以外の処理： ランク０に該当ポリゴンを送信する
    ///
    /// @return POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT
        gather_send_polygons( void );

#endif

protected:
    //=======================================================================
    // クラス変数
    //=======================================================================

    //--------------------------------------
    // 属性
    //--------------------------------------

    /// グループID。
    //  全PolygonGroupに一意のID
    //                              内部処理用識別ID、外部ファイル(.tppファイル）に出力しない
    //                              全体でユニークなIDとなる
    int                 m_internal_id;

    /// 自グループ名。
    //                              .tppファイルに出力
    std::string         m_name;

    /// moveメソッドにより移動するグループか？
    //                              == trueの時に .tppファイルに出力する
    bool                m_movable;


    /// ユーザ定義属性
    //                              tpファイルに出力する
    //                              Polylib2.xの m_label, m_typeはこちらに設定する
    std::vector<UsrAtr> m_atr;


    //--------------------------------------
    // 親子情報
    //--------------------------------------

    /// 親グループのパス名。
    std::string         m_parent_path;

    /// 親グループへのポインタ。
    PolygonGroup            *m_parent;

    /// 子グループへのポインタリスト。
    std::vector<PolygonGroup*>  m_children;


    //--------------------------------------
    // ポリゴン形状関連
    //--------------------------------------

    /// 形状ファイル名とファイル形式
    ///     フィル名とファイル形式の対
    std::map<std::string, std::string>  m_polygon_files;    // // m_file_nameより名前を変更

    //-------------------------------------
    //  旧三角形Polygonsクラスのメンバー変数
    //-------------------------------------
    // Polygons                *m_polygons;

    /// 三角形ポリゴンのリスト
    std::vector<Triangle*>  *m_tri_list;


    /// 全三角形ポリゴンを外包するBoundingBox
    BBox    m_bbox;

    /// KD木クラス
    VTree   *m_vtree;

    /// MAX要素数
    int     m_max_elements;

    //--------------------------------------
    // 移動関数へのポインタ
    //--------------------------------------
    void (*m_move_func)(PolygonGroup*,PolylibMoveParams*);

    void (*m_move_func_c)(PL_GRP_TAG, ::PolylibMoveParamsStruct*);


    //--- 実装用 start ---------

    /// KD木の再構築が必要か？
    bool                    m_need_rebuild;

    /// move()による移動前三角形一時保存リスト。
    //fj><  実装時再検討、このクラスのメンバーから削除可能であれば削除する
    std::vector<Triangle*>      *m_trias_before_move;

    //--- 実装用 end ---------
};

} //namespace PolylibNS

#endif //polylib_polygongroup_h
