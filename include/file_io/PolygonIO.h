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

#ifndef polygon_io_h
#define polygon_io_h

#include <vector>
#include <map>
#include <fstream>
#include "common/PolylibStat.h"
#include "common/PolylibCommon.h"
#include "polygons/Triangle.h"
#include "polygons/NptTriangle.h"

namespace PolylibNS {

using namespace std;

////////////////////////////////////////////////////////////////////////////
///
/// クラス:PolygonIO
/// 三角形ポリゴン入出力管理。
///
////////////////////////////////////////////////////////////////////////////
//fj><  クラス名変更  旧クラス名：TriMeshIO  --> PolygonIO


class PolygonIO {
public:
    ///
    /// STL/NPTファイルを読み込み、tri_listにセットする。
    ///     複数個のポリゴンファイルを読み込む
    ///
    ///  @param[in,out] tri_list    三角形ポリゴンリストの領域。
    ///  @param[in]     fmap        ファイル名、ファイルフォーマットのセット。
    ///                                複数指定可
    ///  @param[in]     scale       スケール
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    static POLYLIB_STAT load(
        std::vector<Triangle*>              *tri_list,
        const std::map<std::string, std::string>    &fmap,
        PL_REAL scale = 1.0
    );

    ///
    /// STL/NPTファイルを読み込み、tri_listにセットする。
    ///     １個のポリゴンファイルを読み込む
    ///
    ///  @param[in,out] tri_list    三角形ポリゴンリストの領域。
    ///  @param[in] fname           ファイル名
    ///  @param[in] fmt             ファイルフォーマット
    ///  @param[in] scale       スケール
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    static POLYLIB_STAT load(
        std::vector<Triangle*>              *tri_list,
        const std::string&                  fname, 
        const std::string&                  fmt,
        PL_REAL scale = 1.0
    );

    ///
    /// STL/NPTファイルのOpen
    ///     ポリゴンを分割してロードする時に使用する
    ///
    ///  @param[out] ifs            入力ファイルストリーム
    ///  @param[in] fname           ファイル名
    ///  @param[in] fmt             ファイルフォーマット
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///  @attention  STLとNPTのバイナリファイルに関しては
    ///              内部でヘッダ部の読み出しも行う
    ///
    static POLYLIB_STAT load_file_open(
        std::ifstream&                      ifs,
        const std::string&                  fname, 
        const std::string&                  fmt
    );

    ///
    /// STL/NPTファイルのClose
    ///     ポリゴンを分割してロードする時に使用する
    ///
    ///  @param[in] ifs            入力ファイルストリーム
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    static POLYLIB_STAT load_file_close( ifstream& ifs );

    ///
    /// STL/NPTファイルのRead（指定個数読み込む）
    ///     ポリゴンを分割してロードする時に使用する
    ///
    ///  @param[in]  ifs          入力ファイルストリーム
    ///  @param[in]  fmt          ファイルフォーマット
    ///  @param[in,out] tri_list     三角形ポリゴンのリスト(出力内容)
    ///  @param[in]  num_read     読み込み指定数
    ///                             EOFに達すれば途中まで読み込まれる
    ///                             -1 の時、読み込み数制限なし（eofまで読む）
    ///  @param[out] num_tri      実際に読み込んだ数
    ///  @param[out] eof          ファイル終了フラグ(end of file)
    ///  @param[in]  scale       スケール
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    static POLYLIB_STAT load_file_read(
        ifstream&                ifs,
        const std::string&       fmt,
        std::vector<Triangle*>&  tri_list,
        int                      num_read,
        int&                     num_tri,
        bool&                    eof,
        PL_REAL scale = 1.0
    );

    ///
    /// tri_listの内容をSTL形式でファイルへ保存。
    ///
    ///  @param[in] tri_list    三角形ポリゴンのリスト(出力内容)
    ///  @param[in] fname       ファイル名
    ///  @param[in] fmt         ファイルフォーマット
    ///  @return    POLYLIB_STATで定義される値が返る
    ///  @attention 長田パッチの場合、std::vector<Triangle*>* に
    ///     キャストすること
    ///     static_castではコンパイルエラーとなる場合、旧キャストを使用
    ///       (std::vector<Triangle*>*)ptr
    ///
    static POLYLIB_STAT save(
        std::vector<Triangle*>  *tri_list,
        const std::string&      fname, 
        const std::string&      fmt
    );

    ///
    /// ファイル名を元に入力ファイルのフォーマットを取得する。
    ///
    ///  @param[in] filename        入力ファイル名。
    ///  @return    判定したファイルフォーマット。
    ///  @attention ファイル拡張子が"stl"の場合、ファイルを読み込んで判定する。
    ///
    static std::string input_file_format(
        const std::string &filename
    );

    ///
    /// ファイルフォーマットより拡張子を求める
    ///
    ///  @param[in] filename        ファイルフォーマット
    ///                                 PolygonIO::FMT_STL_A 等
    ///  @return    ファイル拡張子
    ///
    static std::string get_extension_format(
        const std::string&  fmt
    );

    ///
    /// ファイルフォーマットよりポリゴンタイプを求める
    ///
    ///  @param[in] filename        ファイルフォーマット
    ///                                 PolygonIO::FMT_STL_A 等
    ///  @return    ポリゴンタイプ
    ///                 PL_TYPE_TRIANGLE/PL_TYPE_NPT/PL_TYPE_UNKNOWN
    ///
    static int get_polygon_type(
        const std::string&  fmt
    );

    /// STL/NPTファイルのフォーマット種別
    ///
    ///  @attention STL/NPTファイルの拡張子とは異なるので注意すること。
    ///
    static const std::string FMT_STL_A;     ///< STLアスキーファイル
    static const std::string FMT_STL_AA;    ///< STLアスキーファイル
    static const std::string FMT_STL_B;     ///< STLバイナリファイル
    static const std::string FMT_STL_BB;    ///< STLバイナリファイル
    static const std::string FMT_NPT_A;     ///< 長田パッチアスキーファイル
    static const std::string FMT_NPT_B;     ///< 長田パッチバイナリファイル
    static const std::string DEFAULT_FMT;   ///< PolygonIO.cxxで定義している値
};

} //namespace PolylibNS

#endif  // polygon_io_h
