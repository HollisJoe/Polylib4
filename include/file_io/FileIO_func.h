/*
 * Polylib - Polygon Management Library
 *
 * Copyright (c) 2010-2011 VCAD System Research Program, RIKEN.
 * All rights reserved.
 *
 * Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 *
 */

#ifndef file_io_h
#define file_io_h

#include <vector>
#include "common/PolylibCommon.h"

namespace PolylibNS {


//**************************************************************
//   STLファイル用
//**************************************************************

///
/// ASCIIモードのSTLファイルを読み込み、tri_listに三角形ポリゴン情報を設定する。
///
///  @param[in,out] tri_list    三角形ポリゴンリストの領域。
///                                 出力は追加される
///  @param[in]     fname       STLファイル名。
///  //@param[in,out] total       ポリゴンIDの通番
///  @param[out]    num_tri     STLファイル内のポリゴン数
///  @return    POLYLIB_STATで定義される値が返る。
///  @attention
///      ポリゴンIDはシステム内で自動で採番される
///      num_triはtri_list全体の個数でないことに注意
///
POLYLIB_STAT stl_a_load(
    std::vector<Triangle*>  *tri_list, 
    const std::string&      fname,
    int*                     num_tri,
    PL_REAL                 scale=1.0
);

///
/// ASCIIモードのSTLファイルからポリゴン指定個数分読み込み
///
///  @param[in]  ifs          入力ファイルストリーム
///  @param[in,out] tri_list     三角形ポリゴンのリスト(出力内容)
///  @param[in]  num_read     読み込み指定数
///                             EOFに達すれば途中まで読み込まれる
///                             -1 の時、読み込み数制限なし（eofまで読む）
///  @param[out] num_tri      実際に読み込んだ数
///  @param[out] eof          ファイル終了フラグ(end of file)
///  @param[in]  scale       スケール
///  @return    POLYLIB_STATで定義される値が返る。
///  @attention
///      ポリゴンIDはシステム内で自動で採番される
///
POLYLIB_STAT stl_a_load_read(
    ifstream&                ifs,
    std::vector<Triangle*>&  tri_list, 
    int                      num_read,
    int&                     num_tri,
    bool&                    eof,
    PL_REAL                 scale=1.0
);


///
/// 三角形ポリゴン情報をASCIIモードでSTLファイルに書き出す。
///
///  @param[in] tri_list    三角形ポリゴン情報。
///  @param[in] fname       STLファイル名。
///  @return    POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT stl_a_save(
    std::vector<Triangle*>  *tri_list, 
    const std::string&      fname
);

///
/// バイナリモードのSTLファイルを読み込み、tri_listに三角形ポリゴン情報を設定
/// する。
///
///  @param[in,out] tri_list    三角形ポリゴンリストの領域。
///                                 出力は追加される
///  @param[in]     fname       ファイル名。
///  //@param[in,out] total       ポリゴンIDの通番
///  @param[out]    num_tri     STLファイル内のポリゴン数
///  @return    POLYLIB_STATで定義される値が返る。
///  @attention
///      ポリゴンIDはシステム内で自動で採番される
///      num_triはtri_list全体の個数でないことに注意
///
POLYLIB_STAT stl_b_load(
    std::vector<Triangle*>  *tri_list, 
    const std::string&      fname,
    int                         *num_tri,
    PL_REAL                 scale=1.0
);

///
/// バイナリモードのSTLファイルのヘッダ部を読み飛ばす
///
///  @param[in]  ifs          入力ファイルストリーム
///  @return    POLYLIB_STATで定義される値が返る。
///  @attention
///      ポリゴンを分割して読みだす時はstl_b_load_read()にて
///      順次呼び出して行くが先頭のヘッダ部は
///      stl_b_load_read()では読み出せないため、先に読んでおく
//       
///
POLYLIB_STAT stl_b_load_read_head(
    ifstream&                ifs
);

///
/// バイナリモードのSTLファイルからポリゴン指定個数分読み込み
///
///  @param[in]  ifs          入力ファイルストリーム
///  @param[in,out] tri_list     三角形ポリゴンのリスト(出力内容)
///  @param[in]  num_read     読み込み指定数
///                             EOFに達すれば途中まで読み込まれる
///                             -1 の時、読み込み数制限なし（eofまで読む）
///  @param[out] num_tri      実際に読み込んだ数
///  @param[out] eof          ファイル終了フラグ(end of file)
///  @param[in]  scale       スケール
///  @return    POLYLIB_STATで定義される値が返る。
///  @attention
///      ポリゴンIDはシステム内で自動で採番される
///      ヘッダ部は、先にstl_b_load_read_head()で読みだしておくこと
///
POLYLIB_STAT stl_b_load_read(
    ifstream&                ifs,
    std::vector<Triangle*>&  tri_list, 
    int                      num_read,
    int&                     num_tri,
    bool&                    eof,
    PL_REAL                 scale=1.0
);

///
/// 三角形ポリゴン情報をバイナリモードでSTLファイルに書き出す。
///
///  @param[in] tri_list    三角形ポリゴン情報。
///  @param[in] fname       STLファイル名。
///  @return    POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT stl_b_save(
    std::vector<Triangle*>  *tri_list, 
    const std::string&      fname
);

///
/// STLファイルを読み込みバイナリかアスキーかを判定する。
///
///  @param[in] STLファイルのフルパス名。
///  @return    true:アスキー形式 / false:バイナリ形式。
/// 
bool is_stl_a(
    const std::string&     path
);


//**************************************************************
//   長田パッチファイル用
//**************************************************************

///
/// ASCIIモードの長田パッチファイルを読み込み、tri_listに三角形ポリゴン情報を設定する。
///
///  @param[in,out] tri_list    三角形ポリゴンリストの領域
///                                 出力は追加される
///  @param[in]     fname       NPTファイル名
///  //@param[in,out] total       ポリゴンIDの通番
///  @param[out]    num_tri     STLファイル内のポリゴン数
///  @return    POLYLIB_STATで定義される値が返る
///  @attention
///      ポリゴンIDはシステム内で自動で採番される
///      num_triはtri_list全体の個数でないことに注意
///
POLYLIB_STAT npt_a_load(
    //std::vector<NptTriangle*>   *tri_list, 
    std::vector<Triangle*>   *tri_list, 
    const std::string&          fname,
    int                         *num_tri,
    PL_REAL                     scale=1.0
);

///
/// ASCIIモードの長田パッチファイルのヘッダ部を読み飛ばす
///
///  @param[in]  ifs          入力ファイルストリーム
///  @return    POLYLIB_STATで定義される値が返る。
///  @attention
///      ポリゴンを分割して読みだす時はnpt_a_load_read()にて
///      順次呼び出して行くが先頭のヘッダ部は
///      npt_a_load_read()では読み出せないため、先に読んでおく
//       
///
POLYLIB_STAT npt_a_load_read_head(
    ifstream&                ifs
);

///
/// ASCIIモードの長田パッチファイルからポリゴン指定個数分読み込み
///
///  @param[in]  ifs          入力ファイルストリーム
///  @param[in,out] tri_list     三角形ポリゴンのリスト(出力内容)
///  @param[in]  num_read     読み込み指定数
///                             EOFに達すれば途中まで読み込まれる
///                             -1 の時、読み込み数制限なし（eofまで読む）
///  @param[out] num_tri      実際に読み込んだ数
///  @param[out] eof          ファイル終了フラグ(end of file)
///  @param[in]  scale       スケール
///  @return    POLYLIB_STATで定義される値が返る。
///  @attention
///      ポリゴンIDはシステム内で自動で採番される
///
POLYLIB_STAT npt_a_load_read(
    ifstream&                ifs,
    std::vector<Triangle*>&  tri_list, 
    int                      num_read,
    int&                     num_tri,
    bool&                    eof,
    PL_REAL                 scale=1.0
);

///
/// 三角形ポリゴン情報をASCIIモードで長田パッチファイルに書き出す。
///
///  @param[in] tri_list    三角形ポリゴン情報
///  @param[in] fname       NPTファイル名
///  @return    POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT npt_a_save(
    std::vector<NptTriangle*>   *tri_list, 
    const std::string&          fname
);

///
/// バイナリモードの長田パッチファイルを読み込み、tri_listに三角形ポリゴン情報を設定
/// する。
///
///  @param[in,out] tri_list    三角形ポリゴンリストの領域。
///                                 出力は追加される
///  @param[in]     fname       ファイル名。
///  @param[out]    num_tri     STLファイル内のポリゴン数
///  @return    POLYLIB_STATで定義される値が返る
///  @attention
///      ポリゴンIDはシステム内で自動で採番される
///      num_triはtri_list全体の個数でないことに注意
///
POLYLIB_STAT npt_b_load(
    //std::vector<NptTriangle*>   *tri_list, 
    std::vector<Triangle*>   *tri_list, 
    const std::string&          fname,
    int                         *num_tri,
    PL_REAL                     scale=1.0
);

///
/// バイナリモードの長田パッチファイルのヘッダ部を読み飛ばす
///
///  @param[in]  ifs          入力ファイルストリーム
///  @return    POLYLIB_STATで定義される値が返る。
///  @attention
///      ポリゴンを分割して読みだす時はnpt_b_load_read()にて
///      順次呼び出して行くが先頭のヘッダ部は
///      npt_b_load_read()では読み出せないため、先に読んでおく
//       
///
POLYLIB_STAT npt_b_load_read_head(
    ifstream&                ifs
);

///
/// バイナリモードの長田パッチファイルからポリゴン指定個数分読み込み
///
///  @param[in]  ifs          入力ファイルストリーム
///  @param[in,out] tri_list     三角形ポリゴンのリスト(出力内容)
///  @param[in]  num_read     読み込み指定数
///                             EOFに達すれば途中まで読み込まれる
///                             -1 の時、読み込み数制限なし（eofまで読む）
///  @param[out] num_tri      実際に読み込んだ数
///  @param[out] eof          ファイル終了フラグ(end of file)
///  @param[in]  scale       スケール
///  @return    POLYLIB_STATで定義される値が返る。
///  @attention
///      ポリゴンIDはシステム内で自動で採番される
///
POLYLIB_STAT npt_b_load_read(
    ifstream&                ifs,
    std::vector<Triangle*>&  tri_list, 
    int                      num_read,
    int&                     num_tri,
    bool&                    eof,
    PL_REAL                 scale=1.0
);

///
/// 三角形ポリゴン情報をバイナリモードで長田パッチファイルに書き出す。
///
///  @param[in] tri_list    三角形ポリゴン情報。
///  @param[in] fname       NPTファイル名。
///  @return    POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT npt_b_save(
    std::vector<NptTriangle*>   *tri_list, 
    const std::string&          fname
);

///
/// 長田パッチファイルを読み込みバイナリかアスキーかを判定する。
///
///  @param[in] NPTファイルのフルパス名。
///  @return    true:アスキー形式 / false:バイナリ形式。
/// 
bool is_npt_a(
    const std::string&     path
);

//**************************************************************
//   共通ＩＯ用
//**************************************************************

///
/// ファイルパスから名称(拡張子を除いた部分)を取得する
///
///  @param[in] ファイルパス
///  @return    拡張子を除いた名称
///  @attention 戻り値のchar *は解放不要
///             内部でstaticで持っているため多重処理NG
///
char *get_fname_fr_path(
    const std::string&     path
);

///
/// ファイルパスから拡張子のみを取得する
///
///  @param[in] ファイルパス
///  @return    拡張子
///  @attention 戻り値のchar *は解放不要
///             内部でstaticで持っているため多重処理NG
///
char *get_ext_fr_path(
    const std::string&     path
);

} //namespace PolylibNS

#endif  // file_io_h

