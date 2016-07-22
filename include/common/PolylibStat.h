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

#ifndef polylib_stat_h
#define polylib_stat_h

////////////////////////////////////////////////////////////////////////////
///
/// Polylibで利用するEnumの定義
///
////////////////////////////////////////////////////////////////////////////
typedef enum {
    PLSTAT_OK = 0,                     ///< 処理が成功した。
    PLSTAT_NG = 1,                     ///< 一般的なエラー。
    PLSTAT_INSTANCE_EXISTED     = 2,   ///< Polylibインスタンスがすでに存在している。
    PLSTAT_INSTANCE_NOT_EXIST   = 3,   ///< Polylibインスタンスが存在しない。
    PLSTAT_MPI_ERROR            = 5,   ///< MPI関数がエラーを戻した。
    PLSTAT_ARGUMENT_NULL        = 6,   ///< 引数のメモリ確保が行われていない。
    PLSTAT_MEMORY_NOT_ALLOC     = 7,    ///< メモリ確保に失敗した。
    PLSTAT_LACK_OF_MEMORY       = 8,    ///< メモリ不足
    PLSTAT_CONFIG_ERROR         = 10,  ///< 定義ファイルでエラー発生
    PLSTAT_STL_IO_ERROR         = 11,  ///< STLファイルIOエラー
    PLSTAT_NPT_IO_ERROR         = 12,  ///< 長田パッチファイルIOエラー
    PLSTAT_UNKNOWN_FILE_FORMAT  = 13,  ///< ファイルが.stla、.stlb、.stl、.npta、.nptb、.npt以外。
    PLSTAT_LACK_OF_LOAD_MEMORY  = 14,  ///< ロード処理時のメモリ不足（メモリ削減版使用時）
    PLSTAT_FILE_NOT_SET         = 20,  ///< リーフグループにファイル名が未設定。
    PLSTAT_GROUP_NOT_FOUND      = 21,  ///< グループ名がPolylibに未登録。
    PLSTAT_GROUP_NAME_EMPTY     = 22,  ///< グループ名が空である。
    PLSTAT_GROUP_NAME_DUP       = 23,  ///< グループ名が重複している。
    PLSTAT_POLYGON_NOT_EXIST    = 24,  ///< ポリゴンが存在しない
    PLSTAT_NODE_NOT_FIND        = 26,  ///< KD木生成時に検索点が見つからなかった。
    PLSTAT_ROOT_NODE_NOT_EXIST  = 27,  ///< KD木のルートノードが存在しない。
    PLSTAT_NOT_NPT              = 28,  ///< 長田パッチではない（NptTriangle*へのdynamic cast失敗など）
    PLSTAT_ATR_NOT_EXIST        = 29,  ///< 属性が未設定
// 以下は未使用
//  PLSTAT_GROUP_UNMATCH,       ///< グループ並びがランク0と一致しなかった。
//  PLSTAT_UNkNOWN_ERROR,       ///< 予期せぬエラー。
//  PLSTAT_FILE_NOT_FOUND,      ///< 読み込みファイルが存在しない。
//  PLSTAT_FILE_NOT_OPEN,       ///< ファイルが開けなかった。
//  PLSTAT_FILE_NOT_SCAN,       ///< ファイル読み込みエラー。
//  PLSTAT_PARAMETER_NOT_FOUND, ///< ConfigファイルにParamタグがない。
} POLYLIB_STAT;

//#ifndef C_LANG // C言語版でも本ヘッダを使用しているため#ifndefを追加 2010.11.04
#ifdef __cplusplus

namespace PolylibNS {
////////////////////////////////////////////////////////////////////////////
///
/// PolylibStat文字列出力用クラス
///
////////////////////////////////////////////////////////////////////////////
class PolylibStat2 {
public:
    ///
    /// PolylibStat文字列出力。
    ///
    ///  @param[in] stat    PolylibStat値。
    ///  @return    PolylibStat値を文字列化したもの。
    ///
    static std::string String(
                POLYLIB_STAT    stat
            )
    {
        if (stat == PLSTAT_OK)                          return "PLSTAT_OK";
        else if (stat == PLSTAT_NG)                     return "PLSTAT_NG";
        else if (stat == PLSTAT_INSTANCE_EXISTED)       return "PLSTAT_INSTANCE_EXISTED";
        else if (stat == PLSTAT_INSTANCE_NOT_EXIST)     return "PLSTAT_INSTANCE_NOT_EXIST";
        else if (stat == PLSTAT_MPI_ERROR)              return "PLSTAT_MPI_ERROR";
        else if (stat == PLSTAT_ARGUMENT_NULL)          return "PLSTAT_ARGUMENT_NULL";
        else if (stat == PLSTAT_MEMORY_NOT_ALLOC)       return "PLSTAT_MEMORY_NOT_ALLOC";
        else if (stat == PLSTAT_LACK_OF_MEMORY)         return "PLSTAT_LACK_OF_MEMORY";
        else if (stat == PLSTAT_CONFIG_ERROR)           return "PLSTAT_CONFIG_ERROR";
        else if (stat == PLSTAT_STL_IO_ERROR)           return "PLSTAT_STL_IO_ERROR";
        else if (stat == PLSTAT_NPT_IO_ERROR)           return "PLSTAT_NPT_IO_ERROR";
        else if (stat == PLSTAT_UNKNOWN_FILE_FORMAT)    return "PLSTAT_UNKNOWN_FILE_FORMAT";
        else if (stat == PLSTAT_LACK_OF_LOAD_MEMORY)    return "PLSTAT_LACK_OF_LOAD_MEMORY";
        else if (stat == PLSTAT_FILE_NOT_SET)           return "PLSTAT_FILE_NOT_SET";
        else if (stat == PLSTAT_GROUP_NOT_FOUND)        return "PLSTAT_GROUP_NOT_FOUND";
        else if (stat == PLSTAT_GROUP_NAME_EMPTY)       return "PLSTAT_GROUP_NAME_EMPTY";
        else if (stat == PLSTAT_GROUP_NAME_DUP)         return "PLSTAT_GROUP_NAME_DUP";
        else if (stat == PLSTAT_POLYGON_NOT_EXIST)      return "PLSTAT_POLYGON_NOT_EXIST";
        else if (stat == PLSTAT_NODE_NOT_FIND)          return "PLSTAT_NODE_NOT_FIND";
        else if (stat == PLSTAT_ROOT_NODE_NOT_EXIST)    return "PLSTAT_ROOT_NODE_NOT_EXIST";
        else if (stat == PLSTAT_NOT_NPT)                return "PLSTAT_NO_NPT";
        else if (stat == PLSTAT_ATR_NOT_EXIST)          return "PLSTAT_ATR_NOT_EXIST";
        else                                            return "UNKNOW_STATUS";
    }
};

}
#endif // C_LANG

#endif // polylib_stat_h
