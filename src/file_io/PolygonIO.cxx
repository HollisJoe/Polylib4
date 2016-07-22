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

//  <string> for windows, <string.h> for linux
#include <string>
#include <string.h>

#include "file_io/PolygonIO.h"
#include "file_io/FileIO_func.h"

namespace PolylibNS {

using namespace std;

const string PolygonIO::FMT_STL_A  = "stl_a";
const string PolygonIO::FMT_STL_AA = "stl_aa";
const string PolygonIO::FMT_STL_B  = "stl_b";
const string PolygonIO::FMT_STL_BB = "stl_bb";
const string PolygonIO::FMT_NPT_A  = "npt_a";
const string PolygonIO::FMT_NPT_B  = "npt_b";
const string PolygonIO::DEFAULT_FMT = PolygonIO::FMT_STL_B;

/************************************************************************
 *
 * PolygonIOクラス
 *
 ***********************************************************************/
// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolygonIO::load(
    vector<Triangle*>   *tri_list, 
    const map<string, string>   &fmap,
    PL_REAL scale
) {
    map<string, string>::const_iterator it;
    //int                                   total;
    POLYLIB_STAT                        ret = PLSTAT_OK;

    if (tri_list == NULL) {
        PL_ERROSH << "[ERROR]PolygonIO::load():tri_list is NULL." << endl;
        return PLSTAT_NG;
    }

    //total = 0;    // 通算番号に初期値をセット
    for (it = fmap.begin(); it != fmap.end(); it++) {
        string fname    = it->first;
        string fmt      = it->second;

        // １ファイル読み込み
        ret = load( tri_list, fname, fmt, scale );

        // １ファイルでも読み込みに失敗したら戻る
        if (ret != PLSTAT_OK)       return ret;
    }

    return ret;
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolygonIO::load(
    vector<Triangle*>   *tri_list, 
    const std::string&  fname, 
    const std::string&  fmt,
    PL_REAL scale
)
{
    int                                 num_tri;
    POLYLIB_STAT                        ret = PLSTAT_OK;

    if (tri_list == NULL) {
        PL_ERROSH << "[ERROR]PolygonIO::load():tri_list is NULL." << endl;
        return PLSTAT_NG;
    }

    if (fmt == FMT_STL_A || fmt == FMT_STL_AA) {
        ret = stl_a_load(tri_list, fname, &num_tri, scale);

    } else if (fmt == FMT_STL_B || fmt == FMT_STL_BB) {
        ret = stl_b_load(tri_list, fname, &num_tri, scale);


    } else if (fmt == FMT_NPT_A) {
        //vector<NptTriangle*>* npt_list = static_cast< vector<NptTriangle*>* >(tri_list);  
        //    static_castを使うとコンパイルエラーになる
        //vector<NptTriangle*>*npt_list = ( vector<NptTriangle*>* )tri_list;  
        //ret = npt_a_load(npt_list, fname, &num_tri, scale);
        ret = npt_a_load(tri_list, fname, &num_tri, scale);


    } else if (fmt == FMT_NPT_B) {
        //vector<NptTriangle*>* npt_list = static_cast< vector<NptTriangle*>* >(tri_list);  
        //     static_castを使うとコンパイルエラーになる
        //vector<NptTriangle*>*npt_list = ( vector<NptTriangle*>* )tri_list;  
        //ret = npt_b_load(npt_list, fname, &num_tri, scale);
        ret = npt_b_load(tri_list, fname, &num_tri, scale);

    } else {
        PL_ERROSH << "[ERROR]:PolygonIO::load():Unknown io format." << endl;
        ret = PLSTAT_UNKNOWN_FILE_FORMAT;
    }

    return ret;
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolygonIO::load_file_open(
        ifstream&           ifs, 
        const std::string&  fname, 
        const std::string&  fmt
    )
{
    // Textファイル　STL Open
    if( fmt == FMT_STL_A    || 
        fmt == FMT_STL_AA      )
    {
        ifs.open( fname.c_str() );
        if (ifs.fail()) {
            PL_ERROSH << "[ERROR]PolygonIO::load_file_open():Can't open " << fname << endl;
            return PLSTAT_STL_IO_ERROR;
        }
    }
    // Textファイル　NPT Open
    else if( fmt == FMT_NPT_A )
    {
        ifs.open( fname.c_str() );
        if (ifs.fail()) {
            PL_ERROSH << "[ERROR]PolygonIO::load_file_open():Can't open " << fname << endl;
            return PLSTAT_STL_IO_ERROR;
        }
        // ヘッダ部読み飛ばし
        npt_a_load_read_head( ifs );
    }
    // Binaryファイル　STL Open
    else if( fmt == FMT_STL_B   ||
             fmt == FMT_STL_BB      )
    {
        ifs.open( fname.c_str(), ios::in | ios::binary );
        if (ifs.fail()) {
            PL_ERROSH << "[ERROR]PolygonIO::load_file_open():Can't open " << fname << endl;
            return PLSTAT_STL_IO_ERROR;
        }
        // ヘッダ部読み飛ばし
        stl_b_load_read_head( ifs );
    }
    // Binaryファイル　NPT Open
    else if( fmt == FMT_NPT_B   )
    {
        ifs.open( fname.c_str(), ios::in | ios::binary );
        if (ifs.fail()) {
            PL_ERROSH << "[ERROR]PolygonIO::load_file_open():Can't open " << fname << endl;
            return PLSTAT_STL_IO_ERROR;
        }
        // ヘッダ部読み飛ばし
        npt_b_load_read_head( ifs );
    }
    // Error
    else
    {
        PL_ERROSH << "[ERROR]:PolygonIO::load_file_open():Unknown io format." << endl;
        return PLSTAT_UNKNOWN_FILE_FORMAT;
    }

    return PLSTAT_OK;
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolygonIO::load_file_close( ifstream&  ifs )
{
    ifs.close();
	return PLSTAT_OK;
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolygonIO::load_file_read(
        ifstream&           ifs, 
        const std::string&  fmt,
        vector<Triangle*>&  tri_list, 
        int                 num_read,
        int&                num_tri,
        bool&               eof,
        PL_REAL             scale
    )
{
    POLYLIB_STAT ret;

    if (fmt == FMT_STL_A || fmt == FMT_STL_AA) {
        ret = stl_a_load_read( ifs, tri_list, num_read, num_tri, eof, scale );

    } else if (fmt == FMT_STL_B || fmt == FMT_STL_BB) {
        ret = stl_b_load_read( ifs, tri_list, num_read, num_tri, eof, scale );

    } else if (fmt == FMT_NPT_A) {
        ret = npt_a_load_read( ifs, tri_list, num_read, num_tri, eof, scale );

    } else if (fmt == FMT_NPT_B) {
        ret = npt_b_load_read( ifs, tri_list, num_read, num_tri, eof, scale );

    } else {
        PL_ERROSH << "[ERROR]:PolygonIO::load_file_read():Unknown io format." << endl;
        ret = PLSTAT_UNKNOWN_FILE_FORMAT;
    }

    return ret;
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolygonIO::save(
    vector<Triangle*>   *tri_list, 
    const string&       fname, 
    const string&       fmt
) {
#ifdef DEBUG
    PL_DBGOSH<<__FUNCTION__ << " PolygonIO::save() in."<<endl;
#endif
    if (tri_list == NULL) {
        PL_ERROSH << "[ERROR]:PolygonIO::save():tri_list is NULL." << endl;
        return PLSTAT_NG;
    }

    if (fmt == FMT_STL_A || fmt == FMT_STL_AA) {
        return stl_a_save(tri_list, fname);

    } else if (fmt == FMT_STL_B || fmt == FMT_STL_BB) {
        return stl_b_save(tri_list, fname);

    } else if (fmt == FMT_NPT_A) {
        //vector<NptTriangle*>* npt_list = static_cast< vector<NptTriangle*>* >(tri_list);  
        //    static_castを使うとコンパイルエラーになる
        vector<NptTriangle*>*npt_list = ( vector<NptTriangle*>* )tri_list;  

        return npt_a_save(npt_list, fname);

    } else if (fmt == FMT_NPT_B) {
        //vector<NptTriangle*>* npt_list = static_cast< vector<NptTriangle*>* >(tri_list);  
        //    static_castを使うとコンパイルエラーになる
        vector<NptTriangle*>*npt_list = ( vector<NptTriangle*>* )tri_list;  

        return npt_b_save(npt_list, fname);

    } else {
        return PLSTAT_UNKNOWN_FILE_FORMAT;
    }
}

// public /////////////////////////////////////////////////////////////////////
string PolygonIO::input_file_format(
    const string &filename
)
{
    //PL_DBGOSH <<" input_file_format(): filename="<<filename <<endl;

    //書式の決定
    char    *ext = get_ext_fr_path(filename);
    if (!strcmp(ext, "stla") || !strcmp(ext, "STLA")) {
         return FMT_STL_A;

    } else if (!strcmp(ext, "stlb") || !strcmp(ext, "STLB")) {
         return FMT_STL_B;

    } else if (!strcmp(ext, "stl") || !strcmp(ext, "STL")) {
        //読み込んで書式を判定する
        if(is_stl_a(filename) == true) {
            return FMT_STL_A;
        } else {
            return FMT_STL_B;
        }
    } else if (!strcmp(ext, "npta") || !strcmp(ext, "NPTA")) {
         return FMT_STL_A;

    } else if (!strcmp(ext, "nptb") || !strcmp(ext, "NPTB")) {
         return FMT_STL_B;

    } else if (!strcmp(ext, "npt") || !strcmp(ext, "NPT")) {
        //読み込んで書式を判定する
        if(is_npt_a(filename) == true) {
            //PL_DBGOSH <<" input_file_format(): format=="<<FMT_NPT_A <<endl;
            return FMT_NPT_A;
        } else {
            //PL_DBGOSH <<" input_file_format(): format=="<<FMT_NPT_B <<endl;
            return FMT_NPT_B;
        }
    }

    return "";
}


// public /////////////////////////////////////////////////////////////////////
//  各フォーマットの拡張子取得
std::string PolygonIO::get_extension_format(
    const std::string&  fmt
)
{
    string      extension;

    if (fmt == FMT_STL_A || fmt == FMT_STL_AA) {
        extension = "stla";

    } else if (fmt == FMT_STL_B || fmt == FMT_STL_BB) {
        extension = "stlb";

    } else if (fmt == FMT_NPT_A) {
        extension = "npta";


    } else if (fmt == FMT_NPT_B) {
        extension = "nptb";

    } else {
        PL_ERROSH << "[ERROR]:PolygonIO::get_extension_format():Unknown io format." << endl;
    }

    return extension;
}

// public /////////////////////////////////////////////////////////////////////
//  各フォーマットのポリゴンタイプ取得
int PolygonIO::get_polygon_type(
    const std::string&  fmt
)
{
    if (fmt == FMT_STL_A || fmt == FMT_STL_AA) {
        return PL_TYPE_TRIANGLE;
    } else if (fmt == FMT_STL_B || fmt == FMT_STL_BB) {
        return PL_TYPE_TRIANGLE;
    } else if (fmt == FMT_NPT_A) {
        return PL_TYPE_NPT;
    } else if (fmt == FMT_NPT_B) {
        return PL_TYPE_NPT;
    } else {
        return PL_TYPE_UNKNOWN;
    }
}

} //namespace PolylibNS
