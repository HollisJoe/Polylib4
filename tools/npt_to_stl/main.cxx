
/*
 * STL file to Nagata patch file Convertor
 *
 *
 * Copyright (c) 2015-2016 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 *
 */

////////////////////////////////////////////////////////////////////////////
///
/// 長田パッチファイルを読み込み、STLファイルとして出力する
///     （逐次実行のみ ）
///
////////////////////////////////////////////////////////////////////////////

#include "Polylib.h"
#include <string.h>

using namespace PolylibNS;
using namespace std;

static void Usage ( void )
{
    cerr<<endl;
    cerr<< "Usage: npt_to_stl [options] npt_file" <<endl;
    cerr<< "Options:" <<endl;
    cerr<< "  --help       Display this information" <<endl;
    cerr<< "Arguments:" <<endl;
    cerr<< "  npt_file     nagata patch file name" <<endl;
    cerr<<endl;
}


//----------------------------------------------------
//  メインルーチン
//----------------------------------------------------

int main(int argc, char** argv )
{
    POLYLIB_STAT ret;

    //-------------------------------------------
    //  プログラム引数チェック
    //-------------------------------------------
    if( argc < 2 ) {
        Usage();
        exit(1);
    }

    if( strcmp(argv[1],"--help")==0 ) {
        Usage();
        exit(1);
    }

    //-------------------------------------------
    //  ファイルパス
    //-------------------------------------------

    std::string npt_file_name;   // 入力：nptファイル名
    std::string fmt_in;          // 入力：ファイルフォーマット
    std::string stl_file_name;   // 出力：stlファイル名
    std::string fmt_out;         // 出力：ファイルフォーマット
    
    npt_file_name = argv[1];
    fmt_in = PolygonIO::input_file_format( npt_file_name );

    char* fname = PolylibNS::get_fname_fr_path( npt_file_name );
    stl_file_name = fname;
    stl_file_name += ".stl";
    if( fmt_in == PolygonIO::FMT_NPT_A ) {
        fmt_out = PolygonIO::FMT_STL_A;     // Text
    } else if( fmt_in == PolygonIO::FMT_NPT_B ) {
        fmt_out = PolygonIO::FMT_STL_B;     // Binary
    } else {
        Usage();
        exit(1);
    }

    //-------------------------------------------
    // nptファイル読み込み
    //-------------------------------------------
    std::vector<Triangle*>* npt_list = new std::vector<Triangle*>;

    // nptファイル読み込み
    ret = PolygonIO::load( npt_list, npt_file_name, fmt_in );
    
    if( npt_list->size() == 0 ) {
        PL_ERROSH << "[ERROR] no polygons in npt file. file_name=" << npt_file_name << endl;
    }

    //-------------------------------------------
    // STLファイル保存
    //-------------------------------------------

    //PL_DBGOSH << "save() start" << endl;
    ret = PolygonIO::save( npt_list, stl_file_name, fmt_out );
    //PL_DBGOSH << "save() end" << endl;

    //------------------------------------------
    // 終了化処理
    //-------------------------------------------
    for( int i=0; i< npt_list->size(); i++ ) {
        delete (*npt_list)[i];
    }

    return 0;
}
