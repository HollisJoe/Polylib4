
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
/// 長田パッチの各辺の中点の曲面補間点を追加してポリゴン数を４倍とする
/// STLファイルに出力する
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
    cerr<< "Usage: npt_to_stl4 [options] npt_file" <<endl;
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

    PL_DBGOSH << "#### npt_to_stl4 (start) ####" <<endl;

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
    stl_file_name += "_4.stl";
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

    int num_npt = npt_list->size();
#ifdef DEBUG
    PL_DBGOSH << "num_npt="<<num_npt <<endl;
#endif

    //-------------------------------------------
    // 各辺の中間補間点を加え、ポリゴンを４分割する
    //-------------------------------------------
    std::vector<NptTriangle*>* pNpt = (std::vector<NptTriangle*>*)(npt_list);

    std::vector<Triangle*> tri_list;
    tri_list.reserve( 4*num_npt );
    Vec3<PL_REAL> p12;
    Vec3<PL_REAL> p23;
    Vec3<PL_REAL> p31;
    PL_REAL eta,xi;
    Vec3<PL_REAL> vertexes_new[3];

    Triangle* p_tri;

    for( int i=0; i<num_npt; i++ ) {

        Vec3<PL_REAL>* vertexes = (*npt_list)[i]->get_vertexes();
        //PL_DBGOSH << "NPT i="<<i<<" vertex[0]="<<vertexes[0]<<" vertex[1]="<<vertexes[1]<<" vertex[2]="<<vertexes[2] << endl;

        // 辺１の中点の曲面補間点
        eta=0.5; xi=0.0;
        (*pNpt)[i]->correct( eta, xi, p12 );  
        //PL_DBGOSH << "   p12="<<p12 << endl;

        // 辺２の中点の曲面補間点
        eta=1.0; xi=0.5;
        (*pNpt)[i]->correct( eta, xi, p23 );  
        //PL_DBGOSH << "   p23="<<p23 << endl;

        // 辺３の中点の曲面補間点
        eta=0.5; xi=0.5;
        (*pNpt)[i]->correct( eta, xi, p31 );  
        //PL_DBGOSH << "   p31="<<p31 << endl;

        //ポリゴン１個目
        vertexes_new[0] = vertexes[0];
        vertexes_new[1] = p12;
        vertexes_new[2] = p31;
        p_tri = new Triangle( vertexes_new );
        tri_list.push_back( p_tri );

        //ポリゴン２個目
        vertexes_new[0] = vertexes[1];
        vertexes_new[1] = p23;
        vertexes_new[2] = p12;
        p_tri = new Triangle( vertexes_new );
        tri_list.push_back( p_tri );

        //ポリゴン３個目
        vertexes_new[0] = vertexes[2];
        vertexes_new[1] = p31;
        vertexes_new[2] = p23;
        p_tri = new Triangle( vertexes_new );
        tri_list.push_back( p_tri );

        //ポリゴン４個目
        vertexes_new[0] = p12;
        vertexes_new[1] = p23;
        vertexes_new[2] = p31;
        p_tri = new Triangle( vertexes_new );
        tri_list.push_back( p_tri );
    }


    //-------------------------------------------
    // STLファイル保存
    //-------------------------------------------

    //PL_DBGOSH << "save() start" << endl;
    ret = PolygonIO::save( &tri_list, stl_file_name, fmt_out );
    //PL_DBGOSH << "save() end" << endl;

    //------------------------------------------
    // 終了化処理
    //-------------------------------------------
    for( int i=0; i< npt_list->size(); i++ ) {
        delete (*npt_list)[i];
    }

    for( int i=0; i< tri_list.size(); i++ ) {
        delete tri_list[i];
    }

    PL_DBGOSH << "#### npt_to_stl4 (end) ####" <<endl;

    return 0;
}
