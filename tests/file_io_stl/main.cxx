
/*
 * STL file I/O test
 *
 *
 * Copyright (c) 2015-2016 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 *
 */

////////////////////////////////////////////////////////////////////////////
///
/// STLファイル(Text/Binary)入出力テスト
///
////////////////////////////////////////////////////////////////////////////

#include "Polylib.h"

using namespace PolylibNS;
using namespace std;

#ifdef MPI_PL
struct MyParallelInfo {
  PL_REAL  bpos[3]; //基準座標
  unsigned bbsize[3]; //number of voxel 計算領域
  unsigned gcsize[3]; //number of guidecell voxel
  PL_REAL  dx[3]; //size of voxel
};

//   想定モデル空間範囲
//      x: -50〜50
//      y: -50〜50
//      z: -50〜50

#if 0
// 並列数１：分割なし
static MyParallelInfo myParaInfos[1] = {
  {{-50,-50,-50,}, {10,10,10,}, {1, 1, 1,}, {10,10,10} },
};
#endif
#if 0
// 並列数２：Y方向を 0 の位置で２分割
static MyParallelInfo myParaInfos[2] = {
  {{-50,-50,-50,}, {10, 5,10,}, {1, 1, 1,}, {10,10,10} },
  {{-50,  0,-50,}, {10, 5,10,}, {1, 1, 1,}, {10,10,10} } 
};
#endif
#if 1
// 並列数４：X,Y方向を 0 の位置で２分割
static MyParallelInfo myParaInfos[4] = {
  {{-50,-50,-50,}, { 5, 5,10,}, {1, 1, 1,}, {10,10,10} },
  {{  0,-50,-50,}, { 5, 5,10,}, {1, 1, 1,}, {10,10,10} },
  {{-50,  0,-50,}, { 5, 5,10,}, {1, 1, 1,}, {10,10,10} },
  {{  0,  0,-50,}, { 5, 5,10,}, {1, 1, 1,}, {10,10,10} } 
};
#endif

#endif

//----------------------------------------------------
//  メインルーチン
//----------------------------------------------------

int main(int argc, char** argv )
{
    POLYLIB_STAT ret;
    int     iret;
    int     num_rank;
    int     myrank;
    std::string config_file_name = "polylib_config.tp";   // 入力：初期化ファイル名

    //-------------------------------------------
    //  初期化
    //-------------------------------------------

    // MPI初期化
#ifdef MPI_PL
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &num_rank );
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    PL_DBGOSH << "---- num_rank="<<num_rank<<" myrank="<<myrank <<endl;
#else
    num_rank = 1;
    myrank   = 0;
#endif

    // Polylib初期化
    Polylib* p_polylib = Polylib::get_instance();

    // 並列計算関連情報の設定と初期化
#ifdef MPI_PL
    PL_DBGOSH << "(rk:"<<myrank<<") "<<"init_parallel_info() start" << endl;
    ret = p_polylib->init_parallel_info(
                MPI_COMM_WORLD,
                myParaInfos[myrank].bpos,
                myParaInfos[myrank].bbsize,
                myParaInfos[myrank].gcsize,
                myParaInfos[myrank].dx
            );
    PL_DBGOSH << "init_parallel_info() end" << endl;
    if( ret != PLSTAT_OK ) {
        PL_ERROSH <<"[ERROR] init_parallel_info() ret="<<ret <<endl;
        exit(1);
    }
#endif

    //-------------------------------------------
    //  ロード
    //-------------------------------------------

    // 初期化ファイルを指定してデータロード
    PL_DBGOSH << "load() start" <<endl;
    ret = p_polylib->load( config_file_name );
    PL_DBGOSH << "load() end   ret="<<ret <<endl;
    if( ret != PLSTAT_OK ) {
        PL_ERROSH <<"[ERROR] load() ret="<<ret <<endl;
        exit(1);
    }

    // リーフポリゴングループ取得
    std::vector<PolygonGroup *>* pg_list = p_polylib->get_leaf_groups();
    if( pg_list == NULL ) {
        PL_ERROSH <<"[ERROR] get_leaf_groups()" <<endl;
        exit(1);
    }
    PL_DBGOSH << "pg_list->size()="<<pg_list->size() <<endl;

    //-------------------------
    // load処理の確認
    //-------------------------
    for(int i=0; i<pg_list->size(); i++ ) {
        std::string pg_name = (*pg_list)[i]->get_name();
        std::string pg_path = (*pg_list)[i]->acq_fullpath();
        bool movable = (*pg_list)[i]->get_movable();
        std::vector<Triangle* > *tri_list =  (*pg_list)[i]->get_triangles();
        int num_tri = 0;
        if( tri_list == NULL || tri_list->size()==0 )  {
            PL_ERROSH <<"[ERROR] Can't get Polygon data" <<endl;
            exit(1);
        }
        if( tri_list != NULL )  {
            num_tri = tri_list->size();
        }
        PL_DBGOSH << "polygon group i="<<i<<" name="<<pg_name<<" pg_path="<<pg_path <<endl;
        PL_DBGOSH << "   movable="<<movable<<" num_tri="<<num_tri <<endl;

    }

    //-------------------------------------------
    //  セーブ（バイナリファイルで保存）
    //-------------------------------------------
    std::string  config_name_out;
    std::string  fmt_out;

    fmt_out = PolygonIO::FMT_STL_B;
    std::string  extend  = "out_bin";
    PL_DBGOSH << "save() binary start" <<endl;
    ret = p_polylib->save( config_name_out, fmt_out, extend );
    PL_DBGOSH << "save() binary end   config_name_out="<<config_name_out <<endl;
    if( ret != PLSTAT_OK ) {
        PL_ERROSH <<"[ERROR] save() binary ret="<<ret <<endl;
        exit(1);
    }

    //-------------------------------------------
    //  セーブ（Textファイルで保存）
    //-------------------------------------------

    fmt_out = PolygonIO::FMT_STL_A;
    PL_DBGOSH << "save() text start" <<endl;
    ret = p_polylib->save( config_name_out, fmt_out );
    PL_DBGOSH << "save() text end   config_name_out="<<config_name_out <<endl;
    if( ret != PLSTAT_OK ) {
        PL_ERROSH <<"[ERROR] save() text ret="<<ret <<endl;
        exit(1);
    }

    //-------------------------------------------
    //  妥当性の検証
    //-------------------------------------------

    //  入力したテキストのSTLのポリゴン取得
    std::string pg_sphere_1_path = "sphere_1";
    PolygonGroup* pg_sphere_1 = p_polylib->get_group( pg_sphere_1_path );
    if( pg_sphere_1 == NULL )  {
        PL_ERROSH <<pg_sphere_1_path<<": [ERROR] Can't get PolygonGroup pointer" <<endl;
        exit(1);
    }
#ifdef MPI_PL
    // MPI版の場合
    std::vector<Triangle* > tri_list_sphere_1;
    pg_sphere_1->gather_polygons( tri_list_sphere_1 );
#else
    // 逐次版の場合
    std::vector<Triangle* > *ptri_list_sphere_1 =  pg_sphere_1->get_triangles();
    std::vector<Triangle* > & tri_list_sphere_1 =  *ptri_list_sphere_1;
#endif

    //  入力したバイナリのSTLのポリゴン取得
    std::string pg_sphere_2_path = "sphere_2";
    PolygonGroup* pg_sphere_2 = p_polylib->get_group( pg_sphere_2_path );
    if( pg_sphere_2 == NULL )  {
        PL_ERROSH <<pg_sphere_2_path<<": [ERROR] Can't get PolygonGroup pointer" <<endl;
        exit(1);
    }
#ifdef MPI_PL
    // MPI版の場合
    std::vector<Triangle* > tri_list_sphere_2;
    pg_sphere_1->gather_polygons( tri_list_sphere_2 );
#else
    // 逐次版の場合
    std::vector<Triangle* > *ptri_list_sphere_2 =  pg_sphere_2->get_triangles();
    std::vector<Triangle* > & tri_list_sphere_2 =  *ptri_list_sphere_2;
#endif

    //-------------------------------------------------
    //  以下はランク０のみで値検証
    //-------------------------------------------------

#ifdef MPI_PL
    if( myrank == 0 )  {
#endif
        bool err_flg=false;

        //  バイナリで書きだしたファイルを再度読み込み
        std::string file_name_bin = "sphere_1_out_bin.stlb";
        fmt_out = PolygonIO::FMT_STL_B;
        std::vector<Triangle* > tri_list_sphere_bin;
        PolygonIO::load( &tri_list_sphere_bin, file_name_bin, fmt_out );

        //---- ポリゴン数の確認　-------
        int num_tri_sphere_1   = tri_list_sphere_1.size();
        int num_tri_sphere_2   = tri_list_sphere_2.size();
        int num_tri_sphere_bin = tri_list_sphere_bin.size();
        if( num_tri_sphere_1 != num_tri_sphere_2   ) err_flg = true;
        if( num_tri_sphere_1 != num_tri_sphere_bin ) err_flg = true;
        if( err_flg ) {
            PL_ERROSH << "#### ERROR number of polygons ####" << endl;
            PL_DBGOSH << "#### ERROR number of polygons ####" << endl;
            PL_DBGOSH << "  num_tri_sphere_1   = "<<num_tri_sphere_1   <<endl;
            PL_DBGOSH << "  num_tri_sphere_2   = "<<num_tri_sphere_2   <<endl;
            PL_DBGOSH << "  num_tri_sphere_bin = "<<num_tri_sphere_bin <<endl;
            exit(1);
        }

        //---- ポリゴン座標の確認　-------
        PL_REAL eps = 0.001;
        for(int i=0; i<num_tri_sphere_1; i++ )  {
            Vec3<PL_REAL>* v1 = tri_list_sphere_1  [i]->get_vertexes();
            Vec3<PL_REAL>* v2 = tri_list_sphere_2  [i]->get_vertexes();
            Vec3<PL_REAL>* v3 = tri_list_sphere_bin[i]->get_vertexes();

            for(int j=0; j<3; j++ )  {
                if( fabs(v1[j].x-v2[j].x) > eps ) err_flg=true; 
                if( fabs(v1[j].y-v2[j].y) > eps ) err_flg=true; 
                if( fabs(v1[j].z-v2[j].z) > eps ) err_flg=true; 
                if( fabs(v1[j].x-v3[j].x) > eps ) err_flg=true; 
                if( fabs(v1[j].y-v3[j].y) > eps ) err_flg=true; 
                if( fabs(v1[j].z-v3[j].z) > eps ) err_flg=true; 
            }

            if( err_flg )  {
                PL_ERROSH << "#### ERROR coordinates i="<<i <<endl;
                PL_DBGOSH << "#### ERROR coordinates i="<<i <<endl;
                std::cout <<"    sphere_1 vertexes:" <<endl;
                std::cout <<"       "<<v1[0].x<<" "<<v1[0].y<<" "<<v1[0].z <<endl; 
                std::cout <<"       "<<v1[1].x<<" "<<v1[1].y<<" "<<v1[1].z <<endl; 
                std::cout <<"       "<<v1[2].x<<" "<<v1[2].y<<" "<<v1[3].z <<endl; 
                std::cout <<"    sphere_2 vertexes:" <<endl;
                std::cout <<"       "<<v2[0].x<<" "<<v2[0].y<<" "<<v2[0].z <<endl; 
                std::cout <<"       "<<v2[1].x<<" "<<v2[1].y<<" "<<v2[1].z <<endl; 
                std::cout <<"       "<<v2[2].x<<" "<<v2[2].y<<" "<<v2[3].z <<endl; 
                std::cout <<"    sphere_bin vertexes:" <<endl;
                std::cout <<"       "<<v3[0].x<<" "<<v3[0].y<<" "<<v3[0].z <<endl; 
                std::cout <<"       "<<v3[1].x<<" "<<v3[1].y<<" "<<v3[1].z <<endl; 
                std::cout <<"       "<<v3[2].x<<" "<<v3[2].y<<" "<<v3[3].z <<endl; 
                //break;
                exit(1);
            }
        }

        for(int i=0; i<tri_list_sphere_bin.size(); i++ ) {
            delete tri_list_sphere_bin[i];
        }

#ifdef MPI_PL
    }  //! if( myrank == 0 ) {
#endif

#ifdef MPI_PL
    for(int i=0; i<tri_list_sphere_1.size(); i++ ) {
        delete tri_list_sphere_1[i];
    }
    for(int i=0; i<tri_list_sphere_2.size(); i++ ) {
        delete tri_list_sphere_2[i];
    }
#endif

    //-------------------------------------------
    //  終了化
    //-------------------------------------------
    delete pg_list;

    // MPI終了化
#ifdef MPI_PL
    PL_DBGOSH << "MPI_Finalize() start" << endl;
    MPI_Finalize();
    PL_DBGOSH << "MPI_Finalize() end" << endl;
#endif

    if( myrank == 0 )  {
        std::cerr<<"------------------------------------------" <<endl;
        std::cerr<<"   PASS :  file_io_stl (Normal End)" <<endl;
        std::cerr<<"------------------------------------------" <<endl;
    }

    return 0;
}
