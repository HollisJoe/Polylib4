
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
/// ロード（メモリ削減版）テスト
///
////////////////////////////////////////////////////////////////////////////

#include "Polylib.h"

using namespace PolylibNS;
using namespace std;

#ifdef MPI_PL
typedef struct {
  PL_REAL  bpos[3]; //基準座標
  unsigned bbsize[3]; //number of voxel 計算領域
  unsigned gcsize[3]; //number of guidecell voxel
  PL_REAL  dx[3]; //size of voxel
} MyParallelInfo ;

// モデル空間
//   想定空間 x: -50〜150
//            y: -50〜150
//            z:   0〜150
//
//   plane    x:   0〜100
//            y:   0〜100
//            z:   0〜  0
//   spehre   x: -50〜 50
//            y: -50〜 50
//            z:  50〜150
//

#if 0
// 並列数１：分割なし
static MyParallelInfo myParaInfos[1] = {
  {{-50,-50,-50,}, {20,20,20,}, {2, 2, 2,}, {10,10,10} },
};
#endif
#if 0
// 並列数２：Y方向を +50 の位置で２分割
static MyParallelInfo myParaInfos[2] = {
  {{-50,-50,-50,}, {20,10,20,}, {2, 2, 2,}, {10,10,10} },
  {{-50, 50,-50,}, {20,10,20,}, {2, 2, 2,}, {10,10,10} } 
};
#endif
#if 1
// 並列数４：
//    X方向を +20 の位置で２分割
//    Y方向を 0 の位置で２分割
static MyParallelInfo myParaInfos[4] = {
  {{-50,-50,-50,}, { 7, 5,20,}, {1, 1, 1,}, {10,10,10} },
  {{ 20,-50,-50,}, {13, 5,20,}, {1, 1, 1,}, {10,10,10} },
  {{-50,  0,-50,}, { 7,15,20,}, {1, 1, 1,}, {10,10,10} },
  {{ 20,  0,-50,}, {13,15,20,}, {1, 1, 1,}, {10,10,10} }
};
#endif

#endif

const int num_global_polygon = 192;  // sphereのポリゴン数


//----------------------------------------------------
//  メインルーチン
//----------------------------------------------------

int main(int argc, char** argv )
{
    POLYLIB_STAT ret;
    int     iret;
    int     num_rank;
    int     myrank;
    std::string config_file_name = "polylib_config.tp"; // 入力：初期化ファイル名

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
    num_rank =1;
    myrank   =0;
#endif

    // Polylib初期化
    Polylib* p_polylib = Polylib::get_instance();

    // Polylibが利用する最大メモリサイズ(MB)を設定する
    //   メモリ削減は並列版のみなので逐次では指定しない
#ifdef MPI_PL
    //int max_size_mb = 1000; // F2 30mm 37万ポリゴン 分割ロードされない
    //int max_size_mb = 40; // F2 30mm 37万ポリゴン  分割ロードされる
    //int max_size_mb = 10; // F2 30mm 37万ポリゴン  メモリ不足でエラー
    int max_size_mb =  1;
    p_polylib->set_max_memory_size_mb( max_size_mb );
#endif

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
    PL_DBGOSH << "load() end"   <<endl;
    if( ret != PLSTAT_OK ) {
        PL_ERROSH << PolylibStat2::String(ret) <<endl;
        exit(1);
    }

    // デバッグ用ダンプ
    std::vector<PolygonGroup *>* pg_list = p_polylib->get_root_groups();
    PL_DBGOSH << "pg_list->size()="<<pg_list->size() <<endl;

    for(int i=0; i<pg_list->size(); i++ ) {
        std::string pg_name = (*pg_list)[i]->get_name();
        bool movable = (*pg_list)[i]->get_movable();
        std::vector<Triangle* > *tri_list =  (*pg_list)[i]->get_triangles();
        PL_DBGOSH << "polygon group i="<<i<<" name="<<pg_name<<" movable="<<movable<<" tri_list->size()="<<tri_list->size() <<endl;

    }

    // sphereのPlolygonGroupポインタ取得
    std::string pg_sphere_path = "sphere";
    PolygonGroup* pg_sphere = p_polylib->get_group( pg_sphere_path );
    if( pg_sphere == NULL )  {
        PL_ERROSH <<pg_sphere_path<<": [ERROR] Can't get PolygonGroup pointer" <<endl;
        exit(1);
    }

    // 全体個数確認
#ifdef MPI_PL
        // MPI版の場合
    std::vector<Triangle* > tri_list_sphere;
    pg_sphere->gather_polygons( tri_list_sphere );
#else
        // 逐次版の場合
    std::vector<Triangle* > *ptri_list_sphere =  pg_sphere->get_triangles();
    std::vector<Triangle* > & tri_list_sphere =  *ptri_list_sphere;
#endif
    if( myrank == 0 )  {
        if( tri_list_sphere.size() != num_global_polygon ) {
            PL_ERROSH <<"[ERROR] Load number of polygon  tri_list_sphere.size()="<<tri_list_sphere.size()<<"  num_global_polygon="<<num_global_polygon <<endl;
            exit(1);
        }
    }
        // 後処理
#ifdef MPI_PL
    for(int i=0; i<tri_list_sphere.size(); i++ ) {
        delete tri_list_sphere[i];
    }
#endif


    //-------------------------------------------
    //  セーブ
    //-------------------------------------------
    std::string  config_name_out;
    std::string  fmt_out = PolygonIO::FMT_STL_A;

#ifdef MPI_PL
    //  デバッグのため分割セーブ
    PL_DBGOSH << "save_parallel() start" <<endl;
    ret = p_polylib->save_parallel( config_name_out, fmt_out );
    if( ret != PLSTAT_OK ) {
        PL_ERROSH <<"[ERROR] save_parallel()" <<endl;
        exit(1);
    }
    PL_DBGOSH << "save_parallel() end   config_name_out="<<config_name_out <<endl;
#endif

    //  通常のセーブ
    PL_DBGOSH << "save() start" <<endl;
    ret = p_polylib->save( config_name_out, fmt_out );
    if( ret != PLSTAT_OK ) {
        PL_ERROSH <<"[ERROR] save()" <<endl;
        exit(1);
    }
    PL_DBGOSH << "save() end   config_name_out="<<config_name_out <<endl;


    //-------------------------------------------
    //  終了化
    //-------------------------------------------
    delete pg_list;

#ifdef MPI_PL
    // MPI終了化
    PL_DBGOSH << "MPI_Finalize() start" << endl;
    MPI_Finalize();
    PL_DBGOSH << "MPI_Finalize() end" << endl;
#endif

    if( myrank == 0 )  {
        std::cerr<<"------------------------------------------" <<endl;
        std::cerr<<"   PASS :  load_reduce_mem (Normal End)" <<endl;
        std::cerr<<"------------------------------------------" <<endl;
    }

    return 0;
}
