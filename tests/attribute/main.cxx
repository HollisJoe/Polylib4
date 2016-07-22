
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
/// ポリゴン属性テスト
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
// 並列数４：X,Y方向を +50 の位置で２分割
static MyParallelInfo myParaInfos[4] = {
  {{-50,-50,-50,}, {10,10,20,}, {2, 2, 2,}, {10,10,10} },
  {{ 50,-50,-50,}, {10,10,20,}, {2, 2, 2,}, {10,10,10} },
  {{-50, 50,-50,}, {10,10,20,}, {2, 2, 2,}, {10,10,10} },
  {{ 50, 50,-50,}, {10,10,20,}, {2, 2, 2,}, {10,10,10} }
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
    int     num_global_polygon;
    PL_REAL eps = 0.001;  // 誤差判定用

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

    // 検索モード設定
    Polylib::set_srch_mode( true );  // 設定しているがSTLだと効かない

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
        PL_ERROSH <<"[ERROR] p_polylib->load() ret="<<ret <<endl;
        exit(1);
    }

    std::vector<PolygonGroup *>* pg_list = p_polylib->get_leaf_groups();
    if( pg_list == NULL ) {
        PL_ERROSH <<"[ERROR] p_polylib->get_leaf_groups() ret="<<ret <<endl;
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
        if( tri_list != NULL )  {
            num_tri = tri_list->size();
        }
        PL_DBGOSH << "polygon group i="<<i<<" name="<<pg_name<<" pg_path="<<pg_path <<endl;
        PL_DBGOSH << "   movable="<<movable<<" num_tri="<<num_tri <<endl;
    }

    // デバッグ用　グループ階層構造を標準出力に出力
    p_polylib->show_group_hierarchy();

    //-------------------------
    // ポリゴングループ 初期設定
    //    ポリゴン属性数設定
    //-------------------------
    std::string pg_path;
    PolygonGroup* pg_sphere;

    // PolygonGroupポインタ取得
    pg_path = "test/sphere";   // 複数あるときはパスで指定する
    PL_DBGOSH << "p_polylib->get_group("<<pg_path<<") start" <<endl;
    pg_sphere = p_polylib->get_group( pg_path );
    if( pg_sphere == NULL )  {
        PL_ERROSH <<pg_path<<": [ERROR] Can't get PolygonGroup pointer" <<endl;
        exit(1);
    }

    // 全体のポリゴン数取得
#ifdef MPI_PL
    num_global_polygon = pg_sphere->get_group_num_global_tria();
#else
    num_global_polygon = pg_sphere->get_group_num_tria();
#endif
    if( num_global_polygon == 0 )  {
        PL_ERROSH <<" [ERROR] num_global_polygon="<<num_global_polygon <<endl;
        exit(1);
    }

    // デバッグ用　グループ階層構造を標準出力に出力
    //p_polylib->show_group_info( pg_path );

    // 移動関数登録
    //pg_sphere->set_move_func( move_func );

    // グループの属性値取得
    std::string key;
    std::string sval;

    key = "int_val1";
    ret = pg_sphere->get_atr( key, sval ); 
    if( ret != PLSTAT_OK )  {
        PL_ERROSH <<key<<": [ERROR] Can't get PolygonGroup attribute" <<endl;
        exit(1);
    }
    //int ival1 = std::stoi( sval );  // use c++11
    int ival1 = std::atoi( sval.c_str() );
    key = "int_val2";
    ret = pg_sphere->get_atr( key, sval ); 
    //int ival2 = std::stoi( sval );
    int ival2 = std::atoi( sval.c_str() );
    key = "int_val3";
    ret = pg_sphere->get_atr( key, sval ); 
    //int ival3 = std::stoi( sval );
    int ival3 = std::atoi( sval.c_str() );
        // 妥当性の検証
    if( ival1!=11 || ival2!=22 || ival3!=33 )  {
        PL_ERROSH <<"[ERROR] PolygonGroup attribute integer values" <<endl;
        exit(1);
    }

    key = "real_val1";
    ret = pg_sphere->get_atr( key, sval ); 
    if( ret != PLSTAT_OK )  {
        PL_ERROSH <<key<<": [ERROR] Can't get PolygonGroup attribute" <<endl;
        exit(1);
    }
    //PL_REAL rval1 = std::stod( sval );  // use c++11
    PL_REAL rval1 = std::atof( sval.c_str() );
    key = "real_val2";
    ret = pg_sphere->get_atr( key, sval ); 
    //PL_REAL rval2 = std::stod( sval );
    PL_REAL rval2 = std::atof( sval.c_str() );
        // 妥当性の検証
    if( fabs(rval1-111.0)>eps || fabs(rval2-222.0)>eps )  {
        PL_ERROSH <<"[ERROR] PolygonGroup attribute real values" <<endl;
        exit(1);
    }

    // ポリゴングループの属性追加　テスト
    key  = "test_val";
    sval = "100";
    pg_sphere->set_atr( key, sval ); 
    sval = "999";  // testのため、一時的に適当な値にする
    ret = pg_sphere->get_atr( key, sval ); 
    if( ret != PLSTAT_OK )  {
        PL_ERROSH <<key<<": [ERROR] Can't get PolygonGroup attribute" <<endl;
        exit(1);
    }
    //int ival_test = std::stoi( sval );
    int ival_test = std::atoi( sval.c_str() );
    if( ival_test != 100 )  {
        PL_ERROSH <<"[ERROR] PolygonGroup attribute integer values  ival_test="<<ival_test <<endl;
        exit(1);
    }


    // ポリゴンの属性数設定
    int num_atrI = 3; 
    int num_atrR = 2;
    pg_sphere->set_num_polygon_atr( num_atrI,num_atrR ); 

    // グループ内の全ポリゴンにグループの属性値設定
    //    同じ値が設定される
    std::vector<Triangle* > *tri_list =  pg_sphere->get_triangles();
    for(int i=0; i<tri_list->size(); i++ ) {
        int*     pAtrI =  (*tri_list)[i]->get_pAtrI();
        PL_REAL* pAtrR =  (*tri_list)[i]->get_pAtrR();
        pAtrI[0]=ival1;  pAtrI[1]=ival2;  pAtrI[2]=ival3;
        pAtrR[0]=rval1;  pAtrR[1]=rval2;
    }

 
    //---------------------------------
    // Min/Maxテストのため一部値を書き換える
    //---------------------------------
    long long int id0 = (*tri_list)[0]->get_id();  // id offset
    for(int i=0; i<tri_list->size(); i++ ) {
        long long int id = (*tri_list)[i]->get_id();  // for debug
        int*     pAtrI =  (*tri_list)[i]->get_pAtrI();
        PL_REAL* pAtrR =  (*tri_list)[i]->get_pAtrR();

        if( id == (id0+50) ) {
            // Min値を更新
            pAtrI[0] = -pAtrI[0];
            pAtrR[0] = -pAtrR[0];
        } else if( id == (id0+150) ) {
            // Max値を更新
            pAtrI[1] = 2*pAtrI[1];
            pAtrR[1] = 2.0*pAtrR[1];
        }
    }

    //---------------------------------
    // ポリゴングループ内属性のMin/Max/Sum取得のテスト
    //---------------------------------
    int     ival_min, ival_max, ival_sum, ival_sum_cal; 
    PL_REAL rval_min, rval_max, rval_sum, rval_sum_cal; 
        // 整数１個目の確認
    pg_sphere->get_polygons_reduce_atrI( PL_OP_MIN, 0, ival_min );
    pg_sphere->get_polygons_reduce_atrI( PL_OP_MAX, 0, ival_max );
    pg_sphere->get_polygons_reduce_atrI( PL_OP_SUM, 0, ival_sum );
    ival_sum_cal = ival1*(num_global_polygon-2);
    if( ival_min != -ival1 || ival_max != ival1 || ival_sum != ival_sum_cal ) {
        PL_ERROSH <<"[ERROR] get_polygons_reduce_atrI() 1" <<endl;
        PL_ERROSH <<"    ival_min="<<ival_min<<" ival_max="<<ival_max<<" ival_sum="<<ival_sum <<endl;
        PL_ERROSH <<"    ival1   ="<<ival1   <<" ival_sum_cal="<<ival_sum_cal <<endl;
        exit(1);
    }
        // 整数２個目の確認
    pg_sphere->get_polygons_reduce_atrI( PL_OP_MIN, 1, ival_min );
    pg_sphere->get_polygons_reduce_atrI( PL_OP_MAX, 1, ival_max );
    pg_sphere->get_polygons_reduce_atrI( PL_OP_SUM, 1, ival_sum );
    ival_sum_cal = ival2*(num_global_polygon+1);
    if( ival_min != ival2 || ival_max != (2*ival2) || ival_sum != ival_sum_cal ) {
        PL_ERROSH <<"[ERROR] get_polygons_reduce_atrI() 2" <<endl;
        PL_ERROSH <<"    ival_min="<<ival_min<<" ival_max="<<ival_max<<" ival_sum="<<ival_sum <<endl;
        PL_ERROSH <<"    ival2   ="<<ival2   <<" ival_sum_cal="<<ival_sum_cal <<endl;
        exit(1);
    }
        // 実数１個目の確認
    pg_sphere->get_polygons_reduce_atrR( PL_OP_MIN, 0, rval_min );
    pg_sphere->get_polygons_reduce_atrR( PL_OP_MAX, 0, rval_max );
    pg_sphere->get_polygons_reduce_atrR( PL_OP_SUM, 0, rval_sum );
    rval_sum_cal = rval1*(num_global_polygon-2);
    if( fabs(rval_min-(-rval1))      > eps  || 
        fabs(rval_max-  rval1 )      > eps  ||
        fabs(rval_sum- rval_sum_cal) > eps        ) {
        PL_ERROSH <<"[ERROR] get_polygons_reduce_atrR() 1" <<endl;
        PL_ERROSH <<"    rval_min="<<rval_min<<" rval_max="<<rval_max<<" rval_sum="<<rval_sum <<endl;
        PL_ERROSH <<"    rval1   ="<<rval1   <<" rval_sum_cal="<<rval_sum_cal <<endl;
        exit(1);
    }
        // 実数２個目の確認
    pg_sphere->get_polygons_reduce_atrR( PL_OP_MIN, 1, rval_min );
    pg_sphere->get_polygons_reduce_atrR( PL_OP_MAX, 1, rval_max );
    pg_sphere->get_polygons_reduce_atrR( PL_OP_SUM, 1, rval_sum );
    rval_sum_cal = rval2*(num_global_polygon+1);
    if( fabs(rval_min-  rval2 )      > eps  || 
        fabs(rval_max-(2*rval2))     > eps  ||
        fabs(rval_sum- rval_sum_cal) > eps        ) {
        PL_ERROSH <<"[ERROR] get_polygons_reduce_atrR() 2" <<endl;
        PL_ERROSH <<"    rval_min="<<rval_min<<" rval_max="<<rval_max<<" rval_sum="<<rval_sum <<endl;
        PL_ERROSH <<"    rval2   ="<<rval1   <<" rval_sum_cal="<<rval_sum_cal <<endl;
        exit(1);
    }

    //---------------------------------
    // 保存処理
    //   属性はユーザ独自でファイルに出力する必要がある
    //---------------------------------
#if 0
    std::vector<Triangle*>  tri_list_tmp;
    ret = pg_sphere->gather_polygons( tri_list_tmp );
    std::string file_name_atr = "sphere_atr.txt";

    ofstream ofs( file_name_atr.c_str() );
    if ( ofs.fail() ) {
        PL_ERROSH << "[ERROR] Can't open " << file_name_atr << endl;
        exit(1);
    }

    for(int i=0; i<tri_list_tmp.size(); i++ ) {
        int*     pAtrI =  tri_list_tmp[i]->get_pAtrI();
        PL_REAL* pAtrR =  tri_list_tmp[i]->get_pAtrR();
        long long int id = tri_list_tmp[i]->get_id();  // for debug

        ofs <<"no="<<i<<" id="<<id<<" atrI="<<pAtrI[0]<<" "<<pAtrI[1]<<" "<<pAtrI[2]<<" atrR="<<pAtrR[0]<<" "<<pAtrR[1] <<endl;
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
        std::cerr<<"   PASS :  attribute (Normal End)" <<endl;
        std::cerr<<"------------------------------------------" <<endl;
    }

    return 0;
}
