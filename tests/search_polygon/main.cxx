
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

//------------------------------
// 妥当性検証用リファレンス値
//------------------------------

// 全体ポリゴン数
const int num_polygon_total_reference = 192;  
// rank0検索ポリゴン数
const int num_srch_rank0_reference = 36;  
// rank3検索ポリゴン数
const int num_srch_rank3_reference = 36;  

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

    // ポリゴングループポインタ取得
    std::string   pg_sphere_path = "sphere_1";
    PolygonGroup* pg_sphere = p_polylib->get_group( pg_sphere_path );
    if( pg_sphere == NULL )  {
        PL_ERROSH <<pg_sphere<<": [ERROR] Can't get PolygonGroup pointer" <<endl;
        exit(1);
    }

    // ポリゴン情報取得（ランク毎）
    std::vector<Triangle* > *tri_list =  pg_sphere->get_triangles();
    if( tri_list == NULL || tri_list->size()==0 )  {
        PL_ERROSH <<"[ERROR] Can't get Polygon data" <<endl;
        exit(1);
    }
    PL_DBGOSH << "number of polygon(rank"<<myrank<<"): tri_list->size()="<<tri_list->size() <<endl;
    // 全体ポリゴン数の確認
    int num_polygon_total;
#ifdef MPI_PL
    num_polygon_total = pg_sphere->get_group_num_global_tria();
#else
    num_polygon_total = tri_list->size();
#endif
    	//  妥当性の検証
    if( num_polygon_total != num_polygon_total_reference ) {
        PL_ERROSH <<"[ERROR] num_polygon_total="<<num_polygon_total<<"  num_polygon_total_reference="<<num_polygon_total_reference <<endl;
        PL_DBGOSH <<"[ERROR] num_polygon_total="<<num_polygon_total<<"  num_polygon_total_reference="<<num_polygon_total_reference <<endl;
        exit(1);
    }

    // ポリゴン検索
    if( myrank == 0 ) {
        std::vector<Triangle* > srch_tri_list;
	BBox bbox;
        bbox.min.x=-40.0;  bbox.min.y=-40.0;  bbox.min.z=-40.0;
        bbox.max.x=-10.0;  bbox.max.y=-10.0;  bbox.max.z=+40.0;

        ret = pg_sphere->search( srch_tri_list, bbox, false);
        if( ret != PLSTAT_OK )  {
            PL_ERROSH <<"[ERROR] pg_sphere->search()" <<endl;
            exit(1);
        }
    		//  妥当性の検証
        if( srch_tri_list.size() != num_srch_rank0_reference ) {
            PL_ERROSH <<"[ERROR] srch_tri_list.size()="<<srch_tri_list.size()<<"  num_srch_rank0_reference="<<num_srch_rank0_reference <<endl;
            PL_DBGOSH <<"[ERROR] srch_tri_list.size()="<<srch_tri_list.size()<<"  num_srch_rank0_reference="<<num_srch_rank0_reference <<endl;
            exit(1);
        }
    }

    if( myrank == 3 ) {
        std::vector<Triangle* > srch_tri_list;
	//BBox bbox;
        //bbox.min.x=+10.0;  bbox.min.y=+10.0;  bbox.min.z=-40.0;
        //bbox.max.x= 40.0;  bbox.max.y= 40.0;  bbox.max.z=+40.0;
        //ret = pg_sphere->search( srch_tri_list, bbox, false);

        Vec3<PL_REAL> min_pos;
        min_pos.x=+10.0;  min_pos.y=+10.0;  min_pos.z=-40.0;
        Vec3<PL_REAL> max_pos;
        max_pos.x=+40.0;  max_pos.y=+40.0;  max_pos.z=+40.0;

        ret = p_polylib->search_polygons(
                             srch_tri_list,
                             pg_sphere_path, min_pos, max_pos, false
                          );
        if( ret != PLSTAT_OK )  {
            PL_ERROSH <<"[ERROR] p_polylib->search_polygons()" <<endl;
            exit(1);
        }
    		//  妥当性の検証
        if( srch_tri_list.size() != num_srch_rank3_reference ) {
            PL_ERROSH <<"[ERROR] srch_tri_list.size()="<<srch_tri_list.size()<<"  num_srch_rank3_reference="<<num_srch_rank3_reference <<endl;
            PL_DBGOSH <<"[ERROR] srch_tri_list.size()="<<srch_tri_list.size()<<"  num_srch_rank3_reference="<<num_srch_rank3_reference <<endl;
            exit(1);
        }
    }

    //-------------------------------------------
    //  終了化
    //-------------------------------------------

    // MPI終了化
#ifdef MPI_PL
    PL_DBGOSH << "MPI_Finalize() start" << endl;
    MPI_Finalize();
    PL_DBGOSH << "MPI_Finalize() end" << endl;
#endif

    if( myrank == 0 )  {
        std::cerr<<"------------------------------------------" <<endl;
        std::cerr<<"   PASS :  search_polygon (Normal End)" <<endl;
        std::cerr<<"------------------------------------------" <<endl;
    }

    return 0;
}
