
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
/// 1PE(ランク）あたり複数領域のテスト
///
////////////////////////////////////////////////////////////////////////////

#include "Polylib.h"
#include "pltest_func.h"

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
//  (初期値）
//   plane    x:   0〜100
//            y:   0〜100
//            z:   0〜  0
//   spehre   x: -50〜 50
//            y: -50〜 50
//            z:  50〜150
//  (移動）
//     sphereが移動する
//     移動量  x: +1.0/step  y: +1.0/step  z: 0.0/step
//     50step実行
//

#if 0
// 並列数１：分割なし
static MyParallelInfo myParaInfos[1] = {
  {{-50,-50,-50,}, {20,20,20,}, {2, 2, 2,}, {10,10,10} },
};
#endif
#if 1
// 並列数４
static ParallelBbox myParaInfos[6] = {
//   X側を短くした領域を4分割
//      X,Y方向を +0 の位置で２分割
//      空間 x: -50〜 50  y: -50〜150 z:   0〜150
//   ** 対角上にrank0,rank1に割り振る
  {{-50, -50,-50,}, { 5, 5,20,}, {1, 1, 1,}, {10,10,10} },  // rank0
  {{  0, -50,-50,}, { 5, 5,20,}, {1, 1, 1,}, {10,10,10} },  // rank1
  {{-50,   0,-50,}, { 5,15,20,}, {1, 1, 1,}, {10,10,10} },  // rank1
  {{  0,   0,-50,}, { 5,15,20,}, {1, 1, 1,}, {10,10,10} },  // rank0
//   X側の残りを2分割
//      Y方向を +0 の位置で２分割
//      空間 x:  50〜150  y: -50〜150 z:   0〜150
  {{ 50, -50,-50,}, {10, 5,20,}, {1, 1, 1,}, {10,10,10} },  // rank2 
  {{ 50,   0,-50,}, {10,15,20,}, {1, 1, 1,}, {10,10,10} }   // rank3 
};
#endif

#endif

#define MOVE_X_BY_STEP 1.0
#define MOVE_Y_BY_STEP 1.0
#define MOVE_Z_BY_STEP 0.0
#define EPS  0.01

const int num_global_polygon = 192;  // sphereのポリゴン数

//----------------------------------------------------
//  移動関数
//----------------------------------------------------

void move_func(
          PolygonGroup*      pg,
          PolylibMoveParams* params
       )
{
    //PL_DBGOSH << "----- move_func() start -----" << endl;

    std::string pg_name = pg->get_name();
    std::vector<Triangle* > *tri_list =  pg->get_triangles();
    //PL_DBGOSH << "move polygon group:  name="<<pg_name<<" tri_list->size()="<<tri_list->size() <<endl;
    // 呼ばれるたびに指定値ずつオフセットする
    PL_REAL x_offset =   MOVE_X_BY_STEP;
    PL_REAL y_offset =   MOVE_Y_BY_STEP;
    PL_REAL z_offset =   MOVE_Z_BY_STEP;

    for(int i=0; i<tri_list->size(); i++ ) {
        Vec3<PL_REAL>* vertex = (*tri_list)[i]->get_vertexes();

        for( int j=0; j<3; j++ ) {
            vertex[j].x += x_offset;
            vertex[j].y += y_offset;
            vertex[j].z += z_offset;
        }
    }

    // 頂点座標が移動したことにより、KD木の再構築が必要
    // 再構築フラグを立てる
    pg->set_need_rebuild();

    //PL_DBGOSH << "----- move_func() end   -----" << endl;
}


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
    std::string  config_name_out;
    std::string  fmt_out = PolygonIO::FMT_STL_A;

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

    // 並列計算関連情報の設定と初期化
#ifdef MPI_PL
    PL_DBGOSH << "(rk:"<<myrank<<") "<<"init_parallel_info() start" << endl;
    std::vector<ParallelBbox>  bboxes;

    // 対象領域を対角上に設定する
    if( myrank == 0 )  {
        bboxes.push_back( myParaInfos[0] );
        bboxes.push_back( myParaInfos[3] );
    } else if( myrank == 1 )  {
        bboxes.push_back( myParaInfos[1] );
        bboxes.push_back( myParaInfos[2] );
    } else if( myrank == 2 )  {
        bboxes.push_back( myParaInfos[4] );
    } else if( myrank == 3 )  {
        bboxes.push_back( myParaInfos[5] );
    } else {
        PL_ERROSH <<"[ERROR] rankNo error" <<endl;
        exit(1);
    }

    ret = p_polylib->init_parallel_info(
                MPI_COMM_WORLD,
                bboxes
            );
    PL_DBGOSH << "init_parallel_info() end" << endl;
#endif

    //-------------------------------------------
    //  ロード
    //-------------------------------------------

    // 初期化ファイルを指定してデータロード
    PL_DBGOSH << "load() start" <<endl;
    ret = p_polylib->load( config_file_name );
    PL_DBGOSH << "load() end"   <<endl;

    // sphereのPlolygonGroupポインタ取得
    std::string pg_sphere_path = "sphere";
    PolygonGroup* pg_sphere = p_polylib->get_group( pg_sphere_path );
    if( pg_sphere == NULL )  {
        PL_ERROSH <<pg_sphere_path<<": [ERROR] Can't get PolygonGroup pointer" <<endl;
        exit(1);
    }

#ifdef MPI_PL
    //  デバッグ用：分割Load確認 (ランク毎に出力） 
    PL_DBGOSH << "save_parallel() in start" <<endl;
    std::string extend_in = "in";
    ret = p_polylib->save_parallel( config_name_out, fmt_out, extend_in );
    if( ret != PLSTAT_OK ) {
        PL_ERROSH <<"[ERROR] save_parallel() in ret="<<ret <<endl;
        exit(1);
    }
    PL_DBGOSH << "save_parallel() in end   config_name_out="<<config_name_out <<endl;
#endif

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
    PL_REAL  x_ave_in, y_ave_in, z_ave_in;
    if( myrank == 0 )  {
        if( tri_list_sphere.size() != num_global_polygon ) {
            PL_ERROSH <<"[ERROR] Load number of polygon  tri_list_sphere.size()="<<tri_list_sphere.size()<<"  num_global_polygon="<<num_global_polygon <<endl;
            exit(1);
        }
        // 妥当性検証用：座標平均値取得
        pltest_get_polygon_average( tri_list_sphere, x_ave_in, y_ave_in, z_ave_in );
    }
        // 後処理
#ifdef MPI_PL
    for(int i=0; i<tri_list_sphere.size(); i++ ) {
        delete tri_list_sphere[i];
    }
#endif

    //-------------------------------------------
    // 移動関数登録
    //-------------------------------------------
    pg_sphere->set_move_func( move_func );

    //---------------------------------
    // タイムステップループ
    //---------------------------------

    // move parameter
    PolylibMoveParams params;
    memset( params.m_params, 0x00, 10*sizeof(PL_REAL) ); 

    int nstep = 50;
    for(int istep=0; istep<nstep ; istep++ )  {

        // moveパラメタ設定
        params.m_current_step = istep;
        params.m_next_step    = istep + 1;
        params.m_delta_t      = 1.0;

        // move実行
        ret = p_polylib->move( params ); 
        if( ret != PLSTAT_OK ) {
            PL_ERROSH <<"[ERROR] p_polylib->move() ret="<<ret <<endl;
            exit(1);
        }

        // migrate実行
#ifdef MPI_PL
        ret = p_polylib->migrate(); 
        if( ret != PLSTAT_OK ) {
            PL_ERROSH <<"[ERROR] p_polylib->migrate() ret="<<ret <<endl;
            exit(1);
        }
#endif

    }

    //-------------------------------------------
    //  セーブ
    //-------------------------------------------

    PL_DBGOSH << "save() start" <<endl;
    ret = p_polylib->save( config_name_out, fmt_out );
    if( ret != PLSTAT_OK ) {
        PL_ERROSH <<"[ERROR] save()" <<endl;
        exit(1);
    }
    PL_DBGOSH << "save() end   config_name_out="<<config_name_out <<endl;

#ifdef MPI_PL
    //  デバッグ用：セーブ（各ランクごと）
    PL_DBGOSH << "save_parallel() out start" <<endl;
    std::string extend_out = "out";
    ret = p_polylib->save_parallel( config_name_out, fmt_out, extend_out );
    if( ret != PLSTAT_OK ) {
        PL_ERROSH <<"[ERROR] save_parallel() out ret="<<ret <<endl;
        exit(1);
    }
    PL_DBGOSH << "save_parallel() out end   config_name_out="<<config_name_out <<endl;
#endif

    // 全体個数確認
#ifdef MPI_PL
    // MPI版の場合
    std::vector<Triangle* > tri_list_sphere_out;
    pg_sphere->gather_polygons( tri_list_sphere_out );
#else
    // 逐次版の場合
    std::vector<Triangle* > *ptri_list_sphere_out =  pg_sphere->get_triangles();
    std::vector<Triangle* > & tri_list_sphere_out =  *ptri_list_sphere_out;
#endif
    PL_REAL  x_ave_out, y_ave_out, z_ave_out;
    if( myrank == 0 )  {
        if( tri_list_sphere_out.size() != num_global_polygon ) {
            PL_ERROSH <<"[ERROR] Load number of polygon  tri_list_sphere_out.size()="<<tri_list_sphere_out.size()<<"  num_global_polygon="<<num_global_polygon <<endl;
            exit(1);
        }
        // 平均値取得
        pltest_get_polygon_average( tri_list_sphere_out, x_ave_out, y_ave_out, z_ave_out );
    }
    // 後処理
#ifdef MPI_PL
    for(int i=0; i<tri_list_sphere_out.size(); i++ ) {
        delete tri_list_sphere_out[i];
    }
#endif

    //---------------------------------
    // 妥当性検証
    //---------------------------------
    if( myrank == 0 )  {
        // 移動推定値
        PL_REAL x_ave_ref = x_ave_in + nstep*MOVE_X_BY_STEP;
        PL_REAL y_ave_ref = y_ave_in + nstep*MOVE_Y_BY_STEP;
        PL_REAL z_ave_ref = z_ave_in + nstep*MOVE_Z_BY_STEP;

        // 誤差判定
        if( fabs(x_ave_out-x_ave_ref) > EPS  ||
            fabs(y_ave_out-y_ave_ref) > EPS  ||
            fabs(z_ave_out-z_ave_ref) > EPS      )  {
                PL_ERROSH <<"[ERROR] move distance error" <<endl;
                PL_ERROSH <<"x_ave_in ="<<x_ave_in<<" y_ave_in ="<<y_ave_in<<" z_ave_in ="<<z_ave_in <<endl;
                PL_ERROSH <<"x_ave_out="<<x_ave_out<<" y_ave_out="<<y_ave_out<<" z_ave_out="<<z_ave_out <<endl;
                exit(1);
        }
    }

    //-------------------------------------------
    //  終了化
    //-------------------------------------------

#ifdef MPI_PL
    // MPI終了化
    PL_DBGOSH << "MPI_Finalize() start" << endl;
    MPI_Finalize();
    PL_DBGOSH << "MPI_Finalize() end" << endl;
#endif

    if( myrank == 0 )  {
        std::cerr<<"------------------------------------------" <<endl;
        std::cerr<<"   PASS :  multi_bbox (Normal End)" <<endl;
        std::cerr<<"------------------------------------------" <<endl;
    }

    return 0;
}

