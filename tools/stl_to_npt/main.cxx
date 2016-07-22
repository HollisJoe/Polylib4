
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
/// STL→長田パッチファイル変換プログラム
///
///     高速化のため、PolilibのKD treeによる検索を使用している
///
////////////////////////////////////////////////////////////////////////////

#ifdef MPI_PL
#include "mpi.h"
#endif
#include "CalcGeo.h"
#include "Npt.h"
#include "StlToNpt.h"
#include <string.h>

//  STLが平面で特定方向が0.0になった時に足す距離
//  Polylibの検索がうまく行くようにサイズを拡大するために使用する
#define  EXPAND_SIZE    0.1

//  近傍ポリゴンの検索距離
#define  SRCH_DISTANCE  0.1

//  同一頂点と見なす距離
#define  SAME_DISTANCE  0.001

//  許容誤差
#define  ALW_EPS        0.00001

//  エッジ判定角度（２面の法線ベクトルがこれ以上違うとエッジとみなす）
//      > 180.0  エッジ判定をしない
#define  DEFAULT_EDGE_DEGREE    180.0

const std::string   group_name_in  = "stl_in";

static void Usage(void)
{
    cerr<<endl;
    cerr<< "Usage: stl_to_npt stl_file [ edge_degree ]" <<endl;
    cerr<< "Arguments:" <<endl;
    cerr<< "  stl_file          stl file name" <<endl;
    cerr<< "  edge_degree       edge angle ( default = 180.0 )" <<endl;
    cerr<< "                        if degree between surfeces > edge_degree, it is edge." <<endl;
    cerr<<endl;
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

    std::vector<Triangle*> *tri_list=NULL;
    int num_tri_total;
    PL_REAL min[3] = { +1.0e+10, +1.0e+10, +1.0e+10 };
    PL_REAL max[3] = { -1.0e+10, -1.0e-10, -1.0e-10 };

#ifdef MPI_PL
    // MPI初期化
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &num_rank );
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
#else
    num_rank = 1;
    myrank   = 0;
#endif

    if( myrank == 0 )  {
        PL_DBGOSH << "#### stl_to_ntp (start) ####" <<endl;
    }

    if( argc > 3 || argc < 2 ) {
       Usage();
       exit(1);
    }

    //-------------------------------------------
    //  ファイルパス
    //-------------------------------------------

    std::string stl_file_name;   // 入力：stlファイル名
    std::string fmt_in;          // 入力：ファイルフォーマット
    std::string npt_file_name;   // 出力：nptファイル名
    std::string fmt_out;         // 出力：ファイルフォーマット

    stl_file_name = argv[1];
    fmt_in = PolygonIO::input_file_format( stl_file_name );

    char* fname = PolylibNS::get_fname_fr_path( stl_file_name );
    npt_file_name = fname;
    npt_file_name += ".npt";
    if( fmt_in == PolygonIO::FMT_STL_A ) {
        fmt_out = PolygonIO::FMT_NPT_A;     // Text
    } else if( fmt_in == PolygonIO::FMT_STL_B ) {
        fmt_out = PolygonIO::FMT_NPT_B;     // Binary 
    } else {
        Usage();
        exit(1);
    }

    if( myrank == 0 )  {
        PL_DBGOSH << " stl_file_name="<<stl_file_name <<endl;
        PL_DBGOSH << " fmt_in       ="<<fmt_in <<endl;
        PL_DBGOSH << " npt_file_name="<<npt_file_name <<endl;
        PL_DBGOSH << " fmt_out      ="<<fmt_out <<endl;
    }

    //-------------------------------------------
    //  エッジ許容角度
    //-------------------------------------------
    PL_REAL edge_degree;
    if( argc < 3 ) {
        edge_degree = DEFAULT_EDGE_DEGREE;
    } else {
        double edge_degree_double = atof(argv[2]);
        if( edge_degree_double < ALW_EPS ) {
            Usage();
            exit(1);
        } else {
            edge_degree = edge_degree_double;
        }
    }
    if( myrank == 0 )  {
        PL_DBGOSH << " edge_degree="<<edge_degree <<endl;
    }
    
    // Polylib初期化
    Polylib* p_polylib = Polylib::get_instance();

    // ランク０処理
    if( myrank == 0 ) 
    {
        tri_list = new vector<Triangle*>;

        //-------------------------------------------
        // stlファイル読み込み
        //-------------------------------------------

        // stlファイル読み込み
        ret = PolygonIO::load( tri_list, stl_file_name, fmt_in );
        if( ret != PLSTAT_OK ) {
            PL_ERROSH << "[ERROR] load error stl_file_name=" << stl_file_name << endl;
            exit(1);
        }
    
        if( tri_list->size() == 0 ) {
            PL_ERROSH << "[ERROR] no traiangle in stl file. stl_file_name=" << stl_file_name << endl;
        }

        num_tri_total = tri_list->size();
#ifdef DEBUG
        PL_DBGOSH << "(rk:"<<myrank<<") "<<"num_tri_total="<<num_tri_total <<endl;
#endif

        // MinMax取得
        for( int i=0; i<tri_list->size(); i++ ) {
            Vec3<PL_REAL>* vertex = (*tri_list)[i]->get_vertexes();
            for( int j=0; j<3; j++ ) {
                if( vertex[j].x < min[0] ) min[0] = vertex[j].x;
                if( vertex[j].y < min[1] ) min[1] = vertex[j].y;
                if( vertex[j].z < min[2] ) min[2] = vertex[j].z;
                if( vertex[j].x > max[0] ) max[0] = vertex[j].x;
                if( vertex[j].y > max[1] ) max[1] = vertex[j].y;
                if( vertex[j].z > max[2] ) max[2] = vertex[j].z;
            }
        }
    } 
    // ランク０以外
    else
    {
        tri_list = NULL;
        num_tri_total = 0;
    }

#ifdef MPI_PL
  // 並列処理

    // min,maxのランク０以外へのブロード
    iret = MPI_Bcast( min, 3, PL_MPI_REAL, 0, MPI_COMM_WORLD );
    iret = MPI_Bcast( max, 3, PL_MPI_REAL, 0, MPI_COMM_WORLD );
#ifdef DEBUG
    PL_DBGOSH << "(rk:"<<myrank<<") "<< "min="<<min[0]<<" "<<min[1]<<" "<<min[2] << endl;
    PL_DBGOSH << "(rk:"<<myrank<<") "<< "max="<<max[0]<<" "<<max[1]<<" "<<max[2] << endl;
#endif


    //-------------------------------------------
    // 並列計算関連情報の設定と初期化
    //-------------------------------------------

    // 各ランクの担当領域を決める
    //      一番長い軸をブロック分割する
    PL_REAL      bpos[3];                   // 自PE担当領域の基点座標
    unsigned int bbsize[3];                 // 同、計算領域のボクセル数
    unsigned int gcsize[3] = { 1, 1, 1 };   // 同、ガイドセルのボクセル数
    PL_REAL      dx[3];                     // 同、ボクセル１辺の長さ
    PL_REAL      endpos[3];                 // 自PE担当領域の終点座標 検索時の判定に使用する

    int nvoxel_tot_base  = 1000;    //  分割方向の総ボクセル数(ベース）
    int nvoxel_tot;                 //  分割方向の総ボクセル数
    int nvoxel_rank;                //  分割方向の各欄ンくのボクセル数=nvoxel_tot/num_rank
    int nvoxel2 = 10;               //  分割方向以外のボクセル数（＝総ボクセル数）

    PL_REAL size_area[3];
    size_area[0] = max[0] - min[0];
    size_area[1] = max[1] - min[1];
    size_area[2] = max[2] - min[2];
    
    // ポリゴンが平面の場合、特定方向のサイズが0.0となる可能性があるので、適当に足す
    if( size_area[0] < EXPAND_SIZE ) size_area[0] = EXPAND_SIZE;
    if( size_area[1] < EXPAND_SIZE ) size_area[1] = EXPAND_SIZE;
    if( size_area[2] < EXPAND_SIZE ) size_area[2] = EXPAND_SIZE;

    // 分割方向のボクセル数決定
    if( nvoxel_tot_base%num_rank == 0 ) {
        nvoxel_rank = nvoxel_tot_base / num_rank;
        nvoxel_tot = nvoxel_tot_base;
    } else {
        nvoxel_rank = nvoxel_tot / num_rank + 1;
        nvoxel_tot  = nvoxel_rank*num_rank;
    }

#ifdef DEBUG
    PL_DBGOSH << "(rk:"<<myrank<<") "<<" size_area="<<size_area[0]<<" "<< size_area[1]<<" "<<size_area[2] << endl;
    PL_DBGOSH << "(rk:"<<myrank<<") "<<" nvoxel_rank="<<nvoxel_rank<<" nvoxel_tot="<<nvoxel_tot << endl;
#endif

    // 一番長い方向を分割の方向とする
    if( size_area[0] >= size_area[1] && size_area[0] >= size_area[2] ) {
#ifdef DEBUG
        PL_DBGOSH << "(rk:"<<myrank<<") "<< "x divide" << endl;
#endif
        // X方向を分割する
        bbsize[0] = nvoxel_rank;
        bbsize[1] = nvoxel2;
        bbsize[2] = nvoxel2;
        dx[0] = size_area[0] / nvoxel_tot;
        dx[1] = size_area[1] / nvoxel2;
        dx[2] = size_area[2] / nvoxel2;

        bpos[0] = min[0] + myrank*(size_area[0]/num_rank);
        bpos[1] = min[1];
        bpos[2] = min[2];
        endpos[0] = min[0] + (myrank+1)*(size_area[0]/num_rank);
        endpos[1] = min[1] + size_area[1];
        endpos[2] = min[2] + size_area[2];

    } else if( size_area[1] >= size_area[2] ) {
#ifdef DEBUG
        PL_DBGOSH << "(rk:"<<myrank<<") "<< "y divide" << endl;
#endif
        // Y方向を分割する
        bbsize[0] = nvoxel2;
        bbsize[1] = nvoxel_rank;
        bbsize[2] = nvoxel2;
        dx[0] = size_area[0] / nvoxel2;
        dx[1] = size_area[1] / nvoxel_tot;
        dx[2] = size_area[2] / nvoxel2;

        bpos[0] = min[0];
        bpos[1] = min[1] + myrank*(size_area[1]/num_rank);
        bpos[2] = min[2];
        endpos[0] = min[0] + size_area[0];
        endpos[1] = min[1] + (myrank+1)*(size_area[1]/num_rank);
        endpos[2] = min[2] + size_area[2];

    } else {
#ifdef DEBUG
        PL_DBGOSH << "(rk:"<<myrank<<") "<< "z divide" << endl;
#endif
        // Z方向を分割する
        bbsize[0] = nvoxel2;
        bbsize[1] = nvoxel2;
        bbsize[2] = nvoxel_rank;
        dx[0] = size_area[0] / nvoxel2;
        dx[1] = size_area[1] / nvoxel2;
        dx[2] = size_area[2] / nvoxel_tot;

        bpos[0] = min[0];
        bpos[1] = min[1];
        bpos[2] = min[2] + myrank*(size_area[2]/num_rank);
        endpos[0] = min[0] + size_area[0];
        endpos[1] = min[1] + size_area[1];
        endpos[2] = min[2] + myrank*(size_area[2]/num_rank);

    }
#ifdef DEBUG
    PL_DBGOSH << "(rk:"<<myrank<<") "<<" bbsize="<<bbsize[0]<<" "<<bbsize[1]<<" "<<bbsize[2] << endl;
    PL_DBGOSH << "(rk:"<<myrank<<") "<<" dx    ="<<dx    [0]<<" "<<dx    [1]<<" "<<dx    [2] << endl;
    PL_DBGOSH << "(rk:"<<myrank<<") "<<" bpos  ="<<bpos  [0]<<" "<<bpos  [1]<<" "<<bpos  [2] << endl;
    PL_DBGOSH << "(rk:"<<myrank<<") "<<" endpos="<<endpos[0]<<" "<<endpos[1]<<" "<<endpos[2] << endl;
#endif

    // 並列計算関連情報の設定と初期化
#ifdef DEBUG
    PL_DBGOSH << "(rk:"<<myrank<<") "<<"init_parallel_info() start" << endl;
#endif
    ret = p_polylib->init_parallel_info(
                MPI_COMM_WORLD,
                bpos, bbsize, gcsize, dx 
            );
#ifdef DEBUG
    PL_DBGOSH << "init_parallel_info() end" << endl;
#endif

#else

  // 逐次処理
    PL_REAL      bpos[3];                   // 自PE担当領域の基点座標
    PL_REAL      endpos[3];                 // 自PE担当領域の終点座標 検索時の判定に使用する

    bpos[0]   = min[0] - ALW_EPS;
    bpos[1]   = min[1] - ALW_EPS;
    bpos[2]   = min[2] - ALW_EPS;
    endpos[0] = min[0] + ALW_EPS;
    endpos[1] = max[1] + ALW_EPS;
    endpos[2] = max[2] + ALW_EPS;

#endif
// !MPI_PL

    //-------------------------------------------
    // KD treeを使った検索を行うため、Polylib環境設定
    //-------------------------------------------

    // STL入力用グループ作成
    PolygonGroup* p_stl_polygrp = new PolygonGroup();
    p_stl_polygrp->set_name( group_name_in );
    
    //PL_DBGOSH << "add_pg_list() start" << endl;
    p_polylib->add_pg_list( p_stl_polygrp );    // ルートポリゴンとして登録
    //PL_DBGOSH << "add_pg_list() end" << endl;

    // 3角形の登録
    //      ランク０以外は0が設定されている
    //PL_DBGOSH << "set_triangles_ptr start" << endl;
    p_stl_polygrp->set_triangles_ptr( tri_list );   // 3角形の設定（ポインタのみコピー）
    //PL_DBGOSH << "set_triangles_ptr end" << endl;
    
    // ポリゴンの法線ベクトルの計算、面積の計算、KD木の生成を行う
    //    scatter_polygons()の中でsearchを使う。KD木が必要になる
    if( myrank == 0 ) {
        //PL_DBGOSH << "build_polygon_tree() 1 start" << endl;
        p_stl_polygrp->build_polygon_tree();
        //PL_DBGOSH << "build_polygon_tree() 1 end" << endl;
    }

#ifdef MPI_PL
    // 各ランクにポリゴン情報を分散する

    //PL_DBGOSH << "scatter_polygons() start" << endl;
    p_stl_polygrp->scatter_polygons();
    //PL_DBGOSH << "scatter_polygons() end" << endl;
    
    // ポリゴンの法線ベクトルの計算、面積の計算、KD木の生成を行う
    //    必要なのは近傍のポリゴンを検索するためのKD木の生成

    //PL_DBGOSH << "build_polygon_tree() 2 start" << endl;
    p_stl_polygrp->build_polygon_tree();
    //PL_DBGOSH << "build_polygon_tree() 2 end" << endl;
#endif

    // 自ランクのtri_list取得
    vector<Triangle*> *tri_list_rank = p_stl_polygrp->get_triangles();
    //PL_DBGOSH << " each rank : tri_list_rank->size()="<<tri_list_rank->size() << endl;
    
    int num_tri_rank = tri_list_rank->size();
    PL_REAL** p1_vertex_norm = alloc_array_2d_real( num_tri_rank, 3 );
    PL_REAL** p2_vertex_norm = alloc_array_2d_real( num_tri_rank, 3 );
    PL_REAL** p3_vertex_norm = alloc_array_2d_real( num_tri_rank, 3 );
    int* p1_setted_flag = new int[num_tri_rank];
    int* p2_setted_flag = new int[num_tri_rank];
    int* p3_setted_flag = new int[num_tri_rank];

    //-------------------------------------------
    // 自ランクの頂点ベクトルの決定
    //-------------------------------------------

#ifdef DEBUG
    PL_DBGOSH << "decide vertex normal start" << endl;
#endif
    for( int i=0; i< num_tri_rank; i++ ) {
    
        // ３角形の情報取得
        Vec3<PL_REAL>  *vertex    = (*tri_list_rank)[i]->get_vertexes();
        Vec3<PL_REAL>  norm = (*tri_list_rank)[i]->get_normal();
        long long int id = (*tri_list_rank)[i]->get_id();
        PL_REAL p1[3], p2[3], p3[3];
        
        p1[0] = vertex[0].x;  p1[1] = vertex[0].y;  p1[2] = vertex[0].z;
        p2[0] = vertex[1].x;  p2[1] = vertex[1].y;  p2[2] = vertex[1].z;
        p3[0] = vertex[2].x;  p3[1] = vertex[2].y;  p3[2] = vertex[2].z;

#ifdef DEBUG
        PL_DBGOSH << "===== i="<<i<<" tri id="<<id<<" =====" <<endl;
        PL_DBGOSH << " p1="<<p1[0]<<" "<<p1[1]<<" "<<p1[2]<<"  p2="<<p2[0]<<" "<<p2[1]<<" "<<p2[2]<<"  p3="<<p3[0]<<" "<<p3[1]<<" "<<p3[2] << endl;
#endif
        PL_REAL plane_norm[3];
        plane_norm[0] = norm.x;  plane_norm[1] = norm.y;  plane_norm[2] = norm.z;
#ifdef DEBUG
        PL_DBGOSH << " plane_norm="<<plane_norm[0]<<" "<<plane_norm[1]<<" "<<plane_norm[2] <<endl;
#endif
        

        // 頂点１ベクトル取得
        iret = get_vertex_normal( 
                        *p_stl_polygrp, bpos, endpos, id, plane_norm, p1, edge_degree, 
                        p1_vertex_norm[i], p1_setted_flag[i]
                    );
        if( iret != 0 ) {
            PL_ERROSH << "### Error: get_vertex_normal() 1 i="<<i<<" iret="<<iret << endl;
            exit(997);
        }
#ifdef DEBUG
        PL_DBGOSH << "  p1_setted_flag[i]="<<p1_setted_flag[i] <<endl;
        if( p1_setted_flag[i] == 1 ) {
            PL_DBGOSH << " p1_vertex_norm[i]="<<p1_vertex_norm[i][0]<<" "<<p1_vertex_norm[i][1]<<" "<<p1_vertex_norm[i][2] <<endl;
        }
#endif

        // 頂点２ベクトル取得
        iret = get_vertex_normal( 
                        *p_stl_polygrp, bpos, endpos, id, plane_norm, p2, edge_degree, 
                        p2_vertex_norm[i], p2_setted_flag[i]
                    );
        if( iret != 0 ) {
            PL_ERROSH << "### Error: get_vertex_normal() 2 i="<<i<<" iret="<<iret << endl;
            exit(998);
        }
#ifdef DEBUG
        PL_DBGOSH << "  p2_setted_flag[i]="<<p2_setted_flag[i] <<endl;
        if( p2_setted_flag[i] == 1 ) {
            PL_DBGOSH << " p2_vertex_norm[i]="<<p2_vertex_norm[i][0]<<" "<<p2_vertex_norm[i][1]<<" "<<p2_vertex_norm[i][2] <<endl;
        }
#endif
        
        // 頂点３ベクトル取得
        iret = get_vertex_normal( 
                        *p_stl_polygrp, bpos, endpos, id, plane_norm, p3, edge_degree,
                        p3_vertex_norm[i], p3_setted_flag[i]
                    );
        if( iret != 0 ) {
            PL_ERROSH << "### Error: get_vertex_normal() 3 i="<<i<<" iret="<<iret << endl;
            exit(999);
        }
#ifdef DEBUG
        PL_DBGOSH << "  p3_setted_flag[i]="<<p3_setted_flag[i] <<endl;
        if( p3_setted_flag[i] == 1 ) {
            PL_DBGOSH << " p3_vertex_norm[i]="<<p3_vertex_norm[i][0]<<" "<<p3_vertex_norm[i][1]<<" "<<p3_vertex_norm[i][2] <<endl;
        }
#endif


    }
#ifdef DEBUG
    PL_DBGOSH << "decide vertex normal end" << endl;
#endif

    std::vector<NptTriangle*>   npt_list;
    npt_list.reserve(num_tri_total);     // allocationを高速化するため、サイズを予約する

#ifdef MPI_PL
    // MPI並列

    //-------------------------------------------
    // 各ランクにあるポリゴンを集約する
    //      メモリの消費量を抑えるため、１対１の通信を複数回繰り返すこととする
    //-------------------------------------------
#ifdef DEBUG
    PL_DBGOSH << "gather_polygons_to_npt() start" << endl;
#endif
    iret = gather_polygons_to_npt(
                    num_rank, myrank,
                    num_tri_rank,
                    tri_list_rank,
                    p1_vertex_norm, p1_setted_flag,
                    p2_vertex_norm, p2_setted_flag,
                    p3_vertex_norm, p3_setted_flag,
                    num_tri_total,
                    npt_list
                );
#ifdef DEBUG
    PL_DBGOSH << "gather_polygons_to_npt() end" << endl;
#endif

#else
    // 逐次
    for(int i=0; i<num_tri_rank; i++ )  {
        PL_REAL     p1[3], p2[3], p3[3];
        PL_REAL     cp_side1_1[3], cp_side1_2[3], cp_side2_1[3], cp_side2_2[3];
        PL_REAL     cp_side3_1[3], cp_side3_2[3], cp_center [3];
        NpatchParam npt_param;

        long long int id = (*tri_list_rank)[i]->get_id();
        Vec3<PL_REAL>  *vertex    = (*tri_list_rank)[i]->get_vertexes();
        p1[0] = vertex[0].x;  p1[1] = vertex[0].y;  p1[2] = vertex[0].z;
        p2[0] = vertex[1].x;  p2[1] = vertex[1].y;  p2[2] = vertex[1].z;
        p3[0] = vertex[2].x;  p3[1] = vertex[2].y;  p3[2] = vertex[2].z;
        // 長田パッチパラメータ取得
        iret = npt_param_crt(
                    p1, p1_vertex_norm[i],
                    p2, p2_vertex_norm[i],
                    p3, p3_vertex_norm[i],
                    cp_side1_1, cp_side1_2,
                    cp_side2_1, cp_side2_2,
                    cp_side3_1, cp_side3_2,
                    cp_center
                );
        if( iret != 0 ) {
            PL_ERROSH << "[ERROR]  npt_param_crt() iret=" << iret
                << " p1=" << p1 << " p2=" << p2 << " p3=" << p3 << endl;
            exit(1);
        }
        // 長田パッチ生成のための型変換
        REAL_TO_VEC3(cp_side1_1,npt_param.cp_side1_1);
        REAL_TO_VEC3(cp_side1_2,npt_param.cp_side1_2);
        REAL_TO_VEC3(cp_side2_1,npt_param.cp_side2_1);
        REAL_TO_VEC3(cp_side2_2,npt_param.cp_side2_2);
        REAL_TO_VEC3(cp_side3_1,npt_param.cp_side3_1);
        REAL_TO_VEC3(cp_side3_2,npt_param.cp_side3_2);
        REAL_TO_VEC3(cp_center,npt_param.cp_center);
        // 長田パッチ生成
        NptTriangle *pNpt = new NptTriangle( vertex, npt_param, id );
        npt_list.push_back( pNpt );
    }

#endif

    //-------------------------------------------
    // 長田パッチファイル保存
    //-------------------------------------------

    // NPTファイル出力
    if( myrank == 0 ) {
        std::vector<Triangle*>* pTtri_list = (std::vector<Triangle*>*)(&npt_list);
#ifdef DEBUG
        PL_DBGOSH << "save() npt_list.size()   ="<<npt_list.size() << endl;
        PL_DBGOSH << "save() pTtri_list->size()="<<pTtri_list->size() << endl;

        PL_DBGOSH << "save() start" << endl;
#endif
        ret = PolygonIO::save( pTtri_list, npt_file_name, fmt_out );
#ifdef DEBUG
        PL_DBGOSH << "save() end" << endl;
#endif
    }

    //------------------------------------------
    // 終了化処理
    //-------------------------------------------

    for( int i=0; i< npt_list.size(); i++ ) {
        delete npt_list[i];
    }

    free_array_2d( (void**)p1_vertex_norm );
    free_array_2d( (void**)p2_vertex_norm );
    free_array_2d( (void**)p3_vertex_norm );
    delete[] p1_setted_flag;
    delete[] p2_setted_flag;
    delete[] p3_setted_flag;

#ifdef MPI_PL
#ifdef DEBUG
    PL_DBGOSH << "MPI_Finalize() start" << endl;
#endif
    // MPI終了化
    MPI_Finalize();
#ifdef DEBUG
    PL_DBGOSH << "MPI_Finalize() end" << endl;
#endif
#endif

    if( myrank == 0 )  {
        PL_DBGOSH << "#### stl_to_ntp (end) ####" <<endl;
    }

    return 0;
}


///
/// 頂点ベクトル決定
///
/// @param [in]    pg           ポリゴングループ
/// @param [in]    min          担当領域Min座標
/// @param [in]    max          担当領域Max座標
/// @param [in]    id           自身の３角形のinternal_id
/// @param [in]    norm_self    自身の３角形の法線ベクトル
/// @param [in]    pos          頂点座標（１点）
/// @param [in]    edge_degree  エッジ判定角度
/// @param [out]   vertex_norm  頂点法線ベクトル（単位ベクトル）
/// @param [out]   isetted_flg  頂点設定フラグ
///                                =0 設定なし
///                                =1 設定あり
//cc///                                =0 自身の３角形以外で設定なし
//cc///                                =1 自身の３角形以外で設定あり
/// @return リターンコード   =0 正常  !=0 異常
/// @attention
///     

int get_vertex_normal (
        PolygonGroup& pg,
        PL_REAL  min[3],
        PL_REAL  max[3],
        long long int id,
        PL_REAL  norm_self[3],
        PL_REAL  pos[3],
        PL_REAL  edge_degree,
        PL_REAL  vertex_norm[3],
        int& isetted_flg
    )
{
    isetted_flg = 0;
    vertex_norm[0] = 0.0;
    vertex_norm[1] = 0.0;
    vertex_norm[2] = 0.0;

    // 自身の担当領域に含まれない場合は処理しない
    //     他のランクで法線ベクトルが設定されるはず
    if( pos[0] < (min[0]-ALW_EPS) ) return 0;
    if( pos[1] < (min[1]-ALW_EPS) ) return 0;
    if( pos[2] < (min[2]-ALW_EPS) ) return 0;
    if( pos[0] > (max[0]+ALW_EPS) ) return 0;
    if( pos[1] > (max[1]+ALW_EPS) ) return 0;
    if( pos[2] > (max[2]+ALW_EPS) ) return 0;

    POLYLIB_STAT ret;
    BBox             bbox;
    bbox.init();
    std::vector<Triangle*> srch_list;
    bool every  = false;    // 3頂点の一部でも検索領域と重なるものを抽出
    
    // 抽出矩形領域の設定
    bbox.min.x = pos[0] - SRCH_DISTANCE;
    bbox.min.y = pos[1] - SRCH_DISTANCE;
    bbox.min.z = pos[2] - SRCH_DISTANCE;

    bbox.max.x = pos[0] + SRCH_DISTANCE;
    bbox.max.y = pos[1] + SRCH_DISTANCE;
    bbox.max.z = pos[2] + SRCH_DISTANCE;
        
    // 指定矩形領域に含まれるポリゴンを抽出する。
    ret = pg.search(
                srch_list,
                bbox, every );
    if( ret != PLSTAT_OK ) {
        PL_ERROSH << "[ERROR] pg.search() failed. returns:"<<PolylibStat2::String(ret)<< endl;
        return 1;
    }
#ifdef DEBUG
    PL_DBGOSH << " srch_list.size()="<<srch_list.size() <<endl;
    for( int i=0; i<srch_list.size(); i++ ) {
        long long int  id_wk   = srch_list[i]->get_id();
        PL_DBGOSH << " srch_list["<<i<<"]->get_id()="<<id_wk <<endl;
    }
#endif

    // 法線ベクトル平均化処理
    PL_REAL vec_wk[3] = { 0.0, 0.0, 0.0 };
 
    // 接点に関連する要素数
    int num_tri = 0;        // 頂点決定対象面数
    int isetted_self_flg = 0;   // 自身の３角形の設定フラグ

    for( int i=0; i<srch_list.size(); i++ ) {
        long long int  id_wk   = srch_list[i]->get_id();
        
        // 自身のポリゴンかどうか判定する
        if( id_wk == id ) {
            isetted_self_flg = 1;

            num_tri++;
            vec_wk[0] += norm_self[0];
            vec_wk[1] += norm_self[1];
            vec_wk[2] += norm_self[2];
        } else {
            // 同一点とみなす頂点があるかどうか判定する
            int same_flg = 0;
            PL_REAL len;
            Vec3<PL_REAL>  *vertex    = srch_list[i]->get_vertexes();
            len = sqrt (   (pos[0]-vertex[0].x)*(pos[0]-vertex[0].x)
                         + (pos[1]-vertex[0].y)*(pos[1]-vertex[0].y)
                         + (pos[2]-vertex[0].z)*(pos[2]-vertex[0].z) );
            if( len < SAME_DISTANCE ) {
                same_flg = 1;
#ifdef DEBUG
                PL_DBGOSH << " set same_flg 1" <<endl;
#endif
            } else {
                len = sqrt (   (pos[0]-vertex[1].x)*(pos[0]-vertex[1].x)
                             + (pos[1]-vertex[1].y)*(pos[1]-vertex[1].y)
                             + (pos[2]-vertex[1].z)*(pos[2]-vertex[1].z) );
                if( len < SAME_DISTANCE ) {
                    same_flg = 1;
#ifdef DEBUG
                    PL_DBGOSH << " set same_flg 2" <<endl;
#endif
                } else {
                    len = sqrt (   (pos[0]-vertex[2].x)*(pos[0]-vertex[2].x)
                                 + (pos[1]-vertex[2].y)*(pos[1]-vertex[2].y)
                                 + (pos[2]-vertex[2].z)*(pos[2]-vertex[2].z) );
                    if( len < SAME_DISTANCE ) {
                        same_flg = 1;
#ifdef DEBUG
                        PL_DBGOSH << " set same_flg 3" <<endl;
#endif
                    }
                }
            }
            if( same_flg == 0 ) { // 頂点を共有していない
                continue;
            }

            // ２面の法線ベクトルよりエッジと判定した場合は
            //  頂点の法線ベクトルを決める対象面から除外する
            Vec3<PL_REAL>  norm_wk = srch_list[i]->get_normal();
            PL_REAL norm_wk2[3];
            norm_wk2[0]=norm_wk.x;  norm_wk2[1]=norm_wk.y;  norm_wk2[2]=norm_wk.z;
#ifdef DEBUG
            PL_DBGOSH << "  srch id="<<id_wk <<endl;
            PL_DBGOSH << "  srch normal="<<norm_wk2[0]<<" "<<norm_wk2[1]<<" "<<norm_wk2[2] <<endl;
            PL_DBGOSH << "  srch vertex[0]="<<vertex[0].x<<" "<<vertex[0].y<<" "<<vertex[0].z <<endl;
            PL_DBGOSH << "  srch vertex[1]="<<vertex[1].x<<" "<<vertex[1].y<<" "<<vertex[1].z <<endl;
            PL_DBGOSH << "  srch vertex[2]="<<vertex[2].x<<" "<<vertex[2].y<<" "<<vertex[2].z <<endl;
#endif

#if 1
                // 2面の角度
            PL_REAL degree = CalcVecAngleDegree( norm_self, norm_wk2);
            if( degree > edge_degree ) {  // エッジでない？
#ifdef DEBUG
                PL_DBGOSH << " reject with degree. degree ="<<degree <<endl;
#endif
            /*
                PL_DBGOSH << "#####  polygon id="<<id<<" degree="<<degree <<endl;
                PL_DBGOSH << "    norm_self="<<norm_self[0]<<" "<<norm_self[1]<<" "<<norm_self[2] <<endl;
                PL_DBGOSH << "    pos="<<pos[0]<<" "<<pos[1]<<" "<<pos[2] <<endl;
                PL_DBGOSH << "    aite id="<<id_wk <<endl;
                PL_DBGOSH << "    aite norm="<<norm_wk2[0]<<" "<<norm_wk2[1]<<" "<<norm_wk2[2] <<endl;
                PL_DBGOSH << "    aite p1="<<vertex[0].x<<" "<<vertex[0].y<<" "<<vertex[0].z <<endl;
                PL_DBGOSH << "    aite p2="<<vertex[1].x<<" "<<vertex[1].y<<" "<<vertex[1].z <<endl;
                PL_DBGOSH << "    aite p3="<<vertex[2].x<<" "<<vertex[2].y<<" "<<vertex[2].z <<endl;
            */
            } else {
#else
            {  // for Test エッジ判定なし
#endif
                num_tri++;
                vec_wk[0] += norm_wk2[0];
                vec_wk[1] += norm_wk2[1];
                vec_wk[2] += norm_wk2[2];
            }
        
        }
    }
#ifdef DEBUG
    PL_DBGOSH << " reference triangle suu : num_tri="<<num_tri <<endl;
#endif

    if( isetted_self_flg == 0 ) {   // 自身が探せていない
        PL_ERROSH << "[ERROR] no vertex  pos=" << pos << endl;
        return 1;
    }

    CalcNormalize2( vec_wk, vertex_norm );  // ベクトルの正規化

    isetted_flg = 1;
    
    return 0;
}
