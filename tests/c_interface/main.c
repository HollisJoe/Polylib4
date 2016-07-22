
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
/// Cインターフェーステスト
///
////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "c_lang/CPolylib.h"

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
//  移動関数
//----------------------------------------------------

void move_func_c (
          PL_GRP_TAG      pg_tag,
          PolylibMoveParamsStruct* params
       )
{
    char pg_name[128];
    int  num_tri;
    int  i,j;
    PL_ELM_TAG* tags_tri;
    PL_REAL vertex[9];

#ifdef DEBUG
    printf( "----- User move_func() start -----\n" );

    polylib_group_get_name( pg_tag, pg_name );
    printf( "  pg_name=%s\n",pg_name );
#endif

    polylib_group_get_triangles( pg_tag, &num_tri, &tags_tri );

    for( i=0; i<num_tri; i++ ) {
        polylib_triangle_get_vertexes( tags_tri[i], vertex );

        PL_REAL x_offset =   1.0;
        PL_REAL y_offset =   1.0;
        PL_REAL z_offset =   0.0;

        for( j=0; j<3; j++ ) {
            vertex[3*j+0] += x_offset;
            vertex[3*j+1] += y_offset;
            vertex[3*j+2] += z_offset;
        }

        polylib_triangle_set_vertexes( tags_tri[i], vertex );
    }

    free(tags_tri);

    // 頂点座標が移動したことにより、KD木の再構築が必要
    // 再構築フラグを立てる
    polylib_group_set_need_rebuild( pg_tag );

#ifdef DEBUG
    printf( "----- User move_func() end    -----\n" );
#endif
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
    char*   config_file_name = "polylib_config.tp"; // 入力：初期化ファイル名
    PL_REAL scale = 1.0;
    int     num_pg;
    PL_GRP_TAG *pg_tags;
    char    pg_name[128];
    int     movable;
    int     num_tri;
    PL_ELM_TAG* tags_tri;
    int     pl_type;
    PL_REAL vertex[9];
    PolylibMoveParamsStruct params;
    char*   p_config_name_out;
    char*   fmt_out = FILE_FMT_STL_A;
    char*   extend  = "";
    int     i;

    //-------------------------------------------
    //  初期化
    //-------------------------------------------

    // MPI初期化
#ifdef MPI_PL
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &num_rank );
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    printf( "---- num_rank=%d  myrank=%d\n",num_rank,myrank);
#else
    num_rank =1;
    myrank   =0;
#endif

    // Polylib初期化
    polylib_instance();

    // 並列計算関連情報の設定と初期化
#ifdef MPI_PL
    printf( "(rk:%d) init_parallel_info() start\n",myrank);
    ret = polylib_init_parallel_info(
                MPI_COMM_WORLD,
                myParaInfos[myrank].bpos,
                myParaInfos[myrank].bbsize,
                myParaInfos[myrank].gcsize,
                myParaInfos[myrank].dx
            );
    printf( "(rk:%d) init_parallel_info() end\n",myrank);
    if( ret != PLSTAT_OK ) {
       printf( "# ERROR :init_parallel_info() ret=%d\n",ret);
       exit(1);
    }
#endif

    //-------------------------------------------
    //  ロード
    //-------------------------------------------

    // 初期化ファイルを指定してデータロード
    printf( "(rk:%d) load() start\n",myrank);
    ret = polylib_load( config_file_name, scale );
    printf( "(rk:%d) load() end ret=%d\n",myrank,ret);
    if( ret != PLSTAT_OK ) {
       printf( "# ERROR :polylib_load() ret=%d\n",ret);
       exit(1);
    }


    ret = polylib_get_root_groups_tags( &num_pg, &pg_tags );
    printf( "(rk:%d) polylib_get_root_groups_tags() num_pg=%d\n",myrank,num_pg);
    if( ret != PLSTAT_OK ) {
       printf( "# ERROR :polylib_get_root_groups_tags() ret=%d\n",ret);
       exit(1);
    }

    for( i=0; i<num_pg; i++ ) {

        polylib_group_get_name( pg_tags[i], pg_name );
        polylib_group_get_movable( pg_tags[i], &movable );
        polylib_group_get_triangles( pg_tags[i], &num_tri, &tags_tri );

        printf( "(rk:%d) polylib_group_get_name  i=%d  pg_name=%s  movable=%d  num_tri=%d\n",
                    myrank,i,pg_name,movable,num_tri );

        pl_type = polylib_triangle_get_pl_type( tags_tri[0] );
        printf( "(rk:%d) pl_type=%d\n",myrank,pl_type );

        //polylib_triangle_get_vertexes( tags_tri[0], vertex );
        //printf( "(rk:%d) vertex p1=%f %f %f\n",myrank,vertex[0],vertex[1],vertex[2] );
        //printf( "(rk:%d) vertex p2=%f %f %f\n",myrank,vertex[3],vertex[4],vertex[5] );
        //printf( "(rk:%d) vertex p3=%f %f %f\n",myrank,vertex[6],vertex[7],vertex[8] );


        // 移動関数登録
        if( movable ) {
            polylib_group_set_move_func_c( pg_tags[i], move_func_c );
        }

        free(tags_tri);
    }

    // move parameter
    memset( params.m_params, 0x00, 10*sizeof(PL_REAL) ); 

    //---------------------------------
    // タイムステップループ
    //---------------------------------
    int nstep = 50;
    int istep;
    for( istep=0; istep<nstep ; istep++ )  {

        // moveパラメタ設定
        params.m_current_step = istep;
        params.m_next_step    = istep + 1;
        params.m_delta_t      = 1.0;

        // move実行
        ret = polylib_move( &params ); 

        // migrate実行
#ifdef MPI_PL
        ret = polylib_migrate(); 
#endif
    }

    //-------------------------------------------
    //  セーブ
    //-------------------------------------------

    printf( "save() start\n" );
    ret = polylib_save( &p_config_name_out, fmt_out, extend );
    printf( "save() end   config_name_out=%s\n", p_config_name_out );
    if( ret != PLSTAT_OK ) {
       printf( "# ERROR: polylib_save() ret=%d\n",ret);
       exit(1);
    }

    //-------------------------------------------
    //  終了化
    //-------------------------------------------
    free(pg_tags);

    // MPI終了化
#ifdef MPI_PL
    printf( "(rk:%d) MPI_Finalize() start\n",myrank);
    MPI_Finalize();
    printf( "(rk:%d) MPI_Finalize() end\n",myrank);
#endif

    if( myrank == 0 )  {
        fprintf(stderr,"------------------------------------------\n");
        fprintf(stderr,"   PASS :  c_interface (Normal End)\n");
        fprintf(stderr,"------------------------------------------\n");
    }


    return 0;
}
