
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
/// 例題サンプル
///
////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include "Polylib.h"
//   CalcGeo_Matrix.h : original is Npatch library header
#define GEO_REAL PL_REAL
#ifndef PAI
#define PAI 3.14159265358979323846       /* πの値 */
#endif
#include "CalcGeo_Matrix.h"

using namespace PolylibNS;
using namespace std;

/*
 参考：STLポリゴンの初期MinMax
    min= -1009     -1706.9  -1707.65
    max=   338.78   1706.9   1707.65

  風向により回転することを考慮し全体のモデル空間として以下を想定する
    min= -1800.0   -1800.0  -1800.0 
    max=  1800.0    1800.0   1800.0 
*/

/*
 参考：
  struct ParallelBbox {
    PL_REAL bpos[3];         // 基点座標
    unsigned int bbsize[3];  // 計算領域のボクセル数
    unsigned int gcsize[3];  // ガイドセルのボクセル数
    PL_REAL dx[3];           // ボクセル１辺の長さ
  };
*/

//
// 分割情報
//
#if 0
// ランク数2：Y方向2分割  均等分割
static ParallelBbox myParaBbox[2] = {
    { // rank0
        {-1800.0, -1800.0, -1800.0,},   // 基点座標 
        {     36,      18,      36,},   // 計算領域のボクセル数
        {      1,       1,       1,},   // ガイドセルのボクセル数
        {  100.0,   100.0,   100.0 }    // ボクセル１辺の長さ
    },
    { // rank1
        {-1800.0,     0.0, -1800.0,},   // 基点座標 
        {     36,      18,      36,},   // 計算領域のボクセル数
        {      1,       1,       1,},   // ガイドセルのボクセル数
        {  100.0,   100.0,   100.0 }    // ボクセル１辺の長さ
    } 
};
#endif
#if 0
// ランク数4：Y方向2分割, Z方向2分割  均等分割
static ParallelBbox myParaBbox[4] = {
    { // rank0
        {-1800.0, -1800.0, -1800.0,},   // 基点座標 
        {     36,      18,      18,},   // 計算領域のボクセル数
        {      1,       1,       1,},   // ガイドセルのボクセル数
        {  100.0,   100.0,   100.0 }    // ボクセル１辺の長さ
    },
    { // rank1
        {-1800.0,     0.0, -1800.0,},   // 基点座標 
        {     36,      18,      18,},   // 計算領域のボクセル数
        {      1,       1,       1,},   // ガイドセルのボクセル数
        {  100.0,   100.0,   100.0 }    // ボクセル１辺の長さ
    }, 
    { // rank2
        {-1800.0, -1800.0,     0.0,},   // 基点座標 
        {     36,      18,      18,},   // 計算領域のボクセル数
        {      1,       1,       1,},   // ガイドセルのボクセル数
        {  100.0,   100.0,   100.0 }    // ボクセル１辺の長さ
    },
    { // rank3
        {-1800.0,     0.0,     0.0,},   // 基点座標 
        {     36,      18,      18,},   // 計算領域のボクセル数
        {      1,       1,       1,},   // ガイドセルのボクセル数
        {  100.0,   100.0,   100.0 }    // ボクセル１辺の長さ
    } 
};
#endif
#if 0
// ランク数4：Y方向2分割, Z方向2分割
//     Y方向の境目をずらす
//     Z方向 +400  spinnerのポリゴンがないランクを作る
static ParallelBbox myParaBbox[4] = {
    { // rank0
        {-1800.0, -1800.0, -1800.0,},   // 基点座標 
        {     36,      20,      24,},   // 計算領域のボクセル数
        {      1,       1,       1,},   // ガイドセルのボクセル数
        {  100.0,   100.0,   100.0 }    // ボクセル１辺の長さ
    },
    { // rank1
        {-1800.0,   200.0, -1800.0,},   // 基点座標 
        {     36,      16,      24,},   // 計算領域のボクセル数
        {      1,       1,       1,},   // ガイドセルのボクセル数
        {  100.0,   100.0,   100.0 }    // ボクセル１辺の長さ
    }, 
    { // rank2
        {-1800.0, -1800.0,   400.0,},   // 基点座標 
        {     36,      16,      12,},   // 計算領域のボクセル数
        {      1,       1,       1,},   // ガイドセルのボクセル数
        {  100.0,   100.0,   100.0 }    // ボクセル１辺の長さ
    },
    { // rank3
        {-1800.0,  -200.0,   400.0,},   // 基点座標 
        {     36,      20,      12,},   // 計算領域のボクセル数
        {      1,       1,       1,},   // ガイドセルのボクセル数
        {  100.0,   100.0,   100.0 }    // ボクセル１辺の長さ
    } 
};
#endif
#if 1
// ランク数8：Y方向2分割, Z方向2分割, X方向2分割    YZ均等分割   境目合致
static ParallelBbox myParaBbox[8] = {
    { // rank0
        {-1800.0, -1800.0, -1800.0,},   // 基点座標 
        {     17,      18,      18,},   // 計算領域のボクセル数
        {      1,       1,       1,},   // ガイドセルのボクセル数
        {  100.0,   100.0,   100.0 }    // ボクセル１辺の長さ
    },
    { // rank1
        {-1800.0,     0.0, -1800.0,},   // 基点座標 
        {     17,      18,      18,},   // 計算領域のボクセル数
        {      1,       1,       1,},   // ガイドセルのボクセル数
        {  100.0,   100.0,   100.0 }    // ボクセル１辺の長さ
    }, 
    { // rank2
        {-1800.0, -1800.0,     0.0,},   // 基点座標 
        {     17,      18,      18,},   // 計算領域のボクセル数
        {      1,       1,       1,},   // ガイドセルのボクセル数
        {  100.0,   100.0,   100.0 }    // ボクセル１辺の長さ
    },
    { // rank3
        {-1800.0,     0.0,     0.0,},   // 基点座標 
        {     17,      18,      18,},   // 計算領域のボクセル数
        {      1,       1,       1,},   // ガイドセルのボクセル数
        {  100.0,   100.0,   100.0 }    // ボクセル１辺の長さ
    },
    { // rank4
        { -100.0, -1800.0, -1800.0,},   // 基点座標 
        {     19,      18,      18,},   // 計算領域のボクセル数
        {      1,       1,       1,},   // ガイドセルのボクセル数
        {  100.0,   100.0,   100.0 }    // ボクセル１辺の長さ
    },
    { // rank5
        { -100.0,     0.0, -1800.0,},   // 基点座標 
        {     19,      18,      18,},   // 計算領域のボクセル数
        {      1,       1,       1,},   // ガイドセルのボクセル数
        {  100.0,   100.0,   100.0 }    // ボクセル１辺の長さ
    }, 
    { // rank6
        { -100.0, -1800.0,     0.0,},   // 基点座標 
        {     19,      18,      18,},   // 計算領域のボクセル数
        {      1,       1,       1,},   // ガイドセルのボクセル数
        {  100.0,   100.0,   100.0 }    // ボクセル１辺の長さ
    },
    { // rank7
        { -100.0,     0.0,     0.0,},   // 基点座標 
        {     19,      18,      18,},   // 計算領域のボクセル数
        {      1,       1,       1,},   // ガイドセルのボクセル数
        {  100.0,   100.0,   100.0 }    // ボクセル１辺の長さ
    } 
};
#endif


//----------------------------------------------------
//  風車座標系
//      風見方向の回転のみ
//----------------------------------------------------
static PL_REAL windmill_origin[3];    // 風車原点
static PL_REAL windmill_x_axis[3];    // 風車X軸方向ベクトル
                                      //     ロータの回転軸
static PL_REAL windmill_z_axis[3];    // 風車Z軸方向ベクトル
                                      //     風見方向の回転軸

//----------------------------------------------------
//  移動関数（funnel ）
//      風見方向の回転のみ
//----------------------------------------------------

void move_func_funnel(
          PolygonGroup*      pg,
          PolylibMoveParams* params
       )
{
    //PL_DBGOSH << "----- move_func_funnel() start -----" << endl;

    // 風車の回転角度(rad)を求める
    PL_REAL delta_rad_yaw = params->m_params[0];  // 風車の回転角差分(rad)/step
    PL_REAL rad_yaw = delta_rad_yaw * (params->m_next_step - params->m_current_step);

    // 風車の回転マトリックスを求める
    PL_REAL mat_yaw[4][4];
    Calc_3dMat4Rot2(
            windmill_origin,   // [in]  回転軸上の１点の座標
            windmill_z_axis,   // [in]  回転軸のベクトル
            rad_yaw,           // [in]  回転角(rad) 
            mat_yaw            // [out] 回転用４×４変換マトリクス（行ベクトル系）
         );


    std::vector<Triangle* > *tri_list =  pg->get_triangles();
#if 0
    std::string pg_name = pg->get_name();
    PL_DBGOSH << "move polygon group:  name="<<pg_name<<" tri_list->size()="<<tri_list->size() <<endl;
#endif

    // 三角リスト内の全ての三角形について頂点座標を更新する
    //   回転軸を中心として回転させる
    for(int i=0; i<tri_list->size(); i++ ) {
        Vec3<PL_REAL>* vertex = (*tri_list)[i]->get_vertexes();

        for( int j=0; j<3; j++ ) {
            PL_REAL  pos[4],pos_o[4];
            pos[0]=vertex[j].x; pos[1]=vertex[j].y; pos[2]=vertex[j].z; pos[3]=1.0;
            // 座標を回転させる
            Calc_3dMat4Multi14( pos, mat_yaw, pos_o );
            // 座標更新
            vertex[j].x=pos_o[0]; vertex[j].y=pos_o[1]; vertex[j].z=pos_o[2];
            // 法線ベクトル更新・面積非更新
            (*tri_list)[i]->update( true, false );
        }
    }

    // 頂点座標が移動したことにより、KD木の再構築が必要
    // 再構築フラグを立てる
    pg->set_need_rebuild();

    //PL_DBGOSH << "----- move_func_funnel() end   -----" << endl;
}


//----------------------------------------------------
//  移動関数（blades/spinner ）
//      風見方向に回転した後、ブレードの回転
//----------------------------------------------------

void move_func_blades(
          PolygonGroup*      pg,
          PolylibMoveParams* params
       )
{
    //PL_DBGOSH << "----- move_func_blades() start -----" << endl;

    // 風車の回転角度(rad)を求める
    PL_REAL delta_rad_yaw = params->m_params[0];  // 風車の回転角差分(rad)/step
    PL_REAL rad_yaw = delta_rad_yaw * (params->m_next_step - params->m_current_step);


    // ブレードの回転角度(rad)を求める
    PL_REAL rpm = params->m_params[1];  // 回転数/分
    PL_REAL rad_roll = 2.0*PAI*(rpm/60.0)
       * (params->m_next_step - params->m_current_step)*params->m_delta_t;
#if 0
    PL_DBGOSH << "move_func_blades() yaw_degree="<<((rad_yaw/PAI)*180)<<" roll_degree="<<((rad_roll/PAI)*180) <<endl;
#endif

    // 風車の回転マトリックスを求める
    PL_REAL mat_yaw[4][4];
    Calc_3dMat4Rot2(
            windmill_origin,   // [in]  回転軸上の１点の座標
            windmill_z_axis,   // [in]  回転軸のベクトル
            rad_yaw,           // [in]  回転角(rad) 
            mat_yaw            // [out] 回転用４×４変換マトリクス（行ベクトル系）
         );

    // ブレードのローカルの回転マトリックスを求める
    PL_REAL mat_roll[4][4];
    Calc_3dMat4Rot2(
            windmill_origin,  // [in]  回転軸上の１点の座標
            windmill_x_axis,  // [in]  回転軸のベクトル(global)
            rad_roll,         // [in]  回転角(rad) 
            mat_roll          // [out] 回転用４×４変換マトリクス（行ベクトル系）
         );

    // 最終的な回転マトリックスを求める
    PL_REAL mat[4][4];
    Calc_3dMat4Multi44( mat_roll, mat_yaw, mat );


    std::vector<Triangle* > *tri_list =  pg->get_triangles();
#if 0
    std::string pg_name = pg->get_name();
    PL_DBGOSH << "move polygon group:  name="<<pg_name<<" tri_list->size()="<<tri_list->size() <<endl;
#endif

#ifdef DEBUG
    // 頂点が隣接セルよりも遠くへ移動した三角形情報チェック（前処理）
    //  デバッグ用
    pg->init_check_leaped();
#endif

    // 三角リスト内の全ての三角形について頂点座標を更新する
    //   回転軸を中心として回転させる
    for(int i=0; i<tri_list->size(); i++ ) {
        Vec3<PL_REAL>* vertex = (*tri_list)[i]->get_vertexes();

        for( int j=0; j<3; j++ ) {
            PL_REAL  pos[4],pos_o[4];
            pos[0]=vertex[j].x; pos[1]=vertex[j].y; pos[2]=vertex[j].z; pos[3]=1.0;
            // 座標を回転させる
            Calc_3dMat4Multi14( pos, mat, pos_o );  // 行ベクトル系
            // 座標更新
            vertex[j].x=pos_o[0]; vertex[j].y=pos_o[1]; vertex[j].z=pos_o[2];
            // 法線ベクトル更新・面積非更新
            (*tri_list)[i]->update( true, false );
        }
    }

    // 頂点座標が移動したことにより、KD木の再構築が必要
    // 再構築フラグを立てる
    pg->set_need_rebuild();

#ifdef DEBUG
    // 頂点が隣接セルよりも遠くへ移動した三角形情報チェック（後処理）
    //  デバッグ用
    Polylib* p_polylib = Polylib::get_instance();
    ParallelAreaInfo* area_info = p_polylib->get_myproc_area();
    std::vector< Vec3<PL_REAL> > origin;
    std::vector< Vec3<PL_REAL> > cell_size;
    for(int i=0; i<area_info->m_areas.size(); i++ ) {
        origin.push_back( area_info->m_areas[i].m_bpos );
        cell_size.push_back( area_info->m_areas[i].m_dx );
    }
    
    pg->check_leaped( origin,cell_size );
#endif


    //PL_DBGOSH << "----- move_func_blades() end   -----" << endl;
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

    //-------------------------------------------
    //  ファイルパス
    //-------------------------------------------

    std::string config_file_name = "polylib_config.tp";   // 入力：初期化ファイル名

    //-------------------------------------------
    //  初期化
    //-------------------------------------------

    // MPI初期化
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &num_rank );
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    PL_DBGOSH << "---- num_rank="<<num_rank<<" myrank="<<myrank <<endl;


    // Polylib初期化
    Polylib* p_polylib = Polylib::get_instance();

    // 並列計算関連情報の設定と初期化

    //PL_DBGOSH << "(rk:"<<myrank<<") "<<"init_parallel_info() start" << endl;
    ret = p_polylib->init_parallel_info(
                MPI_COMM_WORLD,
                myParaBbox[myrank].bpos,
                myParaBbox[myrank].bbsize,
                myParaBbox[myrank].gcsize,
                myParaBbox[myrank].dx
            );
    //PL_DBGOSH << "init_parallel_info() end" << endl;

    //-------------------------------------------
    //  ロード
    //-------------------------------------------

    // 初期化ファイルを指定してデータロード
    //PL_DBGOSH << "load() start" <<endl;
    ret = p_polylib->load( config_file_name );
    //PL_DBGOSH << "load() end"   <<endl;

#if 1
    // デバッグ用　読み込んだ情報の確認
    p_polylib->show_all_group_info();

    {
        std::string  config_name_out;
        std::string  fmt_out = PolygonIO::FMT_STL_A;
        std::string  extend = "orig";

        //  セーブ（処理前）
        //PL_DBGOSH << "save() start" <<endl;
        ret = p_polylib->save( config_name_out, fmt_out, extend );
        //PL_DBGOSH << "save() end   config_name_out="<<config_name_out <<endl;

        //  分割セーブ（処理前）
        //PL_DBGOSH << "save_parallel() start" <<endl;
        ret = p_polylib->save_parallel( config_name_out, fmt_out, extend );
        //PL_DBGOSH << "save_parallel() end   config_name_out="<<config_name_out <<endl;
    }
#endif

    // ポリゴングループポインタ取得
    std::string   pg_windmill_path = "windmill";
    PolygonGroup* pg_windmill = p_polylib->get_group( pg_windmill_path );
    if( pg_windmill == NULL )  {
        PL_ERROSH <<pg_windmill_path<<": [ERROR] Can't get PolygonGroup pointer" <<endl;
        exit(1);
    }

    std::string   pg_blades_path = "windmill/blades";
    PolygonGroup* pg_blades = p_polylib->get_group( pg_blades_path );
    if( pg_blades == NULL )  {
        PL_ERROSH <<pg_blades_path<<": [ERROR] Can't get PolygonGroup pointer" <<endl;
        exit(1);
    }

    std::string   pg_spinner_path = "windmill/spinner";
    PolygonGroup* pg_spinner = p_polylib->get_group( pg_spinner_path );
    if( pg_spinner == NULL )  {
        PL_ERROSH <<pg_spinner_path<<": [ERROR] Can't get PolygonGroup pointer" <<endl;
        exit(1);
    }

    std::string   pg_funnel_path = "windmill/funnel";
    PolygonGroup* pg_funnel = p_polylib->get_group( pg_funnel_path );
    if( pg_funnel == NULL )  {
        PL_ERROSH <<pg_funnel_path<<": [ERROR] Can't get PolygonGroup pointer" <<endl;
        exit(1);
    }

    // 風車座標系設定
    std::string   key, val;

    key = "center_x";
    pg_windmill->get_atr( key, val );
    windmill_origin[0] = atof( val.c_str() );
    key = "center_y";
    pg_windmill->get_atr( key, val );
    windmill_origin[1] = atof( val.c_str() );
    key = "center_z";
    pg_windmill->get_atr( key, val );
    windmill_origin[2] = atof( val.c_str() );
    key = "yaw_axis_vec_x";
    pg_windmill->get_atr( key, val );
    windmill_z_axis[0] = atof( val.c_str() );
    key = "yaw_axis_vec_y";
    pg_windmill->get_atr( key, val );
    windmill_z_axis[1] = atof( val.c_str() );
    key = "yaw_axis_vec_z";
    pg_windmill->get_atr( key, val );
    windmill_z_axis[2] = atof( val.c_str() );
    key = "roll_axis_vec_x";
    pg_windmill->get_atr( key, val );
    windmill_x_axis[0] = atof( val.c_str() );
    key = "roll_axis_vec_y";
    pg_windmill->get_atr( key, val );
    windmill_x_axis[1] = atof( val.c_str() );
    key = "roll_axis_vec_z";
    pg_windmill->get_atr( key, val );
    windmill_x_axis[2] = atof( val.c_str() );
#if 1
    if( myrank == 0 )  {
        PL_DBGOSH << "windmill_origin ="<<windmill_origin[0]<<" "<<windmill_origin[1]<<" "<<windmill_origin[2] <<endl;
        PL_DBGOSH << "windmill_z_axis ="<<windmill_z_axis[0]<<" "<<windmill_z_axis[1]<<" "<<windmill_z_axis[2] <<endl;
        PL_DBGOSH << "windmill_x_axis ="<<windmill_x_axis[0]<<" "<<windmill_x_axis[1]<<" "<<windmill_x_axis[2] <<endl;
    }
#endif


    // 移動関数登録
    //   bladesとspinner は同一中心軸に対して同じ回転をさせれば良いので
    //   同じ移動関数を登録する

    pg_blades->set_move_func ( move_func_blades );
    pg_spinner->set_move_func( move_func_blades );
    pg_funnel->set_move_func ( move_func_funnel );


    // move parameter 初期化
    PolylibMoveParams params;
    memset( params.m_params, 0x00, 10*sizeof(PL_REAL) ); 

#if 1
    // テスト用の検証値出力（処理前）
    int ntria_blades_local  = pg_blades->get_group_num_tria(); 
    int ntria_blades_global = pg_blades->get_group_num_global_tria(); 
    PL_REAL area_blades_global = pg_blades->get_group_global_area(); 
    if( myrank == 0 )  {
        PL_DBGOSH << "istep= -1"<<"  ntria_blades_local="<<ntria_blades_local<<" ntria_blades_global="<<ntria_blades_global<<" area_blades_global="<<area_blades_global <<endl;
    }
#endif
        
    //---------------------------------
    // タイムステップループ
    //---------------------------------

    //int nstep = 180;
    int nstep = 51;

    for(int istep=0; istep<nstep ; istep++ )  {

        //---------------------------------------------
        // 現在のステップで計算実行
        //---------------------------------------------
        /*
        vector<Triangle*> tri_list;
        p_polylib->search_polygons( tri_list, / * 検索条件を設定 * / );
            
        //解析処理実行	 

        */
                // 風向きの偏向による風車の角度変更を設定する
        PL_REAL delta_rad_yaw = (0.5*PAI/180.0);  // 1stepで0.5度偏向
                                                  // テスト用:実際の変化量としては大きすぎる
                // ブレードの回転速度を求めたものとする
        PL_REAL rpm = 20;                         // 1分間に20回転 （1回転/3秒）


        //---------------------------------------------
        // 次計算ステップに進むためにポリゴン情報更新
        //---------------------------------------------

        // moveパラメタ設定
        params.m_current_step = istep;
        params.m_next_step    = istep + 1;
        params.m_delta_t      = 0.01;            // delta秒
        params.m_params[0]    = delta_rad_yaw;   // 風車の角度偏向
        params.m_params[1]    = rpm;             // 回転速度


        //--------------------
        // move実行
        //--------------------
        ret = p_polylib->move( params ); 


        // 風車座標系の更新
        //     原点とZ軸は変わらないのでX軸のみ更新する
            // 回転角度(rad)
        PL_REAL rad_yaw = delta_rad_yaw * (params.m_next_step - params.m_current_step);
            // 風車の回転マトリックスを求める
        PL_REAL mat_yaw[4][4];
        Calc_3dMat4Rot2(
            windmill_origin,    // [in]  回転軸上の１点の座標
            windmill_z_axis,   // [in]  回転軸のベクトル
            rad_yaw,       // [in]  回転角(rad) 
            mat_yaw        // [out] 回転用４×４変換マトリクス（行ベクトル系）
        );
             // X軸の更新
        PL_REAL  vec_in[4],vec_out[4];
        vec_in[0]=windmill_x_axis[0];
        vec_in[1]=windmill_x_axis[1];
        vec_in[2]=windmill_x_axis[2];
        vec_in[3]=1.0;
        Calc_3dMat4Multi14( vec_in, mat_yaw, vec_out );
        windmill_x_axis[0]=vec_out[0];
        windmill_x_axis[1]=vec_out[1];
        windmill_x_axis[2]=vec_out[2];

        //--------------------
        // migrate実行
        //--------------------
        ret = p_polylib->migrate(); 

#if 1
        //---------------------------------------------
        // テスト用の検証
        //---------------------------------------------

        int ntria_blades_local  = pg_blades->get_group_num_tria(); 
        int ntria_blades_global = pg_blades->get_group_num_global_tria(); 
        PL_REAL area_blades_global = pg_blades->get_group_global_area(); 
        if( myrank == 0 )  {
            PL_DBGOSH << "istep="<<istep<<"  ntria_blades_local="<<ntria_blades_local<<" ntria_blades_global="<<ntria_blades_global<<" area_blades_global="<<area_blades_global <<endl;
        }
        
        // テスト用にファイル出力
        if( istep!=0 && (istep%10)==0 ) { 
            std::string  config_name_out;
            std::string  fmt_out = PolygonIO::FMT_STL_A;
            char buff[16];
            sprintf( buff,"%d",istep );
            std::string  sstep = buff;
            std::string  extend = "istep" + sstep;
            // 各ランク毎に途中経過を出力
            if( myrank == 0 )  {
                PL_DBGOSH << "save_parallel() start istep="<<istep <<endl;
            }
            ret = p_polylib->save_parallel( config_name_out, fmt_out, extend );
            if( myrank == 0 )  {
                PL_DBGOSH << "save_parallel() end   config_name_out="<<config_name_out <<endl;
            }
            // 途中経過を出力
            if( myrank == 0 )  {
                PL_DBGOSH << "save() start istep="<<istep <<endl;
            }
            ret = p_polylib->save( config_name_out, fmt_out, extend );
            if( myrank == 0 )  {
                PL_DBGOSH << "save() end   config_name_out="<<config_name_out <<endl;
            }
        }
#endif
    }

    //-------------------------------------------
    //  セーブ（処理後）
    //-------------------------------------------
    {
        std::string  config_name_out;
        std::string  fmt_out = PolygonIO::FMT_STL_A;
        if( myrank == 0 )  {
            PL_DBGOSH << "save() start" <<endl;
        }
        ret = p_polylib->save( config_name_out, fmt_out );
        if( myrank == 0 )  {
            PL_DBGOSH << "save() end   config_name_out="<<config_name_out <<endl;
        }
    }

    //-------------------------------------------
    //  終了化
    //-------------------------------------------

    // MPI終了化
    //PL_DBGOSH << "MPI_Finalize() start" << endl;
    MPI_Finalize();
    //PL_DBGOSH << "MPI_Finalize() end" << endl;

    return 0;
}
