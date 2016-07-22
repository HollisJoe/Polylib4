
/*
 * Polylib - Polygon Management Library
 *
 * Copyright (c) 2010-2011 VCAD System Research Program, RIKEN.
 * All rights reserved.
 *
 * Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 *
 */

////////////////////////////////////////////////////////////////////////////
///
/// ポリゴン情報集約、長田パッチポリゴン作成
///
////////////////////////////////////////////////////////////////////////////

#ifdef MPI_PL

#include "mpi.h"
#include "Polylib.h"
#include "Npt.h"
#include "StlToNpt.h"

//using namespace std;
//using namespace PolylibNS;

//#define USE_UNORDERED_MAP
//   unordered_mapの方が高速ですが環境に依存する場合があるので
//   デフォルトは map を使用するものとします。
//     map           : 2分木検索
//     unordered_map : ハッシュ検索


#include <vector>
#include <algorithm>

#if defined(USE_HASH_MAP)
 #include <hash_map>
#elif defined(USE_UNORDERED_MAP)
 #include <unordered_map>
#else
 #include <map>
#endif

// 長田パッチパラメータ（制御点情報）  通信用等に使用する
//      通信で送信するために、クラスはprimitiveな型に変換する
//
struct Npatch_wk {

    PL_REAL     p1[3];      // 頂点1座標
    PL_REAL     p2[3];      // 頂点2座標
    PL_REAL     p3[3];      // 頂点3座標

    PL_REAL     norm1[3];   // 頂点1の法線ベクトル
    PL_REAL     norm2[3];   // 頂点2の法線ベクトル
    PL_REAL     norm3[3];   // 頂点3の法線ベクトル

    int     setted_flag1;   // 頂点1の法線ベクトルの設定フラグ
    int     setted_flag2;   // 頂点2の法線ベクトルの設定フラグ
    int     setted_flag3;   // 頂点3の法線ベクトルの設定フラグ
    
    long long int id;       // ID（ユニーク）
};



// std::sort用ファンクタ   長田パッチのvectorをソートするのに使用
struct NptLess{
    bool operator()( const NptTriangle *l, const NptTriangle *r ) const
    {
        return l->get_id() < r->get_id();
    }
};



///
/// 各ランクにあるポリゴン情報を集約する
///     ポリゴンの重複を削除する
///     長田パッチを生成する
///
/// @param [in]    num_rank         ランク数
/// @param [in]    myrank           自身のランクID
/// @param [in]    num_tri_rank     各ランク内の３角形の数
/// @param [in]    tri_list_rank    各ランク内の３角形情報
///                                    IDの重複あり
/// @param [in]    p1_vertex_norm   頂点1の法線ベクトル
/// @param [in]    p1_setted_flag   頂点1の法線ベクトルの設定フラグ
/// @param [in]    p2_vertex_norm   頂点2の法線ベクトル
/// @param [in]    p2_setted_flag   頂点2の法線ベクトルの設定フラグ
/// @param [in]    p3_vertex_norm   頂点3の法線ベクトル
/// @param [in]    p3_setted_flag   頂点3の法線ベクトルの設定フラグ
/// @param [in]    num_npt_alloc    長田パッチポリゴン領域獲得数
///                                   （出力長田パッチ数と同一のはず）
/// @param [out]   npt_list         長田パッチポリゴン情報
/// @return リターンコード   =0 正常  !=0 異常
/// @attention
///     出力長田パッチ数とnum_npt_allocが合致しない場合はエラーとする
///     NptTriangleオブジェクトは内部で生成される（使用後deleteが必要）

int
gather_polygons_to_npt(
            int num_rank,
            int myrank,
            int num_tri_rank,
            vector<Triangle*> *tri_list_rank,
            PL_REAL** p1_vertex_norm,
            int*      p1_setted_flag,
            PL_REAL** p2_vertex_norm,
            int*      p2_setted_flag,
            PL_REAL** p3_vertex_norm, 
            int*      p3_setted_flag,
            int num_npt_alloc,
            std::vector<NptTriangle*>&  npt_list
        )
{
#ifdef DEBUG
    PL_DBGOSH << "gather_polygons_to_npt() in. " << endl;
#endif
    POLYLIB_STAT ret;
    int iret;

    //      通信サイズをMPI_CHARで指定すると2GBまでしか送信できないため
    //      MPI_DOUBLE単位で送信する。
    //      MPIの引数なのでsize_tではなくint
    int* nsize_com_double = new int[num_rank]; // 通信サイズ(MPI_DOUBLE単位）
    int* num_tri_ranks    = new int[num_rank];        // 各ランクのポリゴン数

#if defined(USE_HASH_MAP)
    hash_map<long long int,Npatch_wk*> map_pl;      // マップ mapを準備する
                                               //    key(first) : ポリゴンID
                                               //    val(second): Npatch_wkへのポインタ
    pair<hash_map<long long int,Npatch_wk*>::iterator,bool> res_pl;
    hash_map<long long int,Npatch_wk*>::iterator it_pl;
#elif defined(USE_UNORDERED_MAP)
    unordered_map<long long int,Npatch_wk*> map_pl; // マップ mapを準備する
                                               //    key(first) : ポリゴンID
                                               //    val(second): Npatch_wkへのポインタ
    pair<unordered_map<long long int,Npatch_wk*>::iterator,bool> res_pl;
    unordered_map<long long int,Npatch_wk*>::iterator it_pl;
#else
    map<long long int,Npatch_wk*> map_pl;         // マップ mapを準備する
                                               //    key(first) : ポリゴンID
                                               //    val(second): Npatch_wkへのポインタ
    pair<map<long long int,Npatch_wk*>::iterator,bool> res_pl;
    map<long long int,Npatch_wk*>::iterator it_pl;
#endif

    // 各ランクのポリゴン数受信
    iret = MPI_Gather( &num_tri_rank, 1, MPI_INT,
                       num_tri_ranks, 1, MPI_INT,
                       0, MPI_COMM_WORLD
                );

    // 通信サイズ(MPI_DOUBLE単位） 自ランク
    int nsize_com_double_myrank;
    size_t isize = (size_t)(num_tri_rank*sizeof(Npatch_wk));
    if( (isize%sizeof(double)) == 0 ) {
        nsize_com_double_myrank = isize/sizeof(double); 
    } else {
        nsize_com_double_myrank = isize/sizeof(double) + 1;
    }

    // ランク０の処理
    //      ランク０のポリゴンデータより長田パッチポリゴンを生成する
    //      ランク０以外のポリゴンデータを受信する。
    //      長田パッチポリゴンが既に存在した場合、以下の処理を行う。
    //           その３頂点が既に設定されていた場合、受信したデータは捨てる
    //           設定されていない頂点があれば、頂点の法線ベクトルデータを更新する。
    //      長田パッチポリゴンが存在していなかった場合、
    //          受信したポリゴンデータで長田パッチポリゴンを生成し、登録する
    //
    //      ※使用メモリを削減するため、
    //        ポリゴンデータ本体の通信は非同期通信でなく同期通信で行う
    if( myrank == 0 )
    {
        // 通信サイズ(MPI_DOUBLE単位）
        int max_size_com_double = 0;
        for( int i=0; i<num_rank; i++ ) {
            size_t isize = (size_t)(num_tri_ranks[i])*sizeof(Npatch_wk);
            if( (isize%sizeof(double)) == 0 ) {
                nsize_com_double[i] = isize/sizeof(double); 
            } else {
                nsize_com_double[i] = isize/sizeof(double) + 1;
            }
            if( nsize_com_double[i] > max_size_com_double ) {
                max_size_com_double = nsize_com_double[i];
            }
        }

        // ランク０の長田パッチポリゴン生成
        for (int i=0; i<num_tri_rank; i++ ) {
            Npatch_wk  *npt_wk = new Npatch_wk;

            // データ取出し
            Vec3<PL_REAL>* vertex = (*tri_list_rank)[i]->get_vertexes();
            long long int id = (*tri_list_rank)[i]->get_id();
            
            // Npatch_wk設定
            VEC3_3_TO_REAL(vertex,npt_wk->p1,npt_wk->p2,npt_wk->p3);
            npt_wk->norm1[0]  = p1_vertex_norm[i][0];
            npt_wk->norm1[1]  = p1_vertex_norm[i][1];
            npt_wk->norm1[2]  = p1_vertex_norm[i][2];
            npt_wk->norm2[0]  = p2_vertex_norm[i][0];
            npt_wk->norm2[1]  = p2_vertex_norm[i][1];
            npt_wk->norm2[2]  = p2_vertex_norm[i][2];
            npt_wk->norm3[0]  = p3_vertex_norm[i][0];
            npt_wk->norm3[1]  = p3_vertex_norm[i][1];
            npt_wk->norm3[2]  = p3_vertex_norm[i][2];
            npt_wk->setted_flag1 = p1_setted_flag[i];
            npt_wk->setted_flag2 = p2_setted_flag[i];
            npt_wk->setted_flag3 = p3_setted_flag[i];
            npt_wk->id        = id;

            // HASH_MAPに登録
            res_pl = map_pl.insert( pair<long long int,Npatch_wk*>( id, npt_wk ) );
//#ifdef DEBUG
            if( !res_pl.second ) {
                PL_ERROSH << "[ERROR] map insert  id=" << id << endl;
                exit(-1);
            }
//#endif
        }

        // 各ランクよりポリゴンデータ受信
        char*  recv_buff = (char*)malloc( max_size_com_double*sizeof(double) );

        // ランク０以外よりポリゴンデータ受信

        for( int i=1; i<num_rank; i++ ) {
            char*  pbuff = recv_buff;
            MPI_Status status;

            // 各ランクからのポリゴンデータ受信
            //iret = MPI_Recv( recv_buff, (num_tri_ranks[i]*sizeof(Npatch_wk)), MPI_CHAR,
            iret = MPI_Recv( recv_buff, nsize_com_double[i], MPI_DOUBLE,
                             i, 1001, MPI_COMM_WORLD, &status );
            
            // ポリゴンデータ毎の処理
            for( int j=0; j<num_tri_ranks[i]; j++ ) {
                Npatch_wk  *npt_recv = (Npatch_wk*)pbuff;
#ifdef DEBUG
                if( npt_recv->id < 1 ) {
                    PL_DBGOSH << "gather_polygons_to_npt() error rank="<<i<<" j="<<j<<" npt_recv->id="<<npt_recv->id<< endl;
                }
#endif
                
                //
                // 既に登録されているポリゴンかどうか調べる
                //     登録されていれば、更新するかどうか確認する
                //     登録されていなければ、新規に登録
                
                // マップ検索
                it_pl = map_pl.find( npt_recv->id );
                
                // 登録済
                if( it_pl != map_pl.end() ) 
                {
                    Npatch_wk  *npt_exist = it_pl->second;  // 登録済のポインタ
                    
                    // 既に3頂点とも登録済であれば何もしない
                        // 頂点1の処理
                    if( npt_exist->setted_flag1==0 && npt_recv->setted_flag1==1 ) {
                        npt_exist->norm1[0] = npt_recv->norm1[0];
                        npt_exist->norm1[1] = npt_recv->norm1[1];
                        npt_exist->norm1[2] = npt_recv->norm1[2];
                        npt_exist->setted_flag1 = 1;    // 登録済に変更
                    }
                        // 頂点2の処理
                    if( npt_exist->setted_flag2==0 && npt_recv->setted_flag2==1 ) {
                        npt_exist->norm2[0] = npt_recv->norm2[0];
                        npt_exist->norm2[1] = npt_recv->norm2[1];
                        npt_exist->norm2[2] = npt_recv->norm2[2];
                        npt_exist->setted_flag2 = 1;    // 登録済に変更
                    }
                        // 頂点3の処理
                    if( npt_exist->setted_flag3==0 && npt_recv->setted_flag3==1 ) {
                        npt_exist->norm3[0] = npt_recv->norm3[0];
                        npt_exist->norm3[1] = npt_recv->norm3[1];
                        npt_exist->norm3[2] = npt_recv->norm3[2];
                        npt_exist->setted_flag3 = 1;    // 登録済に変更
                    }
                }
                // 未登録
                else 
                {
                    // マップに新規登録
                    Npatch_wk  *npt_wk = new Npatch_wk;
                    // Npatch_wk設定
                    (*npt_wk) = (*npt_recv);
#ifdef DEBUG
                    if( npt_recv->id < 1 ) {
                        PL_DBGOSH << "map_pl.insert() error npt_wk->id="<<npt_wk->id<< endl;
                    }
#endif

                    // HASH_MAPに登録
                    res_pl = map_pl.insert( pair<long long int,Npatch_wk*>( npt_wk->id, npt_wk ) );
//#ifdef DEBUG
                    if( !res_pl.second ) {
                        PL_ERROSH << "[ERROR] map insert  id=" << npt_wk->id << endl;
                     exit(-1);
                    }
//#endif
                }
            
                pbuff += sizeof(Npatch_wk);
            }
        }
    }

    // ランク０以外
    //      ランク０に送信する
    //      送信サイズをMPI_CHARで指定すると2GBまでしか送信できないため
    //      MPI_DOUBLEで送信する。
    else
    {
        Npatch_wk npt_wk_send;
        
        char*  send_buff = (char*)malloc( nsize_com_double_myrank*sizeof(double) );
        char*  pbuff = send_buff;
        
        // 各ランクのポリゴンデータを送信バッファに設定
        for( int i=0; i<num_tri_rank; i++ ) {
            // データ取出し
            Vec3<PL_REAL>* vertex = (*tri_list_rank)[i]->get_vertexes();
            long long int id = (*tri_list_rank)[i]->get_id();
            
            // Npatch_wk設定
            VEC3_TO_REAL(vertex[0],npt_wk_send.p1) 
            VEC3_TO_REAL(vertex[1],npt_wk_send.p2) 
            VEC3_TO_REAL(vertex[2],npt_wk_send.p3) 
            npt_wk_send.norm1[0]  = p1_vertex_norm[i][0];
            npt_wk_send.norm1[1]  = p1_vertex_norm[i][1];
            npt_wk_send.norm1[2]  = p1_vertex_norm[i][2];
            npt_wk_send.norm2[0]  = p2_vertex_norm[i][0];
            npt_wk_send.norm2[1]  = p2_vertex_norm[i][1];
            npt_wk_send.norm2[2]  = p2_vertex_norm[i][2];
            npt_wk_send.norm3[0]  = p3_vertex_norm[i][0];
            npt_wk_send.norm3[1]  = p3_vertex_norm[i][1];
            npt_wk_send.norm3[2]  = p3_vertex_norm[i][2];
            npt_wk_send.setted_flag1 = p1_setted_flag[i];
            npt_wk_send.setted_flag2 = p2_setted_flag[i];
            npt_wk_send.setted_flag3 = p3_setted_flag[i];
            npt_wk_send.id        = id;
#ifdef DEBUG
            if( npt_wk_send.id < 1 ) {
                PL_DBGOSH << "gather_polygons_to_npt() error npt_wk_send.id="<<npt_wk_send.id<< endl;
            }
#endif

            
            // 送信バッファに設定
            memcpy( pbuff, &npt_wk_send, sizeof(Npatch_wk) );
            pbuff += sizeof(Npatch_wk);
        }
        
        // 各ランクのポリゴンデータ送信
        //iret = MPI_Send( send_buff, (num_tri_rank*sizeof(Npatch_wk)), MPI_CHAR,
        iret = MPI_Send( send_buff, nsize_com_double_myrank, MPI_DOUBLE,
                         0, 1001, MPI_COMM_WORLD );
        
        // termination
        free( send_buff );
    }

    //-----------------------------------------
    // HASH_MAPよりvectorに変更する
    //   unordered mapでなく、mapを使っていればsortの必要はないが
    //   上では重複チェック時の高速化のためunorderd mapを使用している
    //-----------------------------------------
#ifdef DEBUG
    PL_DBGOSH << "gather_polygons_to_npt() create NPT" << endl;
#endif
    for (it_pl = map_pl.begin(); it_pl != map_pl.end(); it_pl++) {
        //PL_REAL       p1[3];  // 頂点1座標
        //PL_REAL       p2[3];  // 頂点2座標
        //PL_REAL       p3[3];  // 頂点3座標
        PL_REAL     cp_side1_1[3];
        PL_REAL     cp_side1_2[3];
        PL_REAL     cp_side2_1[3];
        PL_REAL     cp_side2_2[3];
        PL_REAL     cp_side3_1[3];
        PL_REAL     cp_side3_2[3];
        PL_REAL     cp_center [3];
        Vec3<PL_REAL>   vertex[3];
        NpatchParam npt_param;

        Npatch_wk  *npt_wk = it_pl->second;  // ポインタ取得
#ifdef DEBUG
        PL_DBGOSH << "gather_polygons_to_npt() polygon id="<<npt_wk->id<< endl;
        PL_DBGOSH << "    vertex_norm1="<<npt_wk->norm1[0]<<" "<<npt_wk->norm1[1]<<" "<<npt_wk->norm1[2] <<endl;
        PL_DBGOSH << "    vertex_norm2="<<npt_wk->norm2[0]<<" "<<npt_wk->norm2[1]<<" "<<npt_wk->norm2[2] <<endl;
        PL_DBGOSH << "    vertex_norm3="<<npt_wk->norm3[0]<<" "<<npt_wk->norm3[1]<<" "<<npt_wk->norm3[2] <<endl;
#endif

        // 長田パッチパラメータ取得
        iret = npt_param_crt(
                    npt_wk->p1, npt_wk->norm1, 
                    npt_wk->p2, npt_wk->norm2,
                    npt_wk->p3, npt_wk->norm3,
                    cp_side1_1, cp_side1_2,
                    cp_side2_1, cp_side2_2, 
                    cp_side3_1, cp_side3_2, 
                    cp_center
                );
        //if( iret != 1 ) {
        if( iret != 0 ) {
            PL_ERROSH << "[ERROR]  npt_param_crt() iret=" << iret 
                << " p1=" << npt_wk->p1 << " p2=" << npt_wk->p2 << " p3=" << npt_wk->p3 << endl;
            return 1;
        }
        
        
        // 長田パッチ生成のための型変換
        REAL_TO_VEC3_3(npt_wk->p1,npt_wk->p2,npt_wk->p3,vertex);

        REAL_TO_VEC3(cp_side1_1,npt_param.cp_side1_1);
        REAL_TO_VEC3(cp_side1_2,npt_param.cp_side1_2);
        REAL_TO_VEC3(cp_side2_1,npt_param.cp_side2_1);
        REAL_TO_VEC3(cp_side2_2,npt_param.cp_side2_2);
        REAL_TO_VEC3(cp_side3_1,npt_param.cp_side3_1);
        REAL_TO_VEC3(cp_side3_2,npt_param.cp_side3_2);
        REAL_TO_VEC3(cp_center,npt_param.cp_center);
        
        // 長田パッチ生成
        NptTriangle *pNpt = new NptTriangle( vertex,
                                            npt_param,
                                            npt_wk->id    // IDは同一とするため指定する
                                        );
        
        // ベクターに登録
        npt_list.push_back( pNpt );
    }
#ifdef DEBUG
    PL_DBGOSH << "  npt_list.size()="<<npt_list.size()<< endl;
#endif

    // 長田パッチのvectorをIDでソートする
    std::sort( npt_list.begin(), npt_list.end(), NptLess() );   // NptLess():比較用
    
    if( npt_list.size() != num_npt_alloc ) {
        PL_ERROSH << "[ERROR]  output npt_list.size() =" << npt_list.size() 
                  << " estimated size =" << num_npt_alloc  << endl;
        //for(int i=0;i<npt_list.size();i++ ) {
        //    PL_DBGOSH << "  i="<<i<<" polygon_id="<<npt_list[i]->get_id()<< endl;
        //}
        return 1;
    }
    
    //-----------------------------------------
    // 終了化処理
    //-----------------------------------------
    if( myrank == 0 ) {
        for (it_pl = map_pl.begin(); it_pl != map_pl.end(); it_pl++) {
            Npatch_wk  *npt_wk = it_pl->second;  // ポインタ取得
            delete npt_wk;
        }
    }
    delete[] nsize_com_double;
    delete[] num_tri_ranks;

    return 0;
}

#endif
