/*
 * Polylib - Polygon Management Library
 *
 * Copyright (c) 2010-2011 VCAD System Research Program, RIKEN.
 * All rights reserved.
 *
 * Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

// MPI環境のみ有効にする
#ifdef MPI_PL

#include "mpi.h"
#include <vector>
#include "Polylib.h"

using namespace std;
using namespace PolylibNS;


////////////////////////////////////////////////////////////////////////////
/// 
/// クラス:PolygonGroup
/// PolygonGroupクラス MPI環境に関係するもの主に記載する
/// 
////////////////////////////////////////////////////////////////////////////

// std::sort用ファンクタ   ポリゴンのvectorをソートするのに使用
struct TriangleLess{
    bool operator()( const Triangle *l, const Triangle *r ) const
    {
        return l->get_id() < r->get_id();
    }
};

// std::equal用ファンクタ   ポリゴンのvectorのID重複を削除するのに使用
struct TriangleEqual{
    bool operator()( const Triangle *l, const Triangle *r ) const
    {
        return l->get_id() == r->get_id();
    }
};


// 重複ポリゴン削除用の構造体
struct TriangleIndex {
    /// キー
    long long int  id;      // ポリゴン
    int            index;   // vectorの何番目かを示す値
};

// std::sort用ファンクタ   ポリゴンのvectorをソートするのに使用
struct TriangleIndexLess{
    bool operator()( const TriangleIndex& l, const TriangleIndex& r ) const
    {
        return l.id < r.id;
    }
};

// std::equal用ファンクタ   ポリゴンのvectorのID重複を削除するのに使用
struct TriangleIndexEqual{
    bool operator()( const TriangleIndex& l, const TriangleIndex& r ) const
    {
        return l.id == r.id;
    }
};


//--- public ------------------------
// グループ内のポリゴンの要素数(global)を返す
//   全プロセスを集約した要素数（重複ポリゴン分は無視される）
int PolygonGroup::get_group_num_global_tria( void )
{
    POLYLIB_STAT ret;
    // 並列環境の場合、重複ポリゴンの削除が必要のため特殊処理
    int num_tri_global = 0;
    std::vector<long long int> ids;
    std::vector<int>           atrs;

        // 全体としてポリゴンが存在する場合でも
        // ランクによってはポリゴンが存在しない場合がある
    if( m_tri_list != NULL ) {
        ids.reserve ( m_tri_list->size() );
        atrs.reserve( m_tri_list->size() );

        for(int i=0; i<m_tri_list->size(); i++ ) {
            ids.push_back ( (*m_tri_list)[i]->get_id() );
            atrs.push_back( 1 );  // 個数をカウントするだけなので１を設定 
        }
    }

    ret = get_polygons_reduce_sum_atrI( ids, atrs, num_tri_global );
    if( ret != PLSTAT_OK ) {
        return 0;
    }

    return num_tri_global;
}


//--- public ------------------------
// グループ内のポリゴンの面積を積算して返す
//   全プロセスを集約した面積（重複ポリゴン分は無視される）
PL_REAL PolygonGroup::get_group_global_area( void )
{
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::get_group_global_area() in. " << endl;
#endif
    POLYLIB_STAT ret;
    // 並列環境の場合、重複ポリゴンの削除が必要のため特殊処理
    PL_REAL area = 0.0;
    std::vector<long long int> ids;
    std::vector<PL_REAL>       atrs;

        // 全体としてポリゴンが存在する場合でも
        // ランクによってはポリゴンが存在しない場合がある
    if( m_tri_list != NULL ) {
        ids.reserve ( m_tri_list->size() );
        atrs.reserve( m_tri_list->size() );

        for(int i=0; i<m_tri_list->size(); i++ ) {
            ids.push_back (  (*m_tri_list)[i]->get_id() );
            atrs.push_back(  (*m_tri_list)[i]->get_area() );
        }
    }

    ret = get_polygons_reduce_sum_atrR( ids, atrs, area );
    if( ret != PLSTAT_OK ) {
        return 0.0;
    }

    return area;
}

//--- public ------------------------
// グループ内のポリゴン属性（整数）の総和を返す
//   並列版
POLYLIB_STAT PolygonGroup::get_polygons_reduce_atrI(
            PL_OP_TYPE op,
            int        atr_no,
            int&       val
        )
{
    POLYLIB_STAT ret;
    val = 0;
    int val_tmp;

    if ( m_tri_list != NULL && m_tri_list->size() > 0 ) {
        int num_atrI = (*m_tri_list)[0]->get_num_atrI();
        int num_atrR = (*m_tri_list)[0]->get_num_atrR();
        if( atr_no >= num_atrI ) {
            return PLSTAT_ATR_NOT_EXIST;
        }
    }

    if( op == PL_OP_MAX ) {
        val_tmp = PL_INT_MIN;
        if ( m_tri_list != NULL  ) {
            for(int i=0; i<m_tri_list->size(); i++ ) {
                int* pAtr = (*m_tri_list)[i]->get_pAtrI();
                if( pAtr[atr_no] > val_tmp ) {
                   val_tmp = pAtr[atr_no];
                }
            }
        }
        // プロセス間MAX取得
        MPI_Comm comm = Polylib::get_instance()->get_MPI_Comm();
        MPI_Allreduce( &val_tmp, &val, 1, MPI_INT, MPI_MAX, comm );
        if( val == PL_INT_MIN ) {
            return PLSTAT_POLYGON_NOT_EXIST;
        }
    } else if( op == PL_OP_MIN ) {
        val_tmp = PL_INT_MAX;
        if ( m_tri_list != NULL  ) {
            for(int i=0; i<m_tri_list->size(); i++ ) {
                int* pAtr = (*m_tri_list)[i]->get_pAtrI();
                if( pAtr[atr_no] < val_tmp ) {
                   val_tmp = pAtr[atr_no];
                }
            }
        }
        // プロセス間MAX取得
        MPI_Comm comm = Polylib::get_instance()->get_MPI_Comm();
        MPI_Allreduce( &val_tmp, &val, 1, MPI_INT, MPI_MIN, comm );
        if( val == PL_INT_MAX ) {
            return PLSTAT_POLYGON_NOT_EXIST;
        }
    } else if( op == PL_OP_SUM ) {
        // 並列環境の場合、重複ポリゴンの削除が必要のため特殊処理
        std::vector<long long int> ids;
        std::vector<int> atrs;
        if( m_tri_list != NULL ) {
            ids.reserve ( m_tri_list->size() );
            atrs.reserve( m_tri_list->size() );

            for(int i=0; i<m_tri_list->size(); i++ ) {
                int* pAtr = (*m_tri_list)[i]->get_pAtrI();
                ids.push_back (  (*m_tri_list)[i]->get_id() );
                atrs.push_back(  pAtr[atr_no] );
            }
        }
        // 重複削除して総和
        ret = get_polygons_reduce_sum_atrI( ids, atrs, val_tmp );
        if( ret != PLSTAT_OK ) {
            return ret;
        }
        val = val_tmp;

    } else {
        return PLSTAT_NG;
    }

    return PLSTAT_OK;
}

//--- public ------------------------
// グループ内のポリゴン属性（実数）の総和を返す
//   並列版
POLYLIB_STAT PolygonGroup::get_polygons_reduce_atrR(
            PL_OP_TYPE op,
            int        atr_no,
            PL_REAL&   val
        )
{
    POLYLIB_STAT ret;
    val = 0.0;
    PL_REAL val_tmp;
    if( m_tri_list != NULL && m_tri_list->size() > 0 ) {
        int num_atrI = (*m_tri_list)[0]->get_num_atrI();
        int num_atrR = (*m_tri_list)[0]->get_num_atrR();
        if( atr_no >= num_atrR ) {
            return PLSTAT_ATR_NOT_EXIST;
        }
    }

    if( op == PL_OP_MAX ) {
        val_tmp = - PL_REAL_MAX;
        if( m_tri_list != NULL ) {
            for(int i=0; i<m_tri_list->size(); i++ ) {
                PL_REAL* pAtr = (*m_tri_list)[i]->get_pAtrR();
                if( pAtr[atr_no] > val_tmp ) {
                   val_tmp = pAtr[atr_no];
                }
            }
        }
        // プロセス間MAX取得
        MPI_Comm comm = Polylib::get_instance()->get_MPI_Comm();
        MPI_Allreduce( &val_tmp, &val, 1, PL_MPI_REAL, MPI_MAX, comm );
        if( val < (-PL_REAL_MAX+0.1) ) {
            return PLSTAT_POLYGON_NOT_EXIST;
        }
    } else if( op == PL_OP_MIN ) {
        val_tmp = PL_REAL_MAX;
        if( m_tri_list != NULL ) {
            for(int i=0; i<m_tri_list->size(); i++ ) {
                PL_REAL* pAtr = (*m_tri_list)[i]->get_pAtrR();
                if( pAtr[atr_no] < val_tmp ) {
                   val_tmp = pAtr[atr_no];
                }
            }
        }
        // プロセス間MIN取得
        MPI_Comm comm = Polylib::get_instance()->get_MPI_Comm();
        MPI_Allreduce( &val_tmp, &val, 1, PL_MPI_REAL, MPI_MIN, comm );
        if( val > (PL_REAL_MAX-0.1) ) {
            return PLSTAT_POLYGON_NOT_EXIST;
        }
    } else if( op == PL_OP_SUM ) {
        // 並列環境の場合、重複ポリゴンの削除が必要のため特殊処理
        std::vector<long long int> ids;
        std::vector<PL_REAL> atrs;
        if( m_tri_list != NULL ) {
            ids.reserve ( m_tri_list->size() );
            atrs.reserve( m_tri_list->size() );

            for(int i=0; i<m_tri_list->size(); i++ ) {
                PL_REAL* pAtr = (*m_tri_list)[i]->get_pAtrR();
                ids.push_back (  (*m_tri_list)[i]->get_id() );
                atrs.push_back(  pAtr[atr_no] );
            }
        }
        // 重複削除して総和
        ret = get_polygons_reduce_sum_atrR( ids, atrs, val_tmp );
        if( ret != PLSTAT_OK ) {
            return ret;
        }
        val = val_tmp;
    } else {
        return PLSTAT_NG;
    }

    return PLSTAT_OK;
}


//--- private ------------------------
// グループ内のポリゴン属性（整数）の総和を返す
POLYLIB_STAT
PolygonGroup::get_polygons_reduce_sum_atrI(
            std::vector<long long int>& ids,
            std::vector<int>&           atrs,
            int&   sum
        )
{
    int numproc,myrank;
    int iret;

    sum = 0;

    Polylib* p_inst = Polylib::get_instance();
    MPI_Comm comm = p_inst->get_MPI_Comm();

    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &numproc);

    int num_tri_myrank = ids.size();
    int num_reqs = 0;
    int* num_tri_ranks = new int[numproc];
    MPI_Request* reqs  = new MPI_Request[numproc*2];
    MPI_Status*  stats = new MPI_Status[numproc*2];

    // 個数の収集
    iret = MPI_Gather(
                 &num_tri_myrank, 1, MPI_INT,
                 num_tri_ranks,   1, MPI_INT, 0, comm
               );

    // ランク０処理        
    if( myrank == 0 )
    {

        //int index[numproc];  // ランク毎の先頭インデックス
        int* index = new int[numproc];
        int num_tri_tot;

        index[0] = 0;
        num_tri_tot = num_tri_ranks[0];
        for(int i=1; i<numproc; i++ ) {
            index[i] = num_tri_tot;
            num_tri_tot += num_tri_ranks[i];
        }

        if( num_tri_tot > 0 ) { 
            long long int *ids_tot   = new long long int[num_tri_tot];
            int           *atrs_tot  = new int          [num_tri_tot];

            // 自身の設定
            int ip = index[0];
            for(int i=0; i<num_tri_myrank; i++ ) {
                ids_tot [ip] = ids [i];
                atrs_tot[ip] = atrs[i];
                ip++;
            }

            //他からの受信
            for(int irank=1; irank<numproc; irank++)  {
                if( num_tri_ranks[irank] > 0 )  {
                    ip = index[irank];
                    // ポリゴンID受信
                    // MPI_LONG/MPI_LONG_LONG  8byte?  -> use MPI_DOUBLE
                    iret = MPI_Irecv( 
                                 &ids_tot[ip], num_tri_ranks[irank], MPI_DOUBLE, 
                                 irank, 1000, comm, &reqs[num_reqs]
                             );
                    num_reqs++;
                    // 属性受信
                    iret = MPI_Irecv( 
                                 &atrs_tot[ip], num_tri_ranks[irank], MPI_INT, 
                                 irank, 1001, comm, &reqs[num_reqs]
                             );
                    num_reqs++;
                }
            }

            iret = MPI_Waitall( num_reqs, reqs, stats );

            //--------------------- 
            // 重複削除処理
            //--------------------- 
            std::vector<TriangleIndex> tri_index_list;
            tri_index_list.reserve( num_tri_tot );
            // 構造体に変更
            for(int i=0; i<num_tri_tot; i++ ) {
                TriangleIndex tri_index;
                tri_index.id    = ids_tot[i];
                tri_index.index = i;
                tri_index_list.push_back( tri_index ); 
            }

            // ポリゴンIDでソート
            //PL_DBGOSH << "  before sort() tri_index_list.size()="<<tri_index_list.size() <<endl;
            std::sort( tri_index_list.begin(), tri_index_list.end(), TriangleIndexLess() );
            // ID重複分を削除
            tri_index_list.erase(
                   std::unique( tri_index_list.begin(), tri_index_list.end(), TriangleIndexEqual() ),
                   tri_index_list.end()
                );
            //PL_DBGOSH << "  afetr erase() tri_index_list.size()="<<tri_index_list.size() <<endl;

            //--------------------- 
            // 総計
            //--------------------- 
            sum = 0;
            for(int i=0; i<tri_index_list.size(); i++ )  {
                sum += atrs_tot[ tri_index_list[i].index ];
            }

            delete[] ids_tot;
            delete[] atrs_tot;
        }

        delete[] index;
    }
    // ランク０以外処理        
    else
    {
        if( num_tri_myrank > 0 )  {
            // ベクター生データのアドレス取得
            long long int *p_ids  = &ids[0];
            int           *p_atrs = &atrs[0];

            // ポリゴンID送信
            // MPI_LONG/MPI_LONG_LONG  8byte?  -> use MPI_DOUBLE
            iret = MPI_Isend( 
                             p_ids, num_tri_myrank, MPI_DOUBLE, 
                             0, 1000, comm, &reqs[num_reqs]
                     );
            num_reqs++;
            // 属性送信
            iret = MPI_Isend( 
                             p_atrs, num_tri_myrank, MPI_INT, 
                             0, 1001, comm, &reqs[num_reqs]
                         );
            num_reqs++;

            iret = MPI_Waitall( num_reqs, reqs, stats );
        }
    }

    //総和値のbcast
    iret = MPI_Bcast( &sum, 1, MPI_INT, 0, comm );

    //終了化処理
    delete[] num_tri_ranks;
    delete[] reqs;
    delete[] stats;

    return PLSTAT_OK;
}


//--- private ------------------------
// グループ内のポリゴン属性（実数）の総和を返す
POLYLIB_STAT
PolygonGroup::get_polygons_reduce_sum_atrR(
            std::vector<long long int>& ids,
            std::vector<PL_REAL>&       atrs,
            PL_REAL&   sum
        )
{
    int numproc,myrank;
    int iret;

    sum = 0;

    Polylib* p_inst = Polylib::get_instance();
    MPI_Comm comm = p_inst->get_MPI_Comm();

    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &numproc);

    int num_tri_myrank = ids.size();
    int num_reqs = 0;
    int* num_tri_ranks = new int[numproc];
    MPI_Request* reqs  = new MPI_Request[numproc*2];
    MPI_Status*  stats = new MPI_Status [numproc*2];

    // 個数の収集
    iret = MPI_Gather(
                 &num_tri_myrank, 1, MPI_INT,
                 num_tri_ranks,   1, MPI_INT, 0, comm
               );

    // ランク０処理        
    if( myrank == 0 )
    {
        int num_tri_tot;
        //int index[numproc];  // ランク毎の先頭インデックス
        int* index = new int[numproc];  // ランク毎の先頭インデックス

        index[0] = 0;
        num_tri_tot = num_tri_ranks[0];
        for(int i=1; i<numproc; i++ ) {
            index[i] = num_tri_tot;
            num_tri_tot += num_tri_ranks[i];
        }

        if( num_tri_tot > 0 ) { 
            long long int *ids_tot   = new long long int[num_tri_tot];
            PL_REAL       *atrs_tot  = new PL_REAL      [num_tri_tot];

            // 自身の設定
            int ip = index[0];
            for(int i=0; i<num_tri_myrank; i++ ) {
                ids_tot [ip] = ids [i];
                atrs_tot[ip] = atrs[i];
                ip++;
            }

            //他からの受信
            for(int irank=1; irank<numproc; irank++)  {
                if( num_tri_ranks[irank] > 0 )  {
                    ip = index[irank];
                    // ポリゴンID受信
                    // MPI_LONG/MPI_LONG_LONG  8byte?  -> use MPI_DOUBLE
                    iret = MPI_Irecv( 
                                 &ids_tot[ip], num_tri_ranks[irank], MPI_DOUBLE, 
                                 irank, 1000, comm, &reqs[num_reqs]
                             );
                    num_reqs++;
                    // 属性受信
                    iret = MPI_Irecv( 
                                 &atrs_tot[ip], num_tri_ranks[irank], PL_MPI_REAL, 
                                 irank, 1001, comm, &reqs[num_reqs]
                             );
                    num_reqs++;
                }
            }

            iret = MPI_Waitall( num_reqs, reqs, stats );

            //--------------------- 
            // 重複削除処理
            //--------------------- 
            std::vector<TriangleIndex> tri_index_list;
            tri_index_list.reserve( num_tri_tot );
            // 構造体に変更
            for(int i=0; i<num_tri_tot; i++ ) {
                TriangleIndex tri_index;
                tri_index.id    = ids_tot[i];
                tri_index.index = i;
                tri_index_list.push_back( tri_index ); 
            }

            // ポリゴンIDでソート
            //PL_DBGOSH << "  before sort() tri_index_list.size()="<<tri_index_list.size() <<endl;
            std::sort( tri_index_list.begin(), tri_index_list.end(), TriangleIndexLess() );
            // ID重複分を削除
            tri_index_list.erase(
                   std::unique( tri_index_list.begin(), tri_index_list.end(), TriangleIndexEqual() ),
                   tri_index_list.end()
                );
            //PL_DBGOSH << "  afetr erase() tri_index_list.size()="<<tri_index_list.size() <<endl;

            //--------------------- 
            // 総計
            //--------------------- 
            sum = 0;
            for(int i=0; i<tri_index_list.size(); i++ )  {
                sum += atrs_tot[ tri_index_list[i].index ];
            }

            delete[] ids_tot;
            delete[] atrs_tot;
        }

        delete[] index;
    }
    // ランク０以外処理        
    else
    {
        if( num_tri_myrank > 0 )  {
            // ベクター生データのアドレス取得
            long long int *p_ids  = &ids[0];
            PL_REAL       *p_atrs = &atrs[0];

            // ポリゴンID送信
            // MPI_LONG/MPI_LONG_LONG  8byte?  -> use MPI_DOUBLE
            iret = MPI_Isend( 
                             p_ids, num_tri_myrank, MPI_DOUBLE, 
                             0, 1000, comm, &reqs[num_reqs]
                     );
            num_reqs++;
            // 属性送信
            iret = MPI_Isend( 
                             p_atrs, num_tri_myrank, PL_MPI_REAL, 
                             0, 1001, comm, &reqs[num_reqs]
                         );
            num_reqs++;

            iret = MPI_Waitall( num_reqs, reqs, stats );
        }
    }

    //総和値のbcast
    iret = MPI_Bcast( &sum, 1, PL_MPI_REAL, 0, comm );

#ifdef DEBUG
    PL_REAL local_sum = 0.0;
    for(int i=0; i<atrs.size(); i++)  {
        local_sum += atrs[i];
    }
    PL_DBGOSH << "PolygonGroup::get_polygons_reduce_sum_atrR()  global sum="<<sum<<" local sum="<<local_sum<<" local num="<<atrs.size() <<endl;
#endif

    //終了化
    delete[] num_tri_ranks;
    delete[] reqs;
    delete[] stats;

    return PLSTAT_OK;
}


// public //////////////////////////////////////////////////////////////////
POLYLIB_STAT
PolygonGroup::get_inbounded_polygons(
        std::vector<Triangle*>  &tri_list
    )
{
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::get_inbounded_polygons() in. " << endl;
#endif
    ParallelAreaInfo* myproc_area = Polylib::get_instance()->get_myproc_area();
    POLYLIB_STAT ret = PLSTAT_OK;

    tri_list.clear();

    // ポリゴン情報を持つグループだけ
    if( m_tri_list != NULL && m_tri_list->size() != 0 ) {
        vector<BBox> bboxes;
        bboxes.clear();

        // 自領域内(複数）に一部でも含まれるポリゴンを検索
        for(int i=0; i<myproc_area->m_areas.size(); i++ ) {
            bboxes.push_back( myproc_area->m_areas[i].m_gcell_bbox );
        }
        bool duplicate = false; // 重複削除

        ret = search( tri_list, bboxes, false, duplicate );
    }

    return ret;
}

// public //////////////////////////////////////////////////////////////////
POLYLIB_STAT
PolygonGroup::erase_outbounded_polygons( void )
{
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::erase_outbounded_polygons() in. " << endl;
#endif
    POLYLIB_STAT ret;

    vector<Triangle*> tri_list;

    ret = get_inbounded_polygons( tri_list );

    //if( tri_list.size()>0 ) {
    if( ret == PLSTAT_OK ) {
        // initの中で最初に既存のポリゴンはポインタを指定して削除される。
        // そのため、別途、外側でnewしないと元のポリゴンを参照できなくなる。
        // initの中でさらにポリゴンのコピーが作成されるため
        // 結局、copy_tri_list内のポリゴンは削除する必要がある。
        //
        vector<Triangle*> copy_tri_list;

        copy_polygons( tri_list, copy_tri_list );

        // 検索結果でポリゴン情報を再構築
        if( (ret = init( &copy_tri_list, true )) != PLSTAT_OK ) {
            PL_ERROSH << "[ERROR]PolygonGroup::erase_outbounded_polygons():pg->init() failed. returns:"
                      << PolylibStat2::String(ret) << endl;
            return ret;
        }

        // copyポリゴン削除
        for(int i=0; i<copy_tri_list.size(); i++ ) {
            delete copy_tri_list[i];
        }
    }

    return ret;
}


// 他ランクにポリゴンデータを分散する
//      リーフのグループしか呼ばれない
//      （非メモリ削減版）
POLYLIB_STAT
PolygonGroup::scatter_polygons( void )
{
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::scatter_polygons() in. " << endl;
#endif
    int myrank = Polylib::get_instance()->get_MPI_myrank();
    POLYLIB_STAT ret;

    if( myrank == 0 ) {
        // ランク０以外に該当ポリゴンを送信する
        ret = scatter_send_polygons();
        
        // 自担当領域分のポリゴン情報とする
        //   自分の担当領域でないポリゴンは削除する
        ret = erase_outbounded_polygons();
    } else {
        std::vector<Triangle*>     tri_list;

        // ポリゴン情報をランク０より受信する
        ret = scatter_receive_polygons( tri_list );

        init( &tri_list, true );
        // pg->init()でcopyされるので元のポリゴン削除
        for(int i=0; i<tri_list.size(); i++) {
            delete tri_list[i];
        }
    }

    return PLSTAT_OK;
}


// 他ランクにポリゴンデータを分散する
//      リーフのグループしか呼ばれない
//      （メモリ削減版）
POLYLIB_STAT
PolygonGroup::scatter_polygons(
            std::vector<Triangle*>  &tri_list_div_all,
            std::vector<Triangle*>  &tri_list_div_local
    )
{
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::scatter_polygons() in. tri_list_div_all.size()="<<tri_list_div_all.size() << endl;
#endif
    int myrank = Polylib::get_instance()->get_MPI_myrank();
    POLYLIB_STAT ret;

    if( tri_list_div_local.size() > 0 ) {
#ifdef DEBUG
        PL_DBGOSH << "   tri_list_div_local is not empty." <<endl;
#endif
        return PLSTAT_NG;
    }

    // ランク０
    if( myrank == 0 ) {
#ifdef DEBUG
        if( m_tri_list == NULL )  {
            PL_DBGOSH << "   m_tri_list is NULL. OK" <<endl;
        } else {
            if( m_tri_list->size() > 0 )  {
               // この処理の前でポリゴンが登録されているのはおかしい
               PL_DBGOSH << "#### ERROR m_tri_list->size()="<<m_tri_list->size() <<endl;
            }
        }
#endif

        // ポリゴンを登録する(一時的に使用するため）
        //     copyを防ぐため、ポインタを登録する
        std::vector<Triangle*>  *tri_list_tmp = new std::vector<Triangle*>;
        *tri_list_tmp = tri_list_div_all;  // 要素のcopy

        // ポリゴンのポインタの登録・KDツリー等の構築
        set_triangles_ptr( tri_list_tmp, true );

        // ランク０以外に該当ポリゴンを送信する
        ret = scatter_send_polygons();
        
        // 自担当領域分のポリゴン情報取得する
        std::vector<Triangle*>  tri_list;
        ret = get_inbounded_polygons( tri_list );
#ifdef DEBUG
        PL_DBGOSH << "  get_inbounded_polygons() tri_list.size()="<<tri_list.size() <<endl;
#endif

        // ポリゴンを出力域にコピーする
        copy_polygons( tri_list, tri_list_div_local );
#ifdef DEBUG
        PL_DBGOSH << "  copy_polygons() tri_list_div_local.size()="<<tri_list_div_local.size() <<endl;
#endif

        // 一時的に登録したポリゴンを削除する
        delete_tri_list();
        // 上記のdelete_tri_list()でtri_list_div_all内の
        // ポリゴンは削除されているため、リストをクリア
        tri_list_div_all.clear();

    }
    // ランク０以外
    else
    {
        // ポリゴン情報をランク０より受信する
        ret = scatter_receive_polygons( tri_list_div_local );
    }

    return PLSTAT_OK;
}


// scatter ランク０の処理 ランク０以外に該当ポリゴンを送信する
//    rank0のみで呼ばれる
//      通信サイズをMPI_CHARで指定すると2GBまでしか送信できないため
//      MPI_DOUBLE単位で送信する。
POLYLIB_STAT
PolygonGroup::scatter_send_polygons( void )
{
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::scatter_send_polygons() in. " << endl;
#endif
    Polylib* p_inst = Polylib::get_instance();
    MPI_Comm comm   = p_inst->get_MPI_Comm();
    std::vector<ParallelAreaInfo>* p_other_procs_area = p_inst->get_other_procs_area();
    std::vector<Triangle*>     tri_list;

    // その他のランク数ループ
#ifdef DEBUG
    PL_DBGOSH << "p_other_procs_area->size()="<<p_other_procs_area->size()<< endl;
    if( m_tri_list == NULL ) {
        PL_DBGOSH << "m_tri_list is NULL"<< endl;
    } else {
        PL_DBGOSH << "m_tri_list->size()="<<m_tri_list->size()<< endl;
    } 
    for( int i=0; i<p_other_procs_area->size(); i++ ) {
        int irank_send = (*p_other_procs_area)[i].m_rank;
        PL_DBGOSH << "m_other_procs_area:irank_send="<<irank_send<< endl;
        for( int j=0; j<(*p_other_procs_area)[i].m_areas.size(); j++ ) {   // 複数担当領域
            BBox* p_bbox = &( (*p_other_procs_area)[i].m_areas[j].m_gcell_bbox);
            Vec3<PL_REAL>* bbox_min = &(p_bbox->min);
            Vec3<PL_REAL>* bbox_max = &(p_bbox->max);
            PL_DBGOSH << "  gcell_bbox min="<<bbox_min->x<<" "<<bbox_min->y<<" "<<bbox_min->z<< endl;
            PL_DBGOSH << "  gcell_bbox max="<<bbox_max->x<<" "<<bbox_max->y<<" "<<bbox_max->z<< endl;
        }
    }
#endif
    for( int i=0; i<p_other_procs_area->size(); i++ ) {
        int irank_send = (*p_other_procs_area)[i].m_rank;
        int num_tri  = 0 ;  // トライアングル数
        int pl_type  = PL_TYPE_UNKNOWN;
        int numAtrI  = 0;
        int numAtrR  = 0;
        int nsize_pl = 0;           // 1ポリゴンのサイズ
        int nsize_com_double = 0;   // 全体の通信サイズ(MPI_DOUBLE単位）
        tri_list.clear();
            
        
        // ポリゴン情報を持つグループだけ
        if( m_tri_list != NULL && m_tri_list->size() > 0 ) {
            for( int j=0; j<(*p_other_procs_area)[i].m_areas.size(); j++ ) {   // 複数担当領域
#ifdef DEBUG
    PL_DBGOSH << "pg->search() start" << endl;
#endif
                search( tri_list, (*p_other_procs_area)[i].m_areas[j].m_gcell_bbox, false );
#ifdef DEBUG
    PL_DBGOSH << "pg->search() end" << endl;
#endif
            }
            // ID重複削除
            //      IDでソート
#ifdef DEBUG
    PL_DBGOSH << "std::sort() start" << endl;
#endif
            std::sort( tri_list.begin(), tri_list.end(), TriangleLess() );  // TriangleLess():比較用
            //      ID重複分を削除
#ifdef DEBUG
    PL_DBGOSH << "tri_list.erase() start" << endl;
#endif
            tri_list.erase(
                            std::unique( tri_list.begin(), tri_list.end(), TriangleEqual() ), // TriangleEqual():比較用
                            tri_list.end()
                        );
            
            if( tri_list.size() > 0 ) {
                num_tri  = tri_list.size();
                pl_type = tri_list[0]->get_pl_type();
                numAtrI = tri_list[0]->get_num_atrI();
                numAtrR = tri_list[0]->get_num_atrR();
                nsize_pl = tri_list[0]->serialized_size();   // 1ポリゴンのバッファサイズ
            }
            // 通信サイズ(MPI_DOUBLE単位）
            size_t isize = (size_t)(num_tri)* (size_t)(nsize_pl);
            if( (isize%sizeof(double)) == 0 ) {
                nsize_com_double = isize/sizeof(double);    
            } else {
                nsize_com_double = isize/sizeof(double) + 1;
            }
        }


        // ヘッダ情報の送信
        //    ポリゴン数、データタイプ(Triangle/NptTriangle)、
        //    ユーザ定義属性(整数型）数、ユーザ定義属性(実数型）数  ポリゴングループ内ではユーザ定義属性数は同じ前提あり
        //    1ポリゴンデータのサイズ(byte)
        int idata[6];
        idata[0] = num_tri;
        idata[1] = pl_type;
        idata[2] = numAtrI;
        idata[3] = numAtrR;
        idata[4] = nsize_pl;
        idata[5] = nsize_com_double;
        int iret;
        
#ifdef DEBUG
    PL_DBGOSH <<"irank_send="<<irank_send<<" num_tri="<<num_tri<<" pl_type="<<pl_type<<" numAtrI="<<numAtrI<<" numAtrR="<<numAtrR<<" nsize_pl="<<nsize_pl<<" nsize_com_double="<<nsize_com_double<<endl;
#endif
        iret = MPI_Send( idata, 6, MPI_INT,
                     irank_send, 1001, comm );
    
        // ポリゴン情報を送信
        if( num_tri > 0 ) {
            // バッファに詰める
            char*  send_buff = (char*)malloc( (size_t)(nsize_com_double)*sizeof(double) );
            char*  pbuff = send_buff;
            for(int j=0; j<num_tri; j++ ) {
                pbuff = tri_list[j]->serialize( pbuff );
            }
            
            // 送信
            iret = MPI_Send( send_buff, nsize_com_double, MPI_DOUBLE,
                     irank_send, 1002, comm );
            
            free( send_buff );
        }
    }

    return PLSTAT_OK;
}


// scatter ランク０以外の処理 ランク０より該当ポリゴンを受信する
//    rank0以外でのみで呼ばれる
POLYLIB_STAT
PolygonGroup::scatter_receive_polygons(
            std::vector<Triangle*>  &tri_list
    )
{
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::scatter_receive_polygons() in. " << endl;
#endif
    Polylib* p_inst = Polylib::get_instance();
    MPI_Comm comm   = p_inst->get_MPI_Comm();
    MPI_Status status;
    int idata[6];
    int iret;

    // ヘッダ受信
    iret = MPI_Recv( idata, 6, MPI_INT, 
                     0, 1001, comm, &status );

    int num_tri  = idata[0];    // トライアングル数
    int pl_type  = idata[1];
    int numAtrI  = idata[2];
    int numAtrR  = idata[3];
    int nsize_pl = idata[4];            // 1ポリゴンのサイズ
    int nsize_com_double = idata[5];    // 全体の通信サイズ(MPI_DOUBLE単位）
#ifdef DEBUG
    PL_DBGOSH << "receive  "<<" num_tri="<<num_tri<<" pl_type="<<pl_type<<" numAtrI="<<numAtrI<<" numAtrR="<<numAtrR<<" nsize_pl="<<nsize_pl<<" nsize_com_double="<<nsize_com_double<<endl;
#endif
    
    if( num_tri == 0 ) {
        return PLSTAT_OK;
    }

    // ランク0より受信
    char*  recv_buff = (char*)malloc( (size_t)(nsize_com_double)*sizeof(double) );
    char*  pbuff = recv_buff;

    iret = MPI_Recv( recv_buff, nsize_com_double, MPI_DOUBLE,
                             0, 1002, comm, &status );

#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::scatter_receive_polygons() deserialize start" << endl;
#endif
    for(int j=0; j<num_tri; j++ ) {
        // deserialize: Triangle/NptTrinangleオブジェクト生成
        Triangle* pTri = deserialize_polygon( pl_type, pbuff );
        pbuff += nsize_pl;
        tri_list.push_back( pTri );
    }
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::scatter_receive_polygons() deserialize end" << endl;
#endif
    
    // ポリゴングループにポリゴンリスト登録、KD木構築
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::scatter_receive_polygons() tri_list.size()="<<tri_list.size()<< endl;
#endif

    // 終了化
    free(recv_buff);

    return PLSTAT_OK;
}


// 分割データをランク０に集約する
//      リーフ（ポリゴンあり）のグループしか呼ばれない
POLYLIB_STAT
PolygonGroup::gather_polygons(
            std::vector<Triangle*>  &tri_list
    )
{
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::gather_polygons() in. " << endl;
    PL_DBGOSH << "   get_name()="<<get_name() <<endl;
#endif
    POLYLIB_STAT ret;
    int myrank = Polylib::get_instance()->get_MPI_myrank();

    if( myrank == 0 ) {
        // 分割ポリゴン情報を他ランクより受信する
        // Rank0は直接設定する
        ret = gather_recv_polygons( tri_list );

    } else {
        // ランク０に分割ポリゴンを送信する
        ret = gather_send_polygons();

    }

    return PLSTAT_OK;
}


// gather ランク０以外の処理 ランク０に該当ポリゴンを送信する
//    rank0以外のみで呼ばれる
//      通信サイズをMPI_CHARで指定すると2GBまでしか送信できないため
//      MPI_DOUBLE単位で送信する。
POLYLIB_STAT
PolygonGroup::gather_send_polygons( void )
{
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::gather_send_polygons() in. " << endl;
#endif
    Polylib* p_inst = Polylib::get_instance();
    MPI_Comm comm   = p_inst->get_MPI_Comm();
    int myrank      = p_inst->get_MPI_myrank();
    int numproc     = p_inst->get_MPI_numproc();
    //PL_DBGOSH << "PolygonGroup::gather_send_polygons() myrank="<<myrank<<" numproc="<<numproc << endl;
    int iret;
    int num_tri  = 0 ;  // トライアングル数
    int pl_type  = PL_TYPE_UNKNOWN;
    int numAtrI  = 0;
    int numAtrR  = 0;
    int nsize_pl = 0;           // 1ポリゴンのサイズ
    int nsize_com_double = 0;   // 全体の通信サイズ(MPI_DOUBLE単位）

    std::vector<Triangle*>* tri_list = m_tri_list;
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::gather_send_polygons() 0" << endl;
    if( tri_list != NULL ) {
        PL_DBGOSH << "PolygonGroup::gather_send_polygons() tri_list is not NULL" << endl;
        PL_DBGOSH << "PolygonGroup::gather_send_polygons() m_tri_list="<<m_tri_list <<endl;
        PL_DBGOSH << "PolygonGroup::gather_send_polygons() tri_list->size()="<<tri_list->size() <<endl;
    }
#endif
    if( tri_list!=NULL && tri_list->size()>0 ) {
        num_tri = tri_list->size();    // トライアングル数
        pl_type = (*tri_list)[0]->get_pl_type();
        numAtrI = (*tri_list)[0]->get_num_atrI();
        numAtrR = (*tri_list)[0]->get_num_atrR();
        // 1ポリゴンのサイズ
        nsize_pl =  (*tri_list)[0]->serialized_size();
        // 通信サイズ(MPI_DOUBLE単位）
        size_t isize = (size_t)(num_tri)* (size_t)(nsize_pl);
        if( (isize%sizeof(double)) == 0 ) {
            nsize_com_double = isize/sizeof(double);    
        } else {
            nsize_com_double = isize/sizeof(double) + 1;
        }
    }
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::gather_send_polygons() 1" << endl;
    PL_DBGOSH << "  num_tri="<<num_tri<<" pl_type="<<pl_type<<" numAtrI="<<numAtrI<<" numAtrR="<<numAtrR<<" nsize_pl="<<nsize_pl<<" nsize_com_double="<<nsize_com_double <<endl;
#endif

    // ポリゴン情報を送信
    //    メモリ使用量を抑えるため、１ランクづつ同期を取って処理を行う
    for( int irank=1; irank<numproc; irank++ ) {
        if( irank == myrank ) {

            // ヘッダ情報の送信
            //    ポリゴン数、データタイプ(Triangle/NptTriangle)、
            //    ユーザ定義属性(整数型）数、ユーザ定義属性(実数型）数  ポリゴングループ内ではユーザ定義属性数は同じ前提あり
            //    1ポリゴンデータのサイズ(byte)
            int idata[6];
            idata[0] = num_tri;
            idata[1] = pl_type;
            idata[2] = numAtrI;
            idata[3] = numAtrR;
            idata[4] = nsize_pl;
            idata[5] = nsize_com_double;
    
            iret = MPI_Send( idata, 6, MPI_INT,
                 0, 1001, comm );

            // ポリゴン情報を送信
            if( num_tri > 0 ) {
                // バッファに詰める
                char*  send_buff = (char*)malloc( (size_t)(nsize_com_double)*sizeof(double) );
                char*  pbuff = send_buff;
                for(int i=0; i<num_tri; i++ ) {
                    pbuff = (*tri_list)[i]->serialize( pbuff );
                }
        
                // 送信
                iret = MPI_Send( send_buff, nsize_com_double, MPI_DOUBLE,
                               0, 1002, comm );
        
                free( send_buff );
            }
        }
        
        MPI_Barrier( comm );  // 1ランクずつ同期を取る
    }

#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::gather_send_polygons() out. " << endl;
#endif

    return PLSTAT_OK;
}


// gather ランク０の処理 ランク０以外より該当ポリゴンを受信する
//    rank0でのみで呼ばれる
POLYLIB_STAT
PolygonGroup::gather_recv_polygons(
            std::vector<Triangle*>  &tri_list
    )
{
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::gather_recv_polygons() in. " << endl;
#endif
    Polylib* p_inst = Polylib::get_instance();
    MPI_Comm comm   = p_inst->get_MPI_Comm();
    int numproc     = p_inst->get_MPI_numproc();
    int iret;
    MPI_Status status;
    int idata[6];
    std::vector<Triangle*>  tri_list_tmp;

    // 自ランク(Rank0)の処理
    std::vector<Triangle*>* tri_list_rank0 = m_tri_list;
    if( tri_list_rank0!=NULL && tri_list_rank0->size()>0 ) {
        // ポリゴンの複製
        copy_polygons( *tri_list_rank0, tri_list_tmp );
    }

    // 他ランクから受信
    for( int irank=1; irank<numproc; irank++ ) {

        // ヘッダ受信
        iret = MPI_Recv( idata, 6, MPI_INT, 
                     irank, 1001, comm, &status );

        int num_tri  = idata[0];    // トライアングル数
        int pl_type  = idata[1];
        int numAtrI  = idata[2];
        int numAtrR  = idata[3];
        int nsize_pl = idata[4];
        int nsize_com_double = idata[5];

#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::gather_recv_polygons() 1 irank="<<irank <<endl;
    PL_DBGOSH << "  num_tri="<<num_tri<<" pl_type="<<pl_type<<" numAtrI="<<numAtrI<<" numAtrR="<<numAtrR<<" nsize_pl="<<nsize_pl<<" nsize_com_double="<<nsize_com_double <<endl;
#endif
    
        if( num_tri >  0 ) {
            // 受信
            char*  recv_buff = (char*)malloc( (size_t)(nsize_com_double)*sizeof(double) );
            char*  pbuff = recv_buff;

            iret = MPI_Recv( recv_buff, nsize_com_double, MPI_DOUBLE,
                             irank, 1002, comm, &status );

            for(int j=0; j<num_tri; j++ ) {
                // deserialize: Triangle/NptTrinangleオブジェクト生成
                Triangle* pTri = deserialize_polygon( pl_type, pbuff );
                pbuff += nsize_pl;
                tri_list_tmp.push_back( pTri );
            }
            
            free(recv_buff);
        }
        
        MPI_Barrier( comm );  // 1ランクずつ同期を取る
    }
    
    if( tri_list_tmp.size() > 0 )  {
        // IDでソート
        //     ポリゴンの順番を入力と同一とする（ユーザ定義属性とのペアで考える必要があるため）
        //PL_DBGOSH << "  before sort() tri_list_tmp.size()="<<tri_list_tmp.size() <<endl;
        std::sort( tri_list_tmp.begin(), tri_list_tmp.end(), TriangleLess() );  // TriangleLess():比較用
        //PL_DBGOSH << "  after  sort() tri_list_tmp.size()="<<tri_list_tmp.size() <<endl;
        // ID重複分を削除
        long long int id_pre = tri_list_tmp[0]->get_id();
        tri_list.push_back( tri_list_tmp[0] );

        for(int i=1; i<tri_list_tmp.size(); i++ ) {
            long long int id = tri_list_tmp[i]->get_id();
            if( id == id_pre ) {
                // IDが同一なので削除
                delete tri_list_tmp[i];
            } else {
                // 新しいIDなのでリストに登録
                tri_list.push_back( tri_list_tmp[i] );
                id_pre = id;
            }
        }

    }

#ifdef DEBUG
    PL_DBGOSH << "  after remove duplicate id() tri_list.size()="<<tri_list.size() <<endl;
    PL_DBGOSH << "PolygonGroup::gather_recv_polygons() out. " << endl;
#endif

    return PLSTAT_OK;
}


// protected ////////////////////////////////////////////////////////////////
// STL/NPTファイルの読み込み(MPI メモリ削減版）
POLYLIB_STAT
PolygonGroup::load_polygons_mem_reduced(
            PL_REAL      scale,
            int          size_mb
    )
{
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::load_polygons_mem_reduced() in. " << endl;
#endif
    Polylib* p_inst = Polylib::get_instance();
    MPI_Comm comm = p_inst->get_MPI_Comm();
    int myrank    = p_inst->get_MPI_myrank();
    POLYLIB_STAT ret;
    int iend = 0;


    // ランク0
    if( myrank == 0 )
    {
        PL_REAL size_mem_init = size_mb*1024*1024; // ロード可能メモリサイズ（実数で表現）
        PL_REAL size_mem_cur;   // 現在の残メモリサイズ（実数で表現）
        std::vector<Triangle*>  *tri_list_regist = new std::vector<Triangle*>;

        // ファイル数ループ
        //   １ポリゴングループに対して複数ファイルがあることあり
        map<string, string>::const_iterator it;
        for (it = m_polygon_files.begin(); it != m_polygon_files.end(); it++) {
            string fname    = it->first;
            string fmt      = it->second;

            // ファイルフォーマットよりTriangle/NptTriangleを判断し、
            //   １ポリゴンのサイズを求める(loadは属性関係なし）
            int pl_type = PolygonIO::get_polygon_type( fmt );
            size_t pl_size;
            if( pl_type == PL_TYPE_TRIANGLE ) {
                Triangle tri_tmp;
                pl_size = tri_tmp.used_memory_size();
            } else if( pl_type == PL_TYPE_NPT ) {
                NptTriangle tri_tmp;
                pl_size = tri_tmp.used_memory_size();
            } else {
                return PLSTAT_NG;
            }
            PL_REAL pl_size_r = pl_size;
#ifdef DEBUG
            PL_DBGOSH << "PolygonGroup::load_polygons_mem_reduced() pl_size="<<pl_size <<endl;
#endif

            // ファイルOpen
            ifstream ifs;
            ret = PolygonIO::load_file_open( ifs, fname, fmt ); 

            bool eof = false;
            // ファイル終了まで読み込む
            while ( !eof )  {
                std::vector<Triangle*>  tri_list_div_all;
                std::vector<Triangle*>  tri_list_div_local;

                // 残メモリサイズ
                size_mem_cur = size_mem_init - pl_size_r*tri_list_regist->size();
                // 読み込み可能数
                int num_read = size_mem_cur / pl_size_r;

#ifdef DEBUG
                PL_DBGOSH << "   tri_list_regist->size()="<<tri_list_regist->size()<<" size_mem_cur="<<size_mem_cur<<" num_read="<<num_read <<endl;
#endif

                if( num_read <= 0 ) {
                    ret = PLSTAT_LACK_OF_LOAD_MEMORY; // メモリ不足によりロード出来ない 
                    PL_ERROSH << PolylibStat2::String(ret) <<endl;
                    return ret;
                }
            
                // 指定ポリゴン数読み込み
                int num_tri; // 実際に読み込んだ数
                ret = PolygonIO::load_file_read( ifs, fmt, 
                          tri_list_div_all, num_read, num_tri, eof, scale );
#ifdef DEBUG
                PL_DBGOSH << "   load_file_read() tri_list_div_all.size()="<<tri_list_div_all.size()<<" num_read="<<num_read<<" num_tri="<<num_tri <<endl;
#endif

                if( tri_list_div_all.size() > 0 )  {
                    // 終了判定フラグ==0(終了なし） 送信
                    iend = 0;
                    MPI_Bcast( &iend, 1, MPI_INT, 0, comm );

                    // ポリゴンのスキャタ
                    ret = scatter_polygons( tri_list_div_all, tri_list_div_local );
#ifdef DEBUG
                PL_DBGOSH << "   scatter_polygons() tri_list_div_all.size()="<<tri_list_div_all.size()<<" tri_list_div_local.size()="<<tri_list_div_local.size() <<endl;
#endif

                    // 登録するポリゴンに追加
                    tri_list_regist->reserve( tri_list_regist->size()+tri_list_div_local.size() );           
                    for(int i=0; i<tri_list_div_local.size();i++ ) {
                        tri_list_regist->push_back( tri_list_div_local[i] );
                    } 
                }

            }   // !eof

            // ファイルClose
            ret = PolygonIO::load_file_close( ifs ); 

        } // !ファイル数
#ifdef DEBUG
        PL_DBGOSH << "   file read end. tri_list_regist->size()="<<tri_list_regist->size() <<endl;
#endif

        // 終了判定フラグ==1 送信
        iend = 1;
        MPI_Bcast( &iend, 1, MPI_INT, 0, comm );

        // ポリゴンの登録(登録したものは呼び出し側でdelete出来ない）
        set_triangles_ptr( tri_list_regist );
    }
    // ランク0以外
    else
    {
        std::vector<Triangle*>  *tri_list_regist = new std::vector<Triangle*>;
        std::vector<Triangle*>  tri_list_div_all_empty;  // 空リスト

        // 終了判定フラグ受信までループ
        while ( true )  {
            std::vector<Triangle*>  tri_list_div_local;

            // 終了判定フラグ受信
            MPI_Bcast( &iend, 1, MPI_INT, 0, comm );
            if( iend == 1 )  {
                break;
            }

            // ポリゴンのスキャタ
            ret = scatter_polygons( tri_list_div_all_empty, tri_list_div_local );

            // 登録するポリゴンに追加
            tri_list_regist->reserve( tri_list_regist->size()+tri_list_div_local.size() );           
            for(int i=0; i<tri_list_div_local.size();i++ ) {
                tri_list_regist->push_back( tri_list_div_local[i] );
            } 
        } 

        // ポリゴンの登録(登録したものは呼び出し側でdelete出来ない）
        set_triangles_ptr( tri_list_regist );
    }


    return PLSTAT_OK;
}


#endif
// eof
