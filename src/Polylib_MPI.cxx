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
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>
#include "Polylib.h"


using namespace std;
using namespace PolylibNS;

////////////////////////////////////////////////////////////////////////////
///
/// クラス:CalcAreaInfo
/// 計算領域情報。
///
////////////////////////////////////////////////////////////////////////////

// allgather_ParallelAreaInfo()処理用　通信用構造体  

struct CalcAreaInfo_com {
    PL_REAL     bpos[3];        // 基点座標
    PL_REAL     bbsize[3];      // 計算領域のボクセル数
    PL_REAL     gcsize[3];      // ガイドセルのボクセル数
    PL_REAL     dx[3];          // ボクセル１辺の長さ
    PL_REAL     gcell_min[3];   // ガイドセルを含めた担当領域の最小位置
    PL_REAL     gcell_max[3];   // ガイドセルを含めた担当領域の最大位置
};

// migrate()処理 内部用
struct ComHeadGroup {
    int     grp_id;     // ポリゴングループID
    int     num_tri;    // ポリゴングループ内のポリゴン数
    int     pl_type;    // ポリゴンタイプ
    int     numAtrI;    // ユーザ定義属性（整数）
    int     numAtrR;    // ユーザ定義属性（実数）
    int     pl_size;    // １ポリゴンの通信バッファサイズ(byte)
};

// migrate()処理 内部用
struct ComHeadRank {
    int             num_grp;            // 1ランクの送受信するポリゴングループ数
    int             nsize_mpi_double;   // MPI_DOUBLE単位としたときのポリゴン個数
    ComHeadGroup*   pGrpHead;           // 1グループ毎のヘッダ情報 ( num_grp数分あり )
};



////////////////////////////////////////////////////////////////////////////
/// 
/// クラス:Polylib
/// Polylibクラス MPI環境に関係するもの主に記載する
/// 
////////////////////////////////////////////////////////////////////////////


// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT
Polylib::init_parallel_info(
        MPI_Comm comm,
        PL_REAL bpos[3], 
        unsigned int bbsize[3], 
        unsigned int gcsize[3], 
        PL_REAL dx[3]
    )
{
    ParallelBbox  bbox;
    vector<ParallelBbox> bboxes;
    
    bbox.bpos[0]   = bpos[0];   bbox.bpos[1]   = bpos[1];    bbox.bpos[2]   = bpos[2];
    bbox.bbsize[0] = bbsize[0]; bbox.bbsize[1] = bbsize[1];  bbox.bbsize[2] = bbsize[2];
    bbox.gcsize[0] = gcsize[0]; bbox.gcsize[1] = gcsize[1];  bbox.gcsize[2] = gcsize[2];
    bbox.dx[0]     = dx[0];     bbox.dx[1]     = dx[1];      bbox.dx[2]     = dx[2];

    bboxes.push_back( bbox );
    
    return init_parallel_info( comm, bboxes );
}


POLYLIB_STAT
Polylib::init_parallel_info(
        MPI_Comm comm,
        const std::vector<ParallelBbox>&  bboxes
    )
{
#ifdef DEBUG
    PL_DBGOSH << "Polylib::init_parallel_info() in. " << endl;
#endif
    POLYLIB_STAT ret;

    // MPI情報の設定
    m_comm = comm;
    MPI_Comm_rank(comm, &m_myrank);
    MPI_Comm_size(comm, &m_numproc);

    // デバッグ出力用ランク番号文字列を設定
    std::ostringstream ostr;
    ostr << m_myrank;
    gs_rankno = "(rk:";
    gs_rankno += ostr.str();
    gs_rankno += ")";

#ifdef DEBUG
    PL_DBGOSH << "m_myrank: " << m_myrank << " m_numproc: " << m_numproc << endl;
#endif

    m_myproc_area.m_rank = m_myrank;
    for( int i=0; i<bboxes.size(); i++ )  {
        PL_REAL bbsize_f[3], gcsize_f[3];
        for (int j= 0; j<3; j++) {
            bbsize_f[j] = (PL_REAL)bboxes[i].bbsize[j];
            gcsize_f[j] = (PL_REAL)bboxes[i].gcsize[j];
        }

        Vec3<PL_REAL> v_bbsize(bbsize_f[0],bbsize_f[1],bbsize_f[2]);
        Vec3<PL_REAL> v_gcsize(gcsize_f[0],gcsize_f[1],gcsize_f[2]);
        Vec3<PL_REAL> v_bpos(bboxes[i].bpos[0],bboxes[i].bpos[1],bboxes[i].bpos[2]);
        Vec3<PL_REAL> v_dx(bboxes[i].dx[0],bboxes[i].dx[1],bboxes[i].dx[2]);
        // 自PE領域情報を設定
        CalcAreaInfo area;
        area.m_gcell_bbox.init();
        area.m_bpos   = v_bpos;
        area.m_bbsize = v_bbsize;
        area.m_gcsize = v_gcsize;
        area.m_dx     = v_dx;
        area.m_gcell_min = v_bpos-( v_gcsize )*v_dx;
        area.m_gcell_max = v_bpos+( v_bbsize+v_gcsize )*v_dx;
        area.m_gcell_bbox.add(area.m_gcell_min);
        area.m_gcell_bbox.add(area.m_gcell_max);

#ifdef DEBUG
        PL_DBGOSH << "(my_rank:" << m_myrank << "):" <<"i =" << i  << endl;
        PL_DBGOSH << "(my_rank:" << m_myrank << "):" <<"bpos      :" << v_bpos  << endl;
        PL_DBGOSH << "(my_rank:" << m_myrank << "):" <<"bbsize    :" << v_bbsize << endl;
        PL_DBGOSH << "(my_rank:" << m_myrank << "):" <<"gcsize    :" << v_gcsize << endl;
        PL_DBGOSH << "(my_rank:" << m_myrank << "):" <<"dx        :" << v_dx << endl;
        PL_DBGOSH << "(my_rank:" << m_myrank << "):" <<"gcell_min :"
            << area.m_gcell_min << endl;
        PL_DBGOSH << "(my_rank:" << m_myrank << "):" <<"gcell_max :"
            << area.m_gcell_max << endl;
#endif

        m_myproc_area.m_areas.push_back( area );
    }

    //  Polylib4.xよりランク内の担当領域数が複数個可能となったため、
    //  MPI_Allgatherは使えない、Rank0に一旦あつめて、その後、broadする
    //  実際はユーザ側は、全体をどう分割すのかは知っているはずであるが．．．
    //
    
    vector<ParallelAreaInfo> all_procs_area;
    // 全ランクの担当領域を設定
    ret = allgather_ParallelAreaInfo ( m_myproc_area, all_procs_area );
    if( ret != PLSTAT_OK ) {
        return ret;
    }
#ifdef DEBUG
    {
        PL_DBGOSH << "m_my_procs_area"<< endl;
        for( int j=0; j<m_myproc_area.m_areas.size(); j++ ) {   // 複数担当領域
            BBox* p_bbox = &(m_myproc_area.m_areas[j].m_gcell_bbox);
            Vec3<PL_REAL>* bbox_min = &(p_bbox->min);
            Vec3<PL_REAL>* bbox_max = &(p_bbox->max);
            PL_DBGOSH << "  gcell_bbox min="<<bbox_min->x<<" "<<bbox_min->y<<" "<<bbox_min->z<< endl;
            PL_DBGOSH << "  gcell_bbox max="<<bbox_max->x<<" "<<bbox_max->y<<" "<<bbox_max->z<< endl;
        }
    }
    for( int i=0; i<m_numproc; i++ ) {
        int irank = all_procs_area[i].m_rank;
        PL_DBGOSH << "all_procs_area:irank="<<irank<< endl;
        for( int j=0; j<all_procs_area[i].m_areas.size(); j++ ) {   // 複数担当領域
            BBox* p_bbox = &(all_procs_area[i].m_areas[j].m_gcell_bbox);
            Vec3<PL_REAL>* bbox_min = &(p_bbox->min);
            Vec3<PL_REAL>* bbox_max = &(p_bbox->max);
            PL_DBGOSH << "  gcell_bbox min="<<bbox_min->x<<" "<<bbox_min->y<<" "<<bbox_min->z<< endl;
            PL_DBGOSH << "  gcell_bbox max="<<bbox_max->x<<" "<<bbox_max->y<<" "<<bbox_max->z<< endl;
        }
    }
#endif

    // 自PEを除く全PE担当領域情報リスト設定
    // 自PE領域と隣接するPE領域情報リスト設定

    for( int irank=0; irank<m_numproc; irank++ ) {
        // 自PE領域情報はスキップ
        if( irank == m_myrank ) continue;

        // 自PEを除く全PE担当領域情報リストに登録
        m_other_procs_area.push_back( all_procs_area[irank] );
        
        // 自PE領域と隣接するPE領域情報はm_neibour_procs_areaにも追加
        bool bcross = false;
            // 自身の担当領域数
        for( int i=0; i<m_myproc_area.m_areas.size(); i++ ) {
            for( int j=0; j<all_procs_area[irank].m_areas.size(); j++ ) {
                if( m_myproc_area.m_areas[i].m_gcell_bbox.crossed( all_procs_area[irank].m_areas[j].m_gcell_bbox ) ) 
                {
                    bcross = true;
                    break;
                }
            }
            if( bcross ) break;
        }
        if( bcross ) {
            // 隣接PE担当領域情報リスト 追加
            m_neibour_procs_area.push_back( all_procs_area[irank] );
            // migrate除外三角形IDマップ の空を追加
            std::map< int, std::vector<long long int> >  exclusion_map_rank;
            m_exclusion_map_procs.push_back( exclusion_map_rank );
        }
    }

    return PLSTAT_OK;
}


// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT
Polylib::load(
    const std::string  config_filename,
    PL_REAL scale
)
{
#ifdef DEBUG
    PL_DBGOSH << m_myrank << ": " << "Polylib::load() in. " << endl;
#endif
    POLYLIB_STAT ret;

    // ポリゴングループtreeの作成：全ランクで同じデータを設定する
    try {
        tp->read(config_filename);

        ret = make_group_tree(tp);
        if( ret != PLSTAT_OK ) return ret;
    }
    catch( POLYLIB_STAT e ){
        return e;
    }
    
    //------------------------------------------------------------------------
    //  メモリ削減のため、ポリゴングループ毎に以下の処理を行う
    //      ・ポリゴンデータのロード
    //      ・ポリゴンデータのscatter
    //      ・ポリゴンIDが重複削除
    //------------------------------------------------------------------------

    // メモリ制限あり
    if( m_max_memory_size_mb > 0 )
    {
#ifdef DEBUG
        PL_DBGOSH << m_myrank << ": " << "Polylib::load() m_max_memory_size_mb="<<m_max_memory_size_mb <<endl;
#endif
        for (int i=0; i<m_pg_list.size(); i++) {
            // リーフの場合
            if ( m_pg_list[i]->get_children().empty() == true) {
                int used_mem_mb = (int)( used_memory_size_mb() );
                int max_mem_mb  = m_max_memory_size_mb - used_mem_mb;
#ifdef DEBUG
                PL_DBGOSH << m_myrank << ": " << "Polylib::load() i="<<i<<" max_mem_mb="<<max_mem_mb <<endl;
#endif
                if ( max_mem_mb <= 0 )  return PLSTAT_LACK_OF_MEMORY;

                ret = m_pg_list[i]->load_polygons_mem_reduced( scale, max_mem_mb );
                if (ret != PLSTAT_OK)       return ret;
            }
        }
    }
    // メモリ制限なし
    else
    {
        for (int i=0; i<m_pg_list.size(); i++) {
            // リーフの場合
            if ( m_pg_list[i]->get_children().empty() == true) {

                // STL/NPTファイルを読み込む
                // KDツリーも作成
                if( m_myrank == 0 )  {
                    ret = m_pg_list[i]->load_polygons_file(scale);
                    if (ret != PLSTAT_OK)       return ret;

                    // ※  MPI版でポリゴン毎の属性を設定するためには
                    //    このあたりに、ユーザ定義属性の処理が必要であるが
                    //    仕様上、各ポリゴン毎に固有の属性を設定する必要は
                    //    ないとのことなので対応しない
                    //    初期値としてはポリゴングループ内の全ポリゴンは同一値
                    //    を設定する仕様である
                }


                // 他ランクにポリゴンデータを分散する
                // 重複IDは削除
                ret = m_pg_list[i]->scatter_polygons();
                if (ret != PLSTAT_OK)       return ret;

            }
        }
    }

    return PLSTAT_OK;
}


// new version 
// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT
Polylib::save(
               std::string&       config_name_out,
               const std::string& file_format,
               std::string        extend
            )
{
#ifdef DEBUG
  PL_DBGOSH << "Polylib::save() in. " << endl;
#endif
    POLYLIB_STAT    ret;
    char    my_extend[128];
    POLYLIB_STAT stat;

    // 拡張文字列がカラであれば、現在時刻から作成
    if (extend == "") {
      time_t timer = time(NULL);
      struct tm *date = localtime(&timer);
      sprintf(my_extend, "%04d%02d%02d%02d%02d%02d",
          date->tm_year+1900,date->tm_mon+1,date->tm_mday,
          date->tm_hour,date->tm_min,date->tm_sec);
    }
    else {
      sprintf(my_extend, "%s", extend.c_str());
    }

    map<string,string> polygon_fname_map;

    for (int i=0; i<m_pg_list.size(); i++) {
        // リーフの場合
        if ( m_pg_list[i]->get_children().empty() == true) {
        
            std::vector<Triangle*>  tri_list;

            // ランク0にポリゴンデータを集約する
            // PolygonGroup内のポリゴンデータは更新しない
            // 集約したポリゴンデータは設定しない
            // 重複IDは削除
            ret = m_pg_list[i]->gather_polygons( tri_list );
            if (ret != PLSTAT_OK)       return ret;

            // STL/NPTファイルの保存
            if( m_myrank == 0 ) {
                string  rank_no = "";

                // 集約したポリゴンデータを出力
                m_pg_list[i]->save_polygons_file(rank_no,my_extend,file_format,
                                    &tri_list, polygon_fname_map);


                // ※MPI版の場合このあたりに、ユーザ定義属性の処理が必要

                // tri_list削除
                for(int j=0; j<tri_list.size(); j++ ) {
                    delete tri_list[j];
                }
            }
        }
    }
    
    // 初期化(tp)ファイルの更新保存
    stat=clearfilepath(tp);
    stat=setfilepath(polygon_fname_map);

    char    *config_name = save_config_file("", my_extend, file_format);
    //  PL_DBGOSH << __FUNCTION__ << " config_name "<< config_name << endl;

    if (config_name == NULL)    return PLSTAT_NG;
    //else    *p_config_filename = string(config_name);
    else    config_name_out = string(config_name);
    
    return PLSTAT_OK;
}


// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT
Polylib::save_parallel(
    //std::string *p_config_filename,
    std::string& config_name_out,
    const std::string& file_format,
    std::string extend
)
{
#ifdef DEBUG
    PL_DBGOSH << "Polylib::save_parallel() in. " << endl;
#endif
    POLYLIB_STAT ret;

    // 各ランク毎に保存
    //if( (ret = Polylib::save_at_rank( p_config_filename, m_myrank, m_numproc-1, extend, file_format)) != PLSTAT_OK ) {
    if( (ret = Polylib::save_at_rank( config_name_out, m_myrank, m_numproc-1, extend, file_format)) != PLSTAT_OK ) {
        PL_ERROSH << "[ERROR]Polylib::save_parallel():Polylib::save_at_rank():failed. returns:" << PolylibStat2::String(ret) << endl;
        return ret;
    }
    return PLSTAT_OK;
}


//TextParser version 
// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT Polylib::save_at_rank(
    //string      *p_config_name,
    string&     config_name_out,
    int         myrank,
    int         maxrank,
    const string&       extend,
    const string&       file_format
)
{
#ifdef DEBUG
    PL_DBGOSH << "Polylib::save_at_rank() in. " << endl;
#endif
    char    my_extend[128];

    // 拡張文字列がカラであれば、現在時刻から作成
    if (extend == "") {
        time_t      timer = time(NULL);
        struct tm   *date = localtime(&timer);
        sprintf(my_extend, "%04d%02d%02d%02d%02d%02d",
            date->tm_year+1900, date->tm_mon+1, date->tm_mday,
            date->tm_hour,      date->tm_min,   date->tm_sec);
    }
    else {
        sprintf(my_extend, "%s", extend.c_str());
    }

    // ランク番号の整形
    char    rank_no[16];
    int     fig = (int)log10((double)maxrank) + 1;
    sprintf(rank_no, "%0*d", fig, myrank);

#ifdef DEBUG
    PL_DBGOSH << "Polylib::save_with_rankno() rank_no" << rank_no<< endl;
#endif

    map<string,string> polygons_fname_map;
    // STL/NPTファイルの保存
    vector<PolygonGroup*>::iterator it;
    for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
        POLYLIB_STAT    stat;
        //リーフのみがポリゴン情報を持っている
        if ((*it)->get_children().empty() == false) continue;

        // ポリゴン数が0ならばファイル出力不要 2010.10.19
        if ((*it)->get_triangles()->size() == 0) continue;

        stat = (*it)->save_polygons_file(rank_no, my_extend, file_format,polygons_fname_map);
        if (stat != PLSTAT_OK)  return stat;
        //stat = (*it)->save_id_file(rank_no, my_extend, id_format);
        //if (stat != PLSTAT_OK)    return stat;
        //string rank_string,my_extend_string;
        //rank_string=rank_no;
        //my_extend_string = my_extend;
        //stat = (*it)->mk_param_tag(tp, "", "", "");
        //if (stat != PLSTAT_OK) return stat;
      
    }
    tp->changeNode("/"); //
    //  cout << "before cleanfilepath" <<endl;
    clearfilepath(tp); //clear whole file path
    //  cout << "after cleanfilepath" <<endl;

    //  cout << "before setfilepath" <<endl;
    setfilepath(polygons_fname_map);
    //cout << "after setfilepath" <<endl;

    // 定義ファイルの保存   
    char    *config_name = save_config_file(rank_no, my_extend, file_format);
#ifdef DEBUG    
    PL_DBGOSH << __FUNCTION__ << " config_name "<< config_name << endl;
#endif 
#ifdef DEBUG
    PL_DBGOSH << "Polylib::save_with_rankno() config_name " << config_name << endl;
#endif
    if (config_name == NULL)    return PLSTAT_NG;
    //else    *p_config_name = string(config_name);
    else    config_name_out = string(config_name);
    return PLSTAT_OK;

}


// public ////////////////////////////////////////////////////////////////////
POLYLIB_STAT
Polylib::move(
    PolylibMoveParams &params
)
{
#ifdef DEBUG
    PL_DBGOSH << "Polylib::move() in. " << endl;
#endif
    POLYLIB_STAT ret;
    vector<PolygonGroup*>::iterator group_itr;
    PolygonGroup *p_pg;

    // move実行前から隣接PE領域に懸かっている三角形をmigrate対象除外リストに載せる
    //     全隣接PE領域、全ポリゴングループ対象
    ret = select_excluded_trias();
    if( ret != PLSTAT_OK ) {
        PL_ERROSH << "[ERROR]Polylib::move():select_exclude_trias() failed. returns:" << PolylibStat2::String(ret) << endl;
        return ret;
    }

    // 各ポリゴングループのmove()を実行
    // 全グループに対して
    for (group_itr = m_pg_list.begin(); group_itr != m_pg_list.end(); group_itr++) {
            p_pg = (*group_itr);

        // 移動する可能性のあるポリゴングループのみ対象
        if( p_pg->get_movable() ) {

            // move実行
            if( (ret = p_pg->move( params )) != PLSTAT_OK ) {
                PL_ERROSH << "[ERROR]Polylib::move():(*group_itr)->move() failed. returns:" << PolylibStat2::String(ret) << endl;
                return ret;
            }

            // KD木を再構築 (三角形同士の位置関係が変化したため、再構築が必要)
            if( (ret = p_pg->rebuild_polygons()) != PLSTAT_OK ) {
                PL_ERROSH << "[ERROR]Polylib::move():(*group_itr)->rebuild_polygons() failed. returns:" << PolylibStat2::String(ret) << endl;
                return ret;
            }
        }
    }
    return PLSTAT_OK;
}


// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT
Polylib::select_excluded_trias( void )
{
#ifdef DEBUG
    PL_DBGOSH << "Polylib::select_excluded_trias() in." <<endl;
#endif
    POLYLIB_STAT ret;

    // migrate除外三角形IDマップ クリア
    m_exclusion_map_procs.clear();

    // 全隣接PEについて
    //     m_neibour_procsとm_exclusion_mapは同一のランク数が設定される
    for(unsigned int i=0; i<m_neibour_procs_area.size(); i++ ) {
        std::map< int, std::vector<long long int> >  exclusion_map_rank;
            //    pg_id            polygon_id
        exclusion_map_rank.clear();

        vector<BBox> bboxes;
        bboxes.clear();
        // 隣接PE内の担当領域数
        for(unsigned int j=0; j<m_neibour_procs_area[i].m_areas.size(); j++) {
            bboxes.push_back( m_neibour_procs_area[i].m_areas[j].m_gcell_bbox );
        }
        
        // 全隣接PEについて
        //     m_neibour_procsとm_exclusion_mapは同一のランク数が設定される
        for(int ig=0; ig<m_pg_list.size(); ig++ ) {
            // 移動する可能性のあるポリゴングループのみ対象
            if( m_pg_list[ig]->get_movable() ) {
            } else {
                continue; // Next Loop
            }

            vector<Triangle*>  tri_list;
            vector<long long int> ids;

            ids.clear();
            tri_list.clear();

            // 隣接PE領域(ガイドセル含)に懸かる三角形IDリストを作成
            bool duplicate = false; // 重複削除なし
            ret = m_pg_list[ig]->search( tri_list, bboxes, false, duplicate );
        
            for(unsigned int k=0; k<tri_list.size(); k++ ) {
                ids.push_back( tri_list[k]->get_id() );
            }
        
            // ポリゴンIDが重複しているが除外リストなので無視する

#ifdef DEBUG
    PL_DBGOSH << "group_name="<<m_pg_list[ig]->get_name()<<" grp internal_id:" << m_pg_list[ig]->get_internal_id() << " neibour_rank:" 
              << m_neibour_procs_area[i].m_rank
              << " 除外三角形数:" << ids.size() << endl;
#endif

            // migrate除外三角形IDマップに追加
            exclusion_map_rank[m_pg_list[ig]->get_internal_id()] = ids;

        }  // !pg loop
        
        m_exclusion_map_procs.push_back( exclusion_map_rank );

    } // ! neibour_procs loop

    return PLSTAT_OK;
}


// public /////////////////////////////////////////////////////////////////////
ParallelAreaInfo*
Polylib::get_proc_area(int rank)
{
    vector<ParallelAreaInfo>::iterator itr;
    itr = m_other_procs_area.begin();
    for (; itr != m_other_procs_area.end(); itr++) {
        if ((*itr).m_rank == rank) {
            return( &(*itr) );
        }
    }
    return NULL;
}


// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT
Polylib::allgather_ParallelAreaInfo(
        ParallelAreaInfo& myproc_area,
        std::vector<ParallelAreaInfo>& all_procs_area
    )
{
    int iret;
    int* num_areas_allrank = (int*)malloc( m_numproc*sizeof(int));  // ランク内の担当領域数テーブル

    all_procs_area.clear();
    
    // 全ランクの担当領域数取集
    int num_areas_myrank = myproc_area.m_areas.size();
    
    int ret = MPI_Allgather( &num_areas_myrank, 1, MPI_INT, num_areas_allrank, 1, MPI_INT, m_comm);
    if( ret  != MPI_SUCCESS) {
        PL_ERROSH << "[ERROR]Polylib::allgather_ParallelAreaInfo():MPI_Allgather "
                  << "faild." << endl;
        free( num_areas_allrank );
        return PLSTAT_MPI_ERROR;
    }

    // 受信する担当領域総数およびMAX数取得
    int num_tot_areas  = 0;

    for( int i=0; i<m_numproc; i++ ) {
        num_tot_areas += num_areas_allrank[i];
    }
    
    CalcAreaInfo_com* all_areas_com = (CalcAreaInfo_com*)malloc( num_tot_areas*sizeof(CalcAreaInfo_com) );
    CalcAreaInfo_com* my_area_com   = (CalcAreaInfo_com*)malloc( num_areas_myrank*sizeof(CalcAreaInfo_com) );
    
    // 自ランクを設定
    for( int i=0; i<num_areas_myrank; i++ ) {
        VEC3_TO_REAL( myproc_area.m_areas[i].m_bpos,   my_area_com[i].bpos   );
        VEC3_TO_REAL( myproc_area.m_areas[i].m_bbsize, my_area_com[i].bbsize );
        VEC3_TO_REAL( myproc_area.m_areas[i].m_gcsize, my_area_com[i].gcsize );
        VEC3_TO_REAL( myproc_area.m_areas[i].m_dx,     my_area_com[i].dx     );
        VEC3_TO_REAL( myproc_area.m_areas[i].m_gcell_min, my_area_com[i].gcell_min );
        VEC3_TO_REAL( myproc_area.m_areas[i].m_gcell_max, my_area_com[i].gcell_max );
    }

    // ランク０に担当領域情報を集約する
        // ランク０の処理
    if( m_myrank == 0 )
    {
        int ip = 0;
        // Rank0の担当領域設定
        memcpy( &all_areas_com[ip], my_area_com, num_areas_myrank*sizeof(CalcAreaInfo_com) );
        ip += num_areas_myrank;

        // Rank0以外の担当領域を受信して設定
        for( int i=1; i<m_numproc; i++ ) {
            MPI_Status status;

            // 各ランクからの担当領域受信
            iret = MPI_Recv( &all_areas_com[ip], (num_areas_allrank[i]*sizeof(CalcAreaInfo_com)), 
                             MPI_CHAR,
                             i, 1001, MPI_COMM_WORLD, &status );

            ip += num_areas_allrank[i];
        }
    
    }
        // ランク０以外の処理
    else
    {
        iret = MPI_Send( my_area_com, (num_areas_myrank*sizeof(CalcAreaInfo_com)), 
                     MPI_CHAR,
                     0, 1001, MPI_COMM_WORLD );
    }

    // ランク０より他ランクに送信する
    iret = MPI_Bcast( all_areas_com, (num_tot_areas*sizeof(CalcAreaInfo_com)), 
                         MPI_CHAR,
                         0, MPI_COMM_WORLD );
    
    
    // ParallelAreaInfoに詰め直す
    int ip=0;
    for( int i=0; i<m_numproc; i++ ) {
        ParallelAreaInfo para_area;
        para_area.m_rank = i;

        for( int j=0; j<num_areas_allrank[i]; j++ ) {
            CalcAreaInfo_com* carea_com = &all_areas_com[ip];
            CalcAreaInfo carea;
            REAL_TO_VEC3( carea_com->bpos,   carea.m_bpos );
            REAL_TO_VEC3( carea_com->bbsize, carea.m_bbsize );
            REAL_TO_VEC3( carea_com->gcsize, carea.m_gcsize );
            REAL_TO_VEC3( carea_com->dx,     carea.m_dx  );
            REAL_TO_VEC3( carea_com->gcell_min, carea.m_gcell_min );
            REAL_TO_VEC3( carea_com->gcell_max, carea.m_gcell_max );
            //REAL_TO_VEC3( carea_com->gcell_bbox_min, carea.m_gcell_bbox.min );
            //REAL_TO_VEC3( carea_com->gcell_bbox_max, carea.m_gcell_bbox.max );
            carea.m_gcell_bbox.add( carea_com->gcell_min );
            carea.m_gcell_bbox.add( carea_com->gcell_max );
            
            para_area.m_areas.push_back( carea );   // ランク内の担当領域追加
            ip += 1;
        }
        all_procs_area.push_back( para_area );  // 各ランクの情報として設定
    }
    
    free( num_areas_allrank );
    free( all_areas_com );
    free( my_area_com );
    
    return PLSTAT_OK;
}


// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT
Polylib::migrate( void )
{
#ifdef DEBUG
    PL_DBGOSH << "Polylib::migrate() in. " << endl;
#endif
    POLYLIB_STAT ret;
    int         iret;
    int num_req_grp_num = 0;
    int num_req_grp     = 0;
    int num_req_tri     = 0;
    MPI_Request* req_grp_num  = new MPI_Request[2*m_neibour_procs_area.size()];  // グループ数
    MPI_Request* req_grp      = new MPI_Request[2*m_neibour_procs_area.size()];  // グループ毎のヘッダ
    MPI_Request* req_tri      = new MPI_Request[2*m_neibour_procs_area.size()];  // ポリゴン情報のリクエスト
    MPI_Status*  stat_grp_num = new MPI_Status [2*m_neibour_procs_area.size()];
    MPI_Status*  stat_grp     = new MPI_Status [2*m_neibour_procs_area.size()];
    MPI_Status*  stat_tri     = new MPI_Status [2*m_neibour_procs_area.size()];

    ComHeadRank* head_rank      = new  ComHeadRank[m_neibour_procs_area.size()]; // 送信用
    ComHeadRank* head_recv_rank = new  ComHeadRank[m_neibour_procs_area.size()]; // 受信用

    char** send_buff_ranks = new char*[m_neibour_procs_area.size()];  // 送信バッファポインタ（ランク毎）
    char** recv_buff_ranks = new char*[m_neibour_procs_area.size()];  // 受信バッファポインタ（ランク毎）

    //-------------------------------------------------------
    // 送信データ準備・設定
    //-------------------------------------------------------
#ifdef DEBUG
    PL_DBGOSH << "Polylib::migrate() prepare send data start" << endl;
#endif

    //  領域確保用個数のカウント、初期設定
    //      隣接PE数ループ
    for (int i=0; i<m_neibour_procs_area.size(); i++) {
        int num_grp_rank_wk =0;
        // ポリゴングループ数
        for(int j=0; j<m_pg_list.size(); j++ ) {
            if( m_pg_list[j]->get_movable() ) { // グループが移動可？
                num_grp_rank_wk += 1;
            }
        }
        
        head_rank[i].num_grp          = 0;
        head_rank[i].nsize_mpi_double = 0;
        head_rank[i].pGrpHead = (ComHeadGroup*) malloc( num_grp_rank_wk*sizeof(ComHeadGroup) );
        head_recv_rank[i].num_grp          = 0;
        head_recv_rank[i].nsize_mpi_double = 0;
        head_recv_rank[i].pGrpHead = NULL;
        send_buff_ranks[i] = NULL;
        recv_buff_ranks[i] = NULL;
    }

    // 隣接PE数ループ
    for (int i=0; i<m_neibour_procs_area.size(); i++) {
        vector<BBox> bboxes;
        bboxes.clear();

        // 隣接PE内の担当領域数
        for(int j=0; j<m_neibour_procs_area[i].m_areas.size(); j++) {
            bboxes.push_back( m_neibour_procs_area[i].m_areas[j].m_gcell_bbox );
        }
        
        vector<Triangle*> tri_list_rank;
        tri_list_rank.clear();
        
        size_t  isize_send_buff = 0;
        
        // ポリゴングループ数
        for(int j=0; j<m_pg_list.size(); j++ ) {
        
            if( m_pg_list[j]->get_movable() ) { // グループが移動可？
            
                vector<Triangle*> tri_list_pg;
                tri_list_pg.clear();
            
                // 当該隣接PE領域への移動除外三角形IDリストを取得
                map< int, vector<long long int> >::iterator const itr =
                    m_exclusion_map_procs[i].find( m_pg_list[j]->get_internal_id() );
                
                // 当該隣接PE領域内にある移動フラグONの三角形を取得
                if( itr != m_exclusion_map_procs[i].end() ) {
                    ret = m_pg_list[j]->search_outbounded( tri_list_pg,
                                              bboxes, ((*itr).second) );
                } else {
#ifdef DEBUG
                    PL_DBGOSH << "Polylib::migrate() m_exclusion_map_procs not matched: i="<<i<<" group name="<<m_pg_list[j]->get_name() <<endl;
#endif
                }

                // 各ランクのヘッダに設定
                int ip= head_rank[i].num_grp;
                head_rank[i].pGrpHead[ip].grp_id  = m_pg_list[j]->get_internal_id();
                head_rank[i].pGrpHead[ip].num_tri = tri_list_pg.size();
                if( tri_list_pg.size() > 0 ) {
                    head_rank[i].pGrpHead[ip].pl_type = tri_list_pg[0]->get_pl_type();
                    head_rank[i].pGrpHead[ip].numAtrI = tri_list_pg[0]->get_num_atrI();
                    head_rank[i].pGrpHead[ip].numAtrR = tri_list_pg[0]->get_num_atrR();
                    head_rank[i].pGrpHead[ip].pl_size = tri_list_pg[0]->serialized_size();
                } else {
                    head_rank[i].pGrpHead[ip].pl_type = PL_TYPE_UNKNOWN;
                    head_rank[i].pGrpHead[ip].numAtrI = 0;
                    head_rank[i].pGrpHead[ip].numAtrR = 0;
                    head_rank[i].pGrpHead[ip].pl_size = 0;
                }

                head_rank[i].num_grp++; // グループ数カウントアップ
                
                // 移動ポリゴンをテーブルに設定
                if( tri_list_pg.size() > 0 ) {
                    tri_list_rank.reserve( tri_list_rank.size()+tri_list_pg.size() );
                    for(int k=0; k<tri_list_pg.size(); k++) {
                        tri_list_rank.push_back( tri_list_pg[k] );
                    }
                    isize_send_buff += head_rank[i].pGrpHead[ip].num_tri*head_rank[i].pGrpHead[ip].pl_size;
                }
            }
        }
        
        // 通信サイズ(MPI_DOUBLE単位）
        if( isize_send_buff > 0 ) {
            if( (isize_send_buff%sizeof(double)) == 0 ) {
                head_rank[i].nsize_mpi_double = isize_send_buff/sizeof(double); 
            } else {
                head_rank[i].nsize_mpi_double = isize_send_buff/sizeof(double) + 1;
            }
        
            send_buff_ranks[i] = (char*)malloc( (size_t)(head_rank[i].nsize_mpi_double)*sizeof(double) );
        } else {
            head_rank[i].nsize_mpi_double = 0;
            send_buff_ranks[i]            = NULL;
        }


        // 送信バッファに設定
        char* pbuff  = send_buff_ranks[i];
        unsigned int ip = 0;
        
        for(int j=0; j<head_rank[i].num_grp; j++ ) {
            int num_tri = head_rank[i].pGrpHead[j].num_tri;
            int pl_size = head_rank[i].pGrpHead[j].pl_size;

            for(int k=0; k<num_tri; k++) {
                // 送信バッファにパック
                pbuff = tri_list_rank[ip]->serialize( pbuff );
                ip++;
            }
        }
        
    }


    //-------------------------------------------------------
    // ランク毎のグループ数・MPI_DOUBLE数の送信／受信
    //-------------------------------------------------------
#ifdef DEBUG
    PL_DBGOSH << "Polylib::migrate() send/recv num of group start" << endl;
#endif
    num_req_grp_num = 0;
    int** idata      = alloc_array_2d_int( m_neibour_procs_area.size(), 2 );
    int** idata_recv = alloc_array_2d_int( m_neibour_procs_area.size(), 2 );

    // 隣接PE数ループ
    for (int i=0; i<m_neibour_procs_area.size(); i++) {
        // 送信
        idata[i][0] = head_rank[i].num_grp;     // 移動可能グループ数
        idata[i][1] = head_rank[i].nsize_mpi_double;    // 通信サイズ(MPI_DOUBLE単位）
        iret = MPI_Isend( idata[i], 2, MPI_INT,
                          m_neibour_procs_area[i].m_rank, 1101,
                          m_comm, &req_grp_num[num_req_grp_num] );
        if( iret != MPI_SUCCESS) {
            PL_ERROSH << "[ERROR]Polylib::migrate():MPI_Isend,"
                      << "num_grp send faild." << endl;
            return PLSTAT_MPI_ERROR;
        }
        num_req_grp_num++;

        // 受信
        iret = MPI_Irecv( idata_recv[i], 2, MPI_INT,
                          m_neibour_procs_area[i].m_rank, 1101,
                          m_comm, &req_grp_num[num_req_grp_num] );
        if( iret != MPI_SUCCESS) {
            PL_ERROSH << "[ERROR]Polylib::migrate():MPI_Irecv,"
                      << "num_grp receive faild." << endl;
            return PLSTAT_MPI_ERROR;
        }
        num_req_grp_num++;
    }

#ifdef DEBUG
PL_DBGOSH << "Polylib::migrate() MPI_Waitall(grp_num) start" << endl;
#endif
    iret = MPI_Waitall( num_req_grp_num, req_grp_num, stat_grp_num );
#ifdef DEBUG
PL_DBGOSH << "Polylib::migrate() MPI_Waitall(grp_num) end" << endl;
#endif


    //-------------------------------------------------------
    // ランク毎のグループ情報ComHeadGroupの送信／受信
    //-------------------------------------------------------
#ifdef DEBUG
    PL_DBGOSH << "Polylib::migrate() send/recv ComHeadGroup start" << endl;
#endif
    
    // 領域確保
    for (int i=0; i<m_neibour_procs_area.size(); i++) {
        head_recv_rank[i].num_grp          = idata_recv[i][0];
        head_recv_rank[i].nsize_mpi_double = idata_recv[i][1];
        if( head_recv_rank[i].num_grp > 0 ) {
            head_recv_rank[i].pGrpHead = (ComHeadGroup*)malloc( (size_t)(head_recv_rank[i].num_grp)*sizeof(ComHeadGroup) );
        } else {
            head_recv_rank[i].pGrpHead = NULL;
        }
        if( head_recv_rank[i].nsize_mpi_double > 0 ) {
            recv_buff_ranks[i] = (char*)malloc( (size_t)(head_recv_rank[i].nsize_mpi_double)*sizeof(double) );
        } else {
            recv_buff_ranks[i] = NULL;
        }
    }

    // ComHeadGroupの送信
    // 隣接PE数ループ
    num_req_grp = 0;
    for (int i=0; i<m_neibour_procs_area.size(); i++) {
        // 送信
        if( head_rank[i].num_grp > 0 ) {
            int nsize_send_char = head_rank[i].num_grp*sizeof(ComHeadGroup);
        
            iret = MPI_Isend( head_rank[i].pGrpHead, nsize_send_char, MPI_CHAR,
                              m_neibour_procs_area[i].m_rank, 1102,
                              m_comm, &req_grp[num_req_grp] );
            if( iret != MPI_SUCCESS) {
                PL_ERROSH << "[ERROR]Polylib::migrate():MPI_Isend,"
                          << "grp send faild." << endl;
                return PLSTAT_MPI_ERROR;
            }
            num_req_grp++;
        }

        // 受信
        if( head_recv_rank[i].num_grp > 0 ) {
            int nsize_recv_char = head_recv_rank[i].num_grp*sizeof(ComHeadGroup);

            iret = MPI_Irecv( head_recv_rank[i].pGrpHead, nsize_recv_char, MPI_CHAR,
                              m_neibour_procs_area[i].m_rank, 1102,
                              m_comm, &req_grp[num_req_grp] );
            if( iret != MPI_SUCCESS) {
                PL_ERROSH << "[ERROR]Polylib::migrate():MPI_Irecv,"
                          << "grp receive faild." << endl;
                return PLSTAT_MPI_ERROR;
            }
            num_req_grp++;
        }
    }

#ifdef DEBUG
PL_DBGOSH << "Polylib::migrate() MPI_Waitall(grp) start" << endl;
#endif
    iret = MPI_Waitall( num_req_grp, req_grp, stat_grp );
#ifdef DEBUG
PL_DBGOSH << "Polylib::migrate() MPI_Waitall(grp) end" << endl;
#endif

    //-------------------------------------------------------
    // ランク毎のポリゴン情報の送信／受信
    //-------------------------------------------------------
#ifdef DEBUG
    PL_DBGOSH << "Polylib::migrate() send/recv polygon data start" << endl;
#endif
    // 隣接PE数ループ
    num_req_tri = 0;
    for (int i=0; i<m_neibour_procs_area.size(); i++) {
        // 送信
        if( head_rank[i].nsize_mpi_double > 0 ) {
            iret = MPI_Isend( send_buff_ranks[i], head_rank[i].nsize_mpi_double, MPI_DOUBLE,
                          m_neibour_procs_area[i].m_rank, 1103,
                          m_comm, &req_tri[num_req_tri] );
            if( iret != MPI_SUCCESS) {
                PL_ERROSH << "[ERROR]Polylib::migrate():MPI_Isend,"
                          << "tri send faild." << endl;
                return PLSTAT_MPI_ERROR;
            }
            num_req_tri++;
        }

        // 受信
        if( head_recv_rank[i].nsize_mpi_double > 0 ) {
            iret = MPI_Irecv( recv_buff_ranks[i], head_recv_rank[i].nsize_mpi_double, MPI_DOUBLE,
                          m_neibour_procs_area[i].m_rank, 1103,
                          m_comm, &req_tri[num_req_tri] ); 
            if( iret != MPI_SUCCESS) {
                PL_ERROSH << "[ERROR]Polylib::migrate():MPI_Irecv,"
                          << "tri receive faild." << endl;
                return PLSTAT_MPI_ERROR;
            }
            num_req_tri++;
        }
    }

#ifdef DEBUG
PL_DBGOSH << "Polylib::migrate() MPI_Waitall(tri) start" << endl;
#endif
    iret = MPI_Waitall( num_req_tri, req_tri, stat_tri );
#ifdef DEBUG
PL_DBGOSH << "Polylib::migrate() MPI_Waitall(tri) end" << endl;
#endif

    //-------------------------------------------------------
    // ポリゴン情報作成（アンパック）およびポリゴン追加
    //-------------------------------------------------------
#ifdef DEBUG
    PL_DBGOSH << "Polylib::migrate() unpack polygon data start" << endl;
#endif

    // 隣接PE数
    for (int i=0; i<m_neibour_procs_area.size(); i++) {
        char* pbuff  = recv_buff_ranks[i];
        //unsigned int ip = 0;
        std::vector<Triangle*> tri_list;

        // グループ数
        int jp = 0;
        for(int j=0; j<m_pg_list.size(); j++ ) {
            if( m_pg_list[j]->get_movable() ) { // グループが移動可？
                int num_tri = head_recv_rank[i].pGrpHead[jp].num_tri;
                int pl_size = head_recv_rank[i].pGrpHead[jp].pl_size;
                int pl_type = head_recv_rank[i].pGrpHead[jp].pl_type;
                tri_list.clear();
                for(int k=0; k<num_tri; k++) {
                    // deserialize: Triangle/NptTrinangleオブジェクト生成
                    Triangle* pTri = deserialize_polygon( pl_type, pbuff );
                    pbuff += pl_size;
                    tri_list.push_back( pTri );
                }
                jp++;
                
                // ポリゴングループに三角形リストを追加
                if( (ret = m_pg_list[j]->add_triangles( tri_list )) != PLSTAT_OK ) {
                    PL_ERROSH << "[ERROR]Polylib::migrate():pg->add_triangles() failed. returns:"
                              << PolylibStat2::String(ret) << endl;
                    return ret;
                }

                // tri_list内のポリゴンは、add_triangles()内でコピーされているので必要なし
                for(int k=0; k<tri_list.size(); k++) {
                    delete tri_list[k];
                }
            }
        }
    }


    //-------------------------------------------------------
    // KD木作成、重複ポリゴン削除
    //-------------------------------------------------------
#ifdef DEBUG
    PL_DBGOSH << "Polylib::migrate() rebuild polygon start" << endl;
#endif

    // 移動してきた三角形を含めたKD木を再構築
    for (int i=0; i<m_pg_list.size(); i++) {

        // 移動可能グループだけ
        if( m_pg_list[i]->get_movable() && m_pg_list[i]->get_triangles() != NULL ) {

            // KD木を再構築
            if( (ret=m_pg_list[i]->rebuild_polygons()) != PLSTAT_OK ) {
                PL_ERROSH << "[ERROR]Polylib::migrate():p_pg->rebuild_polygons() failed. returns:"
                          << PolylibStat2::String(ret) << endl;
                return ret;
            }
            
            // 自PE領域外ポリゴン情報を消去
            if( m_pg_list[i]->erase_outbounded_polygons() != PLSTAT_OK ) {
                PL_ERROSH << "[ERROR]Polylib::migrate():erase_outbounded_polygons() failed." << endl;
            }
        }
    }

    //-------------------------------------------------------
    // 領域解放処理
    //-------------------------------------------------------
    for (int i=0; i<m_neibour_procs_area.size(); i++) {
        if( head_rank[i].pGrpHead != NULL ) free( head_rank[i].pGrpHead );
        if( head_recv_rank[i].pGrpHead != NULL ) free( head_recv_rank[i].pGrpHead );
        if( send_buff_ranks[i] != NULL ) free( send_buff_ranks[i] );
        if( recv_buff_ranks[i] != NULL ) free( recv_buff_ranks[i] );
    }

    delete[] req_grp_num;  delete[] req_grp;  delete[] req_tri;
    delete[] stat_grp_num; delete[] stat_grp; delete[] stat_tri;
    delete[] head_rank; delete[] head_recv_rank;
    delete[] send_buff_ranks; delete[] recv_buff_ranks;

    free_array_2d( (void**)idata ); free_array_2d(  (void**)idata_recv );

#ifdef DEBUG
    PL_DBGOSH << "Polylib::migrate() out normaly." << endl;
#endif
    return PLSTAT_OK;
}

#endif

// eof
