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

#ifdef MPI_PL
#include "mpi.h"
#endif
#include <fstream>
#include <map>
#include "Polylib.h"

using namespace std;
using namespace PolylibNS;

/// デバッグ用ランク番号グローバル文字列
std::string PolylibNS::gs_rankno = "";

/// 検索モード
///    true:
///      曲面補正が可能な場合、曲面補正して検索する  
///      KD木のbboxの大きさが変わる
///    false:
///      3角形平面で検索を行う
///      KD木のbboxの大きさは３頂点のみで決定される
///    デフォルトはfalse
///    設定はPolylib::set_srch_modeで行う
static bool polylib_srch_detail = false;

/************************************************************************
 *
 * Polylibクラス
 *
 ***********************************************************************/
// static & public /////////////////////////////////////////////////////////////////////
Polylib* Polylib::get_instance() {
    // この方法ならば、デストラクタを呼び出さなくてもクラスインスタンスの領域
    // が解放できるし、もし本関数が複数回呼び出されても、クラスインスタンスが
    // 複数作成されることはない(=singletonになる)
    static Polylib m_instance;
    return &m_instance;
}

// static & public /////////////////////////////////////////////////////////////////////
void Polylib::set_srch_mode( bool detail ) {
    polylib_srch_detail = detail;
}

// static & public /////////////////////////////////////////////////////////////////////
bool Polylib::get_srch_mode() {
    return polylib_srch_detail;
}


/// TextParser
// public /////////////////////////////////////////////////////////////////////
// 
#ifndef MPI_PL
// ロード 逐次(非MPI版）
POLYLIB_STAT Polylib::load(
         const string    config_name,
         PL_REAL         scale
    )
{
#ifdef DEBUG
    PL_DBGOSH << "Polylib::load_test() in." << endl;
#endif

    // 設定ファイル読み込み
    try {

      tp->read(config_name);
      //  tp->write("tmp.tpp");
      // グループツリー作成
      POLYLIB_STAT stat = make_group_tree(tp);
      if (stat != PLSTAT_OK)    return stat;

      // STLファイル読み込み (三角形IDファイルは不要なので、第二引数はダミー)
      //return load_polygons(false, ID_BIN, scale);
      return load_polygons( scale );
    }
    catch (POLYLIB_STAT e) {
        return e;
    }
}

// public /////////////////////////////////////////////////////////////////////
//TextParser 版

POLYLIB_STAT Polylib::save(
        string&  config_name_out,
        const string&   file_format,
        string          extend
    ) 
{
#ifdef DEBUG
    PL_DBGOSH << "Polylib::save() in." << endl;
#endif
    char    my_extend[128];
    POLYLIB_STAT stat=PLSTAT_OK;


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
    
    //PL_DBGOSH << __FUNCTION__ << " extend "<< my_extend << endl;

    map<string,string> polygon_fname_map;

    vector<PolygonGroup*>::iterator it;
    for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
      //リーフのみがポリゴン情報を持っている
      if ((*it)->get_children().empty() == false)   continue;

      // ポリゴン数が0ならばファイル出力不要 2010.10.19
      if ((*it)->get_triangles()->size() == 0)  continue;

      // STLファイル保存 (第一引数のランク番号は不要)
      stat = (*it)->save_polygons_file("", my_extend, file_format,
                           polygon_fname_map);
      if (stat != PLSTAT_OK) return stat;

      // Polylib4.xでは使用しない
      //stat = (*it)->mk_param_tag(tp, "", "", "");
      //if (stat != PLSTAT_OK) return stat;
    }


    // update stl filepath 
    // clear file path first
    stat=clearfilepath(tp);
    //tp->write("tmp2.tpp");
    // set filepath
    stat=setfilepath(polygon_fname_map);

    //  string tmp_extend = my_extend;
    //  char    *config_name = save_config_file("", tmp_extend, stl_format);
    char    *config_name = save_config_file("", my_extend, file_format);
    //  PL_DBGOSH << __FUNCTION__ << " config_name "<< config_name << endl;

    if (config_name == NULL)    return PLSTAT_NG;
    //else    *p_config_name = string(config_name);
    else    config_name_out = string(config_name);
    return PLSTAT_OK;
    
}

#endif


#ifndef MPI_PL

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT Polylib::move(
    PolylibMoveParams   &params
) {
#ifdef DEBUG
    PL_DBGOSH << "Polylib::move() in." << endl;
#endif
    POLYLIB_STAT ret;
    vector<PolygonGroup*>::iterator it;

    for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {

        // リーフグループで、movableフラグONのポリゴンを移動
        if ((*it)->get_children().empty() == true && (*it)->get_movable() ) {
            ret = (*it)->move(params);
            if (ret != PLSTAT_OK)   return ret;

            // 座標移動したのでKD木の再構築
            ret = (*it)->rebuild_polygons();
            if( ret != PLSTAT_OK ) return ret;
        }

    }
    return PLSTAT_OK;
}

#endif

// public /////////////////////////////////////////////////////////////////////
vector<PolygonGroup *> *Polylib::get_root_groups() const {
#ifdef DEBUG
    PL_DBGOSH << "Polylib::get_root_groups() in." << endl;
#endif
    vector<PolygonGroup*>                   *root = new vector<PolygonGroup*>;
    vector<PolygonGroup*>::const_iterator   it;

    for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
        if (((*it)->get_parent()) == NULL) {
            root->push_back(*it);
        }
    }
    return root;
}

// public /////////////////////////////////////////////////////////////////////
vector<PolygonGroup *> *Polylib::get_leaf_groups() const {
#ifdef DEBUG
    PL_DBGOSH << "Polylib::get_leaf_groups() in." << endl;
#endif
    vector<PolygonGroup*>                   *leaf = new vector<PolygonGroup*>;
    vector<PolygonGroup*>::const_iterator   it;

    for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
        if (((*it)->get_children()).size() == 0) {
            leaf->push_back(*it);
        }
    }
    return leaf;
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT  Polylib::search_polygons(
        std::vector<Triangle*>& tri_list,
        const string&           group_name, 
        const Vec3<PL_REAL>&    min_pos, 
        const Vec3<PL_REAL>&    max_pos, 
        bool            every
    ) const 
{
#ifdef DEBUG
    //PL_DBGOSH << "Polylib::search_polygons() in." << endl;
#endif
    PolygonGroup* pg = get_group(group_name);

    if (pg == 0) {
        PL_ERROSH << "[ERROR]Polylib::search_polygons():Group not found: " 
                  << group_name << endl;
        return PLSTAT_GROUP_NOT_FOUND;
    }
    vector<PolygonGroup*>* pg_list2 = new vector<PolygonGroup*>;

#ifdef BENCHMARK
    double st1, st2, ut1, ut2, tt1, tt2;
    bool ret1, ret2;
    ret1 = getrusage_sec(&ut1, &st1, &tt1);
#endif

    //子孫を検索
    search_group(pg, pg_list2);

    //自身を追加
    pg_list2->push_back(pg);

    // 検索範囲
    BBox bbox;
    bbox.init();
    bbox.add(min_pos);
    bbox.add(max_pos);

    //全ポリゴングループを検索
    POLYLIB_STAT ret2;
    vector<PolygonGroup*>::iterator it;
    for (it = pg_list2->begin(); it != pg_list2->end(); it++) {

        //リーフ構造からのみ検索を行う
        if ((*it)->get_children().size()==0) {
            ret2 = (*it)->search (tri_list, bbox, every);
            if (ret2 != PLSTAT_OK) {
                delete pg_list2;
                return ret2;
            }
        }
    }
#ifdef BENCHMARK
    ret2 = getrusage_sec(&ut2,&st2,&tt2);
    if (ret1 == false || ret2 == false) {
        PL_DBGOSH << "Search SYS   Time Error" << endl;
        PL_DBGOSH << "Search CPU   Time Error" << endl;
        PL_DBGOSH << "Search Total Time Error" << endl;
    }
    else{
        cout.setf(ios::scientific, ios::floatfield);
        PL_DBGOSH << "Search SYS   Time:" << st2 - st1 << endl;
        PL_DBGOSH << "Search CPU   Time:" << ut2 - ut1 << endl;
        PL_DBGOSH << "Search Total Time:" << tt2 - tt1 << endl;
        std::cout.unsetf(ios::scientific);
    }
#endif

    delete pg_list2;
#ifdef DEBUG
    //PL_DBGOSH << "Polylib::search_polygons() out. tri_list.size()="<< tri_list.size() <<endl;
#endif

    return PLSTAT_OK;
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT  Polylib::search_polygons(
        std::vector<NptTriangle*>& tri_list,
        const string&              group_name, 
        const Vec3<PL_REAL>&       min_pos, 
        const Vec3<PL_REAL>&       max_pos, 
        bool            every
    ) const 
{
#ifdef DEBUG
    PL_DBGOSH << "Polylib::search_polygons() in." << endl;
#endif
    std::vector<Triangle*>  triangle_list;
    
    POLYLIB_STAT ret = search_polygons( triangle_list,group_name,min_pos,max_pos,every );
    if (ret != PLSTAT_OK) {
        return ret;
    }

    // vector<Triangle*> -> vector<NptTriangle*> 変換
    ret = convert_polygons_to_npt( triangle_list, tri_list );

    return ret;
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT Polylib::check_group_name(
    const string    &name, 
    const string    &path
) {
#ifdef DEBUG
    PL_DBGOSH << "Polylib::check_group_name() in." << endl;
#endif
    if (name.empty() == true) {
        PL_ERROSH << "[ERROR]Polylib::check_group_name():Group name is empty." 
                  << endl;
        return PLSTAT_GROUP_NAME_EMPTY;
    }

    vector<PolygonGroup*>::iterator it;
    for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
        if ((*it)->get_name() == name && (*it)->get_parent_path() == path) {
            PL_ERROSH << "[ERROR]Polylib::check_group_name():Group name is "
                << "duplicate:name:" << name << "," << "path:" << path 
                << endl;
            return PLSTAT_GROUP_NAME_DUP;
        }
    }
    return PLSTAT_OK;
}


// public /////////////////////////////////////////////////////////////////////
void Polylib::add_pg_list(PolygonGroup *pg)
{
#ifdef DEBUG
    PL_DBGOSH << "Polylib::add_pg_list() in." << endl;
#endif
    m_pg_list.push_back(pg);
}

// public /////////////////////////////////////////////////////////////////////
/// 2010.10.20 引数FILE *追加。
void Polylib::show_group_hierarchy(FILE *fp)
{
#ifdef DEBUG
    PL_DBGOSH << "Polylib::show_group_hierarchy() in." << endl;
#endif
    string      tab;
    vector<PolygonGroup*>::iterator it;
    for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
        if ((*it)->get_parent() != NULL) {
            // Not Use
        }
        else{
            show_group_name(*it, tab, fp);
        }
    }
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT Polylib::show_group_info(const string& group_name, bool detail)
{
#ifdef DEBUG
    PL_DBGOSH << "Polylib::show_group_info() in." << endl;
#endif
    PolygonGroup *p = get_group(group_name);
    if (p == NULL) {
        PL_ERROSH << "[ERROR]Polylib::show_group_info():Group not found:"
             << group_name << endl;
        return PLSTAT_GROUP_NOT_FOUND;
    }

    return p->show_group_info(-1,detail);
}

// public /////////////////////////////////////////////////////////////////////
void Polylib::show_all_group_info( bool detail )
{
#ifdef DEBUG
    PL_DBGOSH << "Polylib::show_all_group_info() in. " << endl;
    PL_DBGOSH << "Polylib:: number of polygon group ="<<m_pg_list.size() << endl;
#endif
    for(int i=0; i<m_pg_list.size(); i++)  {
        m_pg_list[i]->show_group_info(-1,detail);
    }
}

// public /////////////////////////////////////////////////////////////////////
size_t Polylib::used_memory_size()
{
    size_t  size;
    vector<PolygonGroup*>::iterator     pg;
    
    //-----------------------------------------
    //  逐次部のメモリサイズ
    //-----------------------------------------

    // 自クラス
    size = sizeof(Polylib);

    // ポリゴングループ
#ifdef DEBUG
PL_DBGOSH << "Polylib::used_memory_size:PolygonGroup num=" << m_pg_list.size() << endl;
#endif
    for (pg = m_pg_list.begin(); pg != m_pg_list.end(); pg++) {

        // PolygonGroupクラス
        //      PolygonGroupの属性は対象外としている
        size += sizeof(PolygonGroup);

        // 三角形移動前一時リスト
        size += (*pg)->get_num_of_trias_before_move() * sizeof(vector<void> *);

        // リーフにはポリゴンがある
        if ((*pg)->get_children().empty()) {

            // 三角形ポリゴン
            vector<Triangle*>   *tri_list = (*pg)->get_triangles();
            if( tri_list!=NULL && tri_list->size()>0 ) {
#ifdef DEBUG
PL_DBGOSH << "Polylib::used_memory_size:Triangle num=" << tri_list->size() << endl;
#endif
                if( tri_list->size() > 0 ) {
                    size_t pl_size = (*tri_list)[0]->used_memory_size();
                    size += tri_list->size() * pl_size;
                }
            }

            // KD木
            VTree   *vtree = (*pg)->get_vtree();
            if( vtree != NULL ) {
                size += vtree->memory_size();
            }
        }

    }
    

#ifdef MPI_PL
    //-----------------------------------------
    //  並列処理用のメモリサイズ
    //      MPI内部で使う通信バッファ等は対象外
    //-----------------------------------------

    // 自PE担当領域情報
    size += sizeof(ParallelAreaInfo);
    size += m_myproc_area.m_areas.size()*sizeof(CalcAreaInfo);

    // 自PEを除く全PE担当領域情報リスト
    size += sizeof(vector<ParallelAreaInfo>);
    size += m_other_procs_area.size()*sizeof(ParallelAreaInfo);
    for(int i=0; i<m_other_procs_area.size(); i++) {
        size += m_other_procs_area[i].m_areas.size()*sizeof(CalcAreaInfo);
    }

    // 隣接PE担当領域情報リスト
    size += sizeof(vector<ParallelAreaInfo>);
    size += m_neibour_procs_area.size()*sizeof(ParallelAreaInfo);
    for(int i=0; i<m_neibour_procs_area.size(); i++) {
        size += m_neibour_procs_area[i].m_areas.size()*sizeof(CalcAreaInfo);
    }
       
    // migrate除外三角形IDマップ
    map< int, vector<long long int> >::iterator ex;
    size += sizeof(std::vector< std::map< int, std::vector<long long int> > >);
    size += m_exclusion_map_procs.size()*sizeof(std::map< int, std::vector<long long int> >);
    for(int i=0; i<m_exclusion_map_procs.size(); i++) {
       size += m_exclusion_map_procs[i].size() * (sizeof(int)+sizeof(vector<long long int>)); 
       for (ex = m_exclusion_map_procs[i].begin(); 
                                    ex != m_exclusion_map_procs[i].end(); ex++) {
            size += ex->second.size() * sizeof(long long int);
        }
    }

    // 自プロセスのランク数、全プロセス数
    size += sizeof(int) * 2;

    // 自プロセスが利用するコミュニケーター
    size += sizeof(MPI_Comm);

#endif

    return size;
}

// public /////////////////////////////////////////////////////////////////////
PolygonGroup* Polylib::get_group(const string& name) const
{
#ifdef DEBUG
    //PL_DBGOSH << "Polylib::get_group(" << name << ") in." << endl;
#endif
    vector<PolygonGroup*>::const_iterator it;
    for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
        if (name == (*it)->acq_fullpath()) {
#ifdef DEBUG
/*
            PL_DBGOSH  << "get_group: " << (*it)->get_parent_path()
                       << " name: " << (*it)->get_name()
                       << " size: " << (*it)->get_children().size() << endl;
*/
#endif
            return *it;
        }
    }
    for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
        if (name == (*it)->get_name()) {
            return *it;
        }
    }
#ifdef DEBUG
    //PL_DBGOSH << "Polylib::get_group(" << name << ") returns NULL" << endl;
#endif
    return NULL;
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT  Polylib::search_nearest_polygon(
        Triangle*&              tri,
        const string&           group_name, 
        const Vec3<PL_REAL>&    pos
    ) const
{

    PolygonGroup* pg = get_group(group_name);
    if (pg == 0) {
        PL_ERROSH << "[ERROR]Polylib::search_polygons():Group not found: " 
                  << group_name << endl;
        return PLSTAT_NG;
    }

    vector<PolygonGroup*>* pg_list2 = new vector<PolygonGroup*>;

    //子孫を検索
    search_group(pg, pg_list2);

    //自身を追加
    pg_list2->push_back(pg);

    Triangle* tri_min = 0;
    PL_REAL dist2_min = 0.0;

    //対象ポリゴングループ毎に検索
    vector<PolygonGroup*>::iterator it;
    for (it = pg_list2->begin(); it != pg_list2->end(); it++) {
        //リーフポリゴングループからのみ検索を行う
        if ((*it)->get_children().size()==0) {
            Triangle* tri_near;
            (*it)->search_nearest(tri_near,pos);
            if (tri_near) {
                Vec3<PL_REAL>* v = tri_near->get_vertexes();
                Vec3<PL_REAL>  c((v[0][0]+v[1][0]+v[2][0])/3.0,
                        (v[0][1]+v[1][1]+v[2][1])/3.0,
                        (v[0][2]+v[1][2]+v[2][2])/3.0);
                PL_REAL dist2 = (c - pos).lengthSquared();
                if (tri_min == 0 || dist2 < dist2_min) {
                    tri_min = tri;
                    dist2_min = dist2;
                }
            }
        }
    }

    delete pg_list2;

    tri = tri_min;

    return PLSTAT_OK;
}

// protected //////////////////////////////////////////////////////////////////
Polylib::Polylib()
{
    gs_rankno = "";

    //Polylib にTextParser クラスを持たせる。
    tp = new TextParser;


    m_max_memory_size_mb = 0;

    //PL_DBGOS<< __FUNCTION__ <<" m_factory "<< m_factory << " tp " << tp<<std::endl;
}

// protected //////////////////////////////////////////////////////////////////
Polylib::~Polylib()
{
    vector<PolygonGroup*>::iterator it;
    for (it = m_pg_list.begin(); it != m_pg_list.end();) {
#ifdef DEBUG
PL_DBGOSH << "~Polylib():" << (*it)->get_name() << endl;
#endif
        (*it)->get_children().clear();
        delete *it;
        it = m_pg_list.erase(it);
    }
    if(tp !=0) delete tp;

#ifdef MPI_PL
    //-----------------------------------------
    //  並列処理用の終了化
    //      ※特になし
    //-----------------------------------------
#endif
#ifdef DEBUG
PL_DBGOSH << "~Polylib() normal end." << endl;
#endif

}

// protected //////////////////////////////////////////////////////////////////
// TextParser 版
POLYLIB_STAT Polylib::make_group_tree(
    TextParser* tp
) {
#ifdef DEBUG
    PL_DBGOSH << "Polylib::make_group_tree(TextParser) in." << endl;
#endif
    TextParserError status = TP_NO_ERROR;

    // 念のため階層構造の最上位へ
    string cur; 
    tp->currentNode(cur);
    if(cur != "/Polylib") {
      status=tp->changeNode("/Polylib");
      tp->currentNode(cur);

      if(status!=TP_NO_ERROR){
        PL_ERROSH << 
          "[ERROR]Polylib::make_group_tree(TextParser):Root node not found."
              << endl;
        return PLSTAT_CONFIG_ERROR;
      }

    }
#ifdef DEBUG      
      PL_DBGOSH<<"Top current_node "<<cur <<endl;
#endif // DEBUG   
#if 0 
    if(cur != "/") {
      status=tp->changeNode("/");
      tp->currentNode(cur);

      if(status!=TP_NO_ERROR){
        PL_ERROSH << 
          "[ERROR]Polylib::make_group_tree(TextParser):Root node not found."
              << endl;
        return PLSTAT_CONFIG_ERROR;
      }

    }
#endif

    // ノードとリーフのリストを取る

    //ノードを取得
    vector<string> nodes;
    tp->getNodes(nodes);
    string current_node;


    //  tp->currentNode(current_node);
    // vector<string> leaves;
    // tp->getLabels(leaves);
    // if(nodes.size()==0 && leaves.size()==0){ return PLSTAT_CONFIG_ERROR;}

    // string class_name = "PolygonGroup";
    // if(leaves.size()!=0){
    //   vector<string>::iterator leaf_iter=find(leaves.begin(),
    //                    leaves.end(),
    //                    PolygonGroup::ATT_NAME_CLASS);
    //   if(leaf_iter != leaves.end()){
    //     //       class_name = *leaf_iter;
    //     string value;
    //     status=tp->getValue((*leaf_iter),value);
    //     class_name=value;

    //   }
    // }
#ifdef DEBUG      
      PL_DBGOSH<<"nodes.size()="<<nodes.size() <<endl;
      for(int i=0; i<nodes.size(); i++ ) {
          PL_DBGOSH<<"nodes["<<i<<"]="<<nodes[i] <<endl;
      }
#endif // DEBUG   


    // loop over node recurcively.
    for(vector<string>::iterator nodes_iter = nodes.begin(); 
        nodes_iter != nodes.end();
        nodes_iter++){

      status = tp->changeNode(*nodes_iter);
      tp->currentNode(current_node);
#ifdef DEBUG      
      PL_DBGOSH<<"current_node "<< current_node <<  endl;
#endif // DEBUG   
      vector<string> nodes;
      vector<string> leaves;

      tp->getNodes(nodes);
      tp->getLabels(leaves);
      if(nodes.size()==0 && leaves.size()==0){ return PLSTAT_CONFIG_ERROR;}
/*
      string class_name = "PolygonGroup"; //default

      if(leaves.size()!=0){
        vector<string>::iterator leaf_iter=find(leaves.begin(),
                            leaves.end(),
                            PolygonGroup::ATT_NAME_CLASS);
        if(leaf_iter != leaves.end()){
          //     class_name = *leaf_iter;
          string value;
          status=tp->getValue((*leaf_iter),value);
          class_name=value;
        }
      }
      //pg = m_factory->create_instance(class_name);
*/
      PolygonGroup *pg = new PolygonGroup();
      // ルートポリゴングループとして登録
      add_pg_list(pg);
/*
      if (pg == NULL) {
        PL_ERROSH << "[ERROR]Polylib::make_group_tree():Class name not found."
              << class_name
              << endl;
        return PLSTAT_CONFIG_ERROR;
      }
*/


      // 配下のタグを取得して、PolygonGroupツリーを作成
      //      POLYLIB_STAT res = pg->build_group_tree(this, NULL, elem);
      //  
      POLYLIB_STAT res = pg->build_group_tree(this, NULL, tp);
      if (res != PLSTAT_OK) return res;
      status = tp->changeNode("..");

    } 
    

    return PLSTAT_OK;
}


// protected //////////////////////////////////////////////////////////////////
// 逐次処理用
POLYLIB_STAT Polylib::load_polygons(
    PL_REAL     scale
)
{
#ifdef DEBUG
    PL_DBGOSH << "Polylib::load_polygons() in." << endl;
#endif
    vector<PolygonGroup*>::iterator it;
    for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
        // リーフの場合
        if ((*it)->get_children().empty() == true) {

            //STL/NPTファイルを読み込む
            POLYLIB_STAT ret = (*it)->load_polygons_file(scale);
            if (ret != PLSTAT_OK)       return ret;

        }
    }
    return PLSTAT_OK;
}


// protected //////////////////////////////////////////////////////////////////
//TextPArser 版
char *Polylib::save_config_file(
    const string&   rank_no,
    const string&   extend,
    const string&   format
){
#ifdef DEBUG
  //  PL_DBGOSH << "Polylib::save_config_file() in. " << endl;
#endif

    vector<PolygonGroup*>::iterator it;
    //POLYLIB_STAT                  stat;

    // ファイル出力
    char    *config_name;

    if (rank_no == "") {
        // ランク番号不要
      config_name = polylib_config_save_file("", extend);
      //config_name = PolylibConfig::save_file(doc, "", extend);
    }
    else {
      // 要ランク番号(MPI版)
      //    config_name = PolylibConfig::save_file(doc, rank_no, extend);
      config_name = polylib_config_save_file(rank_no, extend);
    }

#ifdef DEBUG
    //PL_DBGOSH << "save_config_file(): " << config_name << endl;
#endif
//  xmlFreeDoc(doc);

    return config_name;
}

/////////////////////////////////////////　　
// 追加　２０１２ー０８
//////////////////////////////////

char* Polylib::polylib_config_save_file(
                    const string&  rank_no,
                    const string&  extend
                )
{

#define POLYLIB_CONFIG_NAME     "polylib_config"  
#define POLYLIB_CONFIG_EXT      "tpp"

  //  cout <<__FUNCTION__<<" in "<< rank_no <<" "<< extend << endl;

  static char fname[1024];
  if (rank_no == ""){
    sprintf(fname,"%s_%s.%s",POLYLIB_CONFIG_NAME,extend.c_str(),POLYLIB_CONFIG_EXT);
  } else {
    sprintf(fname, "%s_%s_%s.%s", POLYLIB_CONFIG_NAME, rank_no.c_str(),
        extend.c_str(), POLYLIB_CONFIG_EXT);
  }
  string fname_string = fname;

  //config_file の書き込み

  tp->write(fname_string,1);

  return fname;

}



/////////////////////////////////////////　　
// 追加　２０１２ー０８
//////////////////////////////////

POLYLIB_STAT Polylib::setfilepath(map<string,string>& polygons_fname_map){

#ifdef DEBUG
  PL_DBGOS << "stl_map size " <<  polygons_fname_map.size()<<endl;
#endif //DEBUG
  //  tp->changeNode("/");
  tp->changeNode("/Polylib");
  for(map<string,string>::iterator map_iter=polygons_fname_map.begin();
      map_iter != polygons_fname_map.end();
      map_iter++){
#ifdef DEBUG
    PL_DBGOS << "stl_map " <<  map_iter->first <<" "<< map_iter->second<<endl;
#endif // DEBUG
    tp->changeNode(map_iter->first); //
    string cur;
    tp->currentNode(cur);
    
    //    cout << __FUNCTION__ << " " <<  cur << endl;
    string value = "\"" + map_iter->second + "\"";
    //    cout << "before leaf createion" <<endl;
    tp->createLeaf("filepath",value);
    //    cout << "after leaf createion" <<endl;
    //check
    //    vector<string> leaves;
    //   tp->getLabels(leaves);  //label　の取り直し
    //    for(vector<string>::iterator leaf_iter=leaves.begin();
    //  leaf_iter!=leaves.end();
    //  leaf_iter++){
    //      cout << __FUNCTION__ << " "  << *leaf_iter << endl;
    //    }
    //        string value;
    //    tp->getValue("filepath",value);
    //    cout << __FUNCTION__ << " " << cur << " "<< value <<endl;

    //    tp->changeNode("/");

    tp->changeNode("/Polylib");
  }
  return PLSTAT_OK;
}
/////////////////////////////////////////　　
// 追加　２０１２ー０８
//////////////////////////////////
POLYLIB_STAT Polylib::clearfilepath(TextParser* tp_ptr){

  // recursive にするため、TextParserのポインタを引数に取る。

  vector<string> leaves;
  tp_ptr->getLabels(leaves,1);

  vector<string>::iterator leaf_iter=find(leaves.begin(),leaves.end(),"filepath");
  if(leaf_iter!=leaves.end()){ // 見つかったら
    tp_ptr->deleteLeaf("filepath");
  } 
  leaves.clear();
  tp_ptr->getLabels(leaves);


  leaf_iter=find(leaves.begin(),leaves.end(),"filepath[0]");
  if(leaf_iter!=leaves.end()){ // 見つかったら

    int index=0;
    leaf_iter = leaves.begin();
    while(leaf_iter != leaves.end()){
      stringstream ss;
      string tmpstring="filepath";
      ss << tmpstring <<"["<<index<<"]";
      ss >> tmpstring;
      leaf_iter = find(leaf_iter,leaves.end(),tmpstring);
      if(leaf_iter == leaves.end()) break;
      TextParserError tp_error=tp_ptr -> deleteLeaf(tmpstring);
      if (tp_error!=TP_NO_ERROR) {
    PL_ERROSH << "[ERROR]Polylib::save() "
          << "can not remove leaf = " << tmpstring << endl;
    return PLSTAT_NG;
      }
      index++;
      leaf_iter++;
    }
  }

  vector<string> nodes;
  tp_ptr->getNodes(nodes);
  for( vector<string>::iterator node_iter=nodes.begin();
       node_iter!=nodes.end();
       node_iter++){
    tp_ptr->changeNode(*node_iter);
    clearfilepath(tp_ptr);
    tp_ptr->changeNode("..");
  }

  return PLSTAT_OK;

}



// protected //////////////////////////////////////////////////////////////////
// 2010.10.20 引数FILE *追加。
void Polylib::show_group_name(
    PolygonGroup    *pg, 
    string          tab,
    FILE            *fp
){
    vector<PolygonGroup*>::iterator it;

    // ユーザ定義id出力 2010.10.20
    // ユーザ定義ラベル出力 2012.08.31
    // ユーザ定義タイプ出力 2013.07.17
    if (fp == NULL) {
        //PL_DBGOSH << "Polylib::show_group_name: " << tab ;
        //if (pg->get_parent_path().empty() == true)  PL_DBGOS << "+";
        //PL_DBGOS << pg->get_name() << ":" << pg->acq_file_name() << endl;
        std::cout<<gs_rankno << "Polylib::show_group_name: " << tab ;
        if (pg->get_parent_path().empty() == true)  std::cout << "+";
        std::vector<Triangle*>* tri_list = pg->get_triangles();
        if( tri_list == NULL ) {
            std::cout << pg->get_name() << ":" << pg->acq_file_name() << endl;
        } else {
            std::cout << pg->get_name() << ":" << pg->acq_file_name();
            std::cout << " num_tri="<<tri_list->size() << endl;
        }

    }
    else {
        // 出力先が指定されていれば、そちらに出力
        fprintf(fp, "%sPolylib::show_group_name:%s%s%s:%s\n", 
            gs_rankno.c_str(),      tab.c_str(),
            (pg->get_parent_path().empty() == true) ? "+" : "",
            pg->get_name().c_str(), pg->acq_file_name().c_str() );
    }

    tab = tab + "    ";
    for (it = pg->get_children().begin(); it != pg->get_children().end(); it++){
        show_group_name(*it, tab, fp);
    }
}

// protected //////////////////////////////////////////////////////////////////
PolygonGroup* Polylib::get_group(int internal_id) const
{
#ifdef DEBUG
    PL_DBGOSH << "Polylib::get_group(" << internal_id << ") in." << endl;
#endif
    vector<PolygonGroup*>::const_iterator it;
    for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
        if (internal_id == (*it)->get_internal_id()) {
            return *it;
        }
    }
#ifdef DEBUG
    PL_DBGOS << "Polylib::get_group(" << internal_id << ") returns NULL" << endl;
#endif
    return NULL;
}



// private ////////////////////////////////////////////////////////////////////
void Polylib::search_group(
    PolygonGroup            *p, 
    vector<PolygonGroup*>   *pg
) const {
#ifdef DEBUG
    //PL_DBGOSH << "Polylib::search_group() in." << endl;
#endif
    vector<PolygonGroup*>::iterator it;
    for (it = p->get_children().begin(); it != p->get_children().end(); it++) {
        pg->push_back(*it);
        search_group(*it, pg);  // 再帰呼び出し
    }
}

