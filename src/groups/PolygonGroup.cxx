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

#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "Polylib.h"
#include "polygons/Triangle.h"
#include "polygons/NptTriangle.h"
#include "groups/PolygonGroup.h"
#include "file_io/PolygonIO.h"
#include "c_lang/CPolylib.h"

//#define BENCHMARK
//#define DEBUG

using namespace std;
using namespace PolylibNS;

#define M_MAX_ELEMENTS 15   /// VTreeのノードが持つ最大要素数

///
/// 本クラス内でのみ使用するTextParserのタグ
///     ユーザ定義属性のヘッダ（node名として現れる)
///     PolygonGroupとしては登録しない
///
#define ATT_USER_NODE_NAME      "UserAtr"

///
/// 他クラスでも使用するTextParserのタグ
///     Polylib4.xでは使用しないため読み飛ばす
///
//const char *PolygonGroup::ATT_NAME_CLASS = "class_name";

///
/// 本クラス内でのみ使用するTextParserのタグ
///
#define ATT_NAME_PATH       "filepath"
#define ATT_NAME_MOVABLE    "movable"

///
/// 本クラス内でのみ使用するTextParserのタグ
///     Polylib4.xでは使用しないため読み飛ばす
///
//#define ATT_NAME_NAME       "name"
// ユーザ定義ID追加 2010.10.20
//#define ATT_NAME_ID         "id"
// ユーザ定義ラベル追加 2012.08.31
//#define ATT_NAME_LABEL      "label"
// ユーザ定義タイプ追加 2013.07.17
//#define ATT_NAME_TYPE       "type"

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


/************************************************************************
 *
 * PolygonGroupクラス
 *
 ***********************************************************************/
// public /////////////////////////////////////////////////////////////////////
PolygonGroup::PolygonGroup() {
    m_movable      = false;
    //m_parent    = 0;
    m_parent       = NULL;
    m_tri_list     = NULL;
    m_vtree        = NULL;
    m_move_func    = NULL;
    m_move_func_c  = NULL;
    m_need_rebuild = false;
    m_trias_before_move = NULL;
    m_max_elements = M_MAX_ELEMENTS;
}

// public /////////////////////////////////////////////////////////////////////
PolygonGroup::~PolygonGroup()
{
#ifdef DEBUG
    PL_DBGOSH << "~PolygonGroup() in." << endl;
#endif
    if( m_trias_before_move != NULL ) {
        for( unsigned int i=0; i<m_trias_before_move->size(); i++ ) {
            delete m_trias_before_move->at(i);
        }
        delete m_trias_before_move;
    }

    if (m_vtree != NULL) {
        delete m_vtree;
    }

    if (m_tri_list != NULL) {
        for(int i=0; i<m_tri_list->size(); i++ ) {
            delete( (*m_tri_list)[i] );
        }
        delete m_tri_list;
    }

#ifdef DEBUG
    PL_DBGOSH << "~PolygonGroup() normal end." << endl;
#endif
}


// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolygonGroup::init(
            const vector<Triangle*> *tri_list, 
            bool                    clear
        ) 
{
#ifdef DEBUG
    PL_DBGOSH <<"PolygonGroup::init3:clear=" << clear << endl;
#endif
    if (clear == true) {
        init_tri_list();
        // ポリゴンの複製
        copy_polygons( *tri_list, *m_tri_list );
    }
    return build_polygon_tree();
}


// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolygonGroup::init(
            const vector<NptTriangle*>  *tri_list, 
            bool                        clear
    ) 
{
    if (clear == true) {
        init_tri_list();
        m_tri_list->reserve( tri_list->size() );
        for(int i=0; i<tri_list->size(); i++) {
            NptTriangle* pNpt =  (*tri_list)[i];
            m_tri_list->push_back( new NptTriangle(*pNpt) );
        }
    }
    return build_polygon_tree();
}

// public /////////////////////////////////////////////////////////////////////
void PolygonGroup::set_triangles_ptr(
                vector<Triangle*> *tri_list,
                bool build_tree
            )
{
    delete_tri_list();      // 既存データ削除
    m_tri_list = tri_list;     // ポインタ登録
    if( build_tree ) {
        build_polygon_tree();  // KDツリー等の構築
    }
}


//TextParser version
// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolygonGroup::build_group_tree(
    Polylib                 *polylib,
    PolygonGroup            *parent,
    TextParser* tp
) 
{
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::build_group_tree() in." << endl;
#endif

    //current で既に作成されているPolygonGroup:this に属性をつける
    //          ユーザ定義属性を除く
    POLYLIB_STAT    ret = setup_attribute(polylib, parent, tp);
    if (ret != PLSTAT_OK)       return ret;

    // 配下のノード（複数）取得
    //    Polylib4.x系ではユーザ定義属性用のノードが含まれる
    //

    vector<string> nodes;
    nodes.clear();

    TextParserError error=TP_NO_ERROR;
    error=tp->getNodes(nodes);
#ifdef DEBUG
    string node_name;
    tp->currentNode(node_name);
    PL_DBGOSH << "PolygonGroup::build_group_tree() m_name="<<m_name<<" node_name="<<node_name<<" child nodes.size()="<<nodes.size() <<endl;
    for(int i=0; i<nodes.size(); i++ ) {
        PL_DBGOSH << "  nodes["<<i<<"]="<<nodes[i] <<endl;
    }
#endif


    // tp で　情報を読みながら、ループ
    for(vector<string>::iterator nodes_iter=nodes.begin();
        nodes_iter != nodes.end();
        nodes_iter++){

        //-------------------------
        // ユーザ定義属性ノードの処理
        //-------------------------
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::build_group_tree() *nodes_iter="<<*nodes_iter <<endl;
#endif
        // ユーザ定義属性ノード判定
        if( *nodes_iter == ATT_USER_NODE_NAME ) {
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::build_group_tree() setup_user_attribute()" << endl;
#endif
            // ユーザ定義属性設定
            POLYLIB_STAT    ret = setup_user_attribute( *nodes_iter, tp);
            if (ret != PLSTAT_OK) {
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::build_group_tree() setup_user_attribute() ERROR ret="<<ret <<endl;
#endif
                return ret;
            }
            continue;   // 次のノードへ
        }

        //-------------------------
        // 下位ノードの処理
        //-------------------------
        //  下位ノードに移動
        error=tp->changeNode(*nodes_iter);
        if(error!=TP_NO_ERROR){
            PL_ERROSH << "[ERROR]PolygonGroup::build_group_tree():"
                  << " TextParser error " 
                  << tp->TextParserErrorHandler(error,"can not move to ") 
                  << (*nodes_iter) << endl;
            return PLSTAT_CONFIG_ERROR;
        }

        PolygonGroup* pg = new PolygonGroup();
        // ポリゴングループの登録
        polylib->add_pg_list(pg);

        // 再帰呼び出し
        ret = pg->build_group_tree(polylib, this, tp);

        // go up and next
        // 下位ノードから復帰
        tp->changeNode("..");

    }

    return PLSTAT_OK;
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolygonGroup::build_polygon_tree()
{
#ifdef BENCHMARK
    double  st1, st2, ut1, ut2, tt1, tt2;
    bool    ret1, ret2;
    ret1 = getrusage_sec(&ut1, &st1, &tt1);
#endif
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::build_polygon_tree() in.:" << m_name << endl;
#endif

    if( m_tri_list == NULL )  {
        if (m_vtree != NULL) {
            delete m_vtree;
            m_vtree = NULL;
        }
        return PLSTAT_OK;
    }

    //木構造の生成
    //POLYLIB_STAT ret = m_polygons->build();
    {
        // NptTriangleもあるため、各オブジェクトよりBBox取得する
        BBox bbox;
        bbox.init();
        bool detail = Polylib::get_srch_mode();  // Polylib環境より検索モード取得
        //PL_DBGOSH << "PolygonGroup::build_polygon_tree() get_srch_mode() detai="<<detail << endl;

        for( int i=0; i<m_tri_list->size(); i++ ) {
            BBox bbox_tri = (*m_tri_list)[i]->get_bbox( detail );  // 必ずtrueなのか確認要
            bbox.add( bbox_tri.min );
            bbox.add( bbox_tri.max );
        }
        m_bbox = bbox;

#ifdef DEBUG
        Vec3<PL_REAL> min = m_bbox.getPoint(0);
        Vec3<PL_REAL> max = m_bbox.getPoint(7);
        PL_DBGOSH << "Polygons::build:min=(" <<min<< "),max=(" <<max<< ")" << endl;
#endif

        // 木構造作成
        if (m_vtree != NULL) delete m_vtree;
        m_vtree = new VTree(m_max_elements, m_bbox, m_tri_list);
    }

#ifdef BENCHMARK
    ret2 = getrusage_sec(&ut2,&st2,&tt2);
    if (ret1 == false || ret2 == false) {
        PL_DBGOSH << "Reading STL SYS Time Error" << endl;
        PL_DBGOSH << "Reading STL CPU Time Error" << endl;
        PL_DBGOSH << "Reading STL Total Time Error" << endl;
    }
    else{
        cerr.setf(ios::scientific, ios::floatfield);
        PL_DBGOSH << "Reading STL SYS   Time:" << st2 - st1 << endl;
        PL_DBGOSH << "Reading STL CPU   Time:" << ut2 - ut1 << endl;
        PL_DBGOSH << "Reading STL Total Time:" << tt2 - tt1 << endl;
        cerr.unsetf(ios::scientific);
    }
#endif

#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::build_polygon_tree() out." << endl;
#endif
    return PLSTAT_OK;
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolygonGroup::load_polygons_file( PL_REAL scale )
{
#ifdef DEBUG
PL_DBGOSH << "PolygonGroup:load_polygons_file():IN" << endl;
#endif
    //POLYLIB_STAT ret = m_polygons->import(m_polygon_files, scale);
    init_tri_list();
    POLYLIB_STAT ret = PolygonIO::load(m_tri_list, m_polygon_files, scale);
    if (ret != PLSTAT_OK) return ret;

#if 0
// Polylib4.xでは使用しない
    // m_idが指定されていたら、その値で全三角形のm_exidを更新
    // added by tkawanab 2013.06.17
    if( m_id_defined ) {
        ret = m_polygons->set_all_exid( m_id );
    }
    if (ret != PLSTAT_OK) return ret;
#endif

    return build_polygon_tree();
}


// TextParser でのsaveの為、save した stl/Npt ファイルを記憶しておく
/////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolygonGroup::save_polygons_file(
    const string&   rank_no,
    const string&   extend,
    const string&   format,
    map<string,string>& polygon_fname_map
) 
{
    return save_polygons_file( rank_no,extend,format, m_tri_list, polygon_fname_map);
}

// ■新規
/////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolygonGroup::save_polygons_file(
    const string&   rank_no,
    const string&   extend,
    const string&   format,
    std::vector<Triangle*> *tri_list,
    map<string,string>& polygon_fname_map
) 
{
  char  *fname = mk_polygons_fname(rank_no, extend, format,polygon_fname_map);
  return PolygonIO::save(tri_list, fname, format);
}


// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolygonGroup::move(
    PolylibMoveParams   &params
)
{
#ifdef DEBUG
PL_DBGOSH <<  "PolygonGroup::move() in." << endl;
#endif
        // 移動関数は登録されていれば移動関数を呼ぶ出す
            // C++で登録された移動関数の場合
    if( m_move_func !=NULL ) 
    {
#ifdef DEBUG
        PL_DBGOSH <<  "PolygonGroup::move()  m_move_func() start" << endl;
#endif
        m_move_func( this, &params );
    } 
        // Cインターフェースで登録された移動関数の場合
        //    extern "C" int m_move_func_c( PL_TAG, const PolylibMoveParamsStruct* )の宣言が必要
        //    正式にどう書くかはあとで調べる
    else if ( m_move_func_c !=NULL ) 
    {
#ifdef DEBUG
        PL_DBGOSH <<  "PolygonGroup::move()  m_move_func() set params_struct" << endl;
#endif
        //
        //   params_struct <- param の変換要
        //
        ::PolylibMoveParamsStruct params_struct;   // 先頭に::をつけないとPolylibNS::PolylibMoveParamsStruct
                                                   // とされて未定義となる
        params_struct.m_current_step = params.m_current_step;
        params_struct.m_next_step    = params.m_next_step;
        params_struct.m_delta_t      = params.m_delta_t;
        int num = sizeof(params_struct.m_params)/sizeof(PL_REAL);
        for(int i=0; i<num; i++ ) {
            params_struct.m_params[i] = params.m_params[i];
        }
                        
#ifdef DEBUG
        PL_DBGOSH <<  "PolygonGroup::move()  m_move_func_c() start" << endl;
#endif
        PL_GRP_TAG this_tag = reinterpret_cast<PL_GRP_TAG>(this);
        m_move_func_c( this_tag, &params_struct );   // m_internal_id is tagID(handleID)
    }

    return PLSTAT_OK;
}


// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolygonGroup::search(
            vector<Triangle*>&  tri_list,
            const BBox&         bbox, 
            bool                every 
        ) const
{
    return m_vtree->search(tri_list, bbox, every );
}



// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolygonGroup::search(
            vector<Triangle*>&  tri_list,
            const vector<BBox>& bboxes, 
            bool                every, 
            bool                duplicate
        ) const
{
    POLYLIB_STAT ret;
    if( m_vtree == NULL ) {
        return PLSTAT_OK;
    }

    for(int i=0; i<bboxes.size(); i++ ) {
        ret = m_vtree->search(tri_list, bboxes[i], every );
        if( ret != PLSTAT_OK )  return ret;
    }
    
    //重複削除
    if( !duplicate && tri_list.size()>0 ) {
        // IDでソート
        std::sort( tri_list.begin(), tri_list.end(), TriangleLess() );  // TriangleLess():比較用
        // ID重複分を削除
        tri_list.erase(
                        std::unique( tri_list.begin(), tri_list.end(), TriangleEqual() ), // TriangleEqual():比較用
                        tri_list.end()
                    );
    }
    
    return PLSTAT_OK;
}


// public /////////////////////////////////////////////////////////////////////
string PolygonGroup::acq_fullpath() {
    if (m_parent_path.empty() == true)  return m_name;
    else                                return m_parent_path + "/" + m_name;
}

// public /////////////////////////////////////////////////////////////////////
string PolygonGroup::acq_file_name() {
    string                          fnames;
    map<string, string>::iterator   it;
    for (it = m_polygon_files.begin(); it != m_polygon_files.end(); it++) {
        if (it == m_polygon_files.begin()) {
            fnames = it->first;
        }
        else {
            fnames.append(",");
            fnames.append(it->first);
        }
    }
    return fnames;
}

#ifdef MPI_PL

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT
PolygonGroup::search_outbounded(
    std::vector<Triangle*>&     tri_list,
    std::vector<BBox>& neibour_bboxes,
    std::vector<long long int>&  exclude_tria_ids
)
{
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::search_outbounded() in. " << endl;
    PL_DBGOSH << "    group name="<<this->get_name() <<endl;
    PL_DBGOSH << "    exclude_tria_ids.size()="<<exclude_tria_ids.size() <<endl;
#endif
    POLYLIB_STAT ret;

    // 除外IDリストを昇順ソート
    std::sort( exclude_tria_ids.begin(), exclude_tria_ids.end() );

    // 隣接PE領域(ガイドセル含)に懸かる三角形を検索
    for(int i=0; i<neibour_bboxes.size(); i++ ) {
        ret = search( tri_list, neibour_bboxes[i], false );
    }
#ifdef DEBUG
PL_DBGOSH << "p_trias orig tri_list.size()=" << tri_list.size() << endl;
#endif

    // 検索結果から除外対象を除く
    vector<Triangle*>::iterator itr;
    for( itr=tri_list.begin(); itr!=tri_list.end(); ) {
        long long int id = (*itr)->get_id();
        if( std::binary_search(exclude_tria_ids.begin(),
                               exclude_tria_ids.end(), id) ) {
            itr = tri_list.erase(itr);
        }
        else {
            itr++;
        }
    }
#ifdef DEBUG
PL_DBGOSH << "p_trias return tri_list.size()="<<tri_list.size() <<endl;
#endif
    return PLSTAT_OK;
}

// public /////////////////////////////////////////////////////////////////////
// 三角形リストの追加
//   内部IDが重複した三角形は追加しない。KD木の再構築はしない。
POLYLIB_STAT
PolygonGroup::add_triangles(
    const std::vector<Triangle*>& tri_list
)
{
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::add_triangles() in. " << endl;
#endif
    if( tri_list.size()==0 ) {
        return PLSTAT_OK;
    }

    if (m_tri_list == NULL) {
        m_tri_list = new vector<Triangle*>;
    }

    //m_polygons->add( &tri_list );
    int num_reserve = m_tri_list->size() + tri_list.size();
    vector<Triangle*> tri_list_tmp;
    tri_list_tmp.reserve( num_reserve );
    tri_list_tmp = *m_tri_list;
    m_tri_list->clear();
    m_tri_list->reserve( num_reserve );

    // ひとまず全部追加
    copy_polygons( tri_list, tri_list_tmp ); // deep copyして追加

    // 三角形リストをID順にソート
    std::sort( tri_list_tmp.begin(), tri_list_tmp.end(), TriangleLess() );
    // 追加
    long long int id_pre = tri_list_tmp[0]->get_id();
    m_tri_list->push_back( tri_list_tmp[0] );

    for(int i=1; i<tri_list_tmp.size(); i++ ) {
        long long int id = tri_list_tmp[i]->get_id();
        if( id == id_pre ) {
            // IDが同一なので削除
            delete tri_list_tmp[i];
        } else {
            // 新しいIDなのでリストに登録
            m_tri_list->push_back( tri_list_tmp[i] );
            id_pre = id;
        }
    }

    // KD木要再構築フラグを立てる
    m_need_rebuild = true;

    return PLSTAT_OK;
}

#endif

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT
PolygonGroup::rebuild_polygons()
{
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::rebuild_polygons() in. " << endl;
#endif
    // 不要な再構築は行わない
    if( !m_need_rebuild ) {
#ifdef DEBUG
        PL_DBGOSH << "PolygonGroup::rebuild_polygons() didnot need rebuild." << endl;
#endif
        return PLSTAT_OK;
    }

    POLYLIB_STAT ret = build_polygon_tree();
    m_need_rebuild = false;
    return ret;
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolygonGroup::show_group_info(
        int irank,
        bool detail
    )
{
    ostringstream   oss;
    string          rank;

    if (irank < 0) {
        rank = "";
    }
    else{
        oss << setw(3) << setfill('0') << "rank(" << irank << "): ";
        rank = oss.str();
    }
    PL_DBGOSH << "PolygonGroup::show_group_info::rank:" << rank << endl;

    if (m_name.empty() == false) {
        PL_DBGOSH << "  polygon group name: " << m_name << endl;
    }
    else {
        PL_DBGOSH << "  polygon group name: empty." << endl;
    }

    PL_DBGOSH << "  movable: "<<m_movable <<endl;

    if (m_parent_path.empty() == false) {
        PL_DBGOSH << "  parent polygon group name: " << m_parent_path << endl;
    }
    else {
        PL_DBGOSH << "  parent polygon group name: empty." << endl;
    }

    if (m_polygon_files.size() > 0) {
        map<string, string>::iterator it = m_polygon_files.begin();
        for (; it != m_polygon_files.end(); it++) 
            PL_DBGOSH << "  file name: " << (*it).first << endl;
    }
    else {
        PL_DBGOSH << "  file name: empty." << endl;
    }

     PL_DBGOSH << "  group num of attribute: "<<m_atr.size() <<endl;
     for(int i=0; i<m_atr.size(); i++) {
         PL_DBGOSH << "  attribute   key: "<<m_atr[i].key<<"   value: "<<m_atr[i].value <<endl;
     }

    if (m_tri_list == NULL) {
        PL_DBGOSH << "  triangle is nothing." <<endl;
        //return PLSTAT_POLYGON_NOT_EXIST;
        return PLSTAT_OK;
    }

    PL_DBGOSH << "  triangle list size: " << m_tri_list->size() << endl;
    // サマリダンプ
#ifdef MPI_PL
    PL_DBGOSH << "  group local  area: "<<get_group_area() <<endl;
    PL_DBGOSH << "  group global area: "<<get_group_global_area() <<endl;
#else
    PL_DBGOSH << "  group area: "<<get_group_area() <<endl;
#endif

    // 詳細ダンプ
    if( detail )  {
        PL_DBGOSH << "  vertex vector list: " << endl;
        vector<Triangle*>::iterator it;
        for (it = m_tri_list->begin(); it != m_tri_list->end(); it++) {
            Vec3<PL_REAL> *vtx = (*it)->get_vertexes();
            for (int i=0; i<3; i++) {
                PL_DBGOSH << "    id:" << i        << " x:" << vtx[i][0] 
                     << " y:"     << vtx[i][1] << " z:" << vtx[i][2] << endl;
            }
        }

        PL_DBGOSH << "  normal vector list: " << endl;
        for (it = m_tri_list->begin(); it != m_tri_list->end(); it++) {
            Vec3<PL_REAL> vtx = (*it)->get_normal();
            PL_DBGOSH << "    x:" << vtx[0] << " y:" << vtx[1] << " z:" << vtx[2] <<endl;
        }

        PL_DBGOSH << "  triangle area list: " << endl;
        for (it = m_tri_list->begin(); it != m_tri_list->end(); it++) {
            PL_DBGOSH << "    area:" << (*it)->get_area() << endl;
        }
    }


    return PLSTAT_OK;
}

int PolygonGroup::get_group_num_tria( void ) {
  
    return (int)m_tri_list->size();
}

// グループ内のポリゴンの面積を積算して返す
//   ローカルの面積
PL_REAL PolygonGroup::get_group_area( void )
{
    PL_REAL area = 0.0;

    if( m_tri_list == NULL || m_tri_list->size() == 0 ) {
        return 0.0;
    }
    for(int i=0; i<m_tri_list->size(); i++ ) {
        area += (*m_tri_list)[i]->get_area();
    }

    return area;
}

#ifdef MPI_PL
#else

// グループ内のポリゴン属性（整数）の総和を返す
//   逐次版
POLYLIB_STAT PolygonGroup::get_polygons_reduce_atrI(
            PL_OP_TYPE op,
            int        atr_no,
            int&       val
        )
{
    val = 0;
    int val_tmp = 0;
    if( m_tri_list == NULL || m_tri_list->size() == 0 ) {
        return PLSTAT_POLYGON_NOT_EXIST;
    }
    int num_atrI = (*m_tri_list)[0]->get_num_atrI();
    if( atr_no >= num_atrI ) {
        return PLSTAT_ATR_NOT_EXIST;
    }

    if( op == PL_OP_MAX ) {
        int* pAtr = (*m_tri_list)[0]->get_pAtrI();
        val_tmp = pAtr[atr_no];
        for(int i=1; i<m_tri_list->size(); i++ ) {
            pAtr = (*m_tri_list)[i]->get_pAtrI();
            if( pAtr[atr_no] > val_tmp ) {
               val_tmp = pAtr[atr_no];
            }
        }
        val = val_tmp;
    } else if( op == PL_OP_MIN ) {
        int* pAtr = (*m_tri_list)[0]->get_pAtrI();
        val_tmp = pAtr[atr_no];
        for(int i=1; i<m_tri_list->size(); i++ ) {
            pAtr = (*m_tri_list)[i]->get_pAtrI();
            if( pAtr[atr_no] < val_tmp ) {
               val_tmp = pAtr[atr_no];
            }
        }
        val = val_tmp;
    } else if( op == PL_OP_SUM ) {
        for(int i=0; i<m_tri_list->size(); i++ ) {
            int* pAtr = (*m_tri_list)[i]->get_pAtrI();
            val_tmp += pAtr[atr_no];
        }
        val = val_tmp;
    } else {
        return PLSTAT_NG;
    }

    return PLSTAT_OK;
}

// グループ内のポリゴン属性（実数）の総和を返す
//   逐次版
POLYLIB_STAT PolygonGroup::get_polygons_reduce_atrR(
            PL_OP_TYPE op,
            int        atr_no,
            PL_REAL&   val
        )
{
    val = 0.0;
    PL_REAL val_tmp = 0.0;
    if( m_tri_list == NULL || m_tri_list->size() == 0 ) {
        return PLSTAT_POLYGON_NOT_EXIST;
    }
    int num_atrR = (*m_tri_list)[0]->get_num_atrR();
    if( atr_no >= num_atrR ) {
        return PLSTAT_ATR_NOT_EXIST;
    }

    if( op == PL_OP_MAX ) {
        PL_REAL* pAtr = (*m_tri_list)[0]->get_pAtrR();
        val_tmp = pAtr[atr_no];
        for(int i=1; i<m_tri_list->size(); i++ ) {
            pAtr = (*m_tri_list)[i]->get_pAtrR();
            if( pAtr[atr_no] > val_tmp ) {
               val_tmp = pAtr[atr_no];
            }
        }
        val = val_tmp;
    } else if( op == PL_OP_MIN ) {
        PL_REAL* pAtr = (*m_tri_list)[0]->get_pAtrR();
        val_tmp = pAtr[atr_no];
        for(int i=1; i<m_tri_list->size(); i++ ) {
            pAtr = (*m_tri_list)[i]->get_pAtrR();
            if( pAtr[atr_no] < val_tmp ) {
               val_tmp = pAtr[atr_no];
            }
        }
        val = val_tmp;
    } else if( op == PL_OP_SUM ) {
        for(int i=0; i<m_tri_list->size(); i++ ) {
            PL_REAL* pAtr = (*m_tri_list)[i]->get_pAtrR();
            val_tmp += pAtr[atr_no];
        }
        val = val_tmp;
    } else {
        return PLSTAT_NG;
    }

    return PLSTAT_OK;
}

#endif


POLYLIB_STAT PolygonGroup::rescale_polygons( PL_REAL scale )
{
    if( m_tri_list->size() == 0 ) {
        return PLSTAT_OK;
    }

    for( int i=0; i<m_tri_list->size(); i++ ) {
        (*m_tri_list)[0]->rescale( scale );
    }

    m_need_rebuild = true;
    return rebuild_polygons();
}

#if 0
// Polylib4.xでは使用しない
// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolygonGroup::set_all_exid_of_trias( int id )
{
  m_id = id;           // keno 2013-07-20
  m_id_defined = true; // keno 2013-07-20
    return m_polygons->set_all_exid( id );
}
#endif

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolygonGroup::search_nearest(
       Triangle*&              tri,
        const Vec3<PL_REAL>&   pos
    ) const 
{
    //return m_polygons->search_nearest(tri,pos);
    tri = const_cast<Triangle*>(m_vtree->search_nearest(pos)); // constを外す
    return PLSTAT_OK;
}

/// ポリゴングループのユーザ定義属性取得。
///
///  @param[in]     key     キー
///  @param[out]    val     属性値
///  @return OK/NG  NG:キーと属性のペアが登録されていない
///
POLYLIB_STAT PolygonGroup::get_atr( std::string& key, std::string& val ) const
{
    for(int i=0; i<m_atr.size(); i++ ) {
        if( m_atr[i].key == key ) { 
           val = m_atr[i].value;
           return PLSTAT_OK;
        }
    }
    return PLSTAT_NG;
}

///
/// ポリゴングループのユーザ定義属性設定
///
///  @param[in]     key     キー
///  @param[in]     val     属性値
///  @return なし
///  @attention 既に登録されていた場合、上書きする
///
void PolygonGroup::set_atr( std::string& key, std::string& val )
{
    for(int i=0; i<m_atr.size(); i++ ) {
        if( m_atr[i].key == key ) { 
           m_atr[i].value = val;
           return;
        }
    }

    UsrAtr atr;
    atr.key=key;  atr.value=val;
    m_atr.push_back( atr );
}


// TextParser Version
// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT PolygonGroup::setup_attribute (
        Polylib                 *polylib,
        PolygonGroup            *parent, 
        TextParser* tp
    )
{
#ifdef DEBUG
  PL_DBGOS << __FUNCTION__ << " in"  <<endl;
#endif
  TextParserError tp_error = TP_NO_ERROR;
  int ierror;


  //  vector<string> nodes,leaves;
  vector<string> leaves;
  //  tp_error=tp->getNodes(nodes);
  tp_error=tp->getLeaves(leaves,1);

  vector<string>::iterator leaf_iter;
  string movable_string;
  // 移動可能フラグ検索
  leaf_iter = find(leaves.begin(),leaves.end(),ATT_NAME_MOVABLE);

  if(leaf_iter!=leaves.end()) {
    tp_error=tp->getValue((*leaf_iter),movable_string);

    m_movable = tp->convertBool(movable_string,&ierror);      
  }

  // グループ名が重複していないか確認
  // for tp
  string current_node;
  tp->currentNode(current_node);

  //cout << __FUNCTION__ << " current_node = "  << current_node <<endl;
  string pg_name = current_node;

  string    parent_path = "";
  if (parent != NULL)       parent_path = parent->acq_fullpath();
  POLYLIB_STAT ret = polylib->check_group_name(pg_name, parent_path);
  if (ret != PLSTAT_OK)     return ret;

  // ファイルpathの処理（複数ファイルの可能性あり）
  string fname = "";

  leaf_iter = find(leaves.begin(),leaves.end(),"filepath[0]");
  if(leaf_iter != leaves.end()){
    //filepath が複数の場合

    int index=0;
 
    while(leaf_iter != leaves.end()){ //end　にいかなければ。
      stringstream ss;
      string tmpstring=ATT_NAME_PATH;

      ss << tmpstring <<"["<<index<<"]";
      ss >> tmpstring;
#ifdef DEBUG      
      PL_DBGOS << __FUNCTION__<< " multi stl files "<< tmpstring << " "<<*leaf_iter<<endl;
#endif //DEBUG      

      leaf_iter = find(leaf_iter,leaves.end(),tmpstring);
      if(leaf_iter == leaves.end()) break;
      tp_error=tp->getValue((*leaf_iter),fname);

#ifdef DEBUG      
      PL_DBGOS << __FUNCTION__ << " STL/NPTfiles " << index <<" " << fname <<endl;
#endif //DEBUG      

      string format = PolygonIO::input_file_format(fname);
      if (format.empty()) {
        PL_ERROSH << "[ERROR]PolygonGroup::setup_attribute():Unknown"
          << "extention: fname[]=" << fname << endl;
        return PLSTAT_UNKNOWN_FILE_FORMAT;
      }             
    
      m_polygon_files.insert(map<string, string>::value_type(fname, format));
      index++;
      leaf_iter++;
    }
   
  } else { //filepath が単数の場合
    leaf_iter = find(leaves.begin(),leaves.end(),ATT_NAME_PATH);
    if(leaf_iter!=leaves.end()) {
      tp_error=tp->getValue((*leaf_iter),fname);

#ifdef DEBUG
      PL_DBGOS << __FUNCTION__ << " STL/NPTfile "  << fname <<endl;
#endif // DEBUG
      string format = PolygonIO::input_file_format(fname);
      if (format.empty()) {
        PL_ERROSH << "[ERROR]PolygonGroup::setup_attribute():Unknown"
            << "extention: fname=" << fname << endl;
        return PLSTAT_UNKNOWN_FILE_FORMAT;
      }             
    
      m_polygon_files.insert(map<string, string>::value_type(fname, format));
    }
  }


    // 親を設定
  if (parent != NULL) {
    m_parent      = parent;
    m_parent_path = parent->acq_fullpath();
    parent->add_children(this);
  }

  // その他の属性を設定
  // for tp
  m_name = pg_name;

  m_internal_id = create_global_unique_id();

  return PLSTAT_OK;
}


// protected //////////////////////////////////////////////////////////////////
// ユーザ定義属性設定
//     ユーザ定義属性ノード配下の処理
POLYLIB_STAT PolygonGroup::setup_user_attribute (
    const std::string&      node_label,
    TextParser*             tp
) 
{
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::setup_user_attribute in." <<endl;
#endif
    TextParserError tp_error = TP_NO_ERROR;
    int error;

    //  下位ノードに移動
    tp_error=tp->changeNode(node_label);
    if(tp_error!=TP_NO_ERROR){
    PL_ERROSH << "[ERROR]PolygonGroup::setup_user_attribute():"
              << " TextParser error " 
              << tp->TextParserErrorHandler(tp_error,"can not move to ") 
              << node_label << endl;
        return PLSTAT_CONFIG_ERROR;
    }

    // 属性項目名取得
    vector<string> labels;
    tp_error=tp->getLabels(labels,1);
    //PL_DBGOSH << "  labels.size()="<<labels.size() <<endl;

    // 属性項目数ループ
    for( int i=0; i<labels.size(); i++ ) {
        UsrAtr usr_atr;
        usr_atr.key = labels[i];
        tp_error=tp->getValue(usr_atr.key,usr_atr.value);
        //PL_DBGOSH << "    i="<<i<<" usr_atr.key="<<usr_atr.key<<" usr_atr.value="<<usr_atr.value <<endl;
        m_atr.push_back( usr_atr );
    }

    // 下位ノードから復帰
    tp->changeNode("..");

    return PLSTAT_OK;

}


// public //////////////////////////////////////////////////////////////////
POLYLIB_STAT
PolygonGroup::init_check_leaped()
{
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::init_check_leaped() in. " << endl;
#endif

    // 動かないポリゴングループならば何もしないで終了
    if( !m_movable || m_tri_list==NULL || m_tri_list->size()==0 ) {
#ifdef DEBUG
        PL_DBGOSH << "PolygonGroup::init_check_leaped() out. 0 " << endl;
#endif
        return PLSTAT_OK;
    }

    // move後と比較するために三角形ポリゴンリストのディープコピーを保存
    m_trias_before_move = new vector<Triangle*>;
    copy_polygons( *m_tri_list, *m_trias_before_move );

#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::init_check_leaped() out. " << endl;
#endif

    return PLSTAT_OK;
}

// public //////////////////////////////////////////////////////////////////
POLYLIB_STAT
PolygonGroup::check_leaped(
    std::vector< Vec3<PL_REAL> >& origin,
    std::vector< Vec3<PL_REAL> >& cell_size
)
{
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::check_leaped() in. " << endl;
#endif
    unsigned int i, j;
    vector<Triangle*>* p_trias = get_triangles();

    // 動かないポリゴングループならば何もしないで終了
    if( !m_movable || p_trias==NULL || p_trias->size()==0 ) return PLSTAT_OK;

    // move前の三角形と座標を比較。
    for( i=0; i<p_trias->size(); i++ ) {

        for( j=0; j<3; j++ ) {
            // 隣接セルより遠方へ移動してたらcerrにメッセージを出力。

            /*
            if( is_far( origin, cell_size, p_trias->at(i)->get_vertexes()[j],
                        m_trias_before_move->at(i)->get_vertexes()[j]        ) ) {
                PL_ERROSH << "[ERROR]PolygonGroup::check_leaped():Leaped Vertex"
                          << " Detected. GroupID:" << m_internal_id
                          << " TriaID:" << p_trias->at(i)->get_id()
                          << " before:(" << m_trias_before_move->at(i)->get_vertexes()[j]
                          << ") after:(" << p_trias->at(i)->get_vertexes()[j]
                          << ")" << endl;
            }
            */
                // 複数領域に対応する
                //     １領域でも範囲内であればOKとする
            bool err = true;
            for(unsigned int k=0; k<origin.size(); k++ ) {
                bool far = is_far( origin[k], cell_size[k], 
                                   p_trias->at(i)->get_vertexes()[j],
                                   m_trias_before_move->at(i)->get_vertexes()[j] );
                if( !far ) {
                    err = false;
                    break;
                }
            }

            if( err ) {
                PL_ERROSH << "[ERROR]PolygonGroup::check_leaped():Leaped Vertex"
                          << " Detected. GroupID:" << m_internal_id
                          << " TriaID:" << p_trias->at(i)->get_id()
                          << " before:(" << m_trias_before_move->at(i)->get_vertexes()[j]
                          << ") after:(" << p_trias->at(i)->get_vertexes()[j]
                          << ")" << endl;
            }

        }
        
        // 移動前三角形インスタンスはもう不要なので削除
        delete m_trias_before_move->at(i);
    }

    // あとしまつ
    delete m_trias_before_move;
    m_trias_before_move = NULL;
#ifdef DEBUG
    PL_DBGOSH << "PolygonGroup::check_leaped() out. " << endl;
#endif

    return PLSTAT_OK;
}

// protected //////////////////////////////////////////////////////////////////
bool
PolygonGroup::is_far(
    Vec3<PL_REAL>& origin,
    Vec3<PL_REAL>& cell_size,
    Vec3<PL_REAL>& pos1,
    Vec3<PL_REAL>& pos2
)
{
    for( int i=0; i<3; i++ ) {
        // pos1所属ボクセルの起点座標を求める
        PL_REAL p;
        PL_REAL dist = pos1[i] - origin[i];
        if( dist >= 0 ) {
            p = origin[i] + ((int(dist / cell_size[i])) * cell_size[i]);
        }
        else {
            if( fmodf(dist,cell_size[i]) == 0 ) {
                p = origin[i] + ((int(dist / cell_size[i])) * cell_size[i]);
            } else {
                p = origin[i] + ((int(dist / cell_size[i]) - 1) * cell_size[i]);
            }
        }

        // 隣接ボクセルまで含んだ領域のmin/max
        PL_REAL min = p - cell_size[i];
        PL_REAL max = p + cell_size[i] * 2;

        // pos2がmin-max間に含まれなければ真
        if( pos2[i] < min || pos2[i] > max ) return true;
    }
    return false;
}

// private //////////////////////////////////////////////////////////////////
void PolygonGroup::init_tri_list()
{
    if (m_tri_list == NULL) {
        m_tri_list = new vector<Triangle*>;
    }
    else {
        vector<Triangle*>::iterator itr;
        for (itr = m_tri_list->begin(); itr != m_tri_list->end(); itr++) {
            delete *itr;
        }
        m_tri_list->clear();
    }
}

// private //////////////////////////////////////////////////////////////////
void PolygonGroup::delete_tri_list()
{
    if (m_tri_list != NULL) {

        vector<Triangle*>::iterator itr;
        for (itr = m_tri_list->begin(); itr != m_tri_list->end(); itr++) {
            delete *itr;
        }
        delete m_tri_list;
        m_tri_list = NULL;
    }
}

// private //////////////////////////////////////////////////////////////////
char *PolygonGroup::mk_polygons_fname(
    const string&       rank_no,
    const string&       extend,
    const string&       format
) {
    char        fname1[1024];
    //string        prefix;
    string      extension;  // 拡張子
    static char fname2[1024];

    // グループ名のフルパスを取得して、/を_に置き換え
    strcpy(fname1, acq_fullpath().c_str());

    //cout << __FUNCTION__ << " acq_fullpath() " <<acq_fullpath()<<endl;
    
    for (int i = 0; i < (int)strlen(fname1); i++) {
        if (fname1[i] == '/')   fname1[i] = '_';
    }
#ifdef DEBUG
    PL_DBGOS << __FUNCTION__ << " fname1 " <<fname1<<endl;
#endif //  DEBUG

    extension = PolygonIO::get_extension_format( format );

    if (rank_no == "") {
        sprintf(fname2, "%s_%s.%s", fname1, extend.c_str(), extension.c_str());
    }
    else {
        sprintf(fname2, "%s_%s_%s.%s", fname1, rank_no.c_str(), extend.c_str(), 
                                                                extension.c_str());
    }

#ifdef DEBUG
    PL_DBGOS << __FUNCTION__ << " fname2 " <<fname2<<endl;
#endif //DEBUG
    return fname2;
}

// protected //////////////////////////////////////////////////////////////////
char *PolygonGroup::mk_polygons_fname(
    const string&       rank_no,
    const string&       extend,
    const string&       format,
    map<string,string>& polygons_fname_map
) {
    char        fname1[1024];
    string      extension;  // 拡張子
    static char fname2[1024];

    // グループ名のフルパスを取得して、/を_に置き換え
    strcpy(fname1, acq_fullpath().c_str());

    //cout << __FUNCTION__ << " acq_fullpath() " <<acq_fullpath()<<endl;
    
    for (int i = 0; i < (int)strlen(fname1); i++) {
        if (fname1[i] == '/')   fname1[i] = '_';
    }

    //cout << __FUNCTION__ << " fname1 " <<fname1<<endl;

    extension = PolygonIO::get_extension_format( format );

    if (rank_no == "") {
        sprintf(fname2, "%s_%s.%s", fname1, extend.c_str(), extension.c_str());
    }
    else {
        sprintf(fname2, "%s_%s_%s.%s", fname1, rank_no.c_str(), extend.c_str(), 
                                                                extension.c_str());
    }

    //cout << __FUNCTION__ << " fname2 " <<fname2<<endl;

    string tmp_fname = fname2;
    polygons_fname_map.insert(map<string,string>::value_type(acq_fullpath(),tmp_fname));

    return fname2;
}

// protected //////////////////////////////////////////////////////////////////
int PolygonGroup::create_global_unique_id()
{
    static int global_id = 0;
    return global_id++;
}

