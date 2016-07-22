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

#ifndef polylib_vtree_h
#define polylib_vtree_h

#include "common/BBox.h"
#include "common/PolylibStat.h"
#include "common/PolylibCommon.h"
//#include "polygons/Vertex.h"
#include "polygons/Triangle.h"
#include "polygons/NptTriangle.h"

#include <vector>
//#define DEBUG_VTREE
namespace PolylibNS {

//class Vertex;
//class BBox;
//class Triangle;
//class VNode;
//class VElement;

class VElement {
public:
    ///
    /// コンストラクタ。
    ///
    /// @param[in] tri ポリゴン情報のポインタ。
    /// @attention ポインタを格納するが、参照のみ。deleteは行わない。
    ///
    VElement(
        Triangle* tri
        );


    //VElement();

    ///
    /// デストラクタ。
    ///
    virtual ~VElement() {}


    //=======================================================================
    // Setter/Getter
    //=======================================================================
    ///
    /// triangle。
    ///
    Triangle* get_triangle() const
    {
        return m_tri;
    }

    ///
    /// Center position of bbox on triangle.
    /// 
    Vec3<PL_REAL> get_pos() const
    {
        return m_pos;
    }

    ///
    /// Bounding box of this triangle
    ///
    BBox get_bbox() const
    {
        return m_bbox;
    }

private:
    //=======================================================================
    // クラス変数
    //=======================================================================
    /// triangle
    Triangle    *m_tri;

    /// Center position of bbox on triangle.
    Vec3<PL_REAL>           m_pos;  // center of m_bbox

    /// Bounding box of this triangle
    BBox            m_bbox;

};



////////////////////////////////////////////////////////////////////////////
///  
/// VNodeクラス
/// KD木構造のノードクラスです。
///  
////////////////////////////////////////////////////////////////////////////

class VNode {
public:
    ///
    /// コンストラクタ。
    ///
    VNode();

    ///
    /// デストラクタ。
    ///
    ~VNode();

    ///
    /// ノードを２つの子供ノードに分割する。
    ///
    void split(const int& max_elem);

#ifdef USE_DEPTH
    ///
    /// ノードの深さ情報のダンプ。
    ///
    /// @param[in] n タブの出力個数。
    ///
    void dump_depth(int n);
#endif

    //=======================================================================
    // Setter/Getter
    //=======================================================================
    ///
    /// ノードがリーフかどうかの判定結果。
    ///
    /// @return true=リーフ/false=リーフでない。
    ///
    bool is_leaf() const
    {
        return (m_left == 0) ;
    }

    ///
    /// BBoxの値を取得。
    ///
    /// @return bbox。
    ///
    BBox get_bbox() const 
    {
        return m_bbox;
    }

    ///
    /// BBoxの値を設定。
    ///
    /// @param[in] bbox。
    ///
    void set_bbox(const BBox& bbox) 
    {
        m_bbox = bbox;
    }

    ///
    /// 検索用BBoxを取得。
    ///
    /// @return 検索用bbox。
    ///
    BBox get_bbox_search() const
    {
        return m_bbox_search;
    }

    ///
    /// このノードのBounding Boxを引数で与えられる要素を含めた大きさに変更する。
    ///
    /// @param[in] p 要素。
    ///
    void set_bbox_search(const VElement *p)
    {
        m_bbox_search.add(p->get_bbox().min);
        m_bbox_search.add(p->get_bbox().max);
    }

    ///
    /// 左のNodeを取得。
    ///
    /// @return 左のNode。
    ///
    VNode* get_left()
    {
        return m_left;
    }
    
    ///
    /// 右のNodeを取得。
    ///
    /// @return 右のNode。
    ///
    VNode* get_right()
    {
        return m_right;
    }

    ///
    /// Axisを取得
    ///
    /// @return axis
    ///
    AxisEnum get_axis() const
    {
        return m_axis;
    }

    ///
    /// Axisを設定
    ///
    /// @param[in] axis
    ///
    void set_axis(const AxisEnum axis)
    {
        m_axis = axis;
    }

    ///
    /// 要素のリストを取得
    ///
    /// @return 要素のリスト
    ///
    std::vector<VElement*>& get_vlist()
    {
        return m_vlist;
    }

    ///
    /// 木の要素を設定
    ///
    /// @param[in] elm
    ///
    void set_element(VElement* elm)
    {
        m_vlist.push_back(elm);
    }

    ///
    /// ノードが所持する要素の数を取得。
    ///
    /// @return 要素数。
    ///
    int get_elements_num() const
    {
        return m_vlist.size();
    }

#ifdef USE_DEPTH
    ///
    /// ノードの深さ情報を取得。
    ///
    /// @return ノードの深さ。
    ///
    int get_depth() const
    {
        return m_depth;
    }

    ///
    /// ノードの深さ情報の設定。
    ///
    /// @param[in] depth ノードの深さ。
    ///
    void set_depth(int depth)
    {
        m_depth = depth;
    }
#endif

private:
    //=======================================================================
    // クラス変数
    //=======================================================================
    /// 本ノードの左下側ノード。
    VNode                   *m_left;

    /// 本ノードの右下側ノード。
    VNode                   *m_right;

    /// KD木生成用のBouding Box。
    BBox                    m_bbox;

    /// KD木の軸の方向インデックス。
    AxisEnum                m_axis;

    /// ノードの管理する要素リスト。
    std::vector<VElement*>  m_vlist;

    /// KD木検索用のBouding Box。
    BBox                    m_bbox_search;

#ifdef USE_DEPTH
    /// ノードの深さ情報(未使用)。
    int                     m_depth;
#endif
};



////////////////////////////////////////////////////////////////////////////
///
/// クラス:VTree
/// リーフを三角形ポリゴンとするKD木クラスです。
///
////////////////////////////////////////////////////////////////////////////

class VTree {
public:
    ///
    /// コンストラクタ。
    ///
    /// @param[in] max_elem 最大要素数。
    /// @param[in] bbox     VTreeのbox範囲。
    /// @param[in] tri_list 木構造の元になるポリゴンのリスト。
    ///
    VTree(
        int         max_elem, 
        const BBox          bbox, 
        std::vector<Triangle*>  *tri_list
        );

    ///
    /// デストラクタ。
    ///
    ~VTree();

    ///
    /// 木構造を消去する
    ///
    void destroy();

    ///
    /// KD木探索により、指定矩形領域に含まれる三角形ポリゴンを抽出する
    ///
    ///  @param[in] bbox    検索範囲を示す矩形領域
    ///  @param[in] every   true:3頂点が全て検索領域に含まれるものを抽出
    ///                     false:1頂点でも検索領域に含まれるものを抽出
    ///  @return    抽出したポリゴンリストのポインタ
    ///                 vectorのポインタは内部でallocationされているので領域解放が必要
    ///  @attention MPIPolylib用のメソッドなので、ユーザは利用しないで下さい

    std::vector<Triangle*>* search(
        const BBox    *bbox, 
        bool    every
        ) const;

    ///
    /// KD木探索により、指定矩形領域に含まれるポリゴンを抽出する
    ///
    ///  @param[in]     bbox        検索範囲を示す矩形領域
    ///  @param[in]     every       true:3頂点が全て検索領域に含まれるものを抽出
    ///                             false:1頂点でも検索領域に含まれるものを抽出
    ///  @param[in,out] tri_list    抽出した三角形ポリゴンリストへのポインタ
    ///  @return    POLYLIB_STATで定義される値が返る
    ///
    POLYLIB_STAT search(
        std::vector<Triangle*>& tri_list,
        const BBox&             bbox, 
        bool                    every 
        ) const;

    ///
    /// KD木探索により、指定位置に最も近いポリゴンを検索する。
    ///
    ///  @param[in]     pos     指定位置
    ///  @return    検索されたポリゴン
    ///
    const Triangle* search_nearest(
        const Vec3<PL_REAL>&    pos
        ) const;

    ///
    /// KD木探索により、指定位置に最も近いポリゴンを検索する。
    ///
    ///  @param[in]     vn      検索対象のノードへのポインタ。
    ///  @param[in]     pos     指定位置
    ///  @return    検索されたポリゴン
    ///
    const Triangle* search_nearest_recursive(
        VNode       *vn, 
        const Vec3<PL_REAL>&    pos
        ) const;

    ///
    /// KD木クラスが利用しているメモリ量を返す。
    ///
    ///  @return    利用中のメモリ量(byte)
    ///
    size_t memory_size();

private:
    ///
    /// 三角形をKD木構造に組み込む際に、どのノードへ組み込むかを検索する。
    ///
    ///  @param[in]     vn      検索対象のノードへのポインタ。
    ///  @param[in]     elm     組み込む三角形。
    ///  @param[in,out] vnode   検索結果。
    ///
    void traverse(
        VNode       *vn, 
        VElement    *elm, 
        VNode       **vnode
        ) const;

    ///
    /// 三角形ポリゴンをKD木構造から検索する。
    ///
    ///  @param[in]     vn      検索対象のノードへのポインタ。
    ///  @param[in]     bbox    VTreeのbox範囲。
    ///  @param[in]     every   true:ポリゴンの頂点がすべて含まれるNodeを検索。
    ///                         false:それ以外。
    ///  @param[in,out] vlist   検索結果配列へのポインタ。
    ///
    void search_recursive(
        VNode                   *vn, 
        const BBox              &bbox, 
        bool                    every, 
        std::vector<VElement*>  *vlist
        ) const;

    ///
    /// 初期化処理
    ///
    ///  @param[in] max_elem    最大要素数。
    ///  @param[in] bbox        VTreeのbox範囲。
    ///  @param[in] tri_list    木構造の元になるポリゴンのリスト。
    ///  @return    POLYLIB_STATで定義される値が返る。
    ///
    POLYLIB_STAT create(
        int                 max_elem, 
        const BBox              bbox, 
        std::vector<Triangle*>  *tri_list
        );

    ///
    /// KD木の総ノード数と総ポリゴン数を数える。
    ///
    ///  @param[in]  vnode      親ノードへのポインタ。
    ///  @param[out] node_cnt   ノード数。
    ///  @param[out] poly_cnt   ポリゴン数。
    ///
    void node_count(
        VNode           *vnode, 
        unsigned int    *node_cnt, 
        unsigned int    *poly_cnt
        );

    //=======================================================================
    // クラス変数
    //=======================================================================
    /// ルートノードへのポインタ。
    VNode   *m_root;

    /// リーフノードが所持できる最大要素数。
    int     m_max_elements;

#ifdef DEBUG_VTREE
    std::vector<VNode*> m_vnode;
#endif
};

}

#endif  // vtree_h

