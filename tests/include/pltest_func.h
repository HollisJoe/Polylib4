/* -- Mode: c++ --*/
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

#ifndef polylib_test_func_h
#define polylib_test_func_h
#include "Polylib.h"

#define INLINE inline

using namespace Vec3class;

namespace PolylibNS {

///
/// ポリゴンの頂点座標の平均値取得
///   テスト時の妥当性検証に使用する
///
///  @param[in]  tri_list     ポリゴン
///  @param[out] x_ave        X座標の平均値
///  @param[out] y_ave        Y座標の平均値
///  @param[out] z_ave        Z座標の平均値
///
INLINE void  pltest_get_polygon_average(
           const std::vector<Triangle*>& tri_list,
           PL_REAL& x_ave, PL_REAL& y_ave, PL_REAL& z_ave
         )
{
    x_ave=0.0;  y_ave=0.0;  z_ave=0.0;

    if( tri_list.size() == 0 )  {
        return;
    }

    for(int i=0; i<tri_list.size(); i++ ) {
        Vec3<PL_REAL>* vertex = tri_list[i]->get_vertexes();

        for( int j=0; j<3; j++ ) {
            x_ave += vertex[j].x;
            y_ave += vertex[j].y;
            z_ave += vertex[j].z;
        }
    }

    x_ave = x_ave / (3*tri_list.size());
    y_ave = y_ave / (3*tri_list.size());
    z_ave = z_ave / (3*tri_list.size());

}


} //namespace PolylibNS

#endif // polylib_test_func_h
