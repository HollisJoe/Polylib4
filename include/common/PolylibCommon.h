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

#ifndef polylib_common_h
#define polylib_common_h

#include <iostream>
#include "common/PolylibDefine.h"

namespace PolylibNS {

////////////////////////////////////////////////////////////////////////////
/// 
/// 三角形IDファイルフォーマット
///
////////////////////////////////////////////////////////////////////////////
typedef enum {
    ID_BIN,     ///< バイナリ形式で入出力を行う。
    ID_ASCII    ///< アスキー形式で入出力を行う。
} ID_FORMAT;


////////////////////////////////////////////////////////////////////////////
/// 
/// Vec3型とプリミティブ型の型変換
///
////////////////////////////////////////////////////////////////////////////

//  Vec3<PL_REAL> v  -> PL_REAL r[3]
#define VEC3_TO_REAL(v,r)   \
            {   \
                (r)[0]=(v).x;  (r)[1]=(v).y;  (r)[2]=(v).z; \
            }

//  PL_REAL r[3] -> Vec3<PL_REAL> v
#define REAL_TO_VEC3(r,v)   \
            {   \
                (v).x=(r)[0];  (v).y=(r)[1];  (v).z=(r)[2]; \
            }

//  Vec3<PL_REAL> v[3]  -> PL_REAL r1[3],r2[3],r3[3]
#define VEC3_3_TO_REAL(v,r1,r2,r3)  \
            {   \
                (r1)[0]=(v)[0].x;  (r1)[1]=(v)[0].y;  (r1)[2]=(v)[0].z; \
                (r2)[0]=(v)[1].x;  (r2)[1]=(v)[1].y;  (r2)[2]=(v)[1].z; \
                (r3)[0]=(v)[2].x;  (r3)[1]=(v)[2].y;  (r3)[2]=(v)[2].z; \
            }

//  PL_REAL r1[3],r2[3],r3[3] -> Vec3<PL_REAL> v[3]
#define REAL_TO_VEC3_3(r1,r2,r3,v)  \
            {   \
                (v)[0].x=(r1)[0];  (v)[0].y=(r1)[1];  (v)[0].z=(r1)[2]; \
                (v)[1].x=(r2)[0];  (v)[1].y=(r2)[1];  (v)[1].z=(r2)[2]; \
                (v)[2].x=(r3)[0];  (v)[2].y=(r3)[1];  (v)[2].z=(r3)[2]; \
            }

//  Vec3<PL_REAL> v[3]  -> PL_REAL r[9]
#define VEC3_3_TO_REAL9(v,r)    \
            {   \
                (r)[0]=(v)[0].x;  (r)[1]=(v)[0].y;  (r)[2]=(v)[0].z; \
                (r)[3]=(v)[1].x;  (r)[4]=(v)[1].y;  (r)[5]=(v)[1].z; \
                (r)[6]=(v)[2].x;  (r)[7]=(v)[2].y;  (r)[8]=(v)[2].z; \
            }

//  PL_REAL r[9]-> Vec3<PL_REAL> v[3]
#define REAL9_TO_VEC3_3(r,v) \
            {   \
                (v)[0].x=(r)[0];  (v)[0].y=(r)[1];  (v)[0].z=(r)[2]; \
                (v)[1].x=(r)[3];  (v)[1].y=(r)[4];  (v)[1].z=(r)[5]; \
                (v)[2].x=(r)[6];  (v)[2].y=(r)[7];  (v)[2].z=(r)[8]; \
            }


////////////////////////////////////////////////////////////////////////////
///
/// デバッグ出力先、エラー時出力先
///
////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
#define PL_DBGOS    std::cout<<__FILE__<<":"<<__FUNCTION__<<":L"<<__LINE__<<":"
#define PL_DBGOSH   std::cout<<gs_rankno<<"PL:"<<__FILE__<<":"<<__FUNCTION__<<":L"<<__LINE__<<":"
#define PL_ERROS    std::cerr<<__FILE__<<":"<<__FUNCTION__<<":L"<<__LINE__<<":"
#define PL_ERROSH   std::cerr<<gs_rankno<<"PL:"<<__FILE__<<":"<<__FUNCTION__<<":L"<<__LINE__<<":"
#else
#define PL_DBGOS    std::cout
#define PL_DBGOSH   std::cout<<gs_rankno<<"PL:"
#define PL_ERROS    std::cerr
#define PL_ERROSH   std::cerr<<gs_rankno<<"PL:"
#endif

////////////////////////////////////////////////////////////////////////////
///
/// デバッグ出力用ランク番号グローバル文字列
///
////////////////////////////////////////////////////////////////////////////
extern std::string gs_rankno;

} // end of namespace PolylibNS
#endif // polylib_common_h
