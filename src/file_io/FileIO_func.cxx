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

//  <string> for windows, <string.h> for linux
#include <string>
#include <string.h>

#include <fstream>
#include <vector>
#include <iomanip>
#include "common/tt.h"
#include "polygons/Triangle.h"
#include "file_io/PolygonIO.h"
#include "file_io/FileIO_func.h"

namespace PolylibNS {

using namespace std;

#define SCIENTIFIC_OUT      0
#define STL_HEAD            80      // header size for STL binary
#define STL_BUFF_LEN        256
#define TT_OTHER_ENDIAN     1
#define TT_LITTLE_ENDIAN    2
#define TT_BIG_ENDIAN       3

// プライベート（static）関数 プロトタイプ宣言
static void tt_invert_byte_order(void* _mem, int size, int n);
static int  tt_check_machine_endian();
static void tt_read(istream& is, void* _data, int size, int n, int inv);
static void tt_write(ostream& os, const void* _data, int size, int n, int inv);

//////////////////////////////////////////////////////////////////////////////

//**************************************************************
//   STLファイル用
//**************************************************************

// ASCIIモードのSTLファイルを読み込み、tri_listに三角形ポリゴン情報を設定する
POLYLIB_STAT stl_a_load(
    vector<Triangle*>   *tri_list, 
    const string&       fname,
    int                 *num_tri,
    PL_REAL         scale
) 
{
    ifstream is(fname.c_str());
    if (is.fail()) {
        PL_ERROSH << "[ERROR]FileIO_func:stl_a_load():Can't open " << fname << endl;
        return PLSTAT_STL_IO_ERROR;
    }

    int num_read = -1;
    bool eof;
    POLYLIB_STAT ret = stl_a_load_read( is, *tri_list, num_read, *num_tri, eof, scale );

#ifdef DEBUG
PL_DBGOSH << "stl_a_load() ret="<<ret<<" total="<<tri_list->size()<<" eof="<<eof <<endl;
#endif
    return ret;
}


// ASCIIモードのSTLファイルを指定個数分読み込み
POLYLIB_STAT stl_a_load_read(
    ifstream&                is,
    std::vector<Triangle*>&  tri_list,
    int                      num_read,
    int&                     num_tri,
    bool&                    eof,
    PL_REAL                  scale
) 
{
    num_tri   = 0; // 実際に読み込んだ個数
    int n_vtx = 0;    // 頂点の位置( 0-2 )

    eof = false;

    string token;
    Vec3<PL_REAL> nml;
    Vec3<PL_REAL> vtx[3];
    //while (is >> token && !is.eof()) {
    is >> token;
    while ( !is.eof() ) {
        if (token == "solid") {
            string s;
            is >> s;   // solid name  のname部分をread
            // nameは省略可能であり、次のtoken を読んでしまうことがある
            if( s == "facet") {
                token = s;
                continue;
            }
        }
        else if (token == "facet") {
            n_vtx = 0;

            string s;
            is >> s;
            is >> nml;
            nml.normalize();
        }
        else if (token == "vertex") {
            Vec3<PL_REAL> v;
            is >> v;
            if (n_vtx < 3) {
                vtx[n_vtx] = v * scale;
            }
            n_vtx++;
        }
        else if (token == "outer") {
            string s;
            is >> s;
        }
        else if (token == "endloop") {
        }
        else if (token == "endfacet") {
            if (n_vtx == 3) {
                Triangle *tri = new Triangle(vtx, nml); // IDは内部で採番
                tri_list.push_back(tri);
                num_tri++;

                if( num_read > 0 ) {  // -1 の場合、チェックなし 
                    // 指定数読み込んだか否か判定
                    if( num_tri >= num_read ) {
                        break;  // 読み込み終了
                    }
                }
            }
        }
        else if (token == "endsolid") {
            string s;
            is >> s;
        }

        // 次のトークンの読み出し
        is >> token;
    }

    if (!is.eof() && is.fail()) {
        PL_ERROSH << "[ERROR]FileIO_func:stl_a_load_read():Error in loading" <<endl;
        return PLSTAT_STL_IO_ERROR;
    }

    if ( is.eof() ) {
        eof = true;
    }

#ifdef DEBUG
PL_DBGOSH <<"stl_a_load_read() num_read="<<num_read<<" num_tri="<<num_tri << endl;
#endif
    return PLSTAT_OK;
}


//////////////////////////////////////////////////////////////////////////////

// 三角形ポリゴン情報をASCIIモードでSTLファイルに書き出す
POLYLIB_STAT stl_a_save(
    vector<Triangle*>   *tri_list, 
    const string&       fname
) 
{
#ifdef DEBUG
    PL_DBGOSH <<  "stl_a_load() start  fname="<<fname<<" tri_list->size()=" << tri_list->size() <<endl;
#endif
    ofstream os(fname.c_str());

    if (os.fail()) {
        PL_ERROSH << "[ERROR]FileIO_func:stl_a_save():Can't open " << fname << endl;
        return PLSTAT_STL_IO_ERROR;
    }

    os << "solid " << "model1" << endl;

    vector<Triangle*>::iterator itr;
    for (itr = tri_list->begin(); itr != tri_list->end(); itr++) {
#if SCIENTIFIC_OUT
        os  << "  facet " << "normal " << setprecision(6) << scientific 
            << (*itr)->get_normal() << endl;
#else
        os  << "  facet " << "normal " << setprecision(6) 
            << (*itr)->get_normal() << endl;
#endif
        os << "    outer " << "loop" << endl;
        for (int j = 0; j < 3; j++) {
#if SCIENTIFIC_OUT
            os  << "      vertex " << setprecision(6) << scientific 
                << (*itr)->get_vertexes()[j] << endl;
#else
            os  << "      vertex " << setprecision(6) 
                << (*itr)->get_vertexes()[j] << endl;
#endif
        }
        os << "    endloop" << endl;
        os << "  endfacet" << endl;
    }
    os << "endsolid " << "model1" << endl;

    if (!os.eof() && os.fail()) {
        PL_ERROSH << "[ERROR]FileIO_func:stl_a_save():Error in saving: " << fname << endl;
        return PLSTAT_STL_IO_ERROR;
    }

    return PLSTAT_OK;
}

//////////////////////////////////////////////////////////////////////////////

// バイナリモードのSTLファイルを読み込み、tri_listに三角形ポリゴン情報を設定する
//      ファイルからは単精度実数(4byte)で読み込む
POLYLIB_STAT stl_b_load(
    vector<Triangle*>   *tri_list, 
    const string&       fname,
    //int                   *total,
    int                 *num_tri,
    PL_REAL         scale
) 
{
    ifstream ifs(fname.c_str(), ios::in | ios::binary);
    if (ifs.fail()) {
        PL_ERROSH << "[ERROR]FileIO_func:stl_b_load():Can't open " << fname << endl;
        return PLSTAT_STL_IO_ERROR;
    }

    stl_b_load_read_head( ifs );
    int num_read = -1;
    bool eof;
    POLYLIB_STAT ret = stl_b_load_read( ifs, *tri_list, num_read, *num_tri, eof, scale );

#ifdef DEBUG
PL_DBGOSH << "npt_b_load() ret="<<ret<<" total="<<tri_list->size()<<" eof="<<eof <<endl;
#endif
    return ret;
}

// バイナリモードのSTLファイルのヘッダ部を読み飛ばす
POLYLIB_STAT stl_b_load_read_head(
        ifstream&                ifs
    ) 
{
    int inv = tt_check_machine_endian() == TT_LITTLE_ENDIAN ? 0 : 1;

    uint    element = 0;
    char buf[STL_HEAD];
    for (int i = 0; i < STL_HEAD; i++) buf[i] = 0;
    tt_read(ifs, buf, sizeof(char), STL_HEAD, inv);
    tt_read(ifs, &element, sizeof(uint), 1, inv);

    return PLSTAT_OK;
}


// バイナリモードのSTLファイルからポリゴン指定個数分読み込み
//      ファイルからは単精度実数(4byte)で読み込む
POLYLIB_STAT stl_b_load_read(
    ifstream&           ifs,
    vector<Triangle*>&  tri_list, 
    int                 num_read,
    int&                num_tri,
    bool&               eof,
    PL_REAL             scale
) 
{
    int inv = tt_check_machine_endian() == TT_LITTLE_ENDIAN ? 0 : 1;

    num_tri  = 0;
    uint    element = 0;
    ushort padding = 0;

    eof = false;
    while( !eof )  {
        // one plane normal
        float nml[3];
        tt_read(ifs, nml, sizeof(float), 3, inv);

        if ( ifs.eof() ) {
            eof = true;
            break;
        }
        
        Vec3<PL_REAL> normal;
        normal.x = nml[0];
        normal.y = nml[1];
        normal.z = nml[2];

        // three vertices
        Vec3<PL_REAL> vertex[3];

        for (int j = 0; j < 3; j++) {
            float vtx[3];
            tt_read(ifs, vtx, sizeof(float), 3, inv);
            vtx[0] = vtx[0] * scale;
            vtx[1] = vtx[1] * scale;
            vtx[2] = vtx[2] * scale;

            vertex[j].x = vtx[0];
            vertex[j].y = vtx[1];
            vertex[j].z = vtx[2];

        }

        // ２バイト予備領域
        tt_read(ifs, &padding, sizeof(ushort), 1, inv);

        Triangle *tri = new Triangle(vertex, normal);   // IDは内部で採番

        // ２バイト予備領域をユーザ定義IDとして利用(Polylib-2.1より)
        tri->set_exid( (int)padding );

        tri_list.push_back(tri);
        num_tri++;

        if( num_read > 0 ) {  // -1 の場合、チェックなし 
            // 指定数読み込んだか否か判定
            if( num_tri >= num_read ) {
                break;  // 読み込み終了
            }
        }

    }

    if (!ifs.eof() && ifs.fail()) {
        PL_ERROSH << "[ERROR]FileIO_func:stl_b_load_read():Error in loading" <<endl;
        return PLSTAT_STL_IO_ERROR;
    }

    if ( ifs.eof() ) {
        eof = true;
    }

    return PLSTAT_OK;
}


//////////////////////////////////////////////////////////////////////////////

// 三角形ポリゴン情報をバイナリモードでSTLファイルに書き出す
//      ファイルへは単精度実数(4byte)で書き込む
POLYLIB_STAT stl_b_save(
    vector<Triangle*>   *tri_list, 
    const string&       fname
) {
    ofstream ofs(fname.c_str(), ios::out | ios::binary);
    if (ofs.fail()) {
        PL_ERROSH << "[ERROR]FileIO_func:stl_b_save():Can't open " << fname << endl;
        return PLSTAT_STL_IO_ERROR;
    }

    int inv = tt_check_machine_endian() == TT_LITTLE_ENDIAN ? 0 : 1;

    uint element = tri_list->size();

    char buf[STL_HEAD];
    for (int i = 0; i < STL_HEAD; i++) buf[i] = 0;
    strcpy(buf, "default");
    tt_write(ofs, buf, 1, STL_HEAD, inv);
    tt_write(ofs, &element, sizeof(uint), 1, inv);

    for (uint m = 0; m < element; m++) {
        Vec3<PL_REAL>* vertex = (*tri_list)[m]->get_vertexes();
        Vec3<PL_REAL>  norm = (*tri_list)[m]->get_normal();

        // PL_REAL -> float 変換
        float vtx[3][3];
        float nml[3];

        vtx[0][0] = vertex[0].x;  vtx[0][1]= vertex[0].y;  vtx[0][2] = vertex[0].z;
        vtx[1][0] = vertex[1].x;  vtx[1][1]= vertex[1].y;  vtx[1][2] = vertex[1].z;
        vtx[2][0] = vertex[2].x;  vtx[2][1]= vertex[2].y;  vtx[2][2] = vertex[2].z;
        nml[0] = norm.x;  nml[1]= norm.y;  nml[2] = norm.z;

        // write plane normal
        tt_write(ofs, nml, sizeof(float), 3, inv);

        // write vertex
        for(int j=0; j<3; j++ ) {
            tt_write(ofs, vtx[j], sizeof(float), 3, inv);
        }

        // ２バイト予備領域にユーザ定義IDを記録(Polylib-2.1より)
        int exid = (*tri_list)[m]->get_exid();
        ushort exid_2b = exid;
        tt_write(ofs, &exid_2b, sizeof(ushort), 1, inv);
    }

    if (!ofs.eof() && ofs.fail()) {
        PL_ERROSH << "[ERROR]FileIO_func:stl_b_load():Error in saving: " << fname << endl;
        return PLSTAT_STL_IO_ERROR;
    }

    return PLSTAT_OK;
}

//////////////////////////////////////////////////////////////////////////////

// STLファイルを読み込みバイナリかアスキーかを判定する
//     先頭行に"solid"文字列があるかどうかで判断する
bool is_stl_a(const string& path)
{
    const char  dcs[] = "\n\r";
    char        buff[STL_BUFF_LEN];
    size_t      i = 0;
    char        c;

    ifstream ifs(path.c_str());
    if (!ifs)                  return false;

    // ファイル内容の一部を読み込み
    memset(buff, 0, STL_BUFF_LEN);
    while (i < STL_BUFF_LEN - 1 && ifs.get(c) && !strchr(dcs, c)) {
        buff[i++] = c;
    }
    while (ifs.get(c) && strchr(dcs, c))
        ;
    if (ifs.good()) ifs.putback(c);

    if (ifs.eof())              return false;

    // ASCII判定
    if (strstr(buff, "solid"))  return true;
    else                        return false;
}



//**************************************************************
//   長田パッチファイル用
//**************************************************************

// ASCIIモードの長田パッチファイルを読み込み、tri_listに三角形ポリゴン情報を設定する
POLYLIB_STAT npt_a_load(
        std::vector<Triangle*>   *tri_list, 
        const std::string&          fname,
        int          *num_tri,
        PL_REAL      scale
    )
{
    ifstream is(fname.c_str());
    if (is.fail()) {
        PL_ERROSH << "[ERROR]FileIO_func:npt_a_load():Can't open " << fname << endl;
        return PLSTAT_NPT_IO_ERROR;
    }

    npt_a_load_read_head( is );
    int num_read = -1;
    bool eof;
    POLYLIB_STAT ret = npt_a_load_read( is, *tri_list, num_read, *num_tri, eof, scale );

#ifdef DEBUG
PL_DBGOSH << "npt_a_load() ret="<<ret<<" total="<<tri_list->size()<<" eof="<<eof <<endl;
#endif
    return ret;
}


// ASCIIモードの長田パッチファイルのヘッダ部を読み飛ばす
POLYLIB_STAT npt_a_load_read_head(
        ifstream&                is
    )
{
    // ファセット数読み込み
    int num_facets;
    is >> num_facets;
    
    return PLSTAT_OK;
}

// ASCIIモードの長田パッチファイルを読み込み、tri_listに三角形ポリゴン情報を設定する
POLYLIB_STAT npt_a_load_read(
        ifstream&                    is,
        //std::vector<NptTriangle*>&   tri_list, 
        std::vector<Triangle*>&   tri_list, 
        int                          num_read,
        int&                         num_tri,
        bool&                        eof,
        PL_REAL                      scale
    )
{
    num_tri  = 0; // 実際に読み込んだ個数

    string token;
    Vec3<float>   vtx[3];
    Vec3<float>   cp_side1_1;
    Vec3<float>   cp_side1_2;
    Vec3<float>   cp_side2_1;
    Vec3<float>   cp_side2_2;
    Vec3<float>   cp_side3_1;
    Vec3<float>   cp_side3_2;  
    Vec3<float>   cp_center;
    
    Vec3<PL_REAL>   vertex[3];
    NpatchParam     param;

    // ファセット数ループ
    while (is >> token && !is.eof()) {
        string s;

        if (token == "facet") {
        } else {
            PL_ERROSH << "[ERROR]FileIO_func:stl_a_load_read():Error in loading" << endl;
            return PLSTAT_STL_IO_ERROR;
        }
        
        is >> s;
        is >> vtx[0];
        is >> s;
        is >> vtx[1];
        is >> s;
        is >> vtx[2];

        is >> s;
        is >> cp_side1_1;
        is >> s;
        is >> cp_side1_2;
        is >> s;
        is >> cp_side2_1;
        is >> s;
        is >> cp_side2_2;
        is >> s;
        is >> cp_side3_1;
        is >> s;
        is >> cp_side3_2;
        is >> s;
        is >> cp_center;

        // float -> PL_REAL
        vertex[0].x=vtx[0].x;  vertex[0].y=vtx[0].y;  vertex[0].z=vtx[0].z;
        vertex[1].x=vtx[1].x;  vertex[1].y=vtx[1].y;  vertex[1].z=vtx[1].z;
        vertex[2].x=vtx[2].x;  vertex[2].y=vtx[2].y;  vertex[2].z=vtx[2].z;
        
        param.cp_side1_1.x=cp_side1_1.x;  param.cp_side1_1.y=cp_side1_1.y;  param.cp_side1_1.z=cp_side1_1.z;
        param.cp_side1_2.x=cp_side1_2.x;  param.cp_side1_2.y=cp_side1_2.y;  param.cp_side1_2.z=cp_side1_2.z;
        param.cp_side2_1.x=cp_side2_1.x;  param.cp_side2_1.y=cp_side2_1.y;  param.cp_side2_1.z=cp_side2_1.z;
        param.cp_side2_2.x=cp_side2_2.x;  param.cp_side2_2.y=cp_side2_2.y;  param.cp_side2_2.z=cp_side2_2.z;
        param.cp_side3_1.x=cp_side3_1.x;  param.cp_side3_1.y=cp_side3_1.y;  param.cp_side3_1.z=cp_side3_1.z;
        param.cp_side3_2.x=cp_side3_2.x;  param.cp_side3_2.y=cp_side3_2.y;  param.cp_side3_2.z=cp_side3_2.z;
        param.cp_center.x =cp_center.x;   param.cp_center.y =cp_center.y;   param.cp_center.z =cp_center.z;


        NptTriangle *tri = new NptTriangle( vertex, param );
        tri_list.push_back(tri);
        num_tri++;

        if( num_read > 0 ) {  // -1 の場合、チェックなし 
            // 指定数読み込んだか否か判定
            if( num_tri >= num_read ) {
                break;  // 読み込み終了
            }
        }
    }

    if (!is.eof() && is.fail()) {
      PL_ERROSH << "[ERROR]FileIO_func:npt_a_load_read():Error in loading" << endl;
      return PLSTAT_STL_IO_ERROR;
    }

    if ( is.eof() ) {
        eof = true;
    }

#ifdef DEBUG
PL_DBGOSH <<"npt_a_load_read() num_read="<<num_read<<" num_tri="<<num_tri << endl;
#endif
    return PLSTAT_OK;
}

// 三角形ポリゴン情報をASCIIモードで長田パッチファイルに書き出す
POLYLIB_STAT npt_a_save(
    std::vector<NptTriangle*>   *tri_list, 
    const std::string&          fname
)
{
#ifdef DEBUG
    PL_DBGOSH<<__FUNCTION__ << " npt_a_save in."<<endl;
#endif
    ofstream os(fname.c_str());

    if (os.fail()) {
        PL_ERROSH << "[ERROR]FileIO_func:npt_a_save():Can't open " << fname << endl;
        return PLSTAT_STL_IO_ERROR;
    }
    
    int num_tri = tri_list->size();
//#ifdef DEBUG
//    PL_DBGOSH<<__FUNCTION__ << " num_tri="<<num_tri<<endl;
//#endif

    os << num_tri << endl;

    for ( int i=0; i<num_tri; i++ ) {
        Vec3<PL_REAL>* vertex = (*tri_list)[i]->get_vertexes();
        const NpatchParam*   npatch = (*tri_list)[i]->get_npatch_param();
        
        
        os  << "facet" << endl;

#if SCIENTIFIC_OUT
        os  << "  vertex " << setprecision(6) << scientific << vertex[0] << endl;
        os  << "  vertex " << setprecision(6) << scientific << vertex[1] << endl;
        os  << "  vertex " << setprecision(6) << scientific << vertex[2] << endl;

        os  << "  coef1 " << setprecision(6) << scientific << npatch->cp_side1_1 << endl;
        os  << "  coef2 " << setprecision(6) << scientific << npatch->cp_side1_2 << endl;
        os  << "  coef3 " << setprecision(6) << scientific << npatch->cp_side2_1 << endl;
        os  << "  coef4 " << setprecision(6) << scientific << npatch->cp_side2_2 << endl;
        os  << "  coef5 " << setprecision(6) << scientific << npatch->cp_side3_1 << endl;
        os  << "  coef6 " << setprecision(6) << scientific << npatch->cp_side3_2 << endl;
        os  << "  coef7 " << setprecision(6) << scientific << npatch->cp_center  << endl;
#else
        os  << "  vertex " << setprecision(6) << vertex[0] << endl;
        os  << "  vertex " << setprecision(6) << vertex[1] << endl;
        os  << "  vertex " << setprecision(6) << vertex[2] << endl;

        os  << "  coef1 " << setprecision(6) << npatch->cp_side1_1 << endl;
        os  << "  coef2 " << setprecision(6) << npatch->cp_side1_2 << endl;
        os  << "  coef3 " << setprecision(6) << npatch->cp_side2_1 << endl;
        os  << "  coef4 " << setprecision(6) << npatch->cp_side2_2 << endl;
        os  << "  coef5 " << setprecision(6) << npatch->cp_side3_1 << endl;
        os  << "  coef6 " << setprecision(6) << npatch->cp_side3_2 << endl;
        os  << "  coef7 " << setprecision(6) << npatch->cp_center  << endl;
#endif
    }

    if (!os.eof() && os.fail()) {
        PL_ERROSH << "[ERROR]FileIO_func:npt_a_save():Error in saving: " << fname << endl;
        return PLSTAT_STL_IO_ERROR;
    }

    return PLSTAT_OK;
}


// バイナリモードの長田パッチファイルを読み込み、tri_listに三角形ポリゴン情報を設定する
//      ファイルからは単精度実数(4byte)で読み込む
POLYLIB_STAT npt_b_load(
        std::vector<Triangle*>   *tri_list, 
        const std::string&          fname,
        int                         *num_tri,
        PL_REAL                     scale
    )
{
    ifstream ifs(fname.c_str(), ios::in | ios::binary);
    if (ifs.fail()) {
        PL_ERROSH << "[ERROR]FileIO_func:npt_b_load():Can't open " << fname << endl;
        return PLSTAT_NPT_IO_ERROR;
    }

    npt_b_load_read_head( ifs );
    int num_read = -1;
    bool eof;
    POLYLIB_STAT ret = npt_b_load_read( ifs, *tri_list, num_read, *num_tri, eof, scale );

#ifdef DEBUG
PL_DBGOSH << "npt_b_load() ret="<<ret<<" total="<<tri_list->size()<<" eof="<<eof <<endl;
#endif
    return ret;
}


// バイナリモードの長田パッチファイルのヘッダ部を読み飛ばす
POLYLIB_STAT npt_b_load_read_head( ifstream& ifs )
{
    int inv = tt_check_machine_endian() == TT_LITTLE_ENDIAN ? 0 : 1;
    uint    element = 0;
    tt_read(ifs, &element, sizeof(uint), 1, inv);
	return PLSTAT_OK;
}

// バイナリモードの長田パッチファイルからポリゴン指定個数分読み込み
//      ファイルからは単精度実数(4byte)で読み込む
POLYLIB_STAT npt_b_load_read(
        ifstream&                ifs,
        std::vector<Triangle*>&  tri_list, 
        int                      num_read,
        int&                     num_tri,
        bool&                    eof,
        PL_REAL                  scale
    )
{
    int inv = tt_check_machine_endian() == TT_LITTLE_ENDIAN ? 0 : 1;
    num_tri  = 0;

    eof = false;
    while( !eof )  {
        float vtx[3][3];
        float param[7][3];

        // read vertex
        for(int j=0; j<3; j++ ) {
            tt_read(ifs, vtx[j], sizeof(float), 3, inv);
            if ( ifs.eof() ) {
               eof = true;
               break;
            }

            //  scaleをかける
            for(int k=0; k<3; k++ ) {
                vtx[j][k] *= scale;
            }
        }
        if( eof ) {
            break;
        }
        
        // read npach parameter
        for(int j=0; j<7; j++ ) {
            tt_read(ifs, param[j], sizeof(float), 3, inv);
            //  パラメータは制御点なので頂点と同様にscaleをかける
            for(int k=0; k<3; k++ ) {
                param[j][k] *= scale;
            }
        }
        
        // float -> PL_REAL変換
        Vec3<PL_REAL> vertex[3];
        NpatchParam   npatch;
        
        vertex[0].x = vtx[0][0];  vertex[0].y = vtx[0][1];  vertex[0].z = vtx[0][2];
        vertex[1].x = vtx[1][0];  vertex[1].y = vtx[1][1];  vertex[1].z = vtx[1][2];
        vertex[2].x = vtx[2][0];  vertex[2].y = vtx[2][1];  vertex[2].z = vtx[2][2];
        
        npatch.cp_side1_1.x = param[0][0];  npatch.cp_side1_1.y = param[0][1];  npatch.cp_side1_1.z = param[0][2];
        npatch.cp_side1_2.x = param[1][0];  npatch.cp_side1_2.y = param[1][1];  npatch.cp_side1_2.z = param[1][2];
        npatch.cp_side2_1.x = param[2][0];  npatch.cp_side2_1.y = param[2][1];  npatch.cp_side2_1.z = param[2][2];
        npatch.cp_side2_2.x = param[3][0];  npatch.cp_side2_2.y = param[3][1];  npatch.cp_side2_2.z = param[3][2];
        npatch.cp_side3_1.x = param[4][0];  npatch.cp_side3_1.y = param[4][1];  npatch.cp_side3_1.z = param[4][2];
        npatch.cp_side3_2.x = param[5][0];  npatch.cp_side3_2.y = param[5][1];  npatch.cp_side3_2.z = param[5][2];
        npatch.cp_center.x  = param[6][0];  npatch.cp_center.y  = param[6][1];  npatch.cp_center.z  = param[6][2];


        NptTriangle *tri = new NptTriangle(vertex, npatch);
        tri_list.push_back(tri);
        num_tri++;

        if( num_read > 0 ) {  // -1 の場合、チェックなし 
            // 指定数読み込んだか否か判定
            if( num_tri >= num_read ) {
                break;  // 読み込み終了
            }
        }

    }

    if (!ifs.eof() && ifs.fail()) {
        PL_ERROSH << "[ERROR]FileIO_func:stl_b_load_read():Error in loading" <<endl;
        return PLSTAT_STL_IO_ERROR;
    }

    if ( ifs.eof() ) {
        eof = true;
    }

    return PLSTAT_OK;
}


// 三角形ポリゴン情報をバイナリモードで長田パッチファイルに書き出す
//      ファイルへは単精度実数(4byte)で書き込む
POLYLIB_STAT npt_b_save(
    std::vector<NptTriangle*>   *tri_list, 
    const std::string&          fname
)
{
    ofstream ofs(fname.c_str(), ios::out | ios::binary);
    if (ofs.fail()) {
        PL_ERROSH << "[ERROR]FileIO_func:npt_b_save():Can't open " << fname << endl;
        return PLSTAT_NPT_IO_ERROR;
    }

    int inv = tt_check_machine_endian() == TT_LITTLE_ENDIAN ? 0 : 1;

    uint element = tri_list->size();

    tt_write(ofs, &element, sizeof(uint), 1, inv);

    for (uint i = 0; i < element; i++) {
        Vec3<PL_REAL>* vertex = (*tri_list)[i]->get_vertexes();
        const NpatchParam*   npatch = (*tri_list)[i]->get_npatch_param();

        PL_REAL vtx_tmp[9];
        PL_REAL param_tmp[21];

        VEC3_3_TO_REAL9(vertex,vtx_tmp);
        VEC3_TO_REAL( npatch->cp_side1_1, (param_tmp+ 0) );
        VEC3_TO_REAL( npatch->cp_side1_2, (param_tmp+ 3) );
        VEC3_TO_REAL( npatch->cp_side2_1, (param_tmp+ 6) );
        VEC3_TO_REAL( npatch->cp_side2_2, (param_tmp+ 9) );
        VEC3_TO_REAL( npatch->cp_side3_1, (param_tmp+12) );
        VEC3_TO_REAL( npatch->cp_side3_2, (param_tmp+15) );
        VEC3_TO_REAL( npatch->cp_center,  (param_tmp+18) );

        // PL_REAL -> float 変換
        float vtx[9];
        float param[21];
        for(int j=0; j<9; j++) {
            vtx[j] = vtx_tmp[j];
        }
        for(int j=0; j<21; j++) {
            param[j] = param_tmp[j];
        }

        // write vertex
        tt_write(ofs, vtx,   sizeof(float),  9, inv);
        // write npach parameter
        tt_write(ofs, param, sizeof(float), 21, inv);
    }

    if (!ofs.eof() && ofs.fail()) {
        PL_ERROSH << "[ERROR]FileIO_func:npt_b_load():Error in saving: " << fname << endl;
        return PLSTAT_STL_IO_ERROR;
    }

    return PLSTAT_OK;
}

// 長田パッチファイルを読み込みバイナリかアスキーかを判定する
//     //２番目の行に"file name :"文字列があるかどうかで判断する
//     千さんの時からフォーマットを変更する
//          1行目のファセット数は残す
//          fine name出力行は必要なし
//          ファセットIDは出力しない
//     2行目に"facet"文字列があるかどうかで判断する
bool is_npt_a(
    const std::string&     path
)
{
    const char  dcs[] = "\n\r";
    char        buff_tmp[STL_BUFF_LEN];
    char        buff[STL_BUFF_LEN];
    size_t      i = 0;
    char        c;

    ifstream ifs(path.c_str());
    if (!ifs)                  return false;

    // 先頭行を読み飛ばす
    ifs.getline(buff_tmp,STL_BUFF_LEN);
    //PL_DBGOSH << "is_npt_a(): buff_tmp="<<buff_tmp <<endl;
    if (ifs.eof())              return false;
    
    // 2番目行のファイル内容の一部を読み込み
    memset(buff, 0, STL_BUFF_LEN);
    //while (i < STL_BUFF_LEN - 1 && ifs.get(c) && !strchr(dcs, c)) {
    while (i < 8 - 1 && ifs.get(c) && !strchr(dcs, c)) {
        //PL_DBGOSH << "    c="<<c <<endl;
        buff[i++] = c;
    }
    while (ifs.get(c) && strchr(dcs, c))
        ;
    if (ifs.good()) ifs.putback(c);

    if (ifs.eof())              return false;

    // ASCII判定
    if (strstr(buff, "facet"))  return true;
    else                        return false;
}

//**************************************************************
//   共通ＩＯ用
//**************************************************************

//////////////////////////////////////////////////////////////////////////////
char *get_fname_fr_path( const std::string& path )
{
    static char fname[256];

    int pos = path.find_last_of("."); // 拡張子の手前までの位置
    memset(fname, 0, sizeof(fname));
    path.copy(fname, pos, 0);
    return fname;
}

//////////////////////////////////////////////////////////////////////////////
char *get_ext_fr_path( const std::string& path )
{
    static char ext[16];

    int pos = path.find_last_of("."); // 拡張子の手前までの位置
    memset(ext, 0, sizeof(ext));
    path.copy(ext, path.length() - pos, pos + 1);
    return ext;
}

//=======================================================================
// プライベート（static）関数
//=======================================================================
//////////////////////////////////////////////////////////////////////////////
static void tt_invert_byte_order(void* _mem, int size, int n)
{
    char* mem = (char*)_mem;
    char c;
    int i;

    if (size == 1)      return;

    while (n) {
        for (i = 0; i < size/2; i++) {
            c = mem[i];
            mem[i] = mem[size-1-i];
            mem[size-1-i] = c;
        }
        mem += size;
        n--;
    }
}

//////////////////////////////////////////////////////////////////////////////
static int tt_check_machine_endian()
{
    int v = 1;
    char* p = (char*)&v;

    if (p[0])                   return TT_LITTLE_ENDIAN;
    else if (p[sizeof(int)-1])  return TT_BIG_ENDIAN;
    else                        return TT_OTHER_ENDIAN;
}

//////////////////////////////////////////////////////////////////////////////
static void tt_read(istream& is, void* _data, int size, int n, int inv)
{
    char* data = (char*)_data;

    is.read(data, size * n);

    if (inv) {
        tt_invert_byte_order(data, size, n);
    }
}

//////////////////////////////////////////////////////////////////////////////
static void tt_write(ostream& os, const void* _data, int size, int n, int inv)
{
    const char* data = (const char*)_data;
    char* tmp = 0;

    if (inv) {
        int sz = size * n;
        tmp = new char[sz];
        for (int i=0; i<sz; i++) tmp[i] = data[i];
        tt_invert_byte_order(tmp, size, n);
        data = tmp;
    }

    os.write(data, size * n);

    if (inv) {
        delete [] tmp;
    }
}

} //namespace PolylibNS
