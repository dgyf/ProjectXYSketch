#include "CUSketch.h"
#include "SC_CUSketch.h"
#include "CMSketch.h"
#include "Csketch.h"
#include "ASketch.h"
#include <bits/stdc++.h>
#include <ctime>
#include <set>
#include <iostream>
#include <bitset>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <math.h>
#include <string>
#include "Pbsketch_distr_space.h"
#include "XY_CMSketch.h"
#include "FilterCUXY.h"
#include "SpaceSaving_XY.h"
#include "SpaceSaving.h"
#include "SC_SpaceSaving.h"
#include "CMHeap.h"
#include "CUHeap.h"
#include "XYHeap.h"
#include "XY_Heap_ID.h"
#include "Heavy_guardian_heap.h"
#include "Heavy_guardian_top.h"

using namespace std;
const long maxn = 100000000;
//  setting


const int me_cm = 160 * 1024;   // memory size of sketches  (B)
int range = 16;  // range of data items: 2^{range}
const char *filename = "/home/yqliu/data/kosarak.dat";   // Top_1024(1200)
const int ss_k = me_cm / 100;
const int memory_cf = me_cm * 0.2;  // 0.05;0.1;...;0.3    (test number: 6)
const int cf_m = (me_cm - memory_cf) / 100;
//const int thre_cf = 450;
const int thre_cf = 750;
const int bucket_num_cf = 16;
const int xy_memory =  me_cm*0.15;
const int xy_k = (me_cm - xy_memory) / 100;
const int hg_memory=me_cm*0.2;
const int hg_heap=(me_cm-hg_memory)/100;



/*
const int me_cm = 120 * 1024;   // memory size of sketches  (B)
int range = 23;  // range of data items: 2^{range}
const char *filename = "/home/yqliu/data/webdocs.dat"; // Top_1024(58760)
//const char *filename = "/home/yqliu/Desktop/data/webdocs.dat"; // Top_1024(58760)
const int ss_k = me_cm / 100;
const int memory_cf = me_cm * 0.25;
const int cf_m = (me_cm - memory_cf) / 100;
//const int cf_m = 1500;
//const int memory_cf = me_cm - cf_m * 100;
const int thre_cf = 20000;
const int bucket_num_cf = 16;
const int xy_memory = me_cm * 0.15;
const int xy_k = (me_cm - xy_memory) / 100;
const int hg_memory = me_cm * 0.2;
const int hg_heap = (me_cm - hg_memory) / 100;
*/

int a_filter = 1024;     // number of items in A-sketch's filter
int capative = 1024;
int cm_memory = me_cm - capative * 100;
//int cm_memory = 2 * 1024 * 1024;
const int total_memory_in_bytes = me_cm;  //memory size of CF
int me_kb = me_cm / 1024;     // exchange memory to KB
const int filter_memory_percent = 40;  // parameter of CF percent of filter's memory
const int filter_memory_percent_70 = 70;
const int filter_memory_percent_90 = 90;
const int bucket_num = 1;


struct Node {
    int cnt = 0;//count nubmers of every item
    int index;//record index of item in hashtable
};
Node hashtable[maxn];
Node *distribution = new Node[maxn];  // collect data items to figure out the entropy


int cmp(Node a, Node b) {
    return a.cnt > b.cnt;
}


//  Data collecting

int *analysis_data(string filename) { // to calculate the true items frequency

    int *statistics = new int[2];  // statistic n and N
    ifstream fin;
    fin.open(filename, ios::in | ios::out);

    if (!fin)//detect exception
    {
        cout << "Open file failed..." << endl;
        exit(1);
    }

    //var temp: read from file item by item
    int temp, MAX = -1, sum = 0;
    int difference = 0;

    while (!fin.eof() && fin >> temp) {
        if (temp > MAX) {
            MAX = temp;
        }
        if (temp > maxn)//detect exception
        {
            cout << "big item..." << endl;
            exit(1);
        }
        hashtable[temp].cnt++;
        hashtable[temp].index = temp;
        sum++;
    }

    cout << "MAX item number is: " << MAX << endl;
    cout << "The total item numbers is: " << sum << endl;
    sort(hashtable, hashtable + MAX + 1, cmp);

    for (int i = 0; i <= MAX; i++) {
        if (!hashtable[i].cnt)
            break;
        difference += 1;
    }
    fin.close();
    cout << "The number of different item is :" << difference << endl;
/*
    int num_out = 1500;
    for (int i = 0; i < num_out; i++) {
        cout << "NO: " << i << " : " << hashtable[i].cnt << endl;
    }
*/

    statistics[0] = difference;      //n
    statistics[1] = sum;             //N
    return statistics;
}


// Updating

CMSketch *cm_sketch_update(string filename, int me, int dp) {
    ifstream fin;
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir

    if (!fin)//detect exception
    {
        cout << "Open file failed..." << endl;
        exit(1);
    }

    //var temp: read from file item by item
    int temp;
    auto cm = new CMSketch(me, dp);
    while (!fin.eof() && fin >> temp) {
        cm->insert(temp, 1);
    }
    fin.close();
    return cm;
}


Csketch *c_sketch_update(string filename, int memo, int d) {
    ifstream fin;
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir

    if (!fin)//detect exception
    {
        cout << "Open file failed..." << endl;
        exit(1);
    }

    //var temp: read from file item by item
    uint32_t temp;

    auto c_sketch = new Csketch(memo, d);

    while (!fin.eof() && fin >> temp) {
        c_sketch->insert(temp, 1);
    }
    fin.close();
    return c_sketch;
}


CUSketch *cu_sketch_update(string filename, int memo, int d) {
    ifstream fin;
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir

    if (!fin)//detect exception
    {
        cout << "Open file failed..." << endl;
        exit(1);
    }

    //var temp: read from file item by item
    uint32_t temp;

    auto cu_sketch = new CUSketch(memo, d);

    while (!fin.eof() && fin >> temp) {
        cu_sketch->insert(temp, 1);
    }
    fin.close();
    return cu_sketch;
}


ASketch *a_sketch_update(string filename, int memo, int d, int filter_size) {
    ifstream fin;
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir

    if (!fin)//detect exception
    {
        cout << "Open file failed..." << endl;
        exit(1);
    }

    //var temp: read from file item by item
    uint32_t temp;

    auto a_sketch = new ASketch(memo, d, filter_size);

    while (!fin.eof() && fin >> temp) {
        a_sketch->insert(temp, 1);
    }
    fin.close();

    return a_sketch;
}

CMHeap *cm_heap_update(string filename, int capa, int memo, int d) {
    ifstream fin;
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir

    if (!fin)//detect exception
    {
        cout << "Open file failed..." << endl;
        exit(1);
    }

    //var temp: read from file item by item
    uint32_t temp;

    auto cm_heap = new CMHeap(capa, memo, d);

    while (!fin.eof() && fin >> temp) {
        cm_heap->insert(temp);
    }
    fin.close();

    return cm_heap;
}

XYHeap *xy_heap_update(string filename, int capa, int memo, int range) {
    ifstream fin;
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir

    if (!fin)//detect exception
    {
        cout << "Open file failed..." << endl;
        exit(1);
    }

    //var temp: read from file item by item
    uint32_t temp;

    auto xy_heap = new XYHeap(capa, memo, range);

    while (!fin.eof() && fin >> temp) {
        xy_heap->insert(temp, 1);
    }
    fin.close();

    return xy_heap;
}

CUHeap *cu_heap_update(string filename, int capa, int memo, int d) {
    ifstream fin;
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir

    if (!fin)//detect exception
    {
        cout << "Open file failed..." << endl;
        exit(1);
    }

    //var temp: read from file item by item
    uint32_t temp;

    auto cu_heap = new CUHeap(capa, memo, d);

    while (!fin.eof() && fin >> temp) {
        cu_heap->insert(temp);
    }
    fin.close();

    return cu_heap;
}

Heavy_guardian_heap *heavy_g_update(string filename, int memo_heap, int memo_hg) {
    ifstream fin;
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir

    if (!fin)//detect exception
    {
        cout << "Open file failed..." << endl;
        exit(1);
    }

    //var temp: read from file item by item
    uint32_t temp;

    auto heavy_g = new Heavy_guardian_heap(memo_heap, memo_hg);

    while (!fin.eof() && fin >> temp) {
        heavy_g->insert(temp, 1);
    }
    fin.close();

    return heavy_g;
}

Heavy_guardian_tok *heavy_without_update(string filename, int memo) {
    ifstream fin;
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir

    if (!fin)//detect exception
    {
        cout << "Open file failed..." << endl;
        exit(1);
    }

    //var temp: read from file item by item
    uint32_t temp;

    auto hg_without = new Heavy_guardian_tok(memo);

    while (!fin.eof() && fin >> temp) {
        hg_without->insert(temp, 1);
    }
    fin.close();

    return hg_without;
}

SpaceSaving *ss_update(string filename, int k_size) {
    ifstream fin;
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir

    if (!fin)//detect exception
    {
        cout << "Open file failed..." << endl;
        exit(1);
    }

    //var temp: read from file item by item
    uint32_t temp;

    auto ss = new SpaceSaving(k_size);

    while (!fin.eof() && fin >> temp) {
        ss->insert(temp, 1);
    }
    fin.close();
    return ss;
}


SC_SpaceSaving<memory_cf, thre_cf, bucket_num_cf> *cf_ss_update(string filename) {
    ifstream fin;
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir

    if (!fin)//detect exception
    {
        cout << "Open file failed..." << endl;
        exit(1);
    }

    //var temp: read from file item by item
    uint32_t temp;

    auto cf_ss = new SC_SpaceSaving<memory_cf, thre_cf, bucket_num_cf>(cf_m);

    while (!fin.eof() && fin >> temp) {
        cf_ss->insert(temp);
    }
    cf_ss->refre();
    fin.close();

    return cf_ss;
}


Pbsketch_distr_space *pb_distr_space_update(string filename, int me, int ran, int ord[], int sum) {
    ifstream fin;
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir


    if (!fin)//detect exception
    {
        cout << "Open file failed..." << endl;
        exit(1);
    }

    int temp;

    auto Pb_space = new Pbsketch_distr_space(me, ran, ord);
    int count = 0;
    while (!fin.eof() && fin >> temp) {
        if (count >= sum) {
            Pb_space->insert(temp, 1);
        }
        count += 1;
    }
    fin.close();

    return Pb_space;
}


XY_CMSketch *
XY_CM_update(string filename, int me, int ran, int ord[], int sum, double sp_input, double ratio_thre) {
    ifstream fin;
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir


    if (!fin)//detect exception
    {
        cout << "Open file failed..." << endl;
        exit(1);
    }

    int temp;

    auto XY_CM_space = new XY_CMSketch(me, ran, ord, sp_input, ratio_thre);
    int count = 0;
    while (!fin.eof() && fin >> temp) {
        if (count >= sum) {
            XY_CM_space->insert(temp, 1);
        }
        count += 1;
    }
    fin.close();

    return XY_CM_space;
}


FilterCUXY *
FilterCUXY_update(string filename, int me, int ran, int ord[], int sum, double sp_input, double ratio_thre) {
    ifstream fin;
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir

    if (!fin)//detect exception
    {
        cout << "Open file failed..." << endl;
        exit(1);
    }

    int temp;

    auto filter_space = new FilterCUXY(me, ran, ord, sp_input, ratio_thre);
    while (!fin.eof() && fin >> temp) {

        filter_space->insert1(temp, 1, sum);

    }
    fin.close();

    return filter_space;
}

SpaceSaving_XY *
ss_xy_update(string filename, int k_size, int me, int ran) {
    ifstream fin;
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir

    if (!fin)//detect exception
    {
        cout << "Open file failed..." << endl;
        exit(1);
    }

    int temp;

    auto ss_xy_count = new SpaceSaving_XY(k_size, me, ran);
    while (!fin.eof() && fin >> temp) {

        ss_xy_count->insert(temp, 1);

    }
    fin.close();

    return ss_xy_count;
}


template<uint64_t total_memory_in_bytes, int filter_memory_percent, int bucket_num>
CUSketchWithSC<total_memory_in_bytes, filter_memory_percent, bucket_num> *cold_sketch_update(string filename) {
    ifstream fin;
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir

    if (!fin)//detect exception
    {
        cout << "Open file failed..." << endl;
        exit(1);
    }

    //var temp: read from file item by item
    uint32_t temp;

    auto cold_sketch = new CUSketchWithSC<total_memory_in_bytes, filter_memory_percent, bucket_num>();

    while (!fin.eof() && fin >> temp) {
        cold_sketch->insert(temp);
    }
    cold_sketch->synchronize();
    fin.close();
    return cold_sketch;
}


// Querying

double *cold_sketch_query(struct Node raw[],
                          CUSketchWithSC<total_memory_in_bytes, filter_memory_percent, bucket_num> cold_sketch,
                          int number) {
    double *aae = new double[2];
    aae[0] = 0;
    aae[1] = 0;
    for (int i = 0; i < number; i++) {
        aae[0] += (double) abs(raw[i].cnt - cold_sketch.query(raw[i].index)) / number;
        aae[1] += (double) abs(raw[i].cnt - cold_sketch.query(raw[i].index)) / raw[i].cnt / number;
    }

    printf("AAE of the Cold Filter 40 : %lf\n", aae[0]);
    printf("ARE of the Cold Filter 40 : %lf\n", aae[1]);
    printf("***********************************\n");
    return aae;
}


double *cold_sketch_query_70(struct Node raw[],
                             CUSketchWithSC<total_memory_in_bytes, filter_memory_percent_70, bucket_num> cold_sketch,
                             int number) {
    double *aae = new double[2];
    aae[0] = 0;
    aae[1] = 0;
    for (int i = 0; i < number; i++) {
        aae[0] += (double) abs(raw[i].cnt - cold_sketch.query(raw[i].index)) / number;
        aae[1] += (double) abs(raw[i].cnt - cold_sketch.query(raw[i].index)) / raw[i].cnt / number;

    }
    printf("AAE of the Cold Filter 70 : %lf\n", aae[0]);
    printf("ARE of the Cold Filter 70 : %lf\n", aae[1]);
    printf("***********************************\n");
    return aae;
}

double *cold_sketch_query_90(struct Node raw[],
                             CUSketchWithSC<total_memory_in_bytes, filter_memory_percent_90, bucket_num> cold_sketch,
                             int number) {
    double *aae = new double[2];
    aae[0] = 0;
    aae[1] = 0;
    for (int i = 0; i < number; i++) {
        aae[0] += (double) abs(raw[i].cnt - cold_sketch.query(raw[i].index)) / number;
        aae[1] += (double) abs(raw[i].cnt - cold_sketch.query(raw[i].index)) / raw[i].cnt / number;
    }

    printf("AAE of the Cold Filter 90 : %lf\n", aae[0]);
    printf("ARE of the Cold Filter 90 : %lf\n", aae[1]);
    printf("***********************************\n");
    return aae;
}


double *cm_sketch_query(struct Node raw[], CMSketch cm_sketch, int number) {
    double *aae = new double[2];
    aae[0] = 0;
    aae[1] = 0;
    for (int i = 0; i < number; i++) {
        aae[0] += (double) abs(raw[i].cnt - cm_sketch.query(raw[i].index)) / number;
        aae[1] += (double) abs(raw[i].cnt - cm_sketch.query(raw[i].index)) / raw[i].cnt / number;
    }

    printf("AAE of the CMsketch : %lf\n", aae[0]);
    printf("ARE of the CMsketch : %lf\n", aae[1]);
    printf("***********************************\n");

    return aae;
}

double *c_sketch_query(struct Node raw[], Csketch c_sketch, int number) {
    double *aae = new double[2];
    aae[0] = 0;
    aae[1] = 0;
    for (int i = 0; i < number; i++) {
        aae[0] += (double) abs(raw[i].cnt - c_sketch.query(raw[i].index)) / number;
        aae[1] += (double) abs(raw[i].cnt - c_sketch.query(raw[i].index)) / raw[i].cnt / number;
    }

    printf("AAE of the Csketch : %lf\n", aae[0]);
    printf("ARE of the Csketch : %lf\n", aae[1]);
    printf("***********************************\n");
    return aae;
}

double *cu_sketch_query(struct Node raw[], CUSketch cu_sketch, int number) {
    double *aae = new double[2];
    aae[0] = 0;
    aae[1] = 0;
    for (int i = 0; i < number; i++) {
        aae[0] += (double) abs(raw[i].cnt - cu_sketch.query(raw[i].index)) / number;
        aae[1] += (double) abs(raw[i].cnt - cu_sketch.query(raw[i].index)) / raw[i].cnt / number;
    }

    printf("AAE of the CUsketch : %lf\n", aae[0]);
    printf("ARE of the CUsketch : %lf\n", aae[1]);
    printf("***********************************\n");

    return aae;
}

double *a_sketch_query(struct Node raw[], ASketch a_sketch, int number) {
    double *aae = new double[2];
    aae[0] = 0;
    aae[1] = 0;
    for (int i = 0; i < number; i++) {
        aae[0] += (double) abs(raw[i].cnt - a_sketch.query(raw[i].index)) / number;
        aae[1] += (double) abs(raw[i].cnt - a_sketch.query(raw[i].index)) / raw[i].cnt / number;
    }

    printf("AAE of the ASketch : %lf\n", aae[0]);
    printf("ARE of the ASketch : %lf\n", aae[1]);
    printf("***********************************\n");

    return aae;
}


double *pb_dis_space_query(struct Node raw[], Pbsketch_distr_space pb, int number) {
    double *aae = new double[2];
    aae[0] = 0;
    aae[1] = 0;

    for (int i = 0; i < number; i++) {
        aae[0] += (double) abs(raw[i].cnt - pb.query(raw[i].index)) / number;
        aae[1] += (double) abs(raw[i].cnt - pb.query(raw[i].index)) / raw[i].cnt / number;

//        cout<<raw[i].cnt<<" "<<pb.query(raw[i].index)<<endl;
    }

    printf("AAE of the Pb_dis_space : %lf\n", aae[0]);
    printf("ARE of the Pb_dis_space : %lf\n", aae[1]);
    printf("***********************************\n");

    return aae;
}

double *XY_CM_query(struct Node raw[], XY_CMSketch xy_cm, int number) {
    double *aae = new double[6];
    for (int i = 0; i < 6; i++) {
        aae[i] = 0;
    }

    for (int i = 0; i < number; i++) {
        aae[0] += (double) abs(raw[i].cnt - xy_cm.query_ratio_cm(raw[i].index)) / number;
        aae[1] += (double) abs(raw[i].cnt - xy_cm.query_ratio_cm(raw[i].index)) / raw[i].cnt / number;
        /*
         aae[2] += (double) abs(raw[i].cnt - xy_cm.query_value_cm(raw[i].index)) / number;
         aae[3] += (double) abs(raw[i].cnt - xy_cm.query_value_cm(raw[i].index)) / raw[i].cnt / number;
         aae[4] += (double) abs(raw[i].cnt - xy_cm.query_ratio_xy(raw[i].index)) / number;
         aae[5] += (double) abs(raw[i].cnt - xy_cm.query_ratio_xy(raw[i].index)) / raw[i].cnt / number;
 */
//        cout<<raw[i].cnt<<" "<<pb.query(raw[i].index)<<endl;
    }

    printf("AAE of the XY_CM (ratio_cm) : %lf\n", aae[0]);
    printf("ARE of the XY_CM (ratio_cm) : %lf\n", aae[1]);
    printf("***********************************\n");
    /*
    printf("AAE of the XY_CM (value_cm) : %lf\n", aae[2]);
    printf("ARE of the XY_CM (value_cm) : %lf\n", aae[3]);
    printf("***********************************\n");
    printf("AAE of the XY_CM (ratio_xy) : %lf\n", aae[4]);
    printf("ARE of the XY_CM (ratio_xy) : %lf\n", aae[5]);
    printf("***********************************\n");
     */
    printf("over-estimated : %d\n", xy_cm.get_count());

    return aae;
}

double *Filter_query(struct Node raw[], FilterCUXY filter_xycm, int number) {
    double *aae = new double[2];
    for (int i = 0; i < 2; i++) {
        aae[i] = 0;
    }

    for (int i = 0; i < number; i++) {
        aae[0] += (double) abs(raw[i].cnt - filter_xycm.query_ratio_cm(raw[i].index)) / number;
        aae[1] += (double) abs(raw[i].cnt - filter_xycm.query_ratio_cm(raw[i].index)) / raw[i].cnt / number;

//        cout<<raw[i].cnt<<" "<<pb.query(raw[i].index)<<endl;
    }

    printf("AAE of the Filter_XYCU (ratio_cm) : %lf\n", aae[0]);
    printf("ARE of the Filter_XYCU (ratio_cm) : %lf\n", aae[1]);
    printf("***********************************\n");

    return aae;
}

double *ss_xy_query(struct Node raw[], SpaceSaving_XY ss_xy_counter, int number) {
    double *aae = new double[2];
    for (int i = 0; i < 2; i++) {
        aae[i] = 0;
    }

    for (int i = 0; i < number; i++) {
        aae[0] += (double) abs(raw[i].cnt - ss_xy_counter.query_fre(raw[i].index)) / number;
        aae[1] += (double) abs(raw[i].cnt - ss_xy_counter.query_fre(raw[i].index)) / raw[i].cnt / number;

//        cout<<raw[i].cnt<<" "<<pb.query(raw[i].index)<<endl;
    }

    printf("AAE of the SpaceSaving_XY : %lf\n", aae[0]);
    printf("ARE of the SpaceSaving_XY : %lf\n", aae[1]);
    printf("***********************************\n");

    return aae;
}

/*    // Origin AAE
pair<double, double> top_k_accuracy(int k, vector<pair<uint32_t, uint32_t>> &result) {

    double num = 0;
    double aae = 0;
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            if (result[i].first == hashtable[j].index) {
                num = num + 1;
                aae = aae + abs((int) result[i].second - hashtable[j].cnt);
                break;
            }
        }
    }

    return make_pair(num / k, aae / num);
}
 */



// All k AAE
pair<double, double> top_k_accuracy(int k, vector<pair<uint32_t, uint32_t>> &result) {

    double num = 0;
    double aae = 0;
    int open_close[k];
    memset(open_close, 0, sizeof(open_close));
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            if (result[i].first == hashtable[j].index) {
                num = num + 1;
                aae = aae + abs((int) result[i].second - hashtable[j].cnt);
                open_close[j] = 1;
                break;
            }
        }
    }
    for (int i = 0; i < k; i++) {
        if (open_close[i] == 0) {
            aae = aae + hashtable[i].cnt;
//            cout<<" lose: "<<hashtable[i].cnt<<endl;
        }
    }

    return make_pair(num / k, aae / k);
}


// Plus +ARE
double *top_k_accuracy_all(int k, vector<pair<uint32_t, uint32_t>> &result) {

    double num = 0;
    double aae = 0;
    double are = 0;
    double *error_three = new double[3];
    for (int i = 0; i < 3; i++) {
        error_three[i] = 0;
    }
    int open_close[k];

    memset(open_close, 0, sizeof(open_close));
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            if (result[i].first == hashtable[j].index) {
                num = num + 1;
                aae = aae + abs((int) result[i].second - hashtable[j].cnt);
                are = are + (double) abs((int) result[i].second - hashtable[j].cnt) / hashtable[j].cnt;
                open_close[j] = 1;
                break;
            }
        }
    }
    for (int i = 0; i < k; i++) {
        if (open_close[i] == 0) {
            aae = aae + hashtable[i].cnt;
            are = are + 1;
//            cout<<" lose: "<<hashtable[i].cnt<<endl;
        }
    }
    error_three[0] = num / k;
    error_three[1] = aae / k / 1000;  // AAE (10^3)
    error_three[2] = are / k;


    return error_three;
}



int main() {      // plus_ARE      //plus_CUsketch
    clock_t start_time = clock();

    int *as = analysis_data(filename);

    cout << "All skecthes's memory is : " << me_cm << endl;
    //cout << "The filter size of A-sketch is : " << a_filter << endl;
    cout << "***************************************" << endl;


    const int num_k = 6;
    int find_k[num_k] = {32, 64, 128, 256, 512, 1024};

//    const int num_k = 1;
//    int find_k[num_k] = {a_filter};

    SpaceSaving *ss_count = ss_update(filename, ss_k);
    SC_SpaceSaving<memory_cf, thre_cf, bucket_num_cf> *ss_cf_count = cf_ss_update(filename);
//    SpaceSaving_XY *ss_xy_count = ss_xy_update(filename, xy_k, xy_memory,
//                                               range);
    XYHeap *xy_heap_count = xy_heap_update(filename, xy_k, xy_memory, range);

//    ASketch *asketch_count = a_sketch_update(filename, me_cm, 3, a_filter);

    CMHeap *cm_heap_count = cm_heap_update(filename, capative, cm_memory, 3);

    CUHeap *cu_heap_count = cu_heap_update(filename, capative, cm_memory, 3);

    Heavy_guardian_heap *heavy_g = heavy_g_update(filename, hg_heap, hg_memory);


    for (int i = 0; i < num_k; i++) {

        cout << " K value is : " << find_k[i] << endl;

        vector<pair<uint32_t, uint32_t>> result(find_k[i]);


//        asketch_count->get_top_k_with_frequency(find_k[i], result);
//        auto outp_a = top_k_accuracy_all(find_k[i], result);
//        printf("\tA-sketch \n");
//        printf("\t  accuracy %lf, AAE %lf, ARE %lf\n", outp_a[0], outp_a[1], outp_a[2]);


        cm_heap_count->get_top_k_with_frequency(find_k[i], result);
        auto outp_cm_heap = top_k_accuracy_all(find_k[i], result);
        printf("\tCM_Heap \n");
        printf("\t  accuracy %lf, AAE %lf, ARE %lf\n", outp_cm_heap[0], outp_cm_heap[1], outp_cm_heap[2]);

        cu_heap_count->get_top_k_with_frequency(find_k[i], result);
        auto outp_cu_heap = top_k_accuracy_all(find_k[i], result);
        printf("\tCU_Heap \n");
        printf("\t  accuracy %lf, AAE %lf, ARE %lf\n", outp_cu_heap[0], outp_cu_heap[1], outp_cu_heap[2]);

        ss_count->get_top_k_with_frequency(find_k[i], result);
        auto outp_ss = top_k_accuracy_all(find_k[i], result);
        printf("\tSpaceSaving \n");
        printf("\t  accuracy %lf, AAE %lf, ARE %lf\n", outp_ss[0], outp_ss[1], outp_ss[2]);


        ss_cf_count->get_top_k_with_frequency(find_k[i], result);
        auto outp_ss_cf = top_k_accuracy_all(find_k[i], result);
        printf("\tSpaceSaving wich CF \n");
        printf("\t  accuracy %lf, AAE %lf, ARE %lf\n", outp_ss_cf[0], outp_ss_cf[1], outp_ss_cf[2]);

        heavy_g->get_top_k_with_frequency(find_k[i], result);
        auto out_heavy_g = top_k_accuracy_all(find_k[i], result);
        printf("\tHeavy_g_Heap \n");
        printf("\t  accuracy %lf, AAE %lf, ARE %lf\n", out_heavy_g[0], out_heavy_g[1], out_heavy_g[2]);


        xy_heap_count->get_top_k_with_frequency(find_k[i], result);
        auto outp_ss_xy = top_k_accuracy_all(find_k[i], result);
        printf("\tPB_Heap \n");
        printf("\t  accuracy %lf, AAE %lf, ARE %lf\n", outp_ss_xy[0], outp_ss_xy[1], outp_ss_xy[2]);


        printf("***********************************\n");

        cout << " The memory is : " << me_cm << endl;
        cout << " The number of bucket of SS is : " << ss_k << endl;
        cout << " The number of bucket of SS is : " << cf_m << " The rest of memory of filter is : " << memory_cf
             << " The Threshold is : " << thre_cf << endl;
        cout << " The memory of PB is : " << xy_memory << endl;


        ofstream fout("finally_kos.dat", ios::app);


        fout << me_cm << " " << me_kb << " " << find_k[i] << " " << a_filter << " "
             << capative << " " << cm_memory << " "
             << ss_k << " " << cf_m << " " << memory_cf << " " << thre_cf << " " << xy_memory << " "
             //             << outp_a[0] << " " << outp_a[1] << " " << outp_a[2] << " "
             << outp_cm_heap[0] << " " << outp_cm_heap[1] << " " << outp_cm_heap[2] << " "
             << outp_ss[0] << " " << outp_ss[1] << " " << outp_ss[2] << " "
             << outp_ss_cf[0] << " " << outp_ss_cf[1] << " " << outp_ss_cf[2] << " "
             << outp_ss_xy[0] << " " << outp_ss_xy[1] << " " << outp_ss_xy[2] << " "
             << out_heavy_g[0] << " " << out_heavy_g[1] << " " << out_heavy_g[2] << " "
             << outp_cu_heap[0] << " " << outp_cu_heap[1] << " " << outp_cu_heap[2] << endl;

        fout.close();
    }


    clock_t end_time = clock();
    cout << "The total time is :" << (double) (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;
    return 0;
}


/*

int main() {      // only test the top-k querying of heavy_guardian
    clock_t start_time = clock();

    int *as = analysis_data(filename);

    cout << "All skecthes's memory is : " << me_cm << endl;
    cout << "***************************************" << endl;

    clock_t time_log[4];
//    const int num_k = 6;
    const int num_k = 1;
//    int find_k[num_k] = {32, 64, 128, 256, 512, 1024};
    int find_k[num_k] = {512};
//    const int num_k = 1;
//    int find_k[num_k] = {512};

    time_log[0] = clock();
    XYHeap *xy_heap_count = xy_heap_update(filename, xy_k, xy_memory, range);
    time_log[1] = clock();

    CUHeap *cu_heap_count = cu_heap_update(filename, capative, cm_memory, 3);
    time_log[2] = clock();
    Heavy_guardian_heap *heavy_g = heavy_g_update(filename, hg_heap, hg_memory);
    time_log[3] = clock();

//    Heavy_guardian_tok *heavy_without_auxi = heavy_without_update(filename,me_cm);



    double timeinsert[3];
    double eff_in[3];

    for (int i = 0; i < 3; i++) {
        timeinsert[i] = (double) (time_log[i + 1] - time_log[i]) / CLOCKS_PER_SEC;
        eff_in[i] = (double) as[1] / timeinsert[i] / 1000000;

        cout << " Efficiency of insert NO: " << i << " : " << eff_in[i] << endl;
    }

    ofstream ftime("time_kos_topk.dat", ios::app);

    ftime << me_cm << " " << find_k[0] << " "
          << eff_in[0] << " "
          << eff_in[1] << " "
          << eff_in[2] << endl;


    ftime.close();

    cout << " the column of hg+heap: " << heavy_g->get_row() << endl;

    for (int i = 0; i < num_k; i++) {

        cout << " K value is : " << find_k[i] << endl;

        vector<pair<uint32_t, uint32_t>> result(find_k[i]);


//        asketch_count->get_top_k_with_frequency(find_k[i], result);
//        auto outp_a = top_k_accuracy_all(find_k[i], result);
//        printf("\tA-sketch \n");
//        printf("\t  accuracy %lf, AAE %lf, ARE %lf\n", outp_a[0], outp_a[1], outp_a[2]);



        cu_heap_count->get_top_k_with_frequency(find_k[i], result);
        auto outp_cu_heap = top_k_accuracy_all(find_k[i], result);
        printf("\tCU_Heap \n");
        printf("\t  accuracy %lf, AAE %lf, ARE %lf\n", outp_cu_heap[0], outp_cu_heap[1], outp_cu_heap[2]);


        xy_heap_count->get_top_k_with_frequency(find_k[i], result);
        auto outp_ss_xy = top_k_accuracy_all(find_k[i], result);
        printf("\tXY_Heap \n");
        printf("\t  accuracy %lf, AAE %lf, ARE %lf\n", outp_ss_xy[0], outp_ss_xy[1], outp_ss_xy[2]);


        heavy_g->get_top_k_with_frequency(find_k[i], result);
        auto out_heavy_g = top_k_accuracy_all(find_k[i], result);
        printf("\tHeavy_g_Heap \n");
        printf("\t  accuracy %lf, AAE %lf, ARE %lf\n", out_heavy_g[0], out_heavy_g[1], out_heavy_g[2]);


//        heavy_without_auxi->get_top_k_with_frequency(find_k[i], result);
//        for (auto itr:result) {
//            cout << " test: " << itr.first << " " << itr.second << endl;
//        }

//        auto outp_heavy_without = top_k_accuracy_all(find_k[i], result);
//        printf("\tHeavy_without_Heap \n");
//        printf("\t  accuracy %lf, AAE %lf, ARE %lf\n", outp_heavy_without[0], outp_heavy_without[1],
//               outp_heavy_without[2]);

        printf("***********************************\n");

        cout << " The memory is : " << me_cm << endl;
        cout << " The memory of XY is : " << xy_memory << endl;

        ofstream fout("top_k_heavyguardian_space_web.dat", ios::app);


        fout << me_cm << " " << me_kb << " " << find_k[i] << " " << a_filter << " "
             << out_heavy_g[0] << " " << out_heavy_g[1] << " " << out_heavy_g[2] << " "
             //             << outp_heavy_without[0] << " " << outp_heavy_without[1] << " " << outp_heavy_without[2] << " "
             << outp_ss_xy[0] << " " << outp_ss_xy[1] << " " << outp_ss_xy[2] << " "
             << outp_cu_heap[0] << " " << outp_cu_heap[1] << " " << outp_cu_heap[2] << endl;

        fout.close();
    }


    clock_t end_time = clock();
    cout << "The total time is :" << (double) (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;
    return 0;
}

*/
/*
int main() {        // time test
    clock_t start_time = clock();

    int *as = analysis_data(filename);

    cout << "All skecthes's memory is : " << me_cm << endl;
    cout << "The filter size of A-sketch is : " << a_filter << endl;
    cout << "***************************************" << endl;


    const int num_k = 6;
    int find_k[num_k] = {32, 64, 128, 256, 512, 1024};

//    const int num_k = 1;
//    int find_k[num_k] = {a_filter};

    clock_t time_log[6];
    clock_t time_query[10];

    time_log[0] = clock();
    CUHeap *cu_heap_count = cu_heap_update(filename, capative, cm_memory, 3);


//    ASketch *asketch_count = a_sketch_update(filename, me_cm, 3, a_filter);
    time_log[1] = clock();

    CMHeap *cm_heap_count = cm_heap_update(filename, capative, cm_memory, 3);
    time_log[2] = clock();

    SpaceSaving *ss_count = ss_update(filename, ss_k);
    time_log[3] = clock();
    SC_SpaceSaving<memory_cf, thre_cf, bucket_num_cf> *ss_cf_count = cf_ss_update(filename);
    time_log[4] = clock();
//    SpaceSaving_XY *ss_xy_count = ss_xy_update(filename, xy_k, xy_memory,
//                                               range);
    XYHeap *xy_heap_count = xy_heap_update(filename, xy_k, xy_memory, range);
    time_log[5] = clock();


    for (int i = 0; i < num_k; i++) {

        cout << " K value is : " << find_k[i] << endl;

        vector<pair<uint32_t, uint32_t>> result(find_k[i]);

        time_query[0] = clock();
        cu_heap_count->get_top_k_with_frequency(find_k[i], result);
//        asketch_count->get_top_k_with_frequency(find_k[i], result);
        time_query[1] = clock();
        auto outp_a = top_k_accuracy(find_k[i], result);
//        printf("\tA-sketch \n");
        printf("\tCU-sketch \n");
        printf("\t  accuracy %lf, AAE %lf\n", outp_a.first, outp_a.second);


        time_query[2] = clock();
        cm_heap_count->get_top_k_with_frequency(find_k[i], result);
        time_query[3] = clock();
        auto outp_cm_heap = top_k_accuracy(find_k[i], result);
        printf("\tCM_Heap \n");
        printf("\t  accuracy %lf, AAE %lf\n", outp_cm_heap.first, outp_cm_heap.second);

        time_query[4] = clock();
        ss_count->get_top_k_with_frequency(find_k[i], result);
        time_query[5] = clock();
        auto outp_ss = top_k_accuracy(find_k[i], result);
        printf("\tSpaceSaving \n");
        printf("\t  accuracy %lf, AAE %lf\n", outp_ss.first, outp_ss.second);

        time_query[6] = clock();
        ss_cf_count->get_top_k_with_frequency(find_k[i], result);
        time_query[7] = clock();
        auto outp_ss_cf = top_k_accuracy(find_k[i], result);
        printf("\tSpaceSaving wich CF \n");
        printf("\t  accuracy %lf, AAE %lf\n", outp_ss_cf.first, outp_ss_cf.second);

        time_query[8] = clock();
//        ss_xy_count->get_top_k_with_frequency(find_k[i], result);
        xy_heap_count->get_top_k_with_frequency(find_k[i], result);
        time_query[9] = clock();
        auto outp_ss_xy = top_k_accuracy(find_k[i], result);
        printf("\tXY_SS \n");
        printf("\t  accuracy %lf, AAE %lf\n", outp_ss_xy.first, outp_ss_xy.second);


        printf("***********************************\n");

        cout << " The memory is : " << me_cm << endl;
        cout << " The number of bucket of SS is : " << ss_k << endl;
        cout << " The number of bucket of SS is : " << cf_m << " The rest of memory of filter is : " << memory_cf
             << " The Threshold is : " << thre_cf << endl;
        cout << " The memory of XY is : " << xy_memory << endl;


        ofstream fout("time_cuplus_web.dat", ios::app);


        fout << me_cm << " " << me_kb << " " << find_k[i] << " " << a_filter << " "
             << capative << " " << cm_memory << " "
             << ss_k << " " << cf_m << " " << memory_cf << " " << thre_cf << " " << xy_memory << " "
             << outp_a.first << " " << outp_a.second << " "
             << outp_cm_heap.first << " " << outp_cm_heap.second << " "
             << outp_ss.first << " " << outp_ss.second << " "
             << outp_ss_cf.first << " " << outp_ss_cf.second
             << " " << outp_ss_xy.first << " " << outp_ss_xy.second << endl;

        fout.close();

        double timeinsert[5];
        double timequery[5];
        double eff_in[5];
        double eff_query[5];

        for (int i = 0; i < 5; i++) {
            timeinsert[i] = (double) (time_log[i + 1] - time_log[i]) / CLOCKS_PER_SEC;
            timequery[i] = (double) (time_query[2 * i + 1] - time_query[2 * i]) / CLOCKS_PER_SEC;
            eff_in[i] = (double) as[1] / timeinsert[i] / 1000000;
            eff_query[i] = (double) find_k[i] / timequery[i] / 1000000;

            cout << " Efficiency of insert : " << eff_in[i] << " Efficiency of query : " << eff_query[i] << endl;
        }


        ofstream ftime("time_web_pluscu.dat", ios::app);


        ftime << me_cm << " " << me_kb << " " << find_k[i] << " " << a_filter << " "
              << capative << " " << cm_memory << " "
              << ss_k << " " << cf_m << " " << memory_cf << " " << thre_cf << " " << xy_memory << " "
              << eff_in[0] << " " << eff_query[0] << " "
              << eff_in[1] << " " << eff_query[1] << " "
              << eff_in[2] << " " << eff_query[2] << " "
              << eff_in[3] << " " << eff_query[3] << " "
              << eff_in[4] << " " << eff_query[4] << endl;


        ftime.close();

    }


    clock_t end_time = clock();
    cout << "The total time is :" << (double) (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;
    return 0;
}

*/