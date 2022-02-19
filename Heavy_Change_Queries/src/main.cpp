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
#include "InsertableIblt.h"
#include "SC_InsertableIBLT.h"
#include "FlowRadar_without_packetcount.h"
#include "XY_FlowRadar.h"
#include "Heavy_guardian_flowradar.h"
#include "HG_XY_heavy.h"
#include "HG_for_XY.h"
#include "Heavy_guardian_list.h"
#include "HG_for_XY_without_auxiliary.h"
#include "HG(heavy_light).h"

using namespace std;
const long maxn = 100000000;
const long max_web = 300000000;
//  setting
//const char *filename = "/home/yqliu/Desktop/data/webdocs.dat";
//const int window_num = 16 * 1024 * 1024 * 2;

const char *filename = "/home/yqliu/data/kosarak.dat";  // (D=210385)
//const char *filename = "/home/yqliu/data/webdocs.dat";  //(D=24580773)
const int range = 16;  // range of data items: 2^{range}
//const int range = 23;  // range of data items: 2^{range}
const int me_cm = 160 * 1024;   // memory size of sketches  (B)
//const double thre = 0.0004;
//const int threshold_cf = 2159; //(kos) 105 126 147 168 189 210

//XY setting
const int memory_xy_radar = 90 * 1024;
const int memory_xy = me_cm - memory_xy_radar;
//HG setting
const int memory_heavy_radar = 40 * 1024;
const int memory_hg = me_cm - memory_heavy_radar;

//CF setting
//const ingt thre = 1000;
const int threshold_cf = 300;   //kosarak (42 for 0.0004) (200 for 120kb)
//const int threshold_cf = 4500;   //webdocs(7000 for 0.0004 (5000  for space) (4500 for thre)

const int memory_cf = me_cm * 0.45;  // 0.35
const int memory_cold_radar = me_cm - memory_cf;
//const int memory_cold_radar = 50000 * 1024;

/*
const int memory_cf_2 = me_cm * 0.40;  // 0.35   0.40  0.45
const int memory_cold_radar_2 = me_cm - memory_cf_2;
//const int memory_cold_radar_2=50000*1024;

const int memory_cf_3 = me_cm * 0.45;
const int memory_cold_radar_3 = me_cm - memory_cf_3;
//const int memory_cold_radar_3=50000*1024;
*/

const int hg_me_heavy = 90 * 1024;
const int hg_me_light = me_cm - hg_me_heavy;

const int total_memory_in_bytes = me_cm;  //memory size of CF  (note! threshold inside!)
int a_filter = 32;     // number of items in A-sketch's filter
int me_kb = me_cm / 1024;     // exchange memory to KB
const int filter_memory_percent = 40;  // parameter of CF percent of filter's memory
const int filter_memory_percent_70 = 70;
const int filter_memory_percent_90 = 90;
const int bucket_num = 1;
const int memory_radar = me_cm;


struct Node {
    int cnt = 0;//count nubmers of every item
    int index;//record index of item in hashtable
};
Node hashtable[maxn];
Node *distribution = new Node[maxn];  // collect data items to figure out the entropy
//uint32_t insert_data_stream[maxn];
uint32_t *insert_data_stream = new uint32_t[max_web];

int cmp(Node a, Node b) {
    return a.cnt > b.cnt;
}

//  Data collecting

int *analysis_data(string filename, int num_top) { // to calculate the true items frequency

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
        insert_data_stream[sum] = temp;
        sum++;
        if (sum == num_top) {
            cout << " The Max item among the first : " << num_top << " is : " << MAX << endl;
        }
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


typedef unordered_map<uint32_t, int> ItemSet;

ItemSet fr_dump_diff(VirtualIBLT *algo1, VirtualIBLT *algo2) // To find union and their counts
{
    ItemSet a, b, c;
    algo1->dump(a);
    algo2->dump(b);

    // try another method to dump diff
    for (auto itr: a) {
        if (b.find(itr.first) != b.end()) {
            c[itr.first] = a[itr.first] - b[itr.first];
        } else {
            int appr_b = algo2->approximate_query(itr.first);
            if (appr_b != -1) {  // this item is not exist!
                c[itr.first] = a[itr.first] - appr_b;
//                cout<<"test_hg: pass a is: "<<a[itr.first]<<" c is : "<<abs(c[itr.first])<<endl;
            }
        }
    }

    for (auto itr: b) {
        if (c.find(itr.first) != c.end()) {
            continue;
        }
        if (a.find(itr.first) != a.end()) {
            c[itr.first] = a[itr.first] - b[itr.first];
            cout << "error" << endl;
        } else {
            int appr_a = algo1->approximate_query(itr.first);
            if (appr_a != -1) {
                c[itr.first] = appr_a - b[itr.first];
//                cout<<"test_hg: pass b is: "<<b[itr.first]<<" c is : "<<abs(c[itr.first])<<endl;
            }
        }
    }

    return c;
}

ItemSet fr_dump_diff_hg(HG_heavy_light *algo1, HG_heavy_light *algo2, int thre) // To find union and their counts
{
    ItemSet a, b, c;
    algo1->query_find(a, thre);
    algo2->query_find(b, thre);

    // try another method to dump diff
    for (auto itr: a) {
        if (b.find(itr.first) != b.end()) {
            c[itr.first] = a[itr.first] - b[itr.first];
        } else {
            int appr_b = algo2->query_point(itr.first);
            c[itr.first] = a[itr.first] - appr_b;

        }
    }

    for (auto itr: b) {
        if (c.find(itr.first) != c.end()) {
            continue;
        }
        if (a.find(itr.first) != a.end()) {
            c[itr.first] = a[itr.first] - b[itr.first];
            cout << "error" << endl;
        } else {
            int appr_a = algo1->query_point(itr.first);
            c[itr.first] = appr_a - b[itr.first];
        }
    }

    return c;
}

ItemSet fr_dump_diff_xy_hg(HG_XY_heavy *algo1, HG_XY_heavy *algo2) // To find union and their counts
{
    ItemSet a, b, c;
    algo1->query_find(a);
    algo2->query_find(b);

    // try another method to dump diff
    for (auto itr: a) {
        if (b.find(itr.first) != b.end()) {
            c[itr.first] = a[itr.first] - b[itr.first];
        } else {
            int appr_b = algo2->query_xy(itr.first);
//            cout << " frequency of items in XY is : " << appr_b << endl;
            c[itr.first] = a[itr.first] - appr_b;
        }

    }

    for (auto itr: b) {
        if (c.find(itr.first) != c.end()) {
            continue;
        }
        if (a.find(itr.first) != a.end()) {
            cout << "error" << endl;
        } else {
            int appr_a = algo1->query_xy(itr.first);
//            cout << " frequency of items in XY is : " << appr_a << endl;
            c[itr.first] = b[itr.first] - appr_a;
        }
    }

    return c;
}

ItemSet fr_dump_diff_xy(HG_for_XY *algo1, HG_for_XY *algo2, Pbsketch_distr_space *xy_a,
                        Pbsketch_distr_space *xy_b, int thre) // To find union and their counts
{
    ItemSet a, b, c;
    algo1->query_find(a, thre);
    algo2->query_find(b, thre);
    // try another method to dump diff
    for (auto itr: a) {
        if (b.find(itr.first) != b.end()) {
            c[itr.first] = a[itr.first] - b[itr.first];
        } else {
            int appr_b = xy_b->query(itr.first);
//            cout << " frequency of items: " << itr.first << " in XY is : " << appr_b << endl;
            c[itr.first] = a[itr.first] - appr_b;
        }

    }

    for (auto itr: b) {
        if (c.find(itr.first) != c.end()) {
            continue;
        }
        if (a.find(itr.first) != a.end()) {
            cout << "error" << endl;
        } else {
            int appr_a = xy_a->query(itr.first);
//            cout << " frequency of items : " << itr.first << " in XY is : " << appr_a << endl;
            c[itr.first] = b[itr.first] - appr_a;
        }
    }

    return c;
}

ItemSet fr_dump_diff_xy_hg(HG_for_XY *algo1, HG_for_XY *algo2, Pbsketch_distr_space *xy_a,
                           Pbsketch_distr_space *xy_b, int thre) // To find union and their counts
{
    ItemSet a, b, c;
    algo1->query_find(a, thre);
    algo2->query_find(b, thre);
    // try another method to dump diff
    for (auto itr: a) {
        if (b.find(itr.first) != b.end()) {
            c[itr.first] = a[itr.first] - b[itr.first];
        } else {
            int op = algo2->query_point(itr.first);
            if (op == 0) {
                int appr_b = xy_b->query(itr.first);
//                cout << " frequency of items: " << itr.first << " in XY is : " << appr_b << endl;
                c[itr.first] = a[itr.first] - appr_b;
            } else {
                c[itr.first] = a[itr.first] - op;
            }
        }

    }

    for (auto itr: b) {
        if (c.find(itr.first) != c.end()) {
            continue;
        }
        if (a.find(itr.first) != a.end()) {
            cout << "error" << endl;
        } else {
            int op = algo1->query_point(itr.first);
            if (op == 0) {
                int appr_a = xy_a->query(itr.first);
//                cout << " frequency of items : " << itr.first << " in XY is : " << appr_a << endl;
                c[itr.first] = b[itr.first] - appr_a;
            } else {
                c[itr.first] = b[itr.first] - op;
            }
        }
    }

    return c;
}

ItemSet fr_dump_diff_xy_hg_without(HG_for_XY_without_auxiliary *algo1, HG_for_XY_without_auxiliary *algo2,
                                   Pbsketch_distr_space *xy_a,
                                   Pbsketch_distr_space *xy_b, int thre) // To find union and their counts
{
    ItemSet a, b, c;
    algo1->query_find(a, thre);
    algo2->query_find(b, thre);
    // try another method to dump diff
    for (auto itr: a) {
        if (b.find(itr.first) != b.end()) {
            c[itr.first] = a[itr.first] - b[itr.first];
        } else {
            int op = algo2->query_point(itr.first);
            if (op == 0) {
                int appr_b = xy_b->query(itr.first);
//                cout << " frequency of items: " << itr.first << " in XY is : " << appr_b << endl;
                c[itr.first] = a[itr.first] - appr_b;
            } else {
                c[itr.first] = a[itr.first] - op;
            }
        }

    }

    for (auto itr: b) {
        if (c.find(itr.first) != c.end()) {
            continue;
        }
        if (a.find(itr.first) != a.end()) {
            cout << "error" << endl;
        } else {
            int op = algo1->query_point(itr.first);
            if (op == 0) {
                int appr_a = xy_a->query(itr.first);
//                cout << " frequency of items : " << itr.first << " in XY is : " << appr_a << endl;
                c[itr.first] = b[itr.first] - appr_a;
            } else {
                c[itr.first] = b[itr.first] - op;
            }
        }
    }

    return c;
}


ItemSet fr_dump_diff_cu_hg_without(HG_for_XY_without_auxiliary *algo1, HG_for_XY_without_auxiliary *algo2,
                                   CUSketch *xy_a,
                                   CUSketch *xy_b, int thre) // To find union and their counts
{
    ItemSet a, b, c;
    algo1->query_find(a, thre);
    algo2->query_find(b, thre);
    // try another method to dump diff
    for (auto itr: a) {
        if (b.find(itr.first) != b.end()) {
            c[itr.first] = a[itr.first] - b[itr.first];
        } else {
            int op = algo2->query_point(itr.first);
            if (op == 0) {
                int appr_b = xy_b->query(itr.first);
//                cout << " frequency of items: " << itr.first << " in XY is : " << appr_b << endl;
                c[itr.first] = a[itr.first] - appr_b;
            } else {
                c[itr.first] = a[itr.first] - op;
            }
        }

    }

    for (auto itr: b) {
        if (c.find(itr.first) != c.end()) {
            continue;
        }
        if (a.find(itr.first) != a.end()) {
            cout << "error" << endl;
        } else {
            int op = algo1->query_point(itr.first);
            if (op == 0) {
                int appr_a = xy_a->query(itr.first);
//                cout << " frequency of items : " << itr.first << " in XY is : " << appr_a << endl;
                c[itr.first] = b[itr.first] - appr_a;
            } else {
                c[itr.first] = b[itr.first] - op;
            }
        }
    }

    return c;
}

void fr_get_heavy_changer(ItemSet &diff, int threshold, set<uint32_t> &detected) {     // find culprit item
    detected.clear();
    for (auto itr: diff) {
        if (abs(itr.second) >= threshold) {
            detected.insert(itr.first);
        }
    }
}

int fr_get_gt_diff(int package_num, ItemSet &gt_diff) {     // figure out the \sum{i=1}{n}|f_i-\hat{f_{i}}|
    int i;
    gt_diff.clear();
    for (i = 0; i < package_num / 2; i++) {
        gt_diff[insert_data_stream[i]]++;
    }
    for (; i < package_num; i++) {
        gt_diff[insert_data_stream[i]]--;
    }

    int ret = 0;
    for (auto itr: gt_diff) {
        ret += abs(itr.second);
    }

    return ret;
}

double *demo_flow_radar(int num_total, int num_distinct, double thre) {
    double *outp = new double[3];
    printf("\nExp for heavy change detection:\n");
//    const int threshold = 200;
    ItemSet detected_diff, gt_diff;
    set<uint32_t> detected_changer, gt_changer;

    int sum_diff = fr_get_gt_diff(num_total, gt_diff);
    int threshold = floor(thre * sum_diff);
//    int threshold = thre;
    cout << " The sum of exchange frequency of items is : " << sum_diff << endl;
    cout << " The threshold is : " << threshold << endl;
    fr_get_heavy_changer(gt_diff, threshold, gt_changer);

    int num_intersection;
    double precision, recall, f1;

    /*
    // Raw FlowRadar
    auto *fr_a = new InsertableIBLT<memory_radar>();
    auto *fr_b = new InsertableIBLT<memory_radar>();
    fr_a->build(insert_data_stream, num_total / 2);          // window 1
    fr_b->build(insert_data_stream + num_total / 2, num_total / 2);  //window 2
    detected_diff = fr_dump_diff(fr_a, fr_b);
    fr_get_heavy_changer(detected_diff, threshold, detected_changer);
    vector<uint32_t> intersection(num_distinct);
    auto end_itr = std::set_intersection(
            gt_changer.begin(), gt_changer.end(),
            detected_changer.begin(), detected_changer.end(),
            intersection.begin()
    );

    // gt_changer.size() : number of true value of culprit item.
    // detected_changer.size() : estimated number of culprit item.
    // num_intersection : correct reported number of culprit item.

    int repo1 = fr_a->test_complete();
    printf(" The report of the recover result 1: %d \n", repo1);
    int repo2 = fr_b->test_complete();
    printf(" The report of the recover result 2: %d \n", repo2);
    num_intersection = end_itr - intersection.begin();
    printf("number of correctly reported item: %d \n", num_intersection);
    printf("number of total reported item: %d \n", detected_changer.size());
    printf("number of correctly item: %d \n", gt_changer.size());
//    for (auto itr: gt_changer) {
//        printf("set: %d\n", itr);
//    }
    precision = detected_changer.size() ? double(num_intersection) / detected_changer.size() : (gt_changer.size() ? 0
                                                                                                                  : 1);
    recall = gt_changer.size() ? double(num_intersection) / gt_changer.size() : 1;
    f1 = (2 * precision * recall) / (precision + recall);
    printf("\tFlowRadar: %dKB\n", memory_radar / 1024);
    printf("\t  precision: %lf, recall: %lf, f1-score: %lf\n", precision, recall, f1);


     */

    cout << "memory_cf : " << memory_cf << " xuailiary: " << memory_cold_radar << endl;

    // FlowRadar with SC
    auto *sc_fr_a = new SC_InsertableIBLT<memory_cold_radar, memory_cf, bucket_num, threshold_cf>();
    auto *sc_fr_b = new SC_InsertableIBLT<memory_cold_radar, memory_cf, bucket_num, threshold_cf>();


    clock_t time_cf1 = clock();
    sc_fr_a->build(insert_data_stream, num_total / 2);
    clock_t time_cf2 = clock();
    double test_time_cf = (double) (time_cf2 - time_cf1) / CLOCKS_PER_SEC;
    double result_time = (double) num_total / 2 / test_time_cf / 1000000;
    cout << " Efficiency of insert CF(Throughput): " << result_time << endl;
    ofstream ftime("time_heavy.dat", ios::app);
    ftime << me_cm << " " << range << " "
          << result_time << " ";
    ftime.close();
    sc_fr_b->build(insert_data_stream + num_total / 2, num_total / 2);
    detected_diff = fr_dump_diff(sc_fr_a, sc_fr_b);
    fr_get_heavy_changer(detected_diff, threshold, detected_changer);
    vector<uint32_t> intersection(num_distinct);
    auto end_itr = std::set_intersection(
            gt_changer.begin(), gt_changer.end(),
            detected_changer.begin(), detected_changer.end(),
            intersection.begin()
    );
    printf(" The threshold of CF is: %d \n", threshold_cf);
    int repo1 = sc_fr_a->test();
    printf(" The report of the recover result 1: %d \n", repo1);
    int repo2 = sc_fr_b->test();
    printf(" The report of the recover result 2: %d \n", repo2);

    num_intersection = end_itr - intersection.begin();
    printf("number of correctly reported item: %d \n", num_intersection);
    printf("number of total reported item: %d \n", detected_changer.size());
    printf("number of correctly item: %d \n", gt_changer.size());
//    for (auto itr: gt_changer) {
//        printf("set: %d\n", itr);
//    }
    precision = detected_changer.size() ? double(num_intersection) / detected_changer.size() : (gt_changer.size() ? 0
                                                                                                                  : 1);
    recall = gt_changer.size() ? double(num_intersection) / gt_changer.size() : 1;
    f1 = (2 * precision * recall) / (precision + recall);
    printf("\t  precision: %lf, recall: %lf, f1-score: %lf\n", precision, recall, f1);
    outp[0] = precision;
    outp[1] = recall;
    outp[2] = f1;

    return outp;
}

/*
double *demo_flow_radar_2(int num_total, int num_distinct, double thre) {
    double *outp = new double[3];
    printf("\nExp for heavy change detection:\n");
//    const int threshold = 200;
    ItemSet detected_diff, gt_diff;
    set<uint32_t> detected_changer, gt_changer;

    int sum_diff = fr_get_gt_diff(num_total, gt_diff);
    int threshold = floor(thre * sum_diff);
//    int threshold = thre;
    cout << " The sum of exchange frequency of items is : " << sum_diff << endl;
    cout << " The threshold is : " << threshold << endl;
    fr_get_heavy_changer(gt_diff, threshold, gt_changer);

    int num_intersection;
    double precision, recall, f1;


    cout<<"memory_cf : "<<memory_cf<<" xuailiary: "<<memory_cold_radar<<endl;

    // FlowRadar with SC
    auto *sc_fr_a = new SC_InsertableIBLT<memory_cold_radar_2, memory_cf_2, bucket_num, threshold_cf>();
    auto *sc_fr_b = new SC_InsertableIBLT<memory_cold_radar_2, memory_cf_2, bucket_num, threshold_cf>();



    clock_t time_cf1 = clock();
    sc_fr_a->build(insert_data_stream, num_total / 2);
    clock_t time_cf2 = clock();
    double test_time_cf = (double) (time_cf2 - time_cf1) / CLOCKS_PER_SEC;
    double result_time = (double) num_total / 2 / test_time_cf / 1000000;
    cout << " Efficiency of insert CF(time): " << result_time << endl;
    ofstream ftime("time_heavy.dat", ios::app);
    ftime << me_cm << " " << range << " "
          << result_time << " ";
    ftime.close();
    sc_fr_b->build(insert_data_stream + num_total / 2, num_total / 2);
    detected_diff = fr_dump_diff(sc_fr_a, sc_fr_b);
    fr_get_heavy_changer(detected_diff, threshold, detected_changer);
    vector<uint32_t> intersection(num_distinct);
    auto end_itr = std::set_intersection(
            gt_changer.begin(), gt_changer.end(),
            detected_changer.begin(), detected_changer.end(),
            intersection.begin()
    );
    printf(" The threshold of CF is: %d \n", threshold_cf);
    int repo1 = sc_fr_a->test();
    printf(" The report of the recover result 1: %d \n", repo1);
    int repo2 = sc_fr_b->test();
    printf(" The report of the recover result 2: %d \n", repo2);

    num_intersection = end_itr - intersection.begin();
    printf("number of correctly reported item: %d \n", num_intersection);
    printf("number of total reported item: %d \n", detected_changer.size());
    printf("number of correctly item: %d \n", gt_changer.size());
//    for (auto itr: gt_changer) {
//        printf("set: %d\n", itr);
//    }
    precision = detected_changer.size() ? double(num_intersection) / detected_changer.size() : (gt_changer.size() ? 0
                                                                                                                  : 1);
    recall = gt_changer.size() ? double(num_intersection) / gt_changer.size() : 1;
    f1 = (2 * precision * recall) / (precision + recall);
    printf("\tFlowRadar with SC: %dKB FR + %dKB SC (2)\n", memory_cold_radar / 1024, total_memory_in_bytes / 1024);
    printf("\t  precision: %lf, recall: %lf, f1-score: %lf\n", precision, recall, f1);
    outp[0] = precision;
    outp[1] = recall;
    outp[2] = f1;

    return outp;
}

double *demo_flow_radar_3(int num_total, int num_distinct, double thre) {
    double *outp = new double[3];
    printf("\nExp for heavy change detection:\n");
//    const int threshold = 200;
    ItemSet detected_diff, gt_diff;
    set<uint32_t> detected_changer, gt_changer;

    int sum_diff = fr_get_gt_diff(num_total, gt_diff);
    int threshold = floor(thre * sum_diff);
//    int threshold = thre;
    cout << " The sum of exchange frequency of items is : " << sum_diff << endl;
    cout << " The threshold is : " << threshold << endl;
    fr_get_heavy_changer(gt_diff, threshold, gt_changer);

    int num_intersection;
    double precision, recall, f1;

    ItemSet real_set1,real_set2;
    real_set1.clear();
    real_set2.clear();
    for(int i=0;i<num_total / 2;i++){
        real_set1[insert_data_stream[i]]++;
        real_set2[insert_data_stream[i+num_total / 2]]++;
    }
    int cout_set1=0;
    int cout_set2=0;
    for(auto itr:real_set1){
        if(itr.second>=threshold){
            cout_set1++;
        }
    }
    for(auto itr:real_set2){
        if(itr.second>=threshold){
            cout_set2++;
        }
    }
    cout<<"the total number of items whose fre is larger than thre in win1 is : "<<cout_set1<<endl;
    cout<<"the total number of items whose fre is larger than thre in win2 is : "<<cout_set2<<endl;

    cout<<"memory_cf : "<<memory_cf<<" xuailiary: "<<memory_cold_radar<<endl;

    // FlowRadar with SC
    auto *sc_fr_a = new SC_InsertableIBLT<memory_cold_radar_3, memory_cf_3, bucket_num, threshold_cf>();
    auto *sc_fr_b = new SC_InsertableIBLT<memory_cold_radar_3, memory_cf_3, bucket_num, threshold_cf>();

    clock_t time_cf1 = clock();
    sc_fr_a->build(insert_data_stream, num_total / 2);
    clock_t time_cf2 = clock();
    double test_time_cf = (double) (time_cf2 - time_cf1) / CLOCKS_PER_SEC;
    double result_time = (double) num_total / 2 / test_time_cf / 1000000;
    cout << " Efficiency of insert CF(time): " << result_time << endl;
    ofstream ftime("time_heavy.dat", ios::app);
    ftime << me_cm << " " << range << " "
          << result_time << " ";
    ftime.close();
    sc_fr_b->build(insert_data_stream + num_total / 2, num_total / 2);
    detected_diff = fr_dump_diff(sc_fr_a, sc_fr_b);
    fr_get_heavy_changer(detected_diff, threshold, detected_changer);
    vector<uint32_t> intersection(num_distinct);
    auto end_itr = std::set_intersection(
            gt_changer.begin(), gt_changer.end(),
            detected_changer.begin(), detected_changer.end(),
            intersection.begin()
    );
    printf(" The threshold of CF is: %d \n", threshold_cf);
    int repo1 = sc_fr_a->test();
    printf(" The report of the recover result 1: %d \n", repo1);
    int repo2 = sc_fr_b->test();
    printf(" The report of the recover result 2: %d \n", repo2);

    num_intersection = end_itr - intersection.begin();
    printf("number of correctly reported item: %d \n", num_intersection);
    printf("number of total reported item: %d \n", detected_changer.size());
    printf("number of correctly item: %d \n", gt_changer.size());
//    for (auto itr: gt_changer) {
//        printf("set: %d\n", itr);
//    }
    precision = detected_changer.size() ? double(num_intersection) / detected_changer.size() : (gt_changer.size() ? 0
                                                                                                                  : 1);
    recall = gt_changer.size() ? double(num_intersection) / gt_changer.size() : 1;
    f1 = (2 * precision * recall) / (precision + recall);
    printf("\tFlowRadar with SC: %dKB FR + %dKB SC (3)\n", memory_cold_radar / 1024, total_memory_in_bytes / 1024);
    printf("\t  precision: %lf, recall: %lf, f1-score: %lf\n", precision, recall, f1);
    outp[0] = precision;
    outp[1] = recall;
    outp[2] = f1;

    return outp;
}

*/

double *demo_flow_hg(int num_total, int num_distinct, double thre) {
    double *outp = new double[3];
//    printf("\nExp for heavy change detection:\n");
//    const int threshold = 200;
    ItemSet detected_diff, gt_diff;
    set<uint32_t> detected_changer, gt_changer;

    int sum_diff = fr_get_gt_diff(num_total, gt_diff);
    int threshold = floor(thre * sum_diff);
//    int threshold = thre;
    fr_get_heavy_changer(gt_diff, threshold, gt_changer);

    int num_intersection;
    double precision, recall, f1;
    // FlowRadar with Heavy guardian
    printf("\nExp for heavy change detection (heavy_guardian):\n");
//    auto *heavy_a = new Heavy_guardian_flowdadar<memory_heavy_radar>(memory_hg, threshold);
//    auto *heavy_b = new Heavy_guardian_flowdadar<memory_heavy_radar>(memory_hg, threshold);
    auto *heavy_a = new Heavy_guardian_list(memory_hg, memory_heavy_radar, threshold);
    auto *heavy_b = new Heavy_guardian_list(memory_hg, memory_heavy_radar, threshold);


    clock_t time_cf1 = clock();
    heavy_a->build(insert_data_stream, num_total / 2);
    clock_t time_cf2 = clock();
    double test_time_cf = (double) (time_cf2 - time_cf1) / CLOCKS_PER_SEC;
    double result_time = (double) num_total / 2 / test_time_cf / 1000000;
    cout << " Efficiency of insert HG(Throughput): " << result_time << endl;
    ofstream ftime("time_heavy.dat", ios::app);
    ftime << result_time << " ";
    ftime.close();

    heavy_b->build(insert_data_stream + num_total / 2, num_total / 2);

    detected_diff = fr_dump_diff(heavy_a, heavy_b);
    fr_get_heavy_changer(detected_diff, threshold, detected_changer);
    vector<uint32_t> intersection(num_distinct);
    auto end_itr = std::set_intersection(
            gt_changer.begin(), gt_changer.end(),
            detected_changer.begin(), detected_changer.end(),
            intersection.begin()
    );
    int repo1 = heavy_a->test();
    printf(" The report of the recover result 1: %d \n", repo1);
    int repo2 = heavy_b->test();
    printf(" The report of the recover result 2: %d \n", repo2);

    num_intersection = end_itr - intersection.begin();
    printf("number of correctly reported item: %d \n", num_intersection);
    printf("number of total reported item: %d \n", detected_changer.size());
    printf("number of correctly item: %d \n", gt_changer.size());

    precision = detected_changer.size() ? double(num_intersection) / detected_changer.size() : (gt_changer.size() ? 0
                                                                                                                  : 1);
    recall = gt_changer.size() ? double(num_intersection) / gt_changer.size() : 1;
    f1 = (2 * precision * recall) / (precision + recall);
    printf("\t  precision: %lf, recall: %lf, f1-score: %lf\n", precision, recall, f1);
    outp[0] = precision;
    outp[1] = recall;
    outp[2] = f1;
    return outp;
}


double *demo_flow_hg_withlight(int num_total, int num_distinct, double thre) {
    double *outp = new double[3];
    printf("\nExp for heavy change detection:\n");
//    const int threshold = 200;
    ItemSet detected_diff, gt_diff;
    set<uint32_t> detected_changer, gt_changer;

    int sum_diff = fr_get_gt_diff(num_total, gt_diff);
    int threshold = floor(thre * sum_diff);
//    int threshold = thre;
    cout << " The sum of exchange frequency of items is : " << sum_diff << endl;
    cout << " The threshold is : " << threshold << endl;
    fr_get_heavy_changer(gt_diff, threshold, gt_changer);

    int num_intersection;
    double precision, recall, f1;
    // FlowRadar with Heavy guardian
    printf("\nExp for heavy change detection (heavy_guardian_with_light):\n");
//    auto *heavy_a = new Heavy_guardian_flowdadar<memory_heavy_radar>(memory_hg, threshold);
//    auto *heavy_b = new Heavy_guardian_flowdadar<memory_heavy_radar>(memory_hg, threshold);
    auto *heavy_a = new HG_heavy_light(me_cm, threshold);
    auto *heavy_b = new HG_heavy_light(me_cm, threshold);


    clock_t time_cf1 = clock();
    for (int i = 0; i < num_total / 2; i++) {
        heavy_a->insert_hg(insert_data_stream[i], 1);
    }
    clock_t time_cf2 = clock();
    double test_time_cf = (double) (time_cf2 - time_cf1) / CLOCKS_PER_SEC;
    double result_time = (double) num_total / 2 / test_time_cf / 1000000;
    cout << " Efficiency of insert HG(time): " << result_time << endl;
    ofstream ftime("time_heavy.dat", ios::app);
    ftime << result_time << " ";
    ftime.close();

    for (int i = num_total / 2; i < num_total; i++) {
        heavy_b->insert_hg(insert_data_stream[i], 1);
    }

    detected_diff = fr_dump_diff_hg(heavy_a, heavy_b, thre);
    fr_get_heavy_changer(detected_diff, threshold, detected_changer);
    vector<uint32_t> intersection(num_distinct);
    auto end_itr = std::set_intersection(
            gt_changer.begin(), gt_changer.end(),
            detected_changer.begin(), detected_changer.end(),
            intersection.begin()
    );

    num_intersection = end_itr - intersection.begin();
    printf("number of correctly reported item: %d \n", num_intersection);
    printf("number of total reported item: %d \n", detected_changer.size());
    printf("number of correctly item: %d \n", gt_changer.size());

    precision = detected_changer.size() ? double(num_intersection) / detected_changer.size() : (gt_changer.size() ? 0
                                                                                                                  : 1);
    recall = gt_changer.size() ? double(num_intersection) / gt_changer.size() : 1;
    f1 = (2 * precision * recall) / (precision + recall);
    printf("\tFlowRadar with Heavy_guardian_with_light: %dKB heavy + %dKB light\n", hg_me_heavy / 1024,
           hg_me_light / 1024);
    printf("\t  precision: %lf, recall: %lf, f1-score: %lf\n", precision, recall, f1);
    outp[0] = precision;
    outp[1] = recall;
    outp[2] = f1;
    return outp;
}


double *demo_xy_heavy_whole(int num_total, int num_distinct, double thre) {
    double *outp = new double[3];
//    const int threshold = 200;
    ItemSet detected_diff, gt_diff;
    set<uint32_t> detected_changer, gt_changer;

    int sum_diff = fr_get_gt_diff(num_total, gt_diff);
    int threshold = floor(thre * sum_diff);
//    int threshold = thre;
    fr_get_heavy_changer(gt_diff, threshold, gt_changer);

    int num_intersection;
    double precision, recall, f1;

    // FlowRadar with xy-sketch
    printf("\nExp for heavy change detection (PB-sketch):\n");
    int rad[range];
    for (int i = 0; i < range; i++) {
        rad[i] = i;
    }
    Pbsketch_distr_space *xy_count_a = new Pbsketch_distr_space(memory_xy, range, rad);
    Pbsketch_distr_space *xy_count_b = new Pbsketch_distr_space(memory_xy, range, rad);
    HG_for_XY *xy_auxi_a = new HG_for_XY(memory_xy_radar);
    HG_for_XY *xy_auxi_b = new HG_for_XY(memory_xy_radar);

    clock_t time_cf1 = clock();
    for (int i = 0; i < num_total / 2; i++) {
        bool in_hg_a = xy_auxi_a->insert(insert_data_stream[i], 1);
        int in_id = 0, in_fre = 0;
        bool flag_in = false;
        if (!in_hg_a) {
            xy_count_a->insert(insert_data_stream[i], 1);
            int back_fre = xy_count_a->query(insert_data_stream[i]);
            bool exchange = xy_auxi_a->insert_hg(insert_data_stream[i], back_fre, threshold, in_id, in_fre, flag_in);
            if (exchange) {
                xy_count_a->insert(insert_data_stream[i], -back_fre);
                if (flag_in) {
                    xy_count_a->insert(in_id, in_fre);
                }
            }

        }

    }
    clock_t time_cf2 = clock();
    double test_time_cf = (double) (time_cf2 - time_cf1) / CLOCKS_PER_SEC;
    double result_time = (double) num_total / 2 / test_time_cf / 1000000;
    cout << " Efficiency of insert PB(Throughput): " << result_time << endl;
    ofstream ftime("time_heavy.dat", ios::app);
    ftime << result_time << endl;
    ftime.close();

    for (int i = num_total / 2; i < num_total; i++) {
        bool in_hg_b = xy_auxi_b->insert(insert_data_stream[i], 1);
        int in_id = 0, in_fre = 0;
        bool flag_in = false;
        if (!in_hg_b) {

            xy_count_b->insert(insert_data_stream[i], 1);
            int back_fre = xy_count_b->query(insert_data_stream[i]);
            bool exchange = xy_auxi_b->insert_hg(insert_data_stream[i], back_fre, threshold, in_id,
                                                 in_fre, flag_in);
            if (exchange) {
                xy_count_b->insert(insert_data_stream[i], -back_fre);
                if (flag_in) {
                    xy_count_b->insert(in_id, in_fre);
                }
            }

        }
    }

//    ItemSet real_a, real_b;
//    real_a.clear();
//    real_b.clear();
//    for (int i = 0; i < num_total / 2; i++) {
//        real_a[insert_data_stream[i]]++;
//        real_b[insert_data_stream[i + num_total / 2]]++;
//    }

//    detected_diff = fr_dump_diff_xy(xy_auxi_a, xy_auxi_b, xy_count_a, xy_count_b, threshold);
    detected_diff = fr_dump_diff_xy_hg(xy_auxi_a, xy_auxi_b, xy_count_a, xy_count_b, threshold);
    fr_get_heavy_changer(detected_diff, threshold, detected_changer);

/*
    int num_test = 0;
    for (auto itr:gt_changer) {
        if (detected_changer.find(itr) == detected_changer.end()) {
            num_test++;
            cout << " Unreport item is : " << itr << endl;
            int win_a = xy_auxi_a->query_point_test(itr);
            if (win_a == 0) {
                win_a = xy_count_a->query(itr);
                cout << "Finding in XY!" << endl;
            }
            cout << "fre in win_a is : " << win_a << endl;
            int win_b = xy_auxi_b->query_point_test(itr);
            if (win_b == 0) {
                win_b = xy_count_b->query(itr);
                cout << "Finding in XY!" << endl;
            }
            cout << "fre in win_b is : " << win_b << endl;
            cout << " The exact fre in win_a is : " << real_a[itr] << endl;
            cout << " The exact fre in win_b is : " << real_b[itr] << endl;
        }

    }
    cout << " The total number of unreport is : " << num_test << endl;
*/
    vector<uint32_t> intersection(num_distinct);
    auto end_itr = std::set_intersection(
            gt_changer.begin(), gt_changer.end(),
            detected_changer.begin(), detected_changer.end(),
            intersection.begin()
    );

    num_intersection = end_itr - intersection.begin();
    printf("number of correctly reported item: %d \n", num_intersection);
    printf("number of total reported item: %d \n", detected_changer.size());
    printf("number of correctly item: %d \n", gt_changer.size());

    precision = detected_changer.size() ? double(num_intersection) / detected_changer.size() : (gt_changer.size() ? 0
                                                                                                                  : 1);
    recall = gt_changer.size() ? double(num_intersection) / gt_changer.size() : 1;
    f1 = (2 * precision * recall) / (precision + recall);
    printf("\t  precision: %lf, recall: %lf, f1-score: %lf\n", precision, recall, f1);
    outp[0] = precision;
    outp[1] = recall;
    outp[2] = f1;
    return outp;
}


double *demo_xy_heavy_without_auxiliary(int num_total, int num_distinct, double thre) {
    double *outp = new double[3];
//    const int threshold = 200;
    ItemSet detected_diff, gt_diff;
    set<uint32_t> detected_changer, gt_changer;

    int sum_diff = fr_get_gt_diff(num_total, gt_diff);
    int threshold = floor(thre * sum_diff);
//    int threshold = thre;
    fr_get_heavy_changer(gt_diff, threshold, gt_changer);

    int num_intersection;
    double precision, recall, f1;

    // FlowRadar with xy-sketch
    printf("\nExp for heavy change detection (XY_without):\n");
    int rad[range];
    for (int i = 0; i < range; i++) {
        rad[i] = i;
    }
    Pbsketch_distr_space *xy_count_a = new Pbsketch_distr_space(memory_xy, range, rad);
    Pbsketch_distr_space *xy_count_b = new Pbsketch_distr_space(memory_xy, range, rad);
    HG_for_XY_without_auxiliary *xy_auxi_a = new HG_for_XY_without_auxiliary(memory_xy_radar);
    HG_for_XY_without_auxiliary *xy_auxi_b = new HG_for_XY_without_auxiliary(memory_xy_radar);

    clock_t time_cf1 = clock();
    for (int i = 0; i < num_total / 2; i++) {
        bool in_hg_a = xy_auxi_a->insert(insert_data_stream[i], 1);
        int in_id = 0, in_fre = 0;
        if (!in_hg_a) {
            xy_count_a->insert(insert_data_stream[i], 1);
            int back_fre = xy_count_a->query(insert_data_stream[i]);
            bool exchange = xy_auxi_a->insert_hg(insert_data_stream[i], back_fre, threshold, in_id, in_fre);
            if (exchange) {
                xy_count_a->insert(insert_data_stream[i], -back_fre);
                xy_count_a->insert(in_id, in_fre);

            }

        }

    }
    clock_t time_cf2 = clock();
    double test_time_cf = (double) (time_cf2 - time_cf1) / CLOCKS_PER_SEC;
    double result_time = (double) num_total / 2 / test_time_cf / 1000000;
    cout << " Efficiency of insert XY(time): " << result_time << endl;
    ofstream ftime("time_heavy.dat", ios::app);
    ftime << result_time << endl;
    ftime.close();

    for (int i = num_total / 2; i < num_total; i++) {
        bool in_hg_b = xy_auxi_b->insert(insert_data_stream[i], 1);
        int in_id = 0, in_fre = 0;
        if (!in_hg_b) {
            xy_count_b->insert(insert_data_stream[i], 1);
            int back_fre = xy_count_b->query(insert_data_stream[i]);
            bool exchange = xy_auxi_b->insert_hg(insert_data_stream[i], back_fre, threshold, in_id,
                                                 in_fre);
            if (exchange) {
                xy_count_b->insert(insert_data_stream[i], -back_fre);
                xy_count_b->insert(in_id, in_fre);

            }

        }
    }


    detected_diff = fr_dump_diff_xy_hg_without(xy_auxi_a, xy_auxi_b, xy_count_a, xy_count_b, threshold);
    fr_get_heavy_changer(detected_diff, threshold, detected_changer);

    vector<uint32_t> intersection(num_distinct);
    auto end_itr = std::set_intersection(
            gt_changer.begin(), gt_changer.end(),
            detected_changer.begin(), detected_changer.end(),
            intersection.begin()
    );

    num_intersection = end_itr - intersection.begin();
    printf("number of correctly reported item: %d \n", num_intersection);
    printf("number of total reported item: %d \n", detected_changer.size());
    printf("number of correctly item: %d \n", gt_changer.size());

    precision = detected_changer.size() ? double(num_intersection) / detected_changer.size() : (gt_changer.size() ? 0
                                                                                                                  : 1);
    recall = gt_changer.size() ? double(num_intersection) / gt_changer.size() : 1;
    f1 = (2 * precision * recall) / (precision + recall);
    printf("\t  precision: %lf, recall: %lf, f1-score: %lf\n", precision, recall, f1);
    outp[0] = precision;
    outp[1] = recall;
    outp[2] = f1;
    return outp;
}

/*
double *demo_cu_heavy_without_auxiliary(int num_total, int num_distinct, double thre) {
    double *outp = new double[3];
//    const int threshold = 200;
    ItemSet detected_diff, gt_diff;
    set<uint32_t> detected_changer, gt_changer;

    int sum_diff = fr_get_gt_diff(num_total, gt_diff);
    int threshold = floor(thre * sum_diff);
//    int threshold = thre;
    fr_get_heavy_changer(gt_diff, threshold, gt_changer);

    int num_intersection;
    double precision, recall, f1;

    // FlowRadar with xy-sketch
    printf("\nExp for heavy change detection (CU_without):\n");
    int rad[range];
    for (int i = 0; i < range; i++) {
        rad[i] = i;
    }
    CUSketch *cu_count_a = new CUSketch(memory_xy, 3);
    CUSketch *cu_count_b = new CUSketch(memory_xy, 3);

    HG_for_XY_without_auxiliary *xy_auxi_a = new HG_for_XY_without_auxiliary(memory_xy_radar);
    HG_for_XY_without_auxiliary *xy_auxi_b = new HG_for_XY_without_auxiliary(memory_xy_radar);

    clock_t time_cf1 = clock();
    for (int i = 0; i < num_total / 2; i++) {
        bool in_hg_a = xy_auxi_a->insert(insert_data_stream[i], 1);
        int in_id = 0, in_fre = 0;
        if (!in_hg_a) {
            cu_count_a->insert(insert_data_stream[i], 1);
            int back_fre = cu_count_a->query(insert_data_stream[i]);
            bool exchange = xy_auxi_a->insert_hg(insert_data_stream[i], back_fre, threshold, in_id, in_fre);
            if (exchange) {
                cu_count_a->insert(insert_data_stream[i], -back_fre);
                cu_count_a->insert(in_id, in_fre);

            }

        }

    }
    clock_t time_cf2 = clock();
    double test_time_cf = (double) (time_cf2 - time_cf1) / CLOCKS_PER_SEC;
    double result_time = (double) num_total / 2 / test_time_cf / 1000000;
    cout << " Efficiency of insert XY(time): " << result_time << endl;
    ofstream ftime("time_heavy.dat", ios::app);
    ftime << result_time << endl;
    ftime.close();

    for (int i = num_total / 2; i < num_total; i++) {
        bool in_hg_b = xy_auxi_b->insert(insert_data_stream[i], 1);
        int in_id = 0, in_fre = 0;
        if (!in_hg_b) {
            cu_count_b->insert(insert_data_stream[i], 1);
            int back_fre = cu_count_b->query(insert_data_stream[i]);
            bool exchange = xy_auxi_b->insert_hg(insert_data_stream[i], back_fre, threshold, in_id,
                                                 in_fre);
            if (exchange) {
                cu_count_b->insert(insert_data_stream[i], -back_fre);
                cu_count_b->insert(in_id, in_fre);

            }

        }
    }


    detected_diff = fr_dump_diff_cu_hg_without(xy_auxi_a, xy_auxi_b, cu_count_a, cu_count_b, threshold);
    fr_get_heavy_changer(detected_diff, threshold, detected_changer);

    vector<uint32_t> intersection(num_distinct);
    auto end_itr = std::set_intersection(
            gt_changer.begin(), gt_changer.end(),
            detected_changer.begin(), detected_changer.end(),
            intersection.begin()
    );

    num_intersection = end_itr - intersection.begin();
    printf("number of correctly reported item: %d \n", num_intersection);
    printf("number of total reported item: %d \n", detected_changer.size());
    printf("number of correctly item: %d \n", gt_changer.size());

    precision = detected_changer.size() ? double(num_intersection) / detected_changer.size() : (gt_changer.size() ? 0
                                                                                                                  : 1);
    recall = gt_changer.size() ? double(num_intersection) / gt_changer.size() : 1;
    f1 = (2 * precision * recall) / (precision + recall);
    printf("\tHG with cu-sketch: %dKB FR + %dKB CU(without)\n", memory_xy_radar / 1024, me_cm / 1024);
    printf("\t  precision: %lf, recall: %lf, f1-score: %lf\n", precision, recall, f1);
    outp[0] = precision;
    outp[1] = recall;
    outp[2] = f1;
    return outp;
}

*/

/*
// test all (modify)
int main() {

    clock_t start_time = clock();
    int *as = analysis_data(filename, 100);

//    const int num_l = 5;
    const int num_l = 1;
//    const int thre[num_l] = {5000, 10000, 20000, 30000, 40000, 50000};
//    const double thre[num_l] = {0.0002, 0.0004, 0.0006, 0.0008, 0.001};
//    const double thre[num_l] = {0.0004, 0.0005, 0.0006, 0.0007, 0.0008};
    const double thre[num_l] = {0.0005};
    const int num_win = 1;
//    const int win_range[num_win] = {4, 5, 6, 7, 8};
    const int win_range[num_win] = {8};

    for (int ij = 0; ij < num_l; ij++) {
        for (int kk = 0; kk < num_win; kk++) {

//            double *radar_cf_out = demo_flow_radar(as[1] * 2 / win_range[kk], as[0], thre[ij]);
//            double *radar_hg_out = demo_flow_hg(as[1] * 2 / win_range[kk], as[0], thre[ij]);
//            double *heavy_out = demo_xy_heavy_whole(as[1] * 2 / win_range[kk], as[0], thre[ij]);
            double *heavy_xy = demo_xy_heavy_without_auxiliary(as[1] * 2 / win_range[kk], as[0], thre[ij]);
//            double *radar_hg_light = demo_flow_hg_withlight(as[1] * 2 / win_range[kk], as[0], thre[ij]);
//            double *heavy_cu = demo_cu_heavy_without_auxiliary(as[1] * 2 / win_range[kk], as[0], thre[ij]);
            ofstream fout("modify_xy_without.dat", ios::app);
            fout << win_range[kk] << " " << me_kb << " " << memory_xy_radar << " " << threshold_cf << " " << thre[ij]
                 << " "
                 //             << xy_out[0] << " " << xy_out[1] << " " << xy_out[2] << " "
//                 << radar_cf_out[0] << " " << radar_cf_out[1] << " " << radar_cf_out[2] << " "
//                 << radar_hg_out[0] << " " << radar_hg_out[1] << " " << radar_hg_out[2] << " "
//                 << heavy_out[0] << " " << heavy_out[1] << " " << heavy_out[2] << " "
                 << heavy_xy[0] << " " << heavy_xy[1] << " " << heavy_xy[2] << endl;
//                 << radar_hg_light[0] << " " << radar_hg_light[1] << " " << radar_hg_light[2] << endl;
//                 << heavy_cu[0] << " " << heavy_cu[1] << " " << heavy_cu[2] << endl;

            fout.close();
        }
    }
    clock_t end_time = clock();
    cout << "The total time is :" << (double) (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;

    return 0;
}

*/


// test all
int main() {

    clock_t start_time = clock();
    int *as = analysis_data(filename, 100);

//    const int num_l = 5;
    const int num_l = 1;
//    const int thre[num_l] = {5000, 10000, 20000, 30000, 40000, 50000};
//    const double thre[num_l] = {0.0002, 0.0004, 0.0006, 0.0008, 0.001};
//    const double thre[num_l] = {0.0004, 0.0005, 0.0006, 0.0007, 0.0008};
    const double thre[num_l] = {0.0005};
    const int num_win = 1;
//    const int win_range[5] = {2, 3, 4, 5, 6};
    const int win_range[num_win] = {4};

    for (int ij = 0; ij < num_l; ij++) {
        for (int kk = 0; kk < num_win; kk++) {


            double *radar_cf_out = demo_flow_radar(as[1] * 2 / win_range[kk], as[0], thre[ij]);
            double *radar_hg_out = demo_flow_hg(as[1] * 2 / win_range[kk], as[0], thre[ij]);
            double *heavy_out = demo_xy_heavy_whole(as[1] * 2 / win_range[kk], as[0], thre[ij]);

            ofstream fout("new_web_heavy.dat", ios::app);
            fout << win_range[kk] << " " << me_kb << " " << memory_xy_radar << " " << memory_radar << " " << thre[ij]
                 << " "
                 //             << xy_out[0] << " " << xy_out[1] << " " << xy_out[2] << " "
                 << radar_cf_out[0] << " " << radar_cf_out[1] << " " << radar_cf_out[2] << " "
                 << radar_hg_out[0] << " " << radar_hg_out[1] << " " << radar_hg_out[2] << " "
                 << heavy_out[0] << " " << heavy_out[1] << " " << heavy_out[2] << endl;

            fout.close();
        }
    }
    clock_t end_time = clock();
    cout << "The total time is :" << (double) (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;

    return 0;
}


