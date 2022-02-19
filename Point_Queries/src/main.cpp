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

using namespace std;
const long maxn = 100000000;
//  setting

int range = 16;  // range of data items: 2^{range}
int me_cm = 160* 1024;   // memory size of sketches  (B)
const int total_memory_in_bytes =160 * 1024;  //memory size of CF
int a_filter = 32;     // number of items in A-sketch's filter
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

    printf("AAE of the PB-sketch : %lf\n", aae[0]);
    printf("ARE of the PB-sketch : %lf\n", aae[1]);
    printf("***********************************\n");

    return aae;
}



int main() {
    clock_t start_time = clock();

    const char *filename = "/home/yqliu/data/kosarak.dat";
//    const char *filename = "/home/yqliu/data/webdocs.dat";

    int *as = analysis_data(filename);

    cout << "All skecthes's memory is : " << me_cm << endl;
    cout << "The filter size of A-sketch is : " << a_filter << endl;
    cout << "***************************************" << endl;

// undating    input:  filename, memory of sketch, the value of d
    CMSketch *cm_counter = cm_sketch_update(filename, me_cm, 3);
// querying    input:  true frequency, sketch, number of n
    double *cm_aae = cm_sketch_query(hashtable, *cm_counter, as[0]);


    CUSketch *cu_counter = cu_sketch_update(filename, me_cm, 3);
    double *cu_aae = cu_sketch_query(hashtable, *cu_counter, as[0]);

    Csketch *c_counter = c_sketch_update(filename, me_cm, 3);
    double *c_aae = c_sketch_query(hashtable, *c_counter, as[0]);

    // undating    input:  filename, memory of sketch, the number items stored in filter
    ASketch *a_counter = a_sketch_update(filename, me_cm, 3, a_filter);
    double *a_aae = a_sketch_query(hashtable, *a_counter, as[0]);

    // undating    input:  memory of sketch, percent of CF, number of filter (filename)
    CUSketchWithSC<total_memory_in_bytes, filter_memory_percent, bucket_num> *cold_counter = cold_sketch_update<total_memory_in_bytes, filter_memory_percent, bucket_num>(
            filename);
    double *cold_40 = cold_sketch_query(hashtable, *cold_counter, as[0]);

    CUSketchWithSC<total_memory_in_bytes, filter_memory_percent_70, bucket_num> *cold_counter_70 = cold_sketch_update<total_memory_in_bytes, filter_memory_percent_70, bucket_num>(
            filename);
    double *cold_70 = cold_sketch_query_70(hashtable, *cold_counter_70, as[0]);

    CUSketchWithSC<total_memory_in_bytes, filter_memory_percent_90, bucket_num> *cold_counter_90 = cold_sketch_update<total_memory_in_bytes, filter_memory_percent_90, bucket_num>(
            filename);
    double *cold_90 = cold_sketch_query_90(hashtable, *cold_counter_90, as[0]);


    ifstream fin;
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir

    if (!fin)//detect exception
    {
        cout << "Open file failed..." << endl;
        exit(1);
    }
    int temp;
    int count = 0;
    const int thre=50000;      //set N1 to 50000;
    int Max = -1;
    int difference = 0;
    while (!fin.eof() && fin >> temp) {
        if (count < thre) {
            if (temp > Max) {
                Max = temp;
            }
            if (temp > maxn)//detect exception
            {
                cout << "big item..." << endl;
                exit(1);
            }
            distribution[temp].cnt++;
            distribution[temp].index = temp;
        } else {
            break;
        }
        count += 1;
    }
    sort(distribution, distribution + Max + 1, cmp);

    for (int i = 0; i <= Max; i++) {
        if (!distribution[i].cnt)
            break;
        difference += 1;
    }
    fin.close();

    int bin[range];
    int bin_0[range];
    int bin_1[range];
    for (int i = 0; i < range; i++) {
        bin_0[i] = 0;
        bin_1[i] = 0;
    }

    for (int i = 0; i < Max; i++) {
        for (int j = 0; j < range; j++) {
            bin[j] = distribution[i].index & (1 << j);
            bin[j] = bin[j] >> j;
            if (bin[j] > 0) {
                bin_1[j] += distribution[i].cnt;
            } else {
                bin_0[j] += distribution[i].cnt;
            }
        }
    }

    double entr[range];
    for (int i = 0; i < range; i++) {     //function
        double p1 = (double) bin_0[i] / (double) (bin_0[i] + bin_1[i]);
        double p2 = (double) bin_1[i] / (double) (bin_0[i] + bin_1[i]);
        if (p1 == 0 || p2 == 0) {
            entr[i] = 0;
        } else {
            entr[i] = (-p1 * log(p1) / log(2) - p2 * log(p2) / log(2));    //entropy
        }
    }
    int order_pb[range];
    for (int i = 0; i < range; i++) {
        double order_min = -1;
        for (int j = 0; j < range; j++) {        // ording
            if (entr[j] >= order_min) {
                order_min = entr[j];
                order_pb[i] = j;
            }
        }
        entr[order_pb[i]] = -1;
    }

    delete distribution;

    // undating    input:  filename, memory of sketch, the number items stored in filter, range of items, decomposition , N1
    Pbsketch_distr_space *pb_dis_count = pb_distr_space_update(filename, me_cm, range, order_pb, thre);
    // querying    input:  true frequency, sketch, number of n
    double *pb_dis_aae = pb_dis_space_query(hashtable, *pb_dis_count, as[0]);


    ofstream fout("kos_point_error.dat", ios::app);
    fout << me_cm << " " << me_kb << " " << cm_aae[0] << " " << cm_aae[1] << " "
         << cu_aae[0] << " " << cu_aae[1] << " " << c_aae[0] << " " << c_aae[1] << " " << a_aae[0] << " "
         << a_aae[1] << " " << cold_40[0] << " " << cold_40[1] << " " << cold_70[0] << " " << cold_70[1] << " "
         << cold_90[0] << " " << cold_90[1] << " " << pb_dis_aae[0] << " " << pb_dis_aae[1] << endl;

    fout.close();

    clock_t end_time = clock();
    cout << "The total time is :" << (double) (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;
    return 0;
}