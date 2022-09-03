#pragma GCC optimize(2)

#include "benchmark.h"
#include <string>
#include <iostream>
#include <chrono>
#include <climits>
#include <vector>
#include <iostream>
#include <sys/mman.h>
#include <sys/types.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <fstream>
#include "SS.h"
#include "USS.h"
#include "Count_Heap.h"
#include "WavingSketch.h"

const long maxn = 100000000;
using namespace std;
const char *filename = "/home/yqliu/data/kosarak.dat";  // (D=210385)
//const char *filename = "/home/yqliu/data/webdocs.dat";  //(D=24580773)
const int me_cm = 160 * 1024;
const int me_cm_kb=me_cm/1024;
const int maxtopk = 1024;
const int interval = 8019015 * 0.25;
//const int interval=299887139*0.25;
//const int thre=74;
//const int thre=8933;

struct Node1 {
    int cnt = 0;//count nubmers of every item
    int index;//record index of item in hashtable
};
struct nod {
    data_type id;
    count_type cnt;
};
Node1 hashtable[maxn];

int cmp(Node1 a, Node1 b) {
    return a.cnt > b.cnt;
}

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

/*
int main() {
	// =========Dataset=========
	// | fid=0: demo           |
	// | fid=1: synthetic      |
	// | fid=2: CAIDA18        |
	// | fid=3: webpage        |
	// | fid=4: campus         |
	// =========================


	int fid = 1;
	printf("\033[0m\033[1;32m===========================================================================================\n\033[0m");
	printf("\033[0m\033[1;32m|                                Waving Sketch Experiments                                |\n\033[0m");
	printf("\033[0m\033[1;32m===========================================================================================\n\033[0m");
	printf("\033[0m\033[1;32m|                                   1. Test on Accuracy                                   |\n\033[0m");
	printf("\033[0m\033[1;32m===========================================================================================\n\033[0m");
	BenchCmp(fid);
	printf("\033[0m\033[1;32m===========================================================================================\n\033[0m");
	printf("\033[0m\033[1;32m|                                  2. Test on Throughput                                  |\n\033[0m");
	printf("\033[0m\033[1;32m===========================================================================================\n\033[0m");
	BenchThp(fid);
	printf("\033[0m\033[1;32m===========================================================================================\n\033[0m");

}
*/


int main(){   // top-k queries
    clock_t start_time = clock();
    int *as = analysis_data(filename, 100);
    cout << "All skecthes's memory is : " << me_cm << endl;
    //cout << "The filter size of A-sketch is : " << a_filter << endl;
    cout << "***************************************" << endl;

    const int num_k = 6;
    int find_k[num_k] = {32, 64, 128, 256, 512, 1024};

    ifstream fin;
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir

    uint32_t temp;

    clock_t in_start=clock();

    auto wavings = new WavingSketch<8,1>((uint32_t)(me_cm-maxtopk*log2(maxtopk)*8) / (8*8+1*2));

    cout<<"Inserting!!!"<<endl;
    while (!fin.eof() && fin >> temp) {
        wavings->Insert(temp);
    }
    fin.close();
    vector<pair<data_type ,count_type >> result(maxtopk);
    wavings->gettopk(maxtopk,result);
    clock_t in_end=clock();
    double in_time=(double) (in_end - in_start) / CLOCKS_PER_SEC;
            cout << "The total time is :" << in_time << " s" << endl;

    for(int i=0;i<num_k;i++){
        int k=find_k[i];
        double num = 0;
        double aae = 0;
        double are = 0;
        double *error_three = new double[3];
        for (int i = 0; i < 3; i++) {
            error_three[i] = 0;
        }
        vector<int> open_close(k,0);

        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) {
                if (result[i].first == hashtable[j].index) {
                    num = num + 1;
                    aae = aae + abs(abs(result[i].second) - hashtable[j].cnt);
                    //cout<<"ID: "<<result[i].first<<" real: "<<hashtable[j].cnt<<" estimated: "<<result[i].second<<endl;
                    are = are + (double) abs(abs(result[i].second) - hashtable[j].cnt) / hashtable[j].cnt;
                    open_close[j] = 1;
                    break;
                }
            }
        }

        for (int i = 0; i < k; i++) {
            if (open_close[i] != 1) {
                aae = aae + hashtable[i].cnt;
                are = are + 1;
            }
        }

        error_three[0] = num / k;
        error_three[1] = aae / k;  // AAE (10^3)
        error_three[2] = are / k;
        printf("\tWavingSketch \n");
        cout<<"Top-k: "<<k<<endl;
        printf("\t  accuracy %lf, AAE %lf, ARE %lf\n", error_three[0], error_three[1], error_three[2]);

        ofstream fout("ws_web_topk.dat", ios::app);

        fout << me_cm << " " << k << " " << error_three[0] << " " << error_three[1] << " "
             << error_three[2]<< " "<<in_time
            << endl;

        fout.close();

    }

    clock_t end_time = clock();
    cout << "The total time is :" << (double) (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;
    return 0;
}

/*
int main() { // heavy changes
    clock_t start_time = clock();
    cout << "All skecthes's memory is : " << me_cm / 1024 << " KB" << endl;
    //cout << "The filter size of A-sketch is : " << a_filter << endl;
    cout << "The interval items of each window is: " << interval << endl;
    cout << "***************************************" << endl;
    unordered_map<data_type, count_type> mp;

    ifstream fin;
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir
    uint32_t temp;
    int num = 0;
    while (!fin.eof() && fin >> temp) {
        if (num > 2 * interval) {
            break;
        }
        if (num < interval) {
            if (mp.find(temp) == mp.end()) {
                mp[temp] = 1;
            } else {
                mp[temp] += 1;
            }
        } else {
            if (mp.find(temp) == mp.end()) {
                mp[temp] = -1;
            } else {
                mp[temp]--;
            }
        }
        num++;
    }
    fin.close();
    auto wavingsa = new WavingSketch<8, 1>(me_cm / (8 * 8 + 1 * 2));
    auto wavingsb = new WavingSketch<8, 1>(me_cm / (8 * 8 + 1 * 2));
    clock_t in_start = clock();
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir
    num = 0;
    while (!fin.eof() && fin >> temp) {
        if (num > 2 * interval) {
            break;
        }
        if (num < interval) {
            wavingsa->Insert(temp);
        } else {
            wavingsb->Insert(temp);
        }
        num++;
    }
    fin.close();
    clock_t in_end = clock();
    double cost_time = (in_end - in_start) / (double) CLOCKS_PER_SEC;
    cout << "The inserting time is : " <<cost_time<< endl;

    const int thre_num = 5;
    const int thre[5] = {60, 75, 90, 105, 120};
//    const int thre[5]={7146,8933,10719,12506,14292};
    for (int i = 0; i < thre_num; i++) {
        cout << "********************************************" << endl;
        cout << "The threshold is : " << thre[i] << endl;
        cout << "The results of WavingSketch: " << endl;
        int all = 0, size = 0, hit = 0;
        double cr, pr, f1;
        for (auto it = mp.begin(); it != mp.end(); ++it) {
            int value = abs(abs(wavingsa->Query(it->first)) - abs(wavingsb->Query(it->first)));
            if (abs(it->second) > thre[i]) {
                all++;
                if (value > thre[i]) {
                    hit += 1;
                }
            }
            if (value > thre[i])
                size += 1;
        }
        cr = hit / (double) all;
        pr = hit / (double) size;
        f1 = 2 * pr * cr / (pr + cr);
        cout << "Hit : " << hit << "  the exact number of items is : " << all << " The estimated number is " << size
             << endl;
        cout << " The precision is: " << pr << endl;
        cout << " The recall is: " << cr << endl;
        cout << " The f1-score is: " << f1 << endl;


        ofstream fout("ws_web_heavyc.dat", ios::app);
        fout << me_cm_kb << " " << thre[i] << " " << pr << " " << cr << " "
             << f1<< " "<<cost_time
             << endl;

        fout.close();
    }


    clock_t end_time = clock();
    cout << "The total time is :" << (double) (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;
    return 0;
}
*/