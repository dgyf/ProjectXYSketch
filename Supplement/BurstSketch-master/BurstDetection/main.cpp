#include<bits/stdc++.h>
#include "BurstDetector.h"
#include "CorrectBurstDetector.h"
#include "CMBurstDetector.h"

using namespace std;
map<uint64_t, int> D;
vector<pair<int, int> > S[2000000];
int w = 0;
const char *filename = "/home/yqliu/data/kosarak.dat";  // (D=210385)
//const char *filename = "/home/yqliu/data/webdocs.dat";  //(D=24580773)
const int interval = 8019015 * 0.25;
//const int interval=299887139*0.25;
const int me_cm = 560 * 1024;
const int mem = me_cm / 1024;

int main() {
    //int mem = 60; // the size of memory
    double l = 0.3; // the ratio of the Running Track threshold to the burst threshold
    double r12 = 1.75; // the ratio of the size of Stage 1 to the size of Stage 2
    int screen_layer_threshold = l * threshold; // Running Track threshold
    int log_size = mem * 1024 / (8 * r12 + 12) / bucket_size; // number of buckets in Stage 2
//    int log_size = mem * 1024 / (8 * r12 + 12) / bucket_size; // number of buckets in Stage 2
//    int log_size = mem * 1024 / (12 * r12 + 20) / bucket_size; // number of buckets in Stage 2
    int screen_layer_size = log_size * r12 * bucket_size; // number of buckets in Stage 1
    int cm_size = mem * 1024 / (window_num + 2) / 4; // the size of CM Sketch
    int window = 0, cnt = 0;
    clock_t start_time = clock();
    cout << "All skecthes's memory is : " << me_cm / 1024 << " KB" << endl;
    //cout << "The filter size of A-sketch is : " << a_filter << endl;
    cout << "The interval items of each window is: " << interval << endl;
    cout << "***************************************" << endl;
    unordered_map<uint64_t, uint32_t> mp;

    ifstream fin;
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir
//    uint32_t temp;
    uint64_t temp;
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
    BurstDetector A(screen_layer_size, screen_layer_threshold, log_size, threshold);
    BurstDetector B(screen_layer_size, screen_layer_threshold, log_size, threshold);

    clock_t in_start = clock();
    fin.open(filename, ios::in | ios::out);//note: modify path to meet your dir
    num = 0;
    while (!fin.eof() && fin >> temp) {
        if (num > 2 * interval) {
            break;
        }
        if (num < interval) {
            A.insert(temp, 0);
        } else {
            A.insert(temp, 1);
//            B.insert(temp,0);
        }
        num++;
    }
    fin.close();
    clock_t in_end = clock();
    double cost_time = (in_end - in_start) / (double)CLOCKS_PER_SEC;
    cout << "The inserting time is : " << cost_time << endl;

    const int thre_num = 5;
    const int thre[5] = {60, 75, 90, 105, 120};
//    const int thre[5]={7146,8933,10719,12506,14292};
    for (int i = 0; i < thre_num; i++) {
        cout << "********************************************" << endl;
        cout << "The threshold is : " << thre[i] << endl;
        cout << "The results of BurstSketch: " << endl;
        int all = 0, size = 0, hit = 0;
        double cr, pr, f1;
        for (auto it = mp.begin(); it != mp.end(); ++it) {
//            int value = abs((int) A.Query(it->first, 0) - (int) B.Query(it->first, 0));
            int value = abs((int) A.Query(it->first, 0) - (int) A.Query(it->first, 1));
//            int value = abs( (int)A.log.query_log(it->first,0)- (int)A.log.query_log(it->first,1));
            if (abs((int) it->second) >= thre[i]) {
                all++;
                if (value > thre[i]) {
                    hit += 1;
                }
//                cout << "(Real) The id is : " << it->first << " The value is : " << value << " The 0-window is : "
//                     << A.Query(it->first, 0) << " The 1-window is : " << A.Query(it->first, 1) << endl;
            }
            if (value > thre[i]) {
                size += 1;
//                cout << "The id is : " << it->first << " The value is : " << value << " The 0-window is : "
//                     << A.Query(it->first, 0) << " The 1-window is : " << A.Query(it->first, 1) << "The real fre is: "
//                     << abs((int) it->second)  << endl;
            }
        }
        cr = hit / (double) all;
        pr = hit / (double) size;
        f1 = 2 * pr * cr / (pr + cr);
        cout << "Hit : " << hit << "  the exact number of items is : " << all << " The estimated number is " << size
             << endl;
        cout << " The precision is: " << pr << endl;
        cout << " The recall is: " << cr << endl;
        cout << " The f1-score is: " << f1 << endl;

        ofstream fout("bursts_kos_heavyc.dat", ios::app);
        fout << mem << " " << thre[i] << " " << pr << " " << cr << " "
             << f1 << " " << cost_time<<" "<<threshold<<" "<<l<<" "<<r12
             << endl;

        fout.close();
    }

    clock_t end_time = clock();
    cout << "The total time is :" << (double) (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;

    return 0;
}

/*
int main()
{
    int mem = 60; // the size of memory
    double l = 0.3; // the ratio of the Running Track threshold to the burst threshold
    double r12 = 3.75; // the ratio of the size of Stage 1 to the size of Stage 2
    int screen_layer_threshold = l * threshold; // Running Track threshold
    int log_size = mem * 1024 / (12 * r12 + 20) / bucket_size; // number of buckets in Stage 2
    int screen_layer_size = log_size * r12 * bucket_size; // number of buckets in Stage 1
    int cm_size = mem * 1024 / (window_num + 2) / 4; // the size of CM Sketch
    BurstDetector A(screen_layer_size, screen_layer_threshold, log_size, threshold);
    CorrectBurstDetector B(threshold);
    CMBurstDetector C(cm_size, threshold);
    uint64_t tmp;
    int window = 0, cnt = 0;
    for(int i = 0; i < 20000000; i++)
    {

        uint64_t tim, id;
        fin.read((char *)&tim, sizeof(char) * 8);
        fin.read((char *)&id, sizeof(char) * 8);

        cnt++;
        if(cnt > 40000)
        {
            cnt = 0;
            window++;
        }
        A.insert(id, window);
        B.insert(id, window);
        C.insert(id, window);
    }

    printf("Totally %d distinct flows\n", B.w);
    printf("BurstSketch totally reports %d bursts!\n", A.log.Record.size());
    printf("CM Burst detector totally reports %d bursts!\n", C.Record.size());
    printf("Actually there exist %d bursts!\n", B.Record.size());

    D.clear();
    for(int i = 0; i < B.Record.size(); i++)
    {
        if(D.find(B.Record[i].flow_id) == D.end())
            D[B.Record[i].flow_id] = w++;
        S[D[B.Record[i].flow_id]].push_back(make_pair(B.Record[i].start_window, B.Record[i].end_window));
    }
    int bd_correct = 0, bd_overlap = 0, bd_error = 0;
    int cm_correct = 0, cm_overlap = 0, cm_error = 0;
    // calculate experimental results
    for(int i = 0; i < A.log.Record.size(); i++)
    {
        if(D.find(A.log.Record[i].flow_id) == D.end())
            continue;
        int z = D[A.log.Record[i].flow_id];
        for(int j = 0; j < S[z].size(); j++)
        {
            if(A.log.Record[i].start_window == S[z][j].first && A.log.Record[i].end_window == S[z][j].second)
            {
                bd_correct++;
                break;
            }
            if(A.log.Record[i].start_window <= S[z][j].second && A.log.Record[i].end_window >= S[z][j].first)
            {
                bd_overlap++;
                bd_error += abs((S[z][j].second - S[z][j].first) - ((int)A.log.Record[i].end_window - (int)A.log.Record[i].start_window));
                break;
            }
        }
    }

    for(int i = 0; i < C.Record.size(); i++)
    {
        if(D.find(C.Record[i].flow_id) == D.end())
            continue;
        int z = D[C.Record[i].flow_id];
        for(int j = 0; j < S[z].size(); j++)
        {
            if(C.Record[i].start_window == S[z][j].first && C.Record[i].end_window == S[z][j].second)
            {
                cm_correct++;
                break;
            }
            if(C.Record[i].start_window <= S[z][j].second && C.Record[i].end_window >= S[z][j].first)
            {
                cm_overlap++;
                cm_error += abs((S[z][j].second - S[z][j].first) - ((int)C.Record[i].end_window - (int)C.Record[i].start_window));
                break;
            }

        }
    }
    // output BurstSketch results
    printf("BurstSketch:\n");
    printf("Precision: %.5lf\n", 100.0 * (bd_correct + bd_overlap) / A.log.Record.size());
    printf("Recall: %.5lf\n", 100.0 *(bd_correct + bd_overlap) / B.Record.size());
    printf("Average deviation: %.5lf\n", 1.0 * bd_error / (bd_correct + bd_overlap));
    printf("Deviation rate: %.5lf\n\n", 100.0 * bd_overlap / (bd_correct + bd_overlap));

    // output the strawman solution's results
    printf("The strawman solution:\n");
    printf("Precision: %.5lf\n", 100.0 * (cm_correct + cm_overlap) / C.Record.size());
    printf("Recall: %.5lf\n", 100.0 *(cm_correct + cm_overlap) / B.Record.size());
    printf("Average deviation: %.5lf\n", 1.0 * cm_error / (cm_correct + cm_overlap));
    printf("Deviation rate: %.5lf\n\n", 100.0 * cm_overlap / (cm_correct + cm_overlap));
    return 0;
}
*/