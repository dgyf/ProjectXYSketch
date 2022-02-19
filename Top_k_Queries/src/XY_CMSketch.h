//
// Created by yqliu on 5/16/20.
//

#ifndef STREAMCLASSIFIER_XY_CMSKETCH_H
#define STREAMCLASSIFIER_XY_CMSKETCH_H

#include <math.h>
#include "SPA.h"
#include <cstring>
#include <algorithm>
#include <iostream>
#include <vector>
#include "Space_allocation.h"
#include "params.h"
#include "BOBHash32.h"

using namespace std;

class XY_CMSketch : public SPA {
private:
    int d = 0;
    vector<vector<int>> pbcount;     // XY-sketch
    int sum = 0;
    vector<int> bin_rand;            //decomposition
    int range = 0;                   // range of items: 2^{range}
    int *last_sum;
    vector<int> c_head;              // number of range in each lary
    int num = 0;                     //number of model in the last pbcpunt

    int CU_w = 0;
    int *CU_counters;
    BOBHash32 *hash[10];
    int CU_d = 3;
    double ratio_thre = 0;
    int value_thre = 0;
    int count_high_number;


public:   //  memory size,  range of items,  decomposition
//    XY_CMSketch(int space_input, int ran, int ord[], double sp_ratio, int value_threshold) {
        XY_CMSketch(int space_input, int ran, int ord[], double sp_value, int value_threshold) {
        int space = sp_value;
//        ratio_thre = ratio_threshold;
        value_thre = value_threshold;
        CU_w = (space_input - space) / 4;
        CU_counters = new int[(space_input - space) / 4];

        count_high_number = 0;

        for (int i = 0; i < (space_input - space) / 4; i++) {
            CU_counters[i] = 0;
        }
//        memset(counters, 0, sizeof(counters));
        for (int i = 0; i < CU_d; i++)
            hash[i] = new BOBHash32(i + 750);

        range = ran;
        Space_allocation itr;
        int *w_group = itr.get_size(space, ran);  // Space allocation function: the first one is d-value
        d = w_group[0];                           // number of elements : d-value
        c_head.resize(d, 0);
        for (int i = 1; i < d + 1; i++) {
            c_head.at(i - 1) = w_group[i];
        }
        pbcount.resize(d);
        int temp = space / 4;
        for (int i = 0; i < d - 1; i++) {
            pbcount.at(i).resize(pow(2, w_group[i + 1]), 0);
            temp = temp - pow(2, w_group[i + 1]);
        }
        num = temp / pow(2, w_group[d]);        // stratification
        pbcount.at(d - 1).resize(num * pow(2, w_group[d]), 0);  // last layer using the idea of stratification
        last_sum = new int[num];
        for (int i = 0; i < num; i++) {
            last_sum[i] = 0;
        }
        bin_rand.resize(ran);
        for (int i = 0; i < ran; i++) {
            bin_rand.at(i) = ord[i];
        }
    }


    virtual ~XY_CMSketch() {
    }


    void insert(uint32_t key, int f) {
        sum += f;
        if (d > 1) {
            int re_num = 0;
            int re_cl = 0;
            int point = 0;
            for (int i = 0; i < (range - c_head[d - 1]); i++) {
                point = point + (((key & (1 << bin_rand[i])) >> bin_rand[i]) << re_cl);
                re_cl += 1;
                if (re_cl == c_head[re_num]) {
                    pbcount.at(re_num).at(point) += f;
                    re_num += 1;
                    re_cl = 0;
                    if (re_num < d - 1) {
                        point = 0;
                    }
                }
            }
            int local = point % num;      // stratification
            point = 0;
            for (int i = range - c_head[d - 1]; i < range; i++) {
                point += ((key & (1 << bin_rand[i])) >> bin_rand[i]) << (i - range + c_head[d - 1]);
            }
            pbcount.at(d - 1).at(local * (1 << c_head[d - 1]) + point) += f;
            last_sum[local] += 1;
        } else {
            pbcount.at(d - 1).at(key) += f;
        }


        int index[CU_d];
        int value[CU_d];
        int min_val = 1 << 30;

        for (int i = 0; i < CU_d; i++) {
            index[i] = (hash[i]->run((const char *) &key, 4)) % CU_w;
            value[i] = CU_counters[index[i]];
            min_val = min(min_val, value[i]);
        }

        int temp = min_val + f;
        for (int i = 0; i < CU_d; i++) {
            CU_counters[index[i]] = max(CU_counters[index[i]], temp);
        }

    }


    int query_ratio_cm(uint32_t key) {
        double outp = sum;
        int re_num = 0;
        int re_cl = 0;
        int point = 0;
        if (d > 1) {
            for (int i = 0; i < range - c_head[d - 1]; i++) {

                point += ((key & (1 << bin_rand[i])) >> bin_rand[i]) << re_cl;
                re_cl += 1;
                if (re_cl == c_head[re_num]) {
                    re_cl = 0;
                    outp = outp * pbcount.at(re_num).at(point) / sum;
                    re_num += 1;
                    if (re_num < d - 1) {
                        point = 0;
                    }
                }
            }
            int local = point % num;     // stratification
            point = 0;
            for (int i = range - c_head[d - 1]; i < range; i++) {
                point += ((key & (1 << bin_rand[i])) >> bin_rand[i]) << (i - range + c_head[d - 1]);
            }
            outp = outp * pbcount.at(d - 1).at(local * (1 << c_head[d - 1]) + point) /
                   last_sum[local];
        } else {
            outp = pbcount.at(d - 1).at(key);
        }


//CU_Sketch
        int ret = 1 << 30;
        for (int i = 0; i < CU_d; i++) {
            int tmp = CU_counters[(hash[i]->run((const char *) &key, 4)) % CU_w];
            ret = min(ret, tmp);
        }

//judgement
        int final = 0;
        if (ret < ceil(outp)) {
            final = ret;
            count_high_number++;
        } else {
//            if ((double) (outp / ret) < thre) {
            if (ret < value_thre) {
                final = ceil(outp);

            } else {
                final = ret;               //  chose CM-sketch
//                final=ceil(outp);        //  chose XY-sketch
            }
        }

        return final;
    }



    int get_count() {
        return count_high_number;
    }


};


#endif //STREAMCLASSIFIER_XY_CMSKETCH_H
