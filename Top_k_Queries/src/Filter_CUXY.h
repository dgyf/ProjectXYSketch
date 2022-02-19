//
// Created by yqliu on 7/10/20.
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

class Filter_CUXY : public SPA {
private:
    int d = 0;
    vector<vector<int>> pbcount;     // XY-sketch
    int sum = 0;
    int filter_sum = 0;
    vector<int> bin_rand;            //decomposition
    int range = 0;                   // range of items: 2^{range}
    int *last_sum;
    vector<int> c_head;              // number of range in each lary
    int num = 0;                     //number of model in the last pbcpunt
    int filter_count = 0;
    int CU_w = 0;
    int *CU_counters;
    BOBHash32 *hash[10];
    int CU_d = 3;
    double thre = 0;
    int value_thre = 0;
    int count_high_number[4];
    vector<vector<int>> filter_bit;


public:   //  memory size,  range of items,  decomposition
    Filter_CUXY(int space_input1, int ran, int ord[], double sp_ratio, double ratio_thre, int value_threshold) {
        int space_input=space_input1-2*ran*4;
        int space = space_input * sp_ratio;
        filter_bit.resize(2);
        filter_bit.at(0).resize(ran, 0);
        filter_bit.at(1).resize(ran, 0);
        thre = ratio_thre;
        value_thre = value_threshold;
        CU_w = (space_input - space) / 4;
        CU_counters = new int[(space_input - space) / 4];
        for (int i = 0; i < 4; i++) {
            count_high_number[i] = 0;
        }
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


    virtual ~Filter_CUXY() {
    }


    void insert(uint32_t key, int f, int N_1, double fre_threshold) {
        filter_count += 1;
        filter_sum += f;
        for (int i = 0; i < range; i++) {
            filter_bit.at((key & (1 << i)) >> i).at(i) += 1;
        }
        int branch = 1;
        int cu_open = 1;
//        int xy_open=1;
        if (filter_count > N_1) {
            branch = 0;
            double filter_fre = filter_sum;
            for (int i = 0; i < range; i++) {
                filter_fre = filter_fre * filter_bit.at((key & (1 << i)) >> i).at(i) / filter_sum;
            }
            if (filter_fre >= fre_threshold *filter_sum) {
                cu_open = 1;
            } else {
                cu_open = 0;
            }
        }

//        if (branch == 1 || (branch == 0 && cu_open == 0)) {
        if (branch == 0 && cu_open == 0) {      // filst N1 items for split method
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
                last_sum[local] += f;
            } else {
                pbcount.at(d - 1).at(key) += f;
            }
        }
        if (branch == 1 || (branch == 0 && cu_open == 1)) {

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

    }


    int query_ratio_cm(uint32_t key, double fre_threshold) {

        double filter_fre = filter_sum;
        for (int i = 0; i < range; i++) {
            filter_fre = filter_fre * filter_bit.at((key & (1 << i)) >> i).at(i) / filter_sum;
        }

        int ret = 0;
        if (filter_fre < fre_threshold * filter_sum) {

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

            ret = ceil(outp);
        } else {

            ret = 1 << 30;
            for (int i = 0; i < CU_d; i++) {
                int tmp = CU_counters[(hash[i]->run((const char *) &key, 4)) % CU_w];
                ret = min(ret, tmp);
            }

        }

        return ret;
    }

    int query_value_cm(uint32_t key, double fre_threshold) {

        double filter_fre = filter_sum;
        for (int i = 0; i < range; i++) {
            filter_fre = filter_fre * filter_bit.at((key & (1 << i)) >> i).at(i) / filter_sum;
        }

        int ret = 0;
        if (filter_fre < fre_threshold) {

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

            ret = ceil(outp);
        } else {

            ret = 1 << 30;
            for (int i = 0; i < CU_d; i++) {
                int tmp = CU_counters[(hash[i]->run((const char *) &key, 4)) % CU_w];
                ret = min(ret, tmp);
            }

        }

        return ret;
    }

    int get_count() {
        return count_high_number[0];
    }


};


#endif //STREAMCLASSIFIER_FILTER_CUXY_H
