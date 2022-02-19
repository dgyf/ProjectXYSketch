//
// Created by yqliu on 8/26/19.
//
#ifndef STREAMCLASSIFIER_PBSKETCH_DISTR_SPACE_H
#define STREAMCLASSIFIER_PBSKETCH_DISTR_SPACE_H

#include <math.h>
#include "SPA.h"
#include <cstring>
#include <algorithm>
#include <iostream>
#include <vector>
#include "Space_allocation.h"

using namespace std;

class Pbsketch_distr_space : public SPA {
private:
    int d = 0;
    vector<vector<int>> pbcount;     // XY-sketch
    int sum = 0;
    vector<int> bin_rand;            //decomposition
    int range = 0;                   // range of items: 2^{range}
    int *last_sum;
    vector<int> c_head;              // number of range in each lary
    int num = 0;                     //number of model in the last pbcpunt


public:   //  memory size,  range of items,  decomposition
    Pbsketch_distr_space(int space, int ran, int ord[]) {
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


    virtual ~Pbsketch_distr_space() {
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
            last_sum[local] += f;
        } else {
            pbcount.at(d - 1).at(key) += f;
        }
    }


    int min_query(uint32_t key) {
        int outp = sum;
        int outp_min = sum;
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
                    if (outp_min > pbcount.at(re_num).at(point)) {
                        outp_min = pbcount.at(re_num).at(point);
                    }
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
            int medi_xy = pbcount.at(d - 1).at(local * (1 << c_head[d - 1]) + point);
            if (outp_min > medi_xy) {
                outp_min = medi_xy;
            }

        } else {
            outp_min = pbcount.at(d - 1).at(key);
            outp = pbcount.at(d - 1).at(key);
        }
        if (outp_min < outp) {
            return outp_min;
        } else {

            return ceil(outp);
        }
    }


    int query(uint32_t key) {
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

        int final = ceil(outp);         // frequency should be integer
        return final;
    }


    int get_sum() {
        return sum;
    }


    int batch_query(uint32_t *data, int n) {
        int ret = 0;
        for (int i = 0; i < n; ++i) {
            ret += query(data[i]);
        }
        return ret;
    }

    void build(uint32_t *data, int n) {
        for (int i = 0; i < n; ++i) {
            insert(data[i], 1);
        }
    }


};


#endif //STREAMCLASSIFIER_PBSKETCH_DISTR_SPACE_H