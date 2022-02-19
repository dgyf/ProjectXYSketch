//
// Created by yqliu on 4/22/21.
//

#ifndef STREAMCLASSIFIER_XY_FLOWRADAR_H
#define STREAMCLASSIFIER_XY_FLOWRADAR_H


#include "InsertableIblt.h"
#include "Pbsketch_distr_space.h"
#include "FlowRadar_without_packetcount.h"


template<uint64_t iblt_size_in_bytes>
//class XY_FlowRadar : public VirtualIBLT {
class XY_FlowRadar {

    FlowRadar_without_packetcount<iblt_size_in_bytes> iblt;
    Pbsketch_distr_space *xy_count;
//    Pbsketch_distr_space *xy_compared;
    int threshold;


public:
    XY_FlowRadar(int memory_xy, int range, int thre) {
        int rad[range];
        for (int i = 1; i < range; i++) {
            rad[i] = i;
        }
        xy_count = new Pbsketch_distr_space(memory_xy, range, rad);
//        xy_compared = new Pbsketch_distr_space(memory_xy, range, rad);
        threshold = thre;

    }


    void build(uint32_t *items, const int begin, int end) {
        for (int i = begin; i < end; ++i) {
            xy_count->insert(items[i], 1);
//            xy_compared->insert(items[i], 1);
            if (xy_count->query(items[i]) >= threshold) {
                iblt.insert(items[i]);
            }
        }
    }

    void insert(uint32_t key, int fre) {
        xy_count->insert(key, fre);
//        xy_compared->insert(key, fre);
        if (xy_count->query(key) >= threshold) {
            iblt.insert(key);
        }
    }


/*
    inline void build(uint32_t *items, int end) {
        for (int i = 0; i < end; ++i) {
            xy_count->insert(items[i], 1);
            if (xy_count->query(items[i]) >= threshold) {
                iblt.insert(items[i]);
            }
        }
    }
*/

    void dump(unordered_map<uint32_t, int> &result) {
        iblt.dump(result);

        for (auto itr: result) {
            int remain = xy_count->query(itr.first);
            result[itr.first] += remain;
        }
    }

    int query_xy(uint32_t key) {
        return xy_count->query(key);
    }



    int get_n() {
        return xy_count->get_n();
    }


    int test(){
        return iblt.test_complete();
    }

    int approximate_query(uint32_t key) {
        int val_sc = xy_count->query(key);
        int val_iblt = iblt.approximate_query(key);

        if (val_iblt == 0) {
            return val_sc;
        } else if (val_sc <= threshold) {
            return val_sc;
        } else {
            return -1;
        }
    }
};

#endif //STREAMCLASSIFIER_XY_FLOWRADAR_H
