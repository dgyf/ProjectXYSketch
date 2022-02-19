//
// Created by yqliu on 4/29/21.
//

#ifndef STREAMCLASSIFIER_HEAVY_GUARDIAN_FLOWRADAR_H
#define STREAMCLASSIFIER_HEAVY_GUARDIAN_FLOWRADAR_H


#include "InsertableIblt.h"
#include "FlowRadar_without_packetcount.h"
#include "Heavy_guardian_top.h"

template<uint64_t iblt_size_in_bytes>
//class XY_FlowRadar : public VirtualIBLT {
class Heavy_guardian_flowdadar : public VirtualIBLT {

    FlowRadar_without_packetcount<iblt_size_in_bytes> iblt;
    Heavy_guardian_tok *heavy_g;

    //    Pbsketch_distr_space *xy_compared;
    int threshold;


public:
    Heavy_guardian_flowdadar(int memory_xy, int thre) {

        heavy_g=new Heavy_guardian_tok(memory_xy);
        threshold = thre;

    }

    void build(uint32_t *items, const int num) {
        for (int i = 0; i < num; ++i) {
            heavy_g->insert(items[i], 1);
//            xy_compared->insert(items[i], 1);
            if (heavy_g->query(items[i]) >= threshold) {
                iblt.insert(items[i]);
            }
        }
    }

    void insert(uint32_t key, int fre) {
        heavy_g->insert(key, fre);
//        xy_compared->insert(key, fre);
        if (heavy_g->query(key) >= threshold) {
            iblt.insert(key);
        }
    }



    void dump(unordered_map<uint32_t, int> &result) {
        iblt.dump(result);
        for (auto itr: result) {
            int remain = heavy_g->query(itr.first);
            result[itr.first] += remain;
        }
    }




    int test(){
        return iblt.test_complete();
    }

    int query_heavy_g(uint32_t key) {
        return heavy_g->query(key);
    }


    int approximate_query(uint32_t key) {
        int val_sc = heavy_g->query(key);
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

#endif //STREAMCLASSIFIER_HEAVY_GUARDIAN_FLOWRADAR_H
