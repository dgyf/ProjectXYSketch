//
// Created by yqliu on 6/16/21.
//

#ifndef STREAMCLASSIFIER_HEAVY_GUARDIAN_LIST_H
#define STREAMCLASSIFIER_HEAVY_GUARDIAN_LIST_H

#include "InsertableIblt.h"
#include "FlowRadar_without_packetcount.h"
#include "Heavy_guardian_top.h"

//template<uint64_t iblt_size_in_bytes>
//class XY_FlowRadar : public VirtualIBLT {
class Heavy_guardian_list : public VirtualIBLT {

//    FlowRadar_without_packetcount<iblt_size_in_bytes> iblt;
    int *id_ary;
    int id_num = 0;
    int id_fixed_num = 0;
    Heavy_guardian_tok *heavy_g;

    //    Pbsketch_distr_space *xy_compared;
    int threshold;


public:
    Heavy_guardian_list(int memory_xy, int memory_ary, int thre) {
        id_fixed_num = memory_ary / 4;
        id_ary = new int[id_fixed_num];
        for (int i = 0; i < id_fixed_num; i++) {
            id_ary[i] = 0;
        }
        heavy_g = new Heavy_guardian_tok(memory_xy);
        threshold = thre;

    }

    void build(uint32_t *items, const int num) {
        for (int i = 0; i < num; ++i) {
            heavy_g->insert(items[i], 1);
//            xy_compared->insert(items[i], 1);
            if (heavy_g->query(items[i]) == threshold) {
                if (id_num == id_fixed_num) {
                    cout << "HG overflows !!!" << endl;
                } else {
                    id_ary[id_num] = items[i];
                    id_num++;
                }
//                cout<<" N0: (i) is : "<<i<<" : id_num : "<<id_num<<" id_fixed : "<<id_fixed_num<<endl;
            }
        }
    }

/*
    void insert(uint32_t key, int fre) {
        heavy_g->insert(key, fre);
//        xy_compared->insert(key, fre);
        if (heavy_g->query(key) == threshold) {
            if (id_num == id_fixed_num) {
                cout << "HG overflows !!!" << endl;

            } else {
                id_ary[id_num] = key;
                id_num++;
            }
        }
    }
*/

    void dump(unordered_map<uint32_t, int> &result) {
        for (int i = 0; i < id_num; i++) {
            result[id_ary[i]] = heavy_g->query(id_ary[i]);
        }
    }

    int test() {
        if (id_num == id_fixed_num) {
            return -1;
        } else {
            return 0;
        }
    }

    int query_heavy_g(uint32_t key) {
        return heavy_g->query(key);
    }


    int approximate_query(uint32_t key) {
        int val_sc = heavy_g->query(key);
        return val_sc;
    }
};


#endif //STREAMCLASSIFIER_HEAVY_GUARDIAN_LIST_H
