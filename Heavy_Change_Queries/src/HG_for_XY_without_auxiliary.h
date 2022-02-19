//
// Created by yqliu on 8/17/21.
//

#ifndef STREAMCLASSIFIER_HG_FOR_XY_WITHOUT_AUXILIARY_H
#define STREAMCLASSIFIER_HG_FOR_XY_WITHOUT_AUXILIARY_H


#include <math.h>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <vector>
#include "BOBHash32.h"
#include <bits/stdc++.h>
#include "Pbsketch_distr_space.h"

using namespace std;

class HG_for_XY_without_auxiliary {

private:
    struct Node_p {
        int id;
        int cnt;
    };

    vector<vector<Node_p>> heavy_g;
    const int column = 8;
    int row;
    BOBHash32 *hash;

//    pair<int, int> *auxiliary;
//    int auxi_test_num = 0;
//    int num_auxi = 0;


public:
    HG_for_XY_without_auxiliary(int me_all) {

        row = me_all  / (4 * 2) / column;
        heavy_g.resize(row);
        for (int i = 0; i < row; i++) {
            heavy_g.at(i).resize(column);
            for (int j = 0; j < column; j++) {
                heavy_g.at(i).at(j).id = -1;
                heavy_g.at(i).at(j).cnt = 0;
            }
        }
        hash = new BOBHash32(column + 745);

    }

    bool insert(int key, int fre) {
        int temp = hash->run((char *) &key, 4) % row;
        bool flag = false;
        for (int i = 0; i < column; i++) {
            if (heavy_g.at(temp).at(i).id == key) {
                heavy_g.at(temp).at(i).cnt += fre;
                flag = true;
                break;
            } else if (heavy_g.at(temp).at(i).id == -1) {
                heavy_g.at(temp).at(i).id = key;
                heavy_g.at(temp).at(i).cnt = fre;
                flag = true;
                break;
            }
        }

        return flag;
    }


    bool insert_hg(int key, int fre, int threshold, int &back_id, int &back_fre) {
        bool flag = false;
        int temp = hash->run((char *) &key, 4) % row;
        int local_id = 0, local_fre = 100000000;
        for (int i = 0; i < column; i++) {
            if (heavy_g.at(temp).at(i).cnt < local_fre) {
                local_fre = heavy_g.at(temp).at(i).cnt;
                local_id = i;
            }
        }
        if (fre > local_fre) {
            flag = true;
                back_id = heavy_g.at(temp).at(local_id).id;
                back_fre = local_fre;
                heavy_g.at(temp).at(local_id).id = key;
                heavy_g.at(temp).at(local_id).cnt = fre;
        }
        return flag;
    }

    void query_find(unordered_map<uint32_t, int> &result, int threshold) {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                if (heavy_g.at(i).at(j).cnt >= threshold) {
                    result[heavy_g.at(i).at(j).id] = heavy_g.at(i).at(j).cnt;
                }
            }
        }

    }

    int query_point(int key) {
        int temp = hash->run((char *) &key, 4) % row;
        for (int i = 0; i < column; i++) {
            if (heavy_g.at(temp).at(i).id == key)
                return heavy_g.at(temp).at(i).cnt;
        }
        return 0;
    }



    ~HG_for_XY_without_auxiliary() {}
};


#endif //STREAMCLASSIFIER_HG_FOR_XY_WITHOUT_AUXILIARY_H
