//
// Created by yqliu on 8/18/21.
//

#ifndef STREAMCLASSIFIER_HG_HEAVY_LIGHT_H
#define STREAMCLASSIFIER_HG_HEAVY_LIGHT_H


#include <math.h>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <vector>
#include "BOBHash32.h"
#include <bits/stdc++.h>
#include "Pbsketch_distr_space.h"

using namespace std;

class HG_heavy_light {

private:
    struct Node_p {
        int id;
        int cnt;
    };

    vector<vector<Node_p>> heavy_g;
    vector<vector<uint16_t>> heavy_light;
    const int column = 8;
    int row;
    BOBHash32 *hash;
    int light_column = 16;
    double b = 1.08;
    double threshold=0;

//    pair<int, int> *auxiliary;
//    int auxi_test_num = 0;
//    int num_auxi = 0;


public:
    HG_heavy_light(int me_all, double thre) {
        threshold=thre;
        row = me_all / (4 * 2+2) / (column+light_column);
        heavy_g.resize(row);
        heavy_light.resize(row);
        for (int i = 0; i < row; i++) {
            heavy_g.at(i).resize(column);
            for (int j = 0; j < column; j++) {
                heavy_g.at(i).at(j).id = -1;
                heavy_g.at(i).at(j).cnt = 0;
            }
            heavy_light.at(i).resize(light_column);
            for (int j = 0; j < light_column; j++) {
                heavy_light.at(i).at(j) = 0;
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


    void insert_hg(int key, int fre) {
        bool test_flag = insert(key, fre);
        if (!test_flag) {

            int temp = hash->run((char *) &key, 4) % row;
            int local_id = 0, local_fre = 100000000;
            for (int i = 0; i < column; i++) {
                if (heavy_g.at(temp).at(i).cnt < local_fre) {
                    local_fre = heavy_g.at(temp).at(i).cnt;
                    local_id = i;
                }
            }
            srand(1);
            if (!(rand() % int(pow(b, fre)))) {
                heavy_g.at(temp).at(local_id).cnt--;
                if (heavy_g.at(temp).at(local_id).cnt == 0) {
                    heavy_g.at(temp).at(local_id).cnt = 1;
                    heavy_g.at(temp).at(local_id).id = key;
                } else {
                    int temp1 = hash->run((char *) &key, 4) % light_column;
                    heavy_light.at(temp).at(temp1)++;
                }
            } else {
                int temp1 = hash->run((char *) &key, 4) % light_column;
                heavy_light.at(temp).at(temp1)++;
            }
        }


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
        int temp1 = hash->run((char *) &key, 4) % light_column;
        return heavy_light.at(temp).at(temp1);
    }


    ~HG_heavy_light() {}
};


#endif //STREAMCLASSIFIER_HG_HEAVY_LIGHT_H
