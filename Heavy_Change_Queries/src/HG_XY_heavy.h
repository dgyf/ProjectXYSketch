//
// Created by yqliu on 5/24/21.
//

#ifndef STREAMCLASSIFIER_HG_XY_HEAVY_H
#define STREAMCLASSIFIER_HG_XY_HEAVY_H

#include <math.h>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <vector>
#include "BOBHash32.h"
#include <bits/stdc++.h>
#include "Pbsketch_distr_space.h"

using namespace std;

class HG_XY_heavy {

private:
    struct Node_p {
        int id;
        int cnt;
    };

    vector<vector<Node_p>> heavy_g;
    int b_parameter = 1.08;
    const int column = 8;
    int row;
    BOBHash32 *hash;
    Pbsketch_distr_space *xy_count;

    pair<int, int> *auxiliary;
    int auxi_test_num = 0;
    int num_auxi = 0;
    int threshold = 0;


public:
    HG_XY_heavy(int me_all, int me_xy, int range, int thre) {
        threshold = thre;
        int memory = 0.1 * me_all;
        num_auxi = memory / 8;
        auxiliary = new pair<int, int>[num_auxi];
        for (int i = 0; i < num_auxi; i++) {
            auxiliary[i].first = -1;
            auxiliary[i].second = 0;
        }
        row = (me_all - memory) / (4 * 2) / column;
        heavy_g.resize(row);
        for (int i = 0; i < row; i++) {
            heavy_g.at(i).resize(column);
            for (int j = 0; j < column; j++) {
                heavy_g.at(i).at(j).id = -1;
                heavy_g.at(i).at(j).cnt = 0;
            }
        }
        hash = new BOBHash32(column + 745);

        int rad[range];
        for (int i = 0; i < range; i++) {
            rad[i] = i;
        }
        xy_count = new Pbsketch_distr_space(me_xy, range, rad);


    }

    void insert(int key, int fre) {
        int temp = hash->run((char *) &key, 4) % row;
        bool flag = false;
        int local_id = 0;
        int local_num = 10000000;
        for (int i = 0; i < column; i++) {
            if (heavy_g.at(temp).at(i).id == key) {
                heavy_g.at(temp).at(i).cnt += fre;
                flag = true;
                if (heavy_g.at(temp).at(i).cnt < local_num) {
                    local_num = heavy_g.at(temp).at(i).cnt;
                    local_id = i;
                }
                break;
            } else if (heavy_g.at(temp).at(i).id == -1) {
                heavy_g.at(temp).at(i).id = key;
                heavy_g.at(temp).at(i).cnt = fre;
                flag = true;
                if (heavy_g.at(temp).at(i).cnt < local_num) {
                    local_num = heavy_g.at(temp).at(i).cnt;
                    local_id = i;
                }
                break;
            }
        }

        if (!flag) {
            for (int ij = 0; ij <= auxi_test_num; ij++) {
                if (auxiliary[ij].first == key) {
                    auxiliary[ij].second += fre;
                    flag = true;
                    break;
                }
            }
        }


        if (!flag ) {
            xy_count->insert(key, fre);
            int fre_xy = xy_count->query(key);
            if (fre_xy > local_num) {
                if (local_num >= threshold) {
                    if (auxi_test_num < num_auxi) {
                        auxiliary[auxi_test_num].first = key;
                        auxiliary[auxi_test_num].second = fre_xy;
                        auxi_test_num += 1;
                        xy_count->insert(key, -fre_xy);
                    } else {
                        printf("The auxiliary is full! \n");
                    }

                } else {
                    int key_hg = heavy_g.at(temp).at(local_id).id;
                    heavy_g.at(temp).at(local_id).id = key;
                    heavy_g.at(temp).at(local_id).cnt = fre_xy;
                    xy_count->insert(key, -fre_xy);
                    xy_count->insert(key_hg, local_num);
                }
            }
        }
    }


    int query_xy(int key) {
        return xy_count->query(key);
    }


    void query_find(unordered_map<uint32_t, int> &result) {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                if (heavy_g.at(i).at(j).cnt >= threshold) {
                    result[heavy_g.at(i).at(j).id] = heavy_g.at(i).at(j).cnt;
                }
            }
        }
        for (int i = 0; i < auxi_test_num; i++) {
            result[auxiliary[i].first] = auxiliary[i].second;
        }

    }


    ~HG_XY_heavy() {}
};

#endif //STREAMCLASSIFIER_HG_XY_HEAVY_H
