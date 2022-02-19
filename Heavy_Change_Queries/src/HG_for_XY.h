//
// Created by yqliu on 5/24/21.
//

#ifndef STREAMCLASSIFIER_HG_FOR_XY_H
#define STREAMCLASSIFIER_HG_FOR_XY_H


#include <math.h>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <vector>
#include "BOBHash32.h"
#include <bits/stdc++.h>
#include "Pbsketch_distr_space.h"

using namespace std;

class HG_for_XY {

private:
    struct Node_p {
        int id;
        int cnt;
    };

    vector<vector<Node_p>> heavy_g;
    const int column = 8;
    int row;
    BOBHash32 *hash;

    pair<int, int> *auxiliary;
    int auxi_test_num = 0;
    int num_auxi = 0;


public:
    HG_for_XY(int me_all) {
        int memory = 0.1 * me_all;   // the ratio of Ac
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

        if (!flag) {
            for (int ij = 0; ij < auxi_test_num; ij++) {
                if (auxiliary[ij].first == key) {
                    auxiliary[ij].second += fre;
                    flag = true;
                    break;
                }
            }
        }

        return flag;
    }

    void insert_auxi(int key, int fre) {
        auxiliary[auxi_test_num].first = key;
        auxiliary[auxi_test_num].second = fre;
        auxi_test_num += 1;
    }

    bool insert_hg(int key, int fre, int threshold, int &back_id, int &back_fre, bool &flag_out) {
        bool flag = false;
        flag_out=false;
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
            if (local_fre >= threshold) {
                insert_auxi(key, fre);
            } else {
                back_id = heavy_g.at(temp).at(local_id).id;
                back_fre = local_fre;
                heavy_g.at(temp).at(local_id).id = key;
                heavy_g.at(temp).at(local_id).cnt = fre;
                flag_out=true;
            }
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
        for (int i = 0; i < auxi_test_num; i++) {
            result[auxiliary[i].first] = auxiliary[i].second;
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


    int query_point_test(int key) {
        int temp = hash->run((char *) &key, 4) % row;
        for (int i = 0; i < column; i++) {
            if (heavy_g.at(temp).at(i).id == key)
                return heavy_g.at(temp).at(i).cnt;
        }
        for (int i = 0; i < auxi_test_num; i++) {
            if (auxiliary[i].first == key)
                return auxiliary[i].second;
//            cout<<"Finding in auxiliary!"<<endl;
        }
        return 0;
    }

    int get_num_auxi(){
        for(int i=0;i<auxi_test_num;i++){
            cout<<" auxiliary: No: "<<i<<" : "<<auxiliary[i].first<<" fre: "<<auxiliary[i].second<<endl;
        }
        cout<<" The total number of auxiliary is : "<<auxi_test_num<<endl;
        return auxi_test_num;
    }


    ~HG_for_XY() {}
};


#endif //STREAMCLASSIFIER_HG_FOR_XY_H
