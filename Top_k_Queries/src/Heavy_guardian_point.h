//
// Created by yqliu on 4/28/21.
//

#ifndef STREAMCLASSIFIER_HEAVY_GUARDIAN_POINT_H
#define STREAMCLASSIFIER_HEAVY_GUARDIAN_POINT_H

#include <math.h>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <vector>
#include "BOBHash32.h"

using namespace std;

class Heavy_guardian_point {

private:
    struct Node_p {
        int id;
        int cnt;
    };

    vector<vector<Node_p>> heavy_g;
    vector<vector<int>> heavy_cold;
    int b_parameter = 1.08;
    const int column = 8;
    const int column_cold = 16;
    int row;
    BOBHash32 *hash;

public:
    Heavy_guardian_point(int memory) {

        row = memory / 2 / (4 * 2) / column;
        heavy_g.resize(row);
        heavy_cold.resize(row);
        for (int i = 0; i < row; i++) {
            heavy_g.at(i).resize(column);
            heavy_cold.at(i).resize(column_cold, 0);
            for (int j = 0; j < column; j++) {
                heavy_g.at(i).at(j).id = -1;
                heavy_g.at(i).at(j).cnt = 0;
            }

        }

        hash = new BOBHash32(column + 745);

    }

    void insert(int key, int fre) {
        int temp = hash->run((char *) &key, 4) % row;
        bool flag = false;
        bool flag_cold = false;
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
            int local, MIN = 1000000000;
            for (int i = 0; i < column; i++) {         // find the smallest count
                int c = heavy_g.at(temp).at(i).cnt;
                if (c < MIN) {
                    MIN = c;
                    local = i;
                }
            }
            if (!(rand() % int(pow(b_parameter, MIN)))) {  // probability reduce
                heavy_g.at(temp).at(local).cnt--;
                if (heavy_g.at(temp).at(local).cnt <= 0) {
                    heavy_g.at(temp).at(local).id = key;
                    heavy_g.at(temp).at(local).cnt = fre;
                    flag_cold = true;
                }
            }
            if (!flag_cold) {
                int temp2 = hash->run((char *) &key, 4) % column_cold;
                heavy_cold.at(temp).at(temp2) += fre;
            }
        }
    }

    int query(int key) {
        int temp = hash->run((char *) &key, 4) % row;
        for (int i = 0; i < column; i++) {
            if (heavy_g.at(temp).at(i).id == key)
                return heavy_g.at(temp).at(i).cnt;
        }
        return 0;
    }

    ~Heavy_guardian_point() {}
};

#endif //STREAMCLASSIFIER_HEAVY_GUARDIAN_POINT_H
