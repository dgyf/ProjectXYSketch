//
// Created by yqliu on 4/28/21.
//

#ifndef STREAMCLASSIFIER_HEAVY_GUARDIAN_TOP_H
#define STREAMCLASSIFIER_HEAVY_GUARDIAN_TOP_H

#include <math.h>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <vector>
#include "BOBHash32.h"
#include <bits/stdc++.h>

using namespace std;

class Heavy_guardian_tok {

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


public:
    Heavy_guardian_tok(int memory) {

        row = memory / (4 * 2) / column;
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

    void insert(int key, int fre) {
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
                }
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

    int get_row(){
        return row;
    }

    void get_top_k_with_frequency(int k_value, vector<pair<uint32_t, uint32_t>> &result) {
        uint32_t *a = new uint32_t[int(row * column)];
        uint32_t *b = new uint32_t[int(row * column)];
        int temp_cnt=0;
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                a[temp_cnt] = heavy_g.at(i).at(j).cnt;         // fre
                b[temp_cnt] = heavy_g.at(i).at(j).id;        // key
                temp_cnt++;
            }
        }
        result = func(k_value, b, a, int(row * column));

    }

    void hg_print() {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                printf(" id: %d   fre: %d  \n", heavy_g.at(i).at(j).id, heavy_g.at(i).at(j).cnt);
            }
        }
    }

    vector<pair<uint32_t, uint32_t> > func(int k, uint32_t *a, uint32_t *b, int len) {

        struct cmp {
            bool operator()(pair<uint32_t, uint32_t> aa, pair<uint32_t, uint32_t> bb) {
                return aa.second > bb.second;
            }
        };
        priority_queue<pair<uint32_t, uint32_t>, vector<pair<uint32_t, uint32_t> >, cmp> que;
        for (int i = 0; i < len; i++) {
            if (que.size() < k) {
                que.push(pair<uint32_t, uint32_t>{a[i], b[i]});
                continue;
            }
            if (que.top().second < b[i]) {
                que.pop();
                que.push(pair<uint32_t, uint32_t>{a[i], b[i]});
            }
        }
        vector<pair<uint32_t, uint32_t> > ret;
        while (que.size()) {
            ret.push_back(que.top());
            que.pop();
        }
        reverse(ret.begin(), ret.end());
        return ret;
    }

    ~Heavy_guardian_tok() {}
};

#endif //STREAMCLASSIFIER_HEAVY_GUARDIAN_TOP_H
