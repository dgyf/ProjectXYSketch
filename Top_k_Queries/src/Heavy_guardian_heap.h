//
// Created by yqliu on 4/29/21.
//

#ifndef STREAMCLASSIFIER_HEAVY_GUARDIAN_HEAP_H
#define STREAMCLASSIFIER_HEAVY_GUARDIAN_HEAP_H


#include <unordered_map>
#include <algorithm>
#include "BOBHash32.h"
#include "Space_allocation.h"
#include "Heavy_guardian_top.h"
#include <bits/stdc++.h>

using std::min;
using std::swap;

//template<int capacity, int cm_memory_in_bytes, int d>
class Heavy_guardian_heap {
    int capacity;
//    int d;
    int memory_in_bytes;
//    int w;

    typedef pair<uint32_t, uint32_t> KV;
    KV *heap;

//    KV heap[capacity];
    int heap_element_num;

    unordered_map<uint32_t, uint32_t> ht;

    Heavy_guardian_tok *heavy_g;


    // heap
    void heap_adjust_down(int i) {
        while (i < heap_element_num / 2) {
            int l_child = 2 * i + 1;
            int r_child = 2 * i + 2;
            int larger_one = i;
            if (l_child < heap_element_num && heap[l_child] < heap[larger_one]) {
                larger_one = l_child;
            }
            if (r_child < heap_element_num && heap[r_child] < heap[larger_one]) {
                larger_one = r_child;
            }
            if (larger_one != i) {
                swap(heap[i], heap[larger_one]);
                swap(ht[heap[i].second], ht[heap[larger_one].second]);
                heap_adjust_down(larger_one);
            } else {
                break;
            }
        }
    }

    void heap_adjust_up(int i) {
        while (i > 1) {
            int parent = (i - 1) / 2;
            if (heap[parent] <= heap[i]) {
                break;
            }
            swap(heap[i], heap[parent]);
            swap(ht[heap[i].second], ht[heap[parent].second]);
            i = parent;
        }
    }

public:
    Heavy_guardian_heap(int capacity1, int cm_memory_in_bytes1) : heap_element_num(0) {
        capacity = capacity1;
        memory_in_bytes = cm_memory_in_bytes1;
        heap = new KV[capacity1];
        for (int i = 0; i < capacity1; i++) {
            heap->first = 0;
            heap->second = 0;
        }

        heavy_g = new Heavy_guardian_tok(cm_memory_in_bytes1);
    }

    int get_row(){
        return heavy_g->get_row();
    }


    void insert(uint32_t key, int freq) {

        heavy_g->insert(key, freq);
        int est_fre = heavy_g->query(key);

        if (est_fre > heap[0].first) {
            if (ht.find(key) != ht.end()) {
                heap[ht[key]].first++;
                heap_adjust_down(ht[key]);
            } else if (heap_element_num < capacity) {
                heap[heap_element_num].second = key;
                heap[heap_element_num].first = freq;
                ht[key] = heap_element_num++;
                heap_adjust_up(heap_element_num - 1);
            } else {
                KV &kv = heap[0];
                ht.erase(kv.second);
                kv.second = key;
                kv.first = est_fre;
                ht[key] = 0;
                heap_adjust_down(0);
            }
        }
    }



    void get_top_k(uint16_t k, uint32_t *result) {
        KV *a = new KV[capacity];
        memcpy(a, heap, sizeof(heap));
        sort(a, a + capacity);
        int i;
        for (i = 0; i < k && i < capacity; ++i) {
            result[i] = a[capacity - 1 - i].second;
        }
        for (; i < k; ++i) {
            result[i] = 0;
        }
    }


/*
void get_top_k_with_frequency(uint16_t k, vector<pair<uint32_t, uint32_t>> &result) {
    KV *a = new KV[capacity];
    memcpy(a, heap, sizeof(heap));
    sort(a, a + capacity);
    int i;
    for (i = 0; i < k && i < capacity; ++i) {
        result[i].first = a[capacity - 1 - i].second;
        result[i].second = a[capacity - 1 - i].first;
    }
    for (; i < k; ++i) {
        result[i].first = 0;
        result[i].second = 0;
    }
}

 */


    vector<pair<uint32_t, uint32_t> > func(int k, uint32_t *a, uint32_t *b, int len) {

        struct cmp {
            bool operator()(pair<uint32_t, uint32_t> aa, pair<uint32_t, uint32_t> bb) {
                return aa.second > bb.second;
            }
        };
        priority_queue <pair<uint32_t, uint32_t>, vector<pair<uint32_t, uint32_t> >, cmp> que;
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

    void get_top_k_with_frequency(uint16_t k, vector<pair<uint32_t, uint32_t>> &result) {
        uint32_t *a = new uint32_t[capacity];
        uint32_t *b = new uint32_t[capacity];
        for (int i = 0; i < capacity; i++) {
            a[i] = heap[i].first;
            b[i] = heap[i].second;
        }

        result = func(k, b, a, capacity);
    }

    ~Heavy_guardian_heap() {

    }

};

#endif //STREAMCLASSIFIER_HEAVY_GUARDIAN_HEAP_H
