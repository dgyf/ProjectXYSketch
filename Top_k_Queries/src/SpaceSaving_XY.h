//
// Created by yqliu on 10/29/20.
//

#ifndef STREAMCLASSIFIER_SPACESAVING_XY_H
#define STREAMCLASSIFIER_SPACESAVING_XY_H


#include <cstdint>
#include <cstdlib>
#include <unordered_map>
#include <vector>
#include "SpaceSavingUtils.h"
#include "SPA.h"
#include "Space_allocation.h"
#include "Pbsketch_distr_space.h"

using namespace std;

#define tail_node nodes[0].prev

//template<int capacity>
class SpaceSaving_XY : public SPA {
    int now_element;
//    uint16_t tail_idx;
//    Node_ss nodes[capacity + 1];
    Node_ss *nodes;
    int capacity;
    unordered_map<uint32_t, Node_ss *> hash_table;
    Pbsketch_distr_space *xy_counter;

    void append_new_key(uint32_t key, int freq) {
        if (now_element < capacity) {
            uint16_t idx = ++now_element; // we use 0 to represent header
            nodes[idx].key = key;
            nodes[idx].val = 0;
            hash_table[key] = nodes + idx;

            // append to tail
            nodes[idx].prev = tail_node;
            tail_node->next = nodes + idx;
            nodes[idx].next = nodes;
            nodes[idx].parent = nodes + idx;
            tail_node = &nodes[idx];
            add_counter(tail_node, freq);
        } else {
            xy_counter->insert((int) key, freq);
            int feed_b = xy_counter->query((int) key);
            if (feed_b > tail_node->val) {
                uint32_t old_key = tail_node->key;
                int gap_exchange = tail_node->val - xy_counter->query((int) old_key);
                hash_table.erase(old_key);
                tail_node->key = key;
                hash_table[key] = tail_node;
                add_counter(tail_node, freq);
                if (gap_exchange > 0) {
                    xy_counter->insert((int) old_key, gap_exchange);
                }
            }
        }
    }

    void add_counter(Node_ss *my, int freq) {
        //
        if (my->parent == my && my->next->val == my->val) {
//            std::swap(my->key, my->next->key);
//            std::swap(hash_table[my->key], hash_table[my->next->key]);
//            my = my->next;
            Node_ss *p = my->next, *nt = my->next;
            while (p && p->val == my->val) {
                p->parent = nt;
                p = p->next;
            }
        }

        my->val += freq;
        int now_freq = my->val;
        Node_ss *prev_node = my->prev;

        if (prev_node->val > now_freq) {
            return;
        }

        Node_ss *next_node = my->next;

        // make next and prev connect
        prev_node->next = my->next;
        next_node->prev = my->prev;

        while (prev_node->val < now_freq) {
            prev_node = prev_node->parent->prev;
        }

        next_node = prev_node->next;

        my->next = prev_node->next;
        prev_node->next = my;

        my->prev = next_node->prev;
        next_node->prev = my;

        my->parent = (prev_node->val == my->val) ? prev_node->parent : my;
    }

public:
    SpaceSaving_XY(int k_size, int me_size, int num_bit) : now_element(0) {
        capacity = k_size;
        nodes = new Node_ss[k_size + 1];
        memset(nodes, 0, sizeof(nodes));
        now_element = 0;
        nodes[0].val = -1;
        nodes[0].parent = nodes;
        tail_node = nodes;
        hash_table.reserve(100 * capacity);
        int spilt[num_bit];
        for (int i = 0; i < num_bit; i++) {        // simple decomposition
            spilt[i] = i;
        }
        xy_counter = new Pbsketch_distr_space(me_size, num_bit, spilt);
    }

    void insert(uint32_t key, int freq) {
//        Node *& my_node = hash_table[key];
        auto itr = hash_table.find(key);
        if (itr == hash_table.end()) {
            // key not found
            append_new_key(key, freq);
        } else {
            // key found
            add_counter(itr->second, freq);
        }
    }


    int query_fre(int key) {
        int outp = 0;

        auto itr = hash_table.find(key);
        if (itr == hash_table.end()) {
            // key not found
            outp = xy_counter->query(key);
        } else {
            // key found
            outp = itr->second->val;
        }
        return outp;
    }

    void get_top_k(uint16_t k, uint32_t *result) {
        Node_ss *idx = nodes[0].next;

        int i;
        for (i = 0; i < k && i < capacity && i < now_element; ++i) {
            result[i] = idx->key;
            idx = idx->next;
        }

        for (; i < k; ++i) {
            result[i] = 0;
        }

        return;
    }

    void get_top_k_with_frequency(uint16_t k, vector<pair<uint32_t, uint32_t>> &result) {
        Node_ss *idx = nodes[0].next;

        int i;
        for (i = 0; i < k && i < capacity && i < now_element; ++i) {
            result[i].first = idx->key;
            result[i].second = idx->val;
            idx = idx->next;
        }

        for (; i < k; ++i) {
            result[i].first = 0;
            result[i].second = 0;
        }

        return;
    }

    inline void build(uint32_t *items, int n) {
        for (int i = 0; i < n; ++i) {
            insert(items[i], 1);
        }
    }
};


#endif //STREAMCLASSIFIER_SPACESAVING_XY_H
