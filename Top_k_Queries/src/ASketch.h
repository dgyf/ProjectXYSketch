#ifndef _ASKETCH_H
#define _ASKETCH_H

#include "params.h"
#include "BOBHash32.h"
#include "SPA.h"
#include <cstring>
#include <x86intrin.h>
#include <bmiintrin.h>
#include <iostream>
#include<bits/stdc++.h>
using namespace std;

class ASketch
{
private:
    int w = 0;
    int bucket_num = 0;

    int *new_count;
    int *old_count;
    int *items;

    int cur_pos=0;
    int d=0;
    int **counter;
    int filter_size;
    BOBHash32 *bobhash[10];

public:
    ASketch(int tot_memory_in_bytes, int dp, int filter_size1)
    {
        filter_size=filter_size1;
        w=(tot_memory_in_bytes - filter_size1 * 12) * 8 / 32 / dp;
        d=dp;
        bucket_num = filter_size1 / 16;
        new_count=new int[filter_size1];
        old_count=new int[filter_size1];
        items=new int[filter_size1];
        counter=(int **)new int*[d];

        //memset(counter, 0, sizeof(counter));
        for(int i=0;i<dp;i++){
            counter[i]=new int[w];
            for (int j=0;j<w;j++){
                counter[i][j]=0;
            }
        }
        for(int i=0;i<filter_size1;i++){
            items[i]=0;
            new_count[i]=0;
            old_count[i]=0;
        }
        cur_pos = 0;

        for (int i = 0; i < dp; i++)
        {
            bobhash[i] = new BOBHash32(i + 1000);
        }
    }

    int * get_items()
    {
        return items;
    }

    int * get_freq()
    {
        return new_count;
    }

    void insert(uint32_t key, int f=1)
    {
        const __m128i item = _mm_set1_epi32((int)key);

        for (int i = 0; i < bucket_num; i++)
        {
            __m128i *keys_p = (__m128i *)(items + (i << 4));

            __m128i a_comp = _mm_cmpeq_epi32(item, keys_p[0]);
            __m128i b_comp = _mm_cmpeq_epi32(item, keys_p[1]);
            __m128i c_comp = _mm_cmpeq_epi32(item, keys_p[2]);
            __m128i d_comp = _mm_cmpeq_epi32(item, keys_p[3]);

            a_comp = _mm_packs_epi32(a_comp, b_comp);
            c_comp = _mm_packs_epi32(c_comp, d_comp);
            a_comp = _mm_packs_epi32(a_comp, c_comp);

            int matched = _mm_movemask_epi8(a_comp);

            if(matched != 0)
            {
                //return 32 if input is zero;
                int matched_index = _tzcnt_u32((uint32_t)matched) + (i << 4);
                new_count[matched_index] += f;
                return;
            }
        }

        if(cur_pos != filter_size)
        {
            items[cur_pos] = key;
            new_count[cur_pos] = f;
            old_count[cur_pos] = 0;
            cur_pos ++;
            return;
        }


        int estimate_value, min_index, min_value, temp;

        //int index[MAX_HASH_NUM];
        int index[d];
        estimate_value = (1 << 30);
        for(int i = 0; i < d; i++)
        {
            index[i] = (bobhash[i]->run((const char *)&key, 4)) % w;

            // counter[i][hash_value] ++;
            temp = counter[i][index[i]];

            estimate_value = estimate_value < temp ? estimate_value : temp;
        }
        estimate_value += f;
        for(int i = 0; i < d; i++)
        {
            if(counter[i][index[i]] < estimate_value)
            {
                counter[i][index[i]] = 	estimate_value;
            }
        }


        min_index = 0;
        min_value = (1 << 30);
        for(int i = 0; i < filter_size; i++)
        {
            if(items[i] != (uint32_t)(-1) && min_value > new_count[i])
            {
                min_value = new_count[i];
                min_index = i;
            }
        }
        if(estimate_value > min_value)
        {
            temp = new_count[min_index] - old_count[min_index];
            if(temp > 0)
            {
                min_value = (1 << 30);
                for(int i = 0; i < d; i++)
                {
                    index[i] = (bobhash[i]->run((const char *)&items[min_index], 4)) % w;
                    min_value = min_value < counter[i][index[i]] ? min_value : counter[i][index[i]];
                }
                for(int i = 0; i < d; i++)
                {
                    if(counter[i][index[i]] < min_value + temp)
                        counter[i][index[i]] = min_value + temp;
                }
            }
            items[min_index] = key;
            new_count[min_index] = estimate_value;
            old_count[min_index] = estimate_value;
        }
    }

    void build(uint32_t * data, int n) {
        for (int i = 0; i < n; ++i) {
            insert(data[i], 1);
        }
    }

    int query(uint32_t key)
    {
        const __m128i item = _mm_set1_epi32((int)key);
        for (int i = 0; i < bucket_num; i++)
        {
            __m128i *keys_p = (__m128i *)(items + (i << 4));

            __m128i a_comp = _mm_cmpeq_epi32(item, keys_p[0]);
            __m128i b_comp = _mm_cmpeq_epi32(item, keys_p[1]);
            __m128i c_comp = _mm_cmpeq_epi32(item, keys_p[2]);
            __m128i d_comp = _mm_cmpeq_epi32(item, keys_p[3]);

            a_comp = _mm_packs_epi32(a_comp, b_comp);
            c_comp = _mm_packs_epi32(c_comp, d_comp);
            a_comp = _mm_packs_epi32(a_comp, c_comp);

            int matched = _mm_movemask_epi8(a_comp);

            if(matched != 0)
            {
                //return 32 if input is zero;
                int matched_index = _tzcnt_u32((uint32_t)matched) + (i << 4);
                return new_count[matched_index];
            }
        }


        int hash_value, temp;
        int estimate_value = (1 << 30);
        for(int i = 0; i < d; i++)
        {
            hash_value = (bobhash[i]->run((const char *)&key, 4)) % w;
            temp = counter[i][hash_value];
            estimate_value = estimate_value < temp ? estimate_value : temp;
        }
        return estimate_value;
    }

    int batch_query(uint32_t * data, int n)
    {
        int ret = 0;
        for (int i = 0; i < n; ++i) {
            ret += query(data[i]);
        }

        return ret;
    }

    ~ASketch()
    {
        for(int i = 0; i < d; i++)
        {
            delete bobhash[i];
        }
    }

    vector<pair<uint32_t, uint32_t> > func(int k, int *a, int *b, int len) {

        struct cmp {
            bool operator()(pair<int, int> aa, pair<int, int> bb) {
                return aa.second > bb.second;
            }
        };
        priority_queue<pair<int, int>, vector<pair<int, int> >, cmp> que;
        for (int i = 0; i < len; i++) {
            if (que.size() < k) {
                que.push(pair<int, int>{a[i], b[i]});
                continue;
            }
            if (que.top().second < b[i]) {
                que.pop();
                que.push(pair<int, int>{a[i], b[i]});
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

    void get_top_k_with_frequency(int k, vector<pair<uint32_t, uint32_t>> &result) {

//        for(int i=0;i<filter_size;i++){
//            cout<<" ID : "<<items[i]<<"  Frequency : "<<new_count<<endl;
//        }

        result = func(k, items,new_count, filter_size);

//        for(int i=0;i<k;i++){
//            cout<<" ID : "<<result[i].first<<"  Frequency : "<<result[i].second<<endl;
//        }
    }

};

#endif//_ASKETCH_H