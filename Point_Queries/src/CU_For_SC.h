//
// Created by yqliu on 10/20/19.
//

#ifndef STREAMCLASSIFIER_CU_FOR_SC_H
#define STREAMCLASSIFIER_CU_FOR_SC_H
#include "params.h"
#include "BOBHash32.h"
#include "SPA.h"
#include <cstring>
#include <algorithm>

using namespace std;

template <int memory_in_bytes, int d>
class CU_For_SC: public SPA
{
private:
    static constexpr int w = memory_in_bytes * 8 / 32;
    int counters[w];
    BOBHash32 * hash[d];
public:
    CU_For_SC()
    {
        memset(counters, 0, sizeof(counters));
        for (int i = 0; i < d; i++)
            hash[i] = new BOBHash32(i + 750);
    }

    void print_basic_info()
    {
        printf("CU sketch\n");
        printf("\tCounters: %d\n", w);
        printf("\tMemory: %.6lfMB\n", w * 4.0 / 1024 / 1024);
    }

    void prints(){
        printf("CM_CU sketch:  \n");
        for(int i=0;i<w;i++){

            if(i%99==0){
                printf("%d\t\n",counters[i]);
            }else{
                printf("%d\t",counters[i]);
            }

        }
        printf("***********************************\n");

    }

    virtual ~CU_For_SC()
    {
        for (int i = 0; i < d; i++)
            delete hash[i];
    }

    void insert(uint32_t key, int f = 1)
    {
        int index[d];
        int value[d];
        int min_val = 1 << 30;

        for (int i = 0; i < d; i++) {
            index[i] = (hash[i]->run((const char *)&key, 4)) % w;
            value[i] = counters[index[i]];
            min_val = min(min_val, value[i]);
        }

        int temp = min_val + f;
        for (int i = 0; i < d; i++) {
            counters[index[i]] = max(counters[index[i]], temp);
        }
    }

    int query(uint32_t key)
    {
        int ret = 1 << 30;
        for (int i = 0; i < d; i++) {
            int tmp = counters[(hash[i]->run((const char *)&key, 4)) % w];
            ret = min(ret, tmp);
        }
        return ret;
    }

    int batch_query(uint32_t * data, int n) {
        int ret = 0;
        for (int i = 0; i < n; ++i) {
            ret += query(data[i]);
        }
        return ret;
    }

    void build(uint32_t * data, int n)
    {
        for (int i = 0; i < n; ++i) {
            insert(data[i], 1);
        }
    }
};
#endif //STREAMCLASSIFIER_CU_FOR_SC_H
