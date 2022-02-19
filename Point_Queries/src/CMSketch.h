#ifndef _CMSKETCH_H
#define _CMSKETCH_H

#include "params.h"
#include "BOBHash32.h"
#include "SPA.h"
#include <cstring>
#include <algorithm>

using namespace std;


class CMSketch: public SPA
{
private:
    int w=0 ;
	int *counters;
	int d=0;
	BOBHash32 * hash[3];
public:
    CMSketch(int memory_in_bytes, int dp)
    {
        w=memory_in_bytes/4;
        counters=new int[memory_in_bytes/4];
        for(int i=0;i<memory_in_bytes/4;i++){
            counters[i]=0;
        }
        d=dp;
//        memset(counters, 0, sizeof(counters));
        for (int i = 0; i < dp; i++)
            hash[i] = new BOBHash32(i + 750);
    }

    void print_basic_info()
    {
        printf("CM sketch\n");
        printf("\tCounters: %d\n", w);
        printf("\tMemory: %.6lfMB\n", w * 4.0 / 1024 / 1024);
    }


    virtual ~CMSketch()
    {
        for (int i = 0; i < d; i++)
            delete hash[i];
    }

    int * get(){
        return counters;
    }

    void insert(uint32_t key, int f)
    {
        for (int i = 0; i < d; i++) {
            int index = (hash[i]->run((const char *)&key, 4)) % w;
            counters[index] += f;
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

#endif //_CMSKETCH_H