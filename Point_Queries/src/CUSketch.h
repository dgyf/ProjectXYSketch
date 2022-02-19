#ifndef _CUSKETCH_H
#define _CUSKETCH_H

#include "params.h"
#include "BOBHash32.h"
#include "SPA.h"
#include <cstring>
#include <algorithm>

using namespace std;


class CUSketch: public SPA
{
private:
    int w = 0;
	int *counters;
	BOBHash32 * hash[10];
	int d=0;
public:
    CUSketch(int memo,int dp)
    {
        d=dp;
        w=memo/4;
        counters=new int[memo/4];
        for(int i=0;i<memo/4;i++){
            counters[i]=0;
        }
//        memset(counters, 0, sizeof(counters));
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

    virtual ~CUSketch()
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

#endif //_CUSKETCH_H