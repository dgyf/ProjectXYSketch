//
// Created by yqliu on 5/10/19.
//

#ifndef STREAMCLASSIFIER_CSKETCH_H
#define STREAMCLASSIFIER_CSKETCH_H

#include "params.h"
#include "BOBHash32.h"
#include "SPA.h"
#include <cstring>
#include <algorithm>
#include <math.h>
//#include <ctime>
//#include <set>
#include <iostream>
//#include <bitset>
//#include <cstdlib>

using namespace std;


class Csketch: public SPA
{
private:
    int w =0;
    int *counters;
    BOBHash32 * hash[10];
    int *binhas;
    int d=0;
public:
    Csketch(int memo,int dp)
    {  binhas=new int[dp];
       w=memo/4;
       d=dp;
       counters=new int[memo/4];
       for(int i=0;i<memo/4;i++){
           counters[i]=0;
       }
        srand(time(0));
//        memset(counters, 0, sizeof(counters));
        for (int i = 0; i < dp; i++){
            hash[i] = new BOBHash32(i + 750);
            binhas[i]=rand()%1024;
        }
    }


    virtual ~Csketch()
    {
        for (int i = 0; i < d; i++)
            delete hash[i];
    }

    int chas(uint32_t key){
        int out=pow(-1,key%2);
        return out;

    }

    void insert(uint32_t key, int f)
    {
        for (int i = 0; i < d; i++) {
            int index = (hash[i]->run((const char *)&key, 4)) % w;
            counters[index] += pow(-1,binhas[i]+key)*f;
        }
    }

    int query(uint32_t key)
    {
        int ret= 1 << 30;
        int med[d];
        for (int i = 0; i < d; i++) {
            int tmp = counters[(hash[i]->run((const char *)&key, 4)) % w];
            med[i]=pow(-1,binhas[i]+key)*tmp;
        }
        nth_element(med,med,med+d);
        ret=abs(med[d/2]);
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

    void print_basic_info()
    {
        printf("CM sketch\n");
        printf("\tCounters: %d\n", w);
        printf("\tMemory: %.6lfMB\n", w * 4.0 / 1024 / 1024);
    }
};




#endif //STREAMCLASSIFIER_CSKETCH_H
