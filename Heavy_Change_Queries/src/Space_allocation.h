//
// Created by yqliu on 8/12/19.
//

#ifndef STREAMCLASSIFIER_SPACE_ALLOCATION_H
#define STREAMCLASSIFIER_SPACE_ALLOCATION_H

#include <cstring>
#include <cmath>

using namespace std;

class Space_allocation {
private:


public:
    Space_allocation() {
    }

    int *get_size(int space, int ran) {  //output: temp[0]:lenth; begining from temp[1];
        int room = space / 4;
        int range = ran;
//        int temp[ran+1];
        int *temp = new int[ran + 1];
        int count = 1;
        int swit = 0;
        while (swit == 0) {
            int index_c = floor(log(room) / log(2));
            double ram = pow(2, range - index_c);
            if (ram > 1 && ram > (room - pow(2, index_c)) && (room - pow(2, index_c)) > 2 * (range - index_c)) {
                temp[count] = index_c;
                range = range - index_c;
                room = room - pow(2, index_c);
            } else if (ram > 1 && room - pow(2, index_c) <= 2 * (range - index_c)) {
//                else if(room-pow(2,index_c)<=2*(range-index_c)&&(room-pow(2,index_c))!=0){
                temp[count] = index_c - 1;
                range = range - index_c + 1;
                room = room - pow(2, index_c - 1);
            } else if (ram <= 1) {
                swit = 1;
                temp[count] = range;
                temp[0] = count;
            }

            if (ram <= room - pow(2, index_c) && ram > 1) {
                temp[count] = index_c;
                temp[count + 1] = range - index_c;
                swit = 1;
                temp[0] = count + 1;  // length of the temp;

            }
            count += 1;
        }
        return temp;
    }


};


#endif //STREAMCLASSIFIER_SPACE_ALLOCATION_H
