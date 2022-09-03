#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <string>
#include <cstring>
#include <map>
#include <fstream>
#include <vector>
//#include "BOBHASH32.h"
#include "BOBHash32.h"
#include "params.h"
#include "ssummary.h"
#include "heavykeeper.h"
#include "spacesaving.h"
#include "LossyCounting.h"
#include "CSS.h"
using namespace std;
const long maxn = 100000000;
map <string ,int> B,C;
struct node {string x;int y;} p[10000005];

struct Node {
    int cnt = 0;//count nubmers of every item
    int index;//record index of item in hashtable
};
Node hashtable[maxn];
//ifstream fin("u1",ios::in|ios::binary);


const char *filename = "/home/yqliu/data/kosarak.dat"; // Top_1024(58760)
//const char *filename = "/home/yqliu/Desktop/data/webdocs.dat"; // Top_1024(58760)

ifstream fin("/home/yqliu/data/kosarak.dat",ios::in|ios::binary);
char a[105];
string Read()
{
    fin.read(a,13);
    a[13]='\0';
    string tmp=a;
    return tmp;
}

//int cmp(node i,node j) {return i.y>j.y;}
int cmp1(Node a, Node b) {
    return a.cnt > b.cnt;
}

int main()
{
    clock_t time_start=clock();
    int MEM,K;

    cin>>MEM>>K;
    cout<<"MEM="<<MEM<<"KB"<<endl<<"Find top"<<K<<endl<<endl;
    cout<<"preparing all algorithms"<<endl;
    //int m=100000;  // the number of flows
    int m=8019015;  // the number of flows
    // preparing heavykeeper
    int hk_M=(MEM*1024-K*54)/8;
    //for (hk_M=1; 32*hk_M*HK_d+432*K<=MEM*1024*8; hk_M++); if (hk_M%2==0) hk_M--;
    //for (hk_M=1; 8*hk_M*HK_d+8*K<=MEM*1024; hk_M++);
    if (hk_M%2==0) hk_M--;
    heavykeeper *hk; hk=new heavykeeper(hk_M,K); hk->clear();


    // Inserting

    clock_t time_log[3];
    time_log[0]=clock();

    ifstream fin;
    fin.open(filename, ios::in | ios::out);
    if(!fin){
        cout<<"open failed!"<<endl;
        exit(1);
    }
    //int temp, MAXv = -1;
    string temp;
    cout<<"Inserting !"<<endl;
    while (!fin.eof() && fin >> temp) {
        //cout<<"now: "<<temp<<endl;
        //B[temp]++;
        hk->Insert(temp);

    }
    fin.close();

/*
    for (int i=1; i<=m; i++)
	{
	    if (i%(m/10)==0) cout<<"Insert "<<i<<endl;
		string s=Read();
		//cout<<" now: "<<s<<endl;
		B[s]++;
		hk->Insert(s);
		//ss->Insert(s);
		//LC->Insert(s,i/LC_M); if (i%LC_M==0) LC->clear(i/LC_M);
		//css->Insert(s);
	}

 */
    time_log[1]=clock();
    double inserttime=(double) (time_log[1] - time_log[0]) / CLOCKS_PER_SEC;
    cout<<" Time of inserting is : "<<inserttime<<" s "<<endl;
	hk->work();
	//ss->work();
	//LC->work();
	//css->work();
    //time_log[2]=clock();

    cout<<"preparing true flow"<<endl;
    int tmp, MAX = -1;
    fin.open(filename, ios::in | ios::out);
    while (!fin.eof() && fin >> tmp) {
        if (tmp > MAX) {
            MAX = tmp;
        }
        hashtable[tmp].cnt++;
        hashtable[tmp].index = tmp;
    }
    cout << "MAX item number is: " << MAX << endl;
    sort(hashtable, hashtable + MAX + 1, cmp1);
    fin.close();

    /*
	// preparing true flow
	int cnt=0;
    for (map <string,int>::iterator sit=B.begin(); sit!=B.end(); sit++)
    {
        p[++cnt].x=sit->first;
        p[cnt].y=sit->second;
    }
    sort(p+1,p+cnt+1,cmp);

    for (int i=1; i<=K; i++) C[p[i].x]=p[i].y;

     */
    // Calculating PRE, ARE, AAE
    cout<<"Calculating"<<endl;
    const int kn=6;
    int kk[kn]={32,64,128,256,512,1024};
    for(int ik=0;ik<kn;ik++) {
        double hk_AAE=0; double hk_ARE=0;
        int hk_sum=0;
        string hk_string; int hk_num;
        vector<int> isexist(kk[ik],0);
        for (int i = 0; i < kk[ik]; i++) {
            hk_string = (hk->Query(i)).first;
            hk_num = (hk->Query(i)).second;
            int isint = atoi(hk_string.c_str());
            for (int j = 0; j < kk[ik]; j++) {
                if (isint == hashtable[j].index) {
                    //hk_AAE+=abs(B[hk_string]-hk_num); hk_ARE+=abs(B[hk_string]-hk_num)/(B[hk_string]+0.0);
                    hk_AAE += abs(hashtable[j].cnt - hk_num);
                    hk_ARE += abs(hashtable[j].cnt - hk_num) / (hashtable[j].cnt + 0.0);
                    hk_sum++;
                    isexist[j]=1;
                }
            }
            /*
            for (int j = 0; j < kk[ik]; j++) {
                if(isexist[j]==0){
                    hk_AAE+=hashtable[j].cnt;
                    hk_ARE+=1;
                }
            }
             */
        }
            double resum=(hk_sum/(kk[ik]+0.0));
            cout<<" Querying : k= "<<kk[ik]<<endl;
            printf("heavkeeper:\nAccepted: %d/%d  %.10f\nARE: %.10f\nAAE: %.10f\n\n",hk_sum,kk[ik],resum,hk_ARE/kk[ik],hk_AAE/(kk[ik]+0.0));
            ofstream ftime("hk_realkos_topk_54.dat", ios::app);

            ftime << MEM << " " << kk[ik] << " "
                  << resum << " "
                  << hk_ARE/kk[ik]<< " "
                  << hk_AAE/(kk[ik]+0.0) <<" "<<inserttime<< endl;
            ftime.close();
    }

    clock_t time_end=clock();
    //printf("heavkeeper:\nAccepted: %d/%d  %.10f\nARE: %.10f\nAAE: %.10f\n\n",hk_sum,K,(hk_sum/(K+0.0)),hk_ARE/K,hk_AAE/(K+0.0));
    cout<<" Time of inserting is : "<<(double) (time_log[1] - time_log[0]) / CLOCKS_PER_SEC<<" s "<<endl;
    //cout<<" Time of working is : "<<(double) (time_log[2] - time_log[1]) / CLOCKS_PER_SEC<<" s "<<endl;
    cout<<"The total time is : "<<(double) (time_start - time_end) / CLOCKS_PER_SEC<<" s "<<endl;
    return 0;
}

/*

int main()
{
    int MEM,K;
    if(!fin){
        cout<<"open failed!"<<endl;
        exit(1);
    }
    cin>>MEM>>K;
    cout<<"MEM="<<MEM<<"KB"<<endl<<"Find top"<<K<<endl<<endl;
    cout<<"preparing all algorithms"<<endl;
    //int m=100000;  // the number of flows
    int m=8019015;  // the number of flows
    // preparing heavykeeper
    int hk_M;
    for (hk_M=1; 32*hk_M*HK_d+432*K<=MEM*1024*8; hk_M++); if (hk_M%2==0) hk_M--;
    heavykeeper *hk; hk=new heavykeeper(hk_M,K); hk->clear();

    // preparing spacesaving
    int ss_M;
    for (ss_M=1; 432*ss_M<=MEM*1024*8; ss_M++);
    spacesaving *ss; ss=new spacesaving(ss_M,K);

    // preparing LossyCounting
    int LC_M;
    for (LC_M=1; 227*LC_M<=MEM*1024*8; LC_M++);
    LossyCounting *LC; LC=new LossyCounting(K);

    // preparing CSS
    int css_M;
    for (css_M=1; 179*css_M+4*css_M*log(css_M)/log(2)<=MEM*1024*8; css_M++);
    CSS *css; css=new CSS(css_M,K); css->clear();
    // Inserting
    for (int i=1; i<=m; i++)
    {
        if (i%(m/10)==0) cout<<"Insert "<<i<<endl;
        string s=Read();
        B[s]++;
        hk->Insert(s);
        ss->Insert(s);
        LC->Insert(s,i/LC_M); if (i%LC_M==0) LC->clear(i/LC_M);
        css->Insert(s);
    }
    hk->work();
    ss->work();
    LC->work();
    css->work();

    cout<<"preparing true flow"<<endl;
    // preparing true flow
    int cnt=0;
    for (map <string,int>::iterator sit=B.begin(); sit!=B.end(); sit++)
    {
        p[++cnt].x=sit->first;
        p[cnt].y=sit->second;
    }
    sort(p+1,p+cnt+1,cmp);
    for (int i=1; i<=K; i++) C[p[i].x]=p[i].y;

    // Calculating PRE, ARE, AAE
    cout<<"Calculating"<<endl;
    int hk_sum=0,hk_AAE=0; double hk_ARE=0;
    string hk_string; int hk_num;
    for (int i=0; i<K; i++)
    {
        hk_string=(hk->Query(i)).first; hk_num=(hk->Query(i)).second;
        hk_AAE+=abs(B[hk_string]-hk_num); hk_ARE+=abs(B[hk_string]-hk_num)/(B[hk_string]+0.0);
        if (C[hk_string]) hk_sum++;
    }

    int LC_sum=0,LC_AAE=0; double LC_ARE=0;
    string LC_string; int LC_num;
    for (int i=0; i<K; i++)
    {
        LC_string=(LC->Query(i)).first; LC_num=(LC->Query(i)).second;
        LC_AAE+=abs(B[LC_string]-LC_num); LC_ARE+=abs(B[LC_string]-LC_num)/(B[LC_string]+0.0);
        if (C[LC_string]) LC_sum++;
    }

    int ss_sum=0,ss_AAE=0; double ss_ARE=0;
    string ss_string; int ss_num;
    for (int i=0; i<K; i++)
    {
        ss_string=(ss->Query(i)).first; ss_num=(ss->Query(i)).second;
        ss_AAE+=abs(B[ss_string]-ss_num); ss_ARE+=abs(B[ss_string]-ss_num)/(B[ss_string]+0.0);
        if (C[ss_string]) ss_sum++;
    }

    int css_sum=0,css_AAE=0; double css_ARE=0;
    string css_string; int css_num;
    for (int i=0; i<K; i++)
    {
        css_string=(css->Query(i)).first; css_num=(css->Query(i)).second;
        css_AAE+=abs(B[css_string]-css_num); css_ARE+=abs(B[css_string]-css_num)/(B[css_string]+0.0);
        if (C[css_string]) css_sum++;
    }
    printf("heavkeeper:\nAccepted: %d/%d  %.10f\nARE: %.10f\nAAE: %.10f\n\n",hk_sum,K,(hk_sum/(K+0.0)),hk_ARE/K,hk_AAE/(K+0.0));
    printf("LossyCounting:\nAccepted: %d/%d  %.10f\nARE: %.10f\nAAE: %.10f\n\n",LC_sum,K,(LC_sum/(K+0.0)),LC_ARE/K,LC_AAE/(K+0.0));
    printf("spacesaving:\nAccepted: %d/%d  %.10f\nARE: %.10f\nAAE: %.10f\n\n",ss_sum,K,(ss_sum/(K+0.0)),ss_ARE/K,ss_AAE/(K+0.0));
    printf("CSS:\nAccepted: %d/%d  %.10f\nARE: %.10f\nAAE: %.10f\n\n",css_sum,K,(css_sum/(K+0.0)),css_ARE/K,css_AAE/(K+0.0));
    return 0;
}
 */