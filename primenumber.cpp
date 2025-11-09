#include "primenumber.h"
#include <vector>
#include <iostream>
#include <random>
#include <fstream>
#include <cstdint>
#include <cmath>
#include<algorithm>
PrimeNumer::PrimeNumer() {
    prime1.resize(32);
    prime2.resize(32);
}
void PrimeNumer::generateOdd() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, UINT32_MAX);
    std::uniform_int_distribution<> dist1(pow(2,31), UINT32_MAX);
    for(int i=0;i<31;i++) prime1[i]=dist(gen);
    prime1[31]=dist1(gen);
    // std::ofstream outf;
    // outf.open("fixprime.txt");
    // int cnt=0;
    // for(int i=2;cnt<40;i++){
    //     if(judge(i)){
    //         outf<<i<<std::endl;
    //         cnt++;
    //     }
    // }
    // outf.close(); //*生成40个固定质数*//
    std::ifstream inf;
    inf.open("fixprime.txt");
    uint32_t t;
    while(inf>>t){
        fixprime.push_back(t);
    }
    inf.close();
}
std::vector<uint32_t> PrimeNumer::findD(int &s){
    std::vector<uint32_t>d=prime1;
    std::uint64_t t=1;
    for(size_t i=0;i<d.size();i++){
        uint64_t rsu=(uint64_t)d[i]-t; //防止-1溢出
        if((int64_t)rsu<0){
            d[i]=(uint32_t)(rsu+((uint64_t)1<<32));
        }
        else {
            d[i]=(uint32_t)rsu;
            t=0;
        }
    }
    while(d.size()>0&&d.back()==0) d.pop_back();
    while((d[0]&1)==0){
        div(d);
        s++;
        if(d.size()==1&&d[0]==0) break; //防止当出现2的幂次时死循环
    }
    return d;
}
void PrimeNumer::div(std::vector<uint32_t>&d){
    if(d.size()==1&&d[0]==1) d[0]=0;
    uint32_t sheng=0;
    for(int i=(int)d.size()-1;i>=0;i--){ //这样写是当i==0时，i--会变为无穷大，报错
        uint64_t needdiv=((uint64_t)sheng << 32) | (uint64_t)d[i];
        d[i]=(uint32_t)(needdiv>>1);
        sheng=(uint32_t)(needdiv&1);
        //std::cout<<d[i]<<std::endl<<std::endl;
    }
    while(d.size()>1&&d.back()==0) d.pop_back();
}
bool PrimeNumer::miller_rabin(){
    /***test***/
    // generatePrime();
    // for(int i=0;i<fixprime.size();i++) std::cout<<fixprime[i]<<' ';
    /*******/
    while(true){
        generateOdd();
        int s=0;
        std::vector<uint32_t>d=findD(s);
        for(int i=0;i<fixprime.size();i++){

        }
    }
     return true;


}
std::vector<uint32_t> PrimeNumer::quickExp(std::vector<uint32_t> base ,std::vector<uint32_t>&exp,std::vector<uint32_t>&mod1){
    std::vector<uint32_t>rsu={1};
    std::vector<uint32_t>d=exp;
    while(!(d.size()==1&&d[0]==0)){

        if((d[0]&1)==1){
            rsu=karatsuba(rsu,base);

            mod(rsu,mod1);
        }
        div(d);
        base=karatsuba(base,base);
        mod(base,mod1);
    }
    return rsu;

}
int PrimeNumer::compare(std::vector<uint32_t>&a,std::vector<uint32_t>&mod){
    while (a.size()>1&& a.back()==0) a.pop_back();
    while (mod.size()>1&& mod.back()==0) mod.pop_back();
    if(a.size()<mod.size()) return -1;
    if(a.size()>mod.size()) return 1;
    for(int i=mod.size()-1;i>=0;i--){
        if(a[i]<mod[i]) return -1;
        if(a[i]>mod[i]) return 1;
    }
    return 0;
}
std::vector<uint32_t> PrimeNumer::sub(std::vector<uint32_t>&a,std::vector<uint32_t>&b){
    uint64_t tmp;
    uint32_t bor=0;
    int n=std::max(a.size(),b.size());
    std::vector<uint32_t>res(n);
    for(int i=0;i<n;i++){
        uint64_t ta=(i<a.size()?a[i]:0);
        uint64_t tb=(i<b.size()?b[i]:0);
        tmp=ta-tb-bor;
        if(int64_t(tmp)<0){
            bor=1;
            res[i]=(uint32_t)(((uint64_t)1<<32)+tmp);
        }
        else{
            bor=0;
            res[i]=(uint32_t)tmp;
        }
    }
    while(res.size()>0&&res.back()==0) res.pop_back();
    return res;
}
//求余数
void PrimeNumer::mod(std::vector<uint32_t>&a,std::vector<uint32_t>&mod1){

    while(compare(a,mod1) >= 0){
        int shift = (int)a.size() - (int)mod1.size();
        std::vector<uint32_t> tmp = mod1;
        tmp.insert(tmp.begin(), shift, 0);
        if(compare(a,tmp) < 0){
            tmp.erase(tmp.begin());
            shift--;
        }
        a = sub(a,tmp);
    }

}
void mul_mod(std::vector<uint32_t>&rsu,std::vector<uint32_t>&base,std::vector<uint32_t>&mod){

}
std::vector<uint32_t> PrimeNumer::mul(std::vector<uint32_t>&a,std::vector<uint32_t>&b){
    std::vector<uint32_t>res(a.size()+b.size());
    for(int i=0;i<a.size();i++){
        uint32_t jin=0;
        for(int j=0;j<b.size();j++){
            uint64_t t=res[i+j]+((uint64_t)a[i])*b[j]+jin;
            res[i+j]=(uint32_t)t;
            jin=t>>32;
        }
        int k=i+b.size();
        while(jin>0){
            uint64_t t=res[k]+jin;

            res[k]=(uint32_t)t;
            jin=t>>32;
            k++;
        }
    }
    while(res.size()>0&&res.back()==0) res.pop_back();
    return res;

}
//求两个大数乘
std::vector<uint32_t> PrimeNumer::karatsuba(std::vector<uint32_t>&a,std::vector<uint32_t>&b){
    int n=std::max(a.size(),b.size());
    if(n<=32) return mul(a,b);
    std::vector<uint32_t>rsu;
    std::vector<uint32_t>res;
    int mid=n/2;
    std::vector<uint32_t>a0(a.begin(),a.begin()+mid);
    std::vector<uint32_t>a1(a.begin()+mid,a.end());
    std::vector<uint32_t>b0(b.begin(),b.begin()+mid);
    std::vector<uint32_t>b1(b.begin()+mid,b.end());
    std::vector<uint32_t>r1=karatsuba(a0,b0); //ac
    std::vector<uint32_t>r2=karatsuba(a1,b1); //bd
    std::vector<uint32_t>adda=add(a0,a1);
    std::vector<uint32_t>addb=add(b0,b1);
    std::vector<uint32_t>t=karatsuba(adda,addb); //ac+ad+bc+bd
    std::vector<uint32_t>t1=sub(t,r1);
    std::vector<uint32_t>r3=sub(t1,r2); //ad+bc
    karatsuba_add(rsu,r1,0);
    karatsuba_add(rsu,r3,mid);
    karatsuba_add(rsu,r2,2*mid);
    return rsu;
}
void PrimeNumer::karatsuba_add(std::vector<uint32_t>&rsu,std::vector<uint32_t>&a,int offset){
    if(rsu.size()<a.size()+offset){
        rsu.resize(a.size()+offset,0);
    }
    uint32_t jin=0;
    for(int i=0;i<a.size()||jin;i++){
        uint64_t t=((uint64_t)(i<a.size()?a[i]:0))+jin+(uint64_t)rsu[i+offset];
        rsu[i+offset]=(uint32_t)t;
        jin=t>>32;
        if (i + offset + 1 >= rsu.size() && jin) {
            rsu.push_back(0);
        }
    }
}
std::vector<uint32_t> PrimeNumer::add(std::vector<uint32_t>&a,std::vector<uint32_t>&b){
    std::vector<uint32_t>res;
    uint32_t jin=0;
    int n=std::max(a.size(),b.size());
    for(int i=0;i<n;i++){
        uint64_t t=((uint64_t)(i<a.size()?a[i]:0))+(i<b.size()?b[i]:0)+jin;
        res.push_back((uint32_t)t);
        jin=t>>32;
    }
    if(jin>0) res.push_back(jin);
    return res;
}
// bool PrimeNumer::judge(int i){
//     if(i==2||i==3) return true;
//     for(int k=2;k<=sqrt(i);k++){
//         if(i%k==0) return false;
//     }
//     return true;
// }
