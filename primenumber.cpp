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
std::vector<uint32_t> PrimeNumer::quickExp(std::vector<uint32_t> base ,std::vector<uint32_t>&exp,std::vector<uint32_t>&mod1,std::vector<uint32_t>&mu){
    std::vector<uint32_t>rsu={1};
    std::vector<uint32_t>d=exp;
    while(!(d.size()==1&&d[0]==0)){

        if((d[0]&1)==1){
            rsu=karatsuba(rsu,base,mod1,mu);

            barrett_mod(rsu,mod1,mu);
        }
        div(d);
        base=karatsuba(base,base,mod1,mu);
        barrett_mod(base,mod1,mu);
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
std::vector<uint32_t> PrimeNumer::sub(const std::vector<uint32_t>&a,const std::vector<uint32_t>&b){
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
std::vector<uint32_t> PrimeNumer::mul(std::vector<uint32_t>&a,std::vector<uint32_t>&b,std::vector<uint32_t>&mod,std::vector<uint32_t>&mu){
    barrett_mod(a,mod,mu);
    barrett_mod(b,mod,mu);
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
    barrett_mod(res,mod,mu);
    while(res.size()>0&&res.back()==0) res.pop_back();
    return res;

}
std::vector<uint32_t> PrimeNumer::mul_barrett(std::vector<uint32_t>&a,std::vector<uint32_t>&b){

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
std::vector<uint32_t> PrimeNumer::karatsuba(std::vector<uint32_t>&a,std::vector<uint32_t>&b,std::vector<uint32_t>&mod,std::vector<uint32_t>&mu){
    int n=std::max(a.size(),b.size());
    if(n<=20) return mul(a,b,mod,mu);
    barrett_mod(a,mod,mu);
    barrett_mod(b,mod,mu);
    std::vector<uint32_t>rsu;
    std::vector<uint32_t>res;
    int mid=n/2;
    // std::vector<uint32_t> a0(a.begin(), a.begin() + std::min(mid, (int)a.size()));
    // std::vector<uint32_t> a1(a.size() > mid ? std::vector<uint32_t>(a.begin() + mid, a.end()) : std::vector<uint32_t>{});
    // std::vector<uint32_t> b0(b.begin(), b.begin() + std::min(mid, (int)b.size()));
    // std::vector<uint32_t> b1(b.size() > mid ? std::vector<uint32_t>(b.begin() + mid, b.end()) : std::vector<uint32_t>{});
    std::vector<uint32_t>a0(a.begin(),a.begin()+mid);
    std::vector<uint32_t>a1(a.begin()+mid,a.end());
    std::vector<uint32_t>b0(b.begin(),b.begin()+mid);
    std::vector<uint32_t>b1(b.begin()+mid,b.end());
    std::vector<uint32_t>r1=karatsuba(a0,b0,mod,mu); //ac
    std::vector<uint32_t>r2=karatsuba(a1,b1,mod,mu); //bd
    std::vector<uint32_t>adda=add(a0,a1);
    std::vector<uint32_t>addb=add(b0,b1);
    std::vector<uint32_t>t=karatsuba(adda,addb,mod,mu); //ac+ad+bc+bd
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
void PrimeNumer::barrett_mod(std::vector<uint32_t>& a, std::vector<uint32_t>& mod, std::vector<uint32_t> & mu){

    trim(a);
    int k = (int)mod.size();
    if (a.size() < (size_t)k) return;
    // q1 = floor(a / b^(k-1)) => high words starting at index k-1
    std::vector<uint32_t> q1;
    if (a.size() > (size_t)(k - 1))
        q1 = std::vector<uint32_t>(a.begin() + (k - 1), a.end());
    else q1 = {0};
    // q2 = q1 * mu
    std::vector<uint32_t> q2 = mul_barrett(q1, mu);
    // q3 = floor(q2 / b^(k+1)) => take high words starting at index k+1
    std::vector<uint32_t> q3;
    if (q2.size() > (size_t)(k + 1))
        q3 = std::vector<uint32_t>(q2.begin() + (k + 1), q2.end());
    else q3 = {0};
    // r = a - q3 * mod
    std::vector<uint32_t> q3m = mul_barrett(q3, mod);
    std::vector<uint32_t> r;
    if (compare(a, q3m) >= 0) r = sub(a, q3m);
    else {
        // r = a + t where t = mod * ceil((q3m - a) / mod) - q3m; simpler: add mod until >= q3m then subtract r = a;
        while (compare(r, q3m) < 0) r = add(r, mod);
        // add_big 使用者已有 add; 下面给实现调用者可替换为 existing add
        r = sub(r, q3m);
    } while (compare(r, mod) >= 0) r = sub(r, mod);
    // write back
    a.assign(r.begin(), r.end());
    trim(a);
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
std::vector<uint32_t> PrimeNumer::compute_mu(std::vector<uint32_t>& mod){
    // call knuth division on b^(2k) by mod
    int k = (int)mod.size();
    std::vector<uint32_t> b2k(2 * k + 1, 0);
    b2k[2 * k] = 1; // b^(2k)
    if (b2k.size() > 10000) {
        std::cerr << "Mul overflow compumu: a=" << b2k.size()<< std::endl;
        exit(1);
    }
    auto qr = div_mod_knuth(b2k, mod); //   returns {quotient, remainder} // qr.first is mu
    std::vector<uint32_t> mu(qr.first.begin(), qr.first.end());
    return mu;
}
std::pair<std::vector<uint32_t>, std::vector<uint32_t>> PrimeNumer::div_mod_knuth(const std::vector<uint32_t>& a_in, const std::vector<uint32_t>& b_in) {

    std::vector<uint32_t> a = a_in;
    std::vector<uint32_t> b = b_in;
    trim(a); trim(b);
    if (b.size() == 1) {

    // short division by single limb
        uint64_t rem = 0;
    std::vector<uint32_t> q(a.size(), 0);
    for (size_t i = a.size(); i-- > 0;) {
        unsigned __int128 cur = ((unsigned __int128)rem << 32) | a[i];
        q[i] = (uint32_t)(cur / b[0]);
        rem = (uint64_t)(cur % b[0]);
    }
    trim(q);
    return {q, std::vector<uint32_t>{(uint32_t)rem}}; }
    if (compare(a, b) < 0)
        return { std::vector<uint32_t>{0}, a };
    size_t n = b.size();
    size_t m = a.size() - n; // normalization: shift so b.back() high bit set
    unsigned shift = 0; uint32_t bn = b.back();
    while ((bn << shift & 0x80000000u) == 0)
        ++shift;
    std::vector<uint32_t> b_norm = lshift_bits(b, shift);
    std::vector<uint32_t> a_norm = lshift_bits(a, shift);
    if (a_norm.size() == a.size()) a_norm.push_back(0);
    if (a_norm.size() < n + m + 1) a_norm.resize(n + m + 1, 0);
    std::vector<uint32_t> q(m + 1, 0);
    for (size_t j = m + 1; j-- > 0;) {

    // estimate qhat
        unsigned __int128 num = ((unsigned __int128)a_norm[j + n] << 32) | a_norm[j + n - 1];
    unsigned __int128 den = b_norm[n - 1];
        unsigned __int128 qhat = num / den;
    if (qhat > 0xFFFFFFFFu) qhat = 0xFFFFFFFFu;
    // adjust qhat
    while (true) {

        unsigned __int128 left = qhat * (unsigned __int128)b_norm[n - 2];
        unsigned __int128 remhat = num - qhat * den;
        unsigned __int128 right = (remhat << 32) + a_norm[j + n - 2];
        if (left <= right) break;
        --qhat;
    }
    // multiply-subtract  qhat * b_norm from a_norm at offset j
    unsigned __int128 borrow = 0;
    for (size_t i = 0; i < n; ++i) {
        unsigned __int128 prod = qhat * (unsigned __int128)b_norm[i];
        unsigned __int128 sub = (unsigned __int128)a_norm[j + i] - (prod & 0xFFFFFFFFu) - borrow;

        a_norm[j + i] = (uint32_t)sub; borrow = (prod >> 32) + ((sub >> 32) & 1);
    }

    unsigned __int128 sub_last = (unsigned __int128)a_norm[j + n] - borrow; a_norm[j + n] = (uint32_t)sub_last;
    bool negative = ((sub_last >> 63) & 1);

    if (negative) { // add back b_norm
        --qhat;
        uint64_t carry = 0;
        for (size_t i = 0; i < n; ++i) {
            unsigned __int128 sum = (unsigned __int128)a_norm[j + i] + b_norm[i] + carry;
            a_norm[j + i] = (uint32_t)sum;
            carry = (uint64_t)(sum >> 32);
        }
        a_norm[j + n] = (uint32_t)((unsigned __int128)a_norm[j + n] + carry);
    }
    q[j] = (uint32_t)qhat;
    }
    std::vector<uint32_t> r(a_norm.begin(), a_norm.begin() + n);
    if (shift) rshift_bits_inplace(r, shift);
    trim(q); trim(r);
    return {q, r};
}
std::vector<uint32_t> PrimeNumer::lshift_bits(const std::vector<uint32_t>& v, unsigned s) {
    if (s == 0) return v;
    std::vector<uint32_t> out(v.size() + 1, 0);
    uint64_t carry = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        uint64_t cur = ((uint64_t)v[i] << s) | carry;
        out[i] = (uint64_t)cur; carry = cur >> 32;
    }
    out[v.size()] = (uint32_t)carry; trim(out); return out;
}
void PrimeNumer::rshift_bits_inplace(std::vector<uint32_t>& v, unsigned s) {
    if (s == 0) return;
    uint64_t carry = 0;
    for (size_t i = v.size(); i-- > 0;) {
        unsigned __int128 cur = ((unsigned __int128)carry << 32) | v[i]; v[i] = (uint64_t)(cur >> s);
        carry = (uint64_t)(cur & (((unsigned __int128)1 << s) - 1));
    }
    trim(v);
}
