#include "primenumber.h"
#include <vector>
#include <iostream>
#include <random>
#include <fstream>
#include <cstdint>
#include <cmath>
#include<algorithm>
#include <QTranslator>
PrimeNumer::PrimeNumer() {
    prime1={478585197,103485950,2904745936,2847268253,296791683,1193118780,1892927142,2044557110,855224915,2341389853,1166781217,3371591298,1400981518,532120207,1563684065,3426311462};
    prime2={757897211,4128549411,3339880654,130291540,756149348,1356868423,2413761437,3391807815,3679400074,484489568,2746056864,2667321572,1959813599,1976730658,1777281610,4168282924};
    p1_n_1={478585196,103485950,2904745936,2847268253,296791683,1193118780,1892927142,2044557110,855224915,2341389853,1166781217,3371591298,1400981518,532120207,1563684065,3426311462};
    p2_n_1={757897210,4128549411,3339880654,130291540,756149348,1356868423,2413761437,3391807815,3679400074,484489568,2746056864,2667321572,1959813599,1976730658,1777281610,4168282924};
    two.push_back(2);
    std::ifstream inf;
    inf.open("fixprime.txt");
    uint32_t t;
    while(inf>>t){
        fixprime.push_back(t);
    }
    inf.close();
}
std::vector<uint32_t> PrimeNumer::generateOdd() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, UINT32_MAX);
    std::uniform_int_distribution<> dist1(pow(2,31), UINT32_MAX);
    std::vector<uint32_t>odd(16);
    for(int i=0;i<15;i++) {
        odd[i]=dist(gen);
    }
    odd[15]=dist1(gen);
    odd[0]|=1;
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
    return odd;
}
std::vector<uint32_t> PrimeNumer::findD(std::vector<uint32_t>&odd,int &s){
    std::vector<uint32_t>d=odd;
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
    n_1=d;
    while(d.size()>0&&d.back()==0) d.pop_back();
    while((d[0]&1)==0){
        div(d);
        s++;
        if(d.size()==1&&d[0]==0) break; //防止当出现2的幂次时死循环
    }
    return d;
}
void PrimeNumer::div(std::vector<uint32_t>&v){
    uint64_t carry = 0;
    for (size_t i = v.size(); i-- > 0; ){
        uint64_t cur = (carry << 32) | v[i];
        v[i] = cur >> 1;
        carry = cur & 1;
    }
    trim(v);
    ensure_nonempty_zero(v);
}
bool PrimeNumer::miller_rabin(){
    int tries = 0;
    bool p1=false,p2=false;
    while(!p1||!p2) {
        std::vector<uint32_t>odd=generateOdd();
        int s=0;
        std::vector<uint32_t>d=findD(odd,s);
        //std::cout<<"d"<<d[0]<<std::endl;
        bool isprime = true;
        std::vector<uint32_t> mu = compute_mu(odd);
        for(size_t i = 0; i < 15; i++) {
            std::vector<uint32_t>a={fixprime[i]};
            std::vector<uint32_t>x=quickExp(a, d, odd, mu);
            // for(auto t:x) std::cout<<t<<' ';
            // std::cout<<std::endl;
            if(x==std::vector<uint32_t>{1}||x==n_1) continue;
            bool passed = false;
            for(int r=0;r<s;r++){
                x=karatsuba(x,x);
                barrett_mod(x,odd,mu);
                if(x==n_1) {
                    passed = true;
                    break;
                }
                if(x == std::vector<uint32_t>{1}) {
                    isprime=false;
                    break;
                }
                // for(auto t:x) std::cout<<t<<' ';
                // std::cout<<std::endl<<std::endl;
            }

            if(!passed) {
                isprime=false;
                break;
            }
        }
        tries++;
        if(isprime) {
            if(!p1){
                prime1=odd;

                p1_n_1=n_1;
                p1=true;
            }
            else if(!p2){
                prime2=odd;
                p2_n_1=n_1;
                p2=true;
            }

        }
        //prime1=add(prime1,two);
    }
    std::cout << "Found prime after " << tries << " tries." << std::endl;
    return true;


}
std::vector<uint32_t> PrimeNumer::quickExp(std::vector<uint32_t> &base ,std::vector<uint32_t>&exp,std::vector<uint32_t>&mod1,std::vector<uint32_t>&mu){
    std::vector<uint32_t>rsu={1};
    std::vector<uint32_t>d=exp;
    while(!(d.size()==1&&d[0]==0)){

        if((d[0]&1)==1){
            rsu=karatsuba(rsu,base);

            barrett_mod(rsu,mod1,mu);
        }
        div(d);
        base=karatsuba(base,base);
        barrett_mod(base,mod1,mu);
    }
    return rsu;

}
int PrimeNumer::compare(const std::vector<uint32_t>&a,const std::vector<uint32_t>&mod){
    int as=a.size();
    int bs=mod.size();
    int i=as,j=bs;
    while(i>1&&a[i-1]==0) --i;
    while(j>1&&a[j-1]==0) --j;
    if(i<j) return -1;
    if(i>j) return 1;
    for(int k=i-1;k>=0;k--){
        //std::cout<<"k"<<a[k]<<std::endl;
        if(a[k]<mod[k]) return -1;
        if(a[k]>mod[k]) return 1;
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
        if(ta>=tb+bor){
            tmp=ta-tb-bor;
            bor=0;
        }
        else{
            tmp=((uint64_t)1<<32)+ta-tb-bor;
            bor=1;
        }
        res[i]=(uint32_t)tmp;
    }
    trim(res);
    return res;
}
//求余数
std::vector<uint32_t> PrimeNumer::mod(const std::vector<uint32_t>&need, const std::vector<uint32_t>&mod1){
    std::vector<uint32_t>a=need;
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
    return a;
}
void mul_mod(std::vector<uint32_t>&rsu,std::vector<uint32_t>&base,std::vector<uint32_t>&mod){

}
std::vector<uint32_t> PrimeNumer::mul(const std::vector<uint32_t>&a,const std::vector<uint32_t>&b){
    if(a.size()==1&&a[0]==0) return std::vector<uint32_t>{0};
    if(b.size()==1&&b[0]==0) return std::vector<uint32_t>{0};
    std::vector<uint32_t>res(a.size()+b.size(),0);
    for(int i=0;i<a.size();i++){
        uint64_t jin=0;
        for(int j=0;j<b.size();j++){
            uint64_t t=(uint64_t)res[i+j]+((uint64_t)a[i])*(uint64_t)b[j]+jin;
            res[i+j]=(uint32_t)t;
            jin=t>>32;
        }
        int k=i+b.size();
        while(jin>0){
            uint64_t t=(uint64_t)res[k]+jin;

            res[k]=(uint32_t)t;
            jin=t>>32;
            k++;
        }

    }
    trim(res);
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

    trim(res);
    return res;

}
//求两个大数乘
std::vector<uint32_t> PrimeNumer::karatsuba(const std::vector<uint32_t>&a,const std::vector<uint32_t>&b){
    int n=std::max(a.size(),b.size());
    if(n==0) return std::vector<uint32_t>{0};
    if(n<=64) return mul(a,b);
    std::vector<uint32_t>res;
    int mid=n/2;
    std::vector<uint32_t> a0(std::min(mid,(int)a.size()));
    std::copy(a.begin(),a.begin()+std::min((size_t)mid,a.size()),a0.begin());
    std::vector<uint32_t> a1;
    if(a.size()>mid) a1.assign(a.begin()+mid,a.end());
    else a1.clear();

    std::vector<uint32_t> b0(std::min(mid,(int)b.size()));
    std::copy(b.begin(),b.begin()+std::min((size_t)mid,b.size()),b0.begin());
    std::vector<uint32_t> b1;
    if(b.size()>mid) b1.assign(b.begin()+mid,b.end());
    else b1.clear();
    // std::vector<uint32_t>a0(a.begin(),a.begin()+mid);
    // std::vector<uint32_t>a1(a.begin()+mid,a.end());
    // std::vector<uint32_t>b0(b.begin(),b.begin()+mid);
    // std::vector<uint32_t>b1(b.begin()+mid,b.end());
    std::vector<uint32_t>r1=karatsuba(a0,b0); //ac
    std::vector<uint32_t>r2=karatsuba(a1,b1); //bd
    std::vector<uint32_t>adda=add(a0,a1);
    std::vector<uint32_t>addb=add(b0,b1);
    std::vector<uint32_t>t=karatsuba(adda,addb); //ac+ad+bc+bd
    std::vector<uint32_t>t1=sub(t,r1);
    std::vector<uint32_t>r3=sub(t1,r2); //ad+bc
    karatsuba_add(res,r1,0);
    karatsuba_add(res,r3,mid);
    karatsuba_add(res,r2,2*mid);
    trim(res);
    return res;
}

void PrimeNumer::karatsuba_add(std::vector<uint32_t>&rsu,const std::vector<uint32_t>&a,int offset){
    if(a.empty()||a.size()==1&&a[0]==0) return ;
    if(rsu.size()<a.size()+offset+1){
        rsu.resize(a.size()+offset+1,0);
    }
    uint64_t jin=0;
    int i=0;
    for(;i<a.size();i++){
        uint64_t t=(uint64_t)a[i]+jin+(uint64_t)rsu[i+offset];
        rsu[i+offset]=(uint32_t)t;
        jin=t>>32;
    }
    int pos=i+offset;
    while(jin){
        uint64_t t=(uint64_t)rsu[pos]+jin;
        rsu[pos]=(uint32_t)t;
        jin=t>>32;
        pos++;
        if(pos>=rsu.size()) rsu.push_back(0);
    }
    trim(rsu);
}
void PrimeNumer::barrett_mod(std::vector<uint32_t>& a,const std::vector<uint32_t>& mod, std::vector<uint32_t> & mu){

    const size_t k = mod.size();            // number of limbs in modulus
    if (k == 0) return;                     // defensive
    if (compare(a, mod) < 0) return;        // already reduced
    if (a.size() > 2 * k) {
        auto qr = div_mod_knuth(a, mod);
        a = qr.second;
        trim(a);
        if (compare(a, mod) < 0) return;
    }
    std::vector<uint32_t> q1;
    if (a.size() > (size_t)(k - 1)) {
        q1.assign(a.begin() + (k - 1), a.end());
    } else {
        q1 = {0};
    }
    trim(q1);
    std::vector<uint32_t> q2 = karatsuba(q1, mu);

    std::vector<uint32_t> q3;
    if (q2.size() > (size_t)(k + 1)) {
        q3.assign(q2.begin() + (k + 1), q2.end());
    } else {
        q3 = {0};
    }
    trim(q3);

    std::vector<uint32_t> q3m = karatsuba(q3, mod);
    trim(q3m);

    std::vector<uint32_t> r;
    if (compare(a, q3m) >= 0) {
        r = sub(a, q3m);
        trim(r);
    } else {
        auto qr = div_mod_knuth(a, mod);
        a = qr.second;
        trim(a);
        return;
    }

    if (compare(r, mod) >= 0) {
        r = sub(r, mod);
        trim(r);
        if (compare(r, mod) >= 0) {
            r = sub(r, mod);
            trim(r);
        }
    }
    a.swap(r);
    trim(a);
}
std::vector<uint32_t> PrimeNumer::add(const std::vector<uint32_t>&a,const std::vector<uint32_t>&b){
    std::vector<uint32_t>res;
    uint32_t jin=0;
    int n=std::max(a.size(),b.size());
    res.resize(n+1);
    for(int i=0;i<n;i++){
        uint64_t t=((uint64_t)(i<a.size()?a[i]:0))+(i<b.size()?b[i]:0)+jin;
        res[i]=((uint32_t)t);
        jin=t>>32;
    }
    if(jin>0) res[n]=(jin);
    trim(res);
    return res;
}
// bool PrimeNumer::judge(int i){
//     if(i==2||i==3) return true;
//     for(int k=2;k<=sqrt(i);k++){
//         if(i%k==0) return false;
//     }
//     return true;
// }
std::vector<uint32_t> PrimeNumer::compute_mu(const std::vector<uint32_t>& m){
    size_t k = m.size();

    std::vector<uint32_t> b2k(2 * k + 1, 0);
    b2k[2 * k] = 1;  // b^(2k)

    auto qr = div_mod_knuth(b2k, m);
    auto mu = qr.first;

    trim(mu);
    ensure_nonempty_zero(mu);
    return mu;
}
std::pair<std::vector<uint32_t>, std::vector<uint32_t>> PrimeNumer::div_mod_knuth(const std::vector<uint32_t>& a_in, const std::vector<uint32_t>& b_in) {
    std::vector<uint32_t> a = a_in;
    std::vector<uint32_t> b = b_in;
    trim(a); trim(b);
    ensure_nonempty_zero(a); ensure_nonempty_zero(b);


    if (b.size() == 1) {
        uint64_t rem = 0;
        std::vector<uint32_t> q(a.size(), 0);
        for (size_t ii = a.size(); ii-- > 0;) {
            unsigned __int128 cur = ((unsigned __int128)rem << 32) | a[ii];
            uint64_t qdigit = (uint64_t)(cur / b[0]);
            rem = (uint64_t)(cur % b[0]);
            q[ii] = (uint32_t)qdigit;
            if (ii==0) break;
        }
        trim(q);
        std::vector<uint32_t> r(1, (uint32_t)rem);
        trim(r);
        ensure_nonempty_zero(q);
        ensure_nonempty_zero(r);
        return {q, r};
    }

    if (compare(a, b) < 0) {

        std::vector<uint32_t> q(1, 0);
        trim(a);
        ensure_nonempty_zero(a);
        return {q, a};
    }

    size_t n = b.size();
    size_t m = (a.size() > n) ? (a.size() - n) : 0;


    uint32_t bn = b.back();
    unsigned shift = 0;
    if (bn == 0) {
        trim(b);
        bn = b.back();
    }

    while (shift < 32 && (((uint32_t)(bn << shift)) & 0x80000000u) == 0u) ++shift;

    std::vector<uint32_t> b_norm = lshift_bits(b, shift);
    std::vector<uint32_t> a_norm = lshift_bits(a, shift);

    if (a_norm.size() < n + m + 1) a_norm.resize(n + m + 1, 0u);

    std::vector<uint32_t> q(m + 1, 0u);


    for (size_t j = m + 1; j-- > 0;) {

        unsigned __int128 top = ((unsigned __int128)a_norm[j + n] << 32) | a_norm[j + n - 1];
        unsigned __int128 den = b_norm[n - 1];
        unsigned __int128 qhat = top / den;
        if (qhat > 0xFFFFFFFFull) qhat = 0xFFFFFFFFull;


        while (true) {
            unsigned __int128 left = qhat * (unsigned __int128)b_norm[n - 2];
            unsigned __int128 remhat = top - qhat * den;
            unsigned __int128 right = (remhat << 32) + a_norm[j + n - 2];
            if (left <= right) break;
            --qhat;
        }

        unsigned __int128 borrow = 0;
        for (size_t i = 0; i < n; ++i) {
            unsigned __int128 prod = qhat * (unsigned __int128)b_norm[i];
            unsigned __int128 ai = (unsigned __int128)a_norm[j + i];
            unsigned __int128 sub = ai - (prod & 0xFFFFFFFFull) - borrow;
            a_norm[j + i] = (uint32_t)sub;
            unsigned __int128 prod_high = (prod >> 32);
            unsigned __int128 prod_low_plus_borrow = (prod & 0xFFFFFFFFull) + borrow;
            borrow = prod_high + ((ai < prod_low_plus_borrow) ? 1 : 0);
        }
        unsigned __int128 ai_last = (unsigned __int128)a_norm[j + n];
        unsigned __int128 sub_last = ai_last - borrow;
        a_norm[j + n] = (uint32_t)sub_last;
        bool negative = ((sub_last >> 63) & 1);

        if (negative) {
            --qhat;
            uint64_t carry = 0;
            for (size_t i = 0; i < n; ++i) {
                unsigned __int128 sum = (unsigned __int128)a_norm[j + i] + (unsigned __int128)b_norm[i] + carry;
                a_norm[j + i] = (uint32_t)sum;
                carry = (uint64_t)(sum >> 32);
            }
            unsigned __int128 new_high = (unsigned __int128)a_norm[j + n] + carry;
            a_norm[j + n] = (uint32_t)new_high;
        }
        q[j] = (uint32_t)qhat;
    }
    std::vector<uint32_t> r(a_norm.begin(), a_norm.begin() + n);
    if (shift) rshift_bits_inplace(r, shift);
    trim(q); trim(r);
    ensure_nonempty_zero(q);
    ensure_nonempty_zero(r);
    return {q, r};
}
std::vector<uint32_t> PrimeNumer::lshift_bits(const std::vector<uint32_t>& v, unsigned s) {
    if (s == 0) {
        std::vector<uint32_t> out = v;
        trim(out);
        ensure_nonempty_zero(out);
        return out;
    }
    if (v.empty()) return std::vector<uint32_t>{0};
    // s must be in 1..31
    s &= 31u;
    std::vector<uint32_t> out;
    out.resize(v.size() + 1);
    uint64_t carry = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        uint64_t cur = (((uint64_t)v[i]) << s) | carry;
        out[i] = (uint32_t)cur;
        carry = cur >> 32;
    }
    out[v.size()] = (uint32_t)carry;
    trim(out);
    ensure_nonempty_zero(out);
    return out;
}
void PrimeNumer::rshift_bits_inplace(std::vector<uint32_t>& v, unsigned s) {
    if (s == 0) { trim(v); ensure_nonempty_zero(v); return; }
    s &= 31u;
    if (v.empty()) { v.push_back(0); return; }
    uint64_t carry = 0;
    for (size_t idx = v.size(); idx-- > 0;) {
        unsigned __int128 cur = ((unsigned __int128)carry << 32) | v[idx];
        uint64_t newv = (uint64_t)(cur >> s);
        uint64_t newcarry = (uint64_t)(cur & ((( (unsigned __int128)1 << s) ) - 1));
        v[idx] = (uint32_t)newv;
        carry = newcarry;
        if (idx == 0) break;
    }
    trim(v);
    ensure_nonempty_zero(v);
}
