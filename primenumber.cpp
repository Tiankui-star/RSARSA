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
int PrimeNumer::compare(std::vector<uint32_t>&a,std::vector<uint32_t>&mod){
    int as=a.size();
    int bs=mod.size();
    int i=as,j=bs;
    while(i>1&&a[i-1]==0) --i;
    while(j>1&&a[j-1]==0) --j;
    if(i<j) return -1;
    if(i>j) return 1;
    for(int k=i;k>0;k--){
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
std::vector<uint32_t> PrimeNumer::karatsuba(std::vector<uint32_t>&a,std::vector<uint32_t>&b){
    // int n=std::max(a.size(),b.size());
    // if(n==0) return std::vector<uint32_t>{0};
    // if(n<=64) return mul(a,b);
    // std::vector<uint32_t>res;
    // int mid=n/2;
    // std::vector<uint32_t> a0(std::min(mid,(int)a.size()));
    // std::copy(a.begin(),a.begin()+std::min((size_t)mid,a.size()),a0.begin());
    // std::vector<uint32_t> a1;
    // if(a.size()>mid) a1.assign(a.begin()+mid,a.end());
    // else a1.clear();

    // std::vector<uint32_t> b0(std::min(mid,(int)b.size()));
    // std::copy(b.begin(),b.begin()+std::min((size_t)mid,b.size()),b0.begin());
    // std::vector<uint32_t> b1;
    // if(b.size()>mid) b1.assign(b.begin()+mid,b.end());
    // else b1.clear();
    // // std::vector<uint32_t>a0(a.begin(),a.begin()+mid);
    // // std::vector<uint32_t>a1(a.begin()+mid,a.end());
    // // std::vector<uint32_t>b0(b.begin(),b.begin()+mid);
    // // std::vector<uint32_t>b1(b.begin()+mid,b.end());
    // std::vector<uint32_t>r1=karatsuba(a0,b0); //ac
    // std::vector<uint32_t>r2=karatsuba(a1,b1); //bd
    // std::vector<uint32_t>adda=add(a0,a1);
    // std::vector<uint32_t>addb=add(b0,b1);
    // std::vector<uint32_t>t=karatsuba(adda,addb); //ac+ad+bc+bd
    // std::vector<uint32_t>t1=sub(t,r1);
    // std::vector<uint32_t>r3=sub(t1,r2); //ad+bc
    // karatsuba_add(res,r1,0);
    // karatsuba_add(res,r3,mid);
    // karatsuba_add(res,r2,2*mid);
    // trim(res);
    // return res;
    return karatsuba_rec(a.data(),a.size(),b.data(),b.size());
}
std::vector<uint32_t> PrimeNumer::karatsuba_rec(const uint32_t *a, int an, const uint32_t *b, int bn){
    while (an > 1 && a[an-1] == 0) --an;
    while (bn > 1 && b[bn-1] == 0) --bn;
    if (an == 0 || bn == 0) return std::vector<uint32_t>{0};
    size_t n = std::max(an, bn);
    if(n<=64){
        std::vector<uint32_t> va(a, a + an);
        std::vector<uint32_t> vb(b, b + bn);
        return mul(va, vb);
    }
    int mid=n/2;
    size_t a0_len = (an < mid) ? an : mid;
    size_t a1_len = (an > mid) ? (an - mid) : 0;
    const uint32_t *a0_ptr = a;
    const uint32_t *a1_ptr = (a1_len ? (a + mid) : nullptr);
    size_t b0_len = (bn < mid) ? bn : mid;
    size_t b1_len = (bn > mid) ? (bn - mid) : 0;
    const uint32_t *b0_ptr = b;
    const uint32_t *b1_ptr = (b1_len ? (b + mid) : nullptr);
    // r1 = a0 * b0
    std::vector<uint32_t> r1 = karatsuba_rec(a0_ptr, a0_len, b0_ptr, b0_len);
    // r2 = a1 * b1
    std::vector<uint32_t> r2 = karatsuba_rec(a1_ptr ? a1_ptr : (const uint32_t*)"\0", a1_len, b1_ptr ? b1_ptr : (const uint32_t*)"\0", b1_len);
    // adda = a0 + a1
    std::vector<uint32_t> adda;
    if (a1_len == 0) adda.assign(a0_ptr, a0_ptr + a0_len);
    else {
        std::vector<uint32_t> va0(a0_ptr, a0_ptr + a0_len);
        std::vector<uint32_t> va1(a1_ptr, a1_ptr + a1_len);
        adda = add(va0, va1);
    }
    // addb = b0 + b1
    std::vector<uint32_t> addb;
    if (b1_len == 0) addb.assign(b0_ptr, b0_ptr + b0_len);
    else {
        std::vector<uint32_t> vb0(b0_ptr, b0_ptr + b0_len);
        std::vector<uint32_t> vb1(b1_ptr, b1_ptr + b1_len);
        addb = add(vb0, vb1);
    }
    // t = (a0+a1)*(b0+b1)
    std::vector<uint32_t> t = karatsuba_rec(adda.data(), adda.size(), addb.data(), addb.size());
    // r3 = t - r1 - r2  => r3 must be non-negative and equals ad+bc
    std::vector<uint32_t> t_minus_r1 = sub(t, r1);
    std::vector<uint32_t> r3 = sub(t_minus_r1, r2);

    // assemble result: r1 + (r3 << mid) + (r2 << 2*mid)
    std::vector<uint32_t> res;
    karatsuba_add(res, r1, 0);
    karatsuba_add(res, r3, mid);
    karatsuba_add(res, r2, 2 * mid);
    trim(res);
    return res;

}
void PrimeNumer::karatsuba_add(std::vector<uint32_t>&rsu,std::vector<uint32_t>&a,int offset){
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
void PrimeNumer::barrett_mod(std::vector<uint32_t>& a, std::vector<uint32_t>& mod, std::vector<uint32_t> & mu){

    trim(a);
    trim(mod);
    const size_t k = mod.size();            // number of limbs in modulus
    if (k == 0) return;                     // defensive
    if (compare(a, mod) < 0) return;        // already reduced

    // If a is very large (>> 2k) we can first truncate to at most 2k limbs:
    // Barrett assumes a < b^(2k). If a >= b^(2k), we reduce top part using division or iterate.
    // A simple safe pre-trim: if a has > 2k limbs, compute remainder of high limbs by dividing:
    if (a.size() > 2 * k) {
        // Efficient strategy: split a = a_low + a_high * b^(2k)
        // We can reduce a_high mod mod via div_mod_knuth of the top part.
        // For simplicity & correctness: call div_mod_knuth on entire 'a' once (rare if a oversized).
        // This ensures a < b^(2k) afterwards.
        auto qr = div_mod_knuth(a, mod);
        a = qr.second;
        trim(a);
        if (compare(a, mod) < 0) return;
    }

    // Now ensure a.size() <= 2k (or at least not much larger).
    // Step 1: q1 = floor(a / b^(k-1)) -> take high words starting at index (k-1)
    std::vector<uint32_t> q1;
    if (a.size() > (size_t)(k - 1)) {
        q1.assign(a.begin() + (k - 1), a.end());
    } else {
        q1 = {0};
    }
    trim(q1);

    // Step 2: q2 = q1 * mu (use karatsuba/mul)
    std::vector<uint32_t> q2 = karatsuba(q1, mu); // or mul(q1, mu)

    // Step 3: q3 = floor(q2 / b^(k+1)) -> take high words starting at index (k+1)
    std::vector<uint32_t> q3;
    if (q2.size() > (size_t)(k + 1)) {
        q3.assign(q2.begin() + (k + 1), q2.end());
    } else {
        q3 = {0};
    }
    trim(q3);

    // Step 4: r = a - q3 * mod
    std::vector<uint32_t> q3m = karatsuba(q3, mod); // don't use barrett here
    trim(q3m);

    // If q3m <= a then r = a - q3m, else approximation failed -> fallback to div_mod_knuth
    std::vector<uint32_t> r;
    if (compare(a, q3m) >= 0) {
        r = sub(a, q3m);
        trim(r);
    } else {
        // Very rare: approximation produced q3 too large; fallback to exact division for safety.
        auto qr = div_mod_knuth(a, mod);
        a = qr.second;
        trim(a);
        return;
    }

    // Step 5: At most two subtractions of mod are necessary (theory: q3 <= floor(a/m)+2)
    // Use at most two if-checks (not while) for speed and safety.
    if (compare(r, mod) >= 0) {
        r = sub(r, mod);
        trim(r);
        if (compare(r, mod) >= 0) {
            r = sub(r, mod);
            trim(r);
        }
    }

    // Final write-back
    a.swap(r);
    trim(a);
}
std::vector<uint32_t> PrimeNumer::add(std::vector<uint32_t>&a,std::vector<uint32_t>&b){
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
std::vector<uint32_t> PrimeNumer::compute_mu(std::vector<uint32_t>& m){
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

    // copy inputs and normalize
    std::vector<uint32_t> a = a_in;
    std::vector<uint32_t> b = b_in;
    trim(a); trim(b);
    ensure_nonempty_zero(a); ensure_nonempty_zero(b);

    // trivial cases
    if (b.size() == 1) {
        // short division
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
        // quotient = 0, remainder = a
        std::vector<uint32_t> q(1, 0);
        trim(a);
        ensure_nonempty_zero(a);
        return {q, a};
    }

    size_t n = b.size();
    size_t m = (a.size() > n) ? (a.size() - n) : 0;

    // normalization: shift so highest bit of b.back() is 1
    uint32_t bn = b.back();
    unsigned shift = 0;
    if (bn == 0) { // safety - shouldn't happen after trim
        trim(b);
        bn = b.back();
    }
    // find shift in 0..31 such that (bn << shift) has high bit set
    while (shift < 32 && (((uint32_t)(bn << shift)) & 0x80000000u) == 0u) ++shift;
    // normalize b and a
    std::vector<uint32_t> b_norm = lshift_bits(b, shift);
    std::vector<uint32_t> a_norm = lshift_bits(a, shift);
    // ensure a_norm has length >= n + m + 1
    if (a_norm.size() < n + m + 1) a_norm.resize(n + m + 1, 0u);

    std::vector<uint32_t> q(m + 1, 0u);

    // main loop: j = m .. 0
    for (size_t j = m + 1; j-- > 0;) {
        // estimate qhat = (a_norm[j+n]*b + a_norm[j+n-1]) / b_norm[n-1]
        unsigned __int128 top = ((unsigned __int128)a_norm[j + n] << 32) | a_norm[j + n - 1];
        unsigned __int128 den = b_norm[n - 1];
        unsigned __int128 qhat = top / den;
        if (qhat > 0xFFFFFFFFull) qhat = 0xFFFFFFFFull;

        // correct qhat: ensure qhat*b_norm[n-2] <= (top - qhat*den)*b + a_norm[j+n-2]
        while (true) {
            // compute left = qhat * b_norm[n-2]
            unsigned __int128 left = qhat * (unsigned __int128)b_norm[n - 2];
            unsigned __int128 remhat = top - qhat * den; // remainder of dividing top by den
            unsigned __int128 right = (remhat << 32) + a_norm[j + n - 2];
            if (left <= right) break;
            --qhat;
        }

        // multiply-subtract: a_norm[j .. j+n] -= qhat * b_norm[0..n-1]
        unsigned __int128 borrow = 0;
        for (size_t i = 0; i < n; ++i) {
            unsigned __int128 prod = qhat * (unsigned __int128)b_norm[i];
            unsigned __int128 ai = (unsigned __int128)a_norm[j + i];
            unsigned __int128 sub = ai - (prod & 0xFFFFFFFFull) - borrow;
            a_norm[j + i] = (uint32_t)sub;
            // new borrow: prod_high + ((prod_low + borrow) > ai ? 1 : 0)
            unsigned __int128 prod_high = (prod >> 32);
            // Determine if (ai < (prod_low + borrow))
            unsigned __int128 prod_low_plus_borrow = (prod & 0xFFFFFFFFull) + borrow;
            borrow = prod_high + ((ai < prod_low_plus_borrow) ? 1 : 0);
        }
        // last subtraction for the high word
        unsigned __int128 ai_last = (unsigned __int128)a_norm[j + n];
        unsigned __int128 sub_last = ai_last - borrow;
        a_norm[j + n] = (uint32_t)sub_last;
        bool negative = ((sub_last >> 63) & 1); // if high bit set, it wrapped (negative in two's complement sense)

        if (negative) {
            // qhat was too big, add back b_norm to a_norm[j..j+n-1]
            --qhat;
            uint64_t carry = 0;
            for (size_t i = 0; i < n; ++i) {
                unsigned __int128 sum = (unsigned __int128)a_norm[j + i] + (unsigned __int128)b_norm[i] + carry;
                a_norm[j + i] = (uint32_t)sum;
                carry = (uint64_t)(sum >> 32);
            }
            // add carry to a_norm[j + n]
            unsigned __int128 new_high = (unsigned __int128)a_norm[j + n] + carry;
            a_norm[j + n] = (uint32_t)new_high;
        }
        q[j] = (uint32_t)qhat;
    }

    // remainder r = a_norm[0..n-1] >> shift
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
