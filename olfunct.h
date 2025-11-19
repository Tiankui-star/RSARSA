#ifndef OLFUNCT_H
#define OLFUNCT_H
#include<vector>
#include<cstdint>
#include<iostream>
#include"primenumber.h"
class olfunct {
private:
    PrimeNumer *fun;

public:
    std::vector<uint32_t>d;
    struct BigInt{
        std::vector<uint32_t>nums;
        int flag;
        BigInt(){}
        BigInt(std::vector<uint32_t>t):nums(t),flag(1){}
    };
    std::vector<uint32_t>mu;
    std::vector<uint32_t>mumu;
    std::vector<uint32_t> pube;
    std::vector<uint32_t> olfunction;
    olfunct( PrimeNumer* prime){
        // fun.prime1={952204601,110};
        // fun.prime2={4511491};

        // fun.p1_n_1={952204600,110};
        // fun.p2_n_1={4511490};
        fun=prime;
        mu=fun->karatsuba(fun->prime1,fun->prime2);
        olfunction=fun->karatsuba(fun->p1_n_1,fun->p2_n_1);


        pube={17};
        mumu=fun->compute_mu(mu);
    }
    void init(PrimeNumer* prime){
        fun=prime;
        mu=fun->karatsuba(fun->prime1,fun->prime2);
        olfunction=fun->karatsuba(fun->p1_n_1,fun->p2_n_1);
        for(auto &t:fun->prime1) std::cout<<t<<' ';
        std::cout<<std::endl;
        for(auto &t:fun->prime2) std::cout<<t<<' ';
        std::cout<<std::endl;

        pube={17};
        mumu=fun->compute_mu(mu);
    }
    olfunct(){};
    void solve(){
        BigInt x;
        x.flag=1;
        BigInt y;
        y.flag=1;

        private_e(pube,olfunction,x,y);
        std::vector<uint32_t> barretmu = fun->compute_mu(olfunction);
        if(x.flag==1){

            fun->barrett_mod(x.nums,olfunction,barretmu);
            d=x.nums;
        }
        else{

            fun->barrett_mod(x.nums,olfunction,barretmu);
            d=fun->sub(olfunction,x.nums);
        }
    }
    void private_e(const std::vector<uint32_t>&a,const std::vector<uint32_t>&b,BigInt &x,BigInt&y){
        if(b==std::vector<uint32_t>{0}) {
            x.nums={1};
            x.flag=1;
            y.nums={0};
            y.flag=1;
            return;
        }
        BigInt x1;
        BigInt y1;
        std::vector<uint32_t>tmod;
        tmod=mod(a,b);
        private_e(b,tmod,x1,y1);
        x=y1;
        BigInt tmp;
        tmp.nums=fun->mul(div(a,b),y1.nums);
        tmp.flag=y1.flag;
        y=sub(x1,tmp);

    }

    std::vector<uint32_t>div(const std::vector<uint32_t>&a,const std::vector<uint32_t>&b){
        if(fun->compare(a,b)==-1) return {0};
        std::vector<uint32_t>rsu=a;
        int n=a.size();
        uint32_t chu=b[0];
        uint64_t carry=0;
        for(int i=n-1;i>=0;i--){
            uint64_t tmp = (carry<<32)|a[i];
            rsu[i]=tmp/chu;
            carry=tmp%chu;
        }
        fun->trim(rsu);
        return rsu;

    }
    std::vector<uint32_t> mod(const std::vector<uint32_t>&a,const std::vector<uint32_t>&b){
        if(b.size()==1){
            uint32_t t=b[0];
            uint64_t remind=0;
            for(int i=a.size()-1;i>=0;i--){
                remind=(remind<<32)|a[i];
                remind%=t;
            }
            return std::vector<uint32_t>{(uint32_t)remind};
        }
        else{
            return a;
        }


    }
    BigInt sub(const BigInt&a,const BigInt&b){

        BigInt res;
        if(a.flag!=b.flag){
            res.flag=a.flag;
            res.nums=fun->add(a.nums,b.nums);
        }
        else{
            if(fun->compare(a.nums,b.nums)==0){
                res.flag=1;
                res.nums={0};
            }
            else if(fun->compare(a.nums,b.nums)==1){
                res.flag=a.flag;
                res.nums=fun->sub(a.nums,b.nums);
            }
            else{
                res.flag=!b.flag;
                res.nums=fun->sub(b.nums,a.nums);
            }
        }
        fun->trim(res.nums);
        return res;
    }


};

#endif // OLFUNCT_H
