#ifndef PRIMENUMER_H
#define PRIMENUMER_H
#include <vector>
#include<cstdint>
class PrimeNumer
{
private:

    std::vector<uint32_t> prime2;
    std::vector<uint32_t>fixprime;
    std::vector<uint32_t>n_1;
    std::vector<uint32_t>yu;
public:
    PrimeNumer();
    void generateOdd();
    bool miller_rabin();
    std::vector<uint32_t> prime1;
    std::pair<std::vector<uint32_t>, std::vector<uint32_t>> div_mod_knuth(const std::vector<uint32_t>& a_in, const std::vector<uint32_t>& b_in);
    std::vector<uint32_t> findD(int &s);
    void div(std::vector<uint32_t>&d);
    std::vector<uint32_t> quickExp(std::vector<uint32_t> base,std::vector<uint32_t>&d,std::vector<uint32_t>&mod,std::vector<uint32_t>&mu);
    int compare(std::vector<uint32_t>&a,std::vector<uint32_t>&mod);
    std::vector<uint32_t> sub(const std::vector<uint32_t>&a,const std::vector<uint32_t>&mod);
    void mod(std::vector<uint32_t>&a,std::vector<uint32_t>&mod);
    void mul_mod(std::vector<uint32_t>&rsu,std::vector<uint32_t>&base,std::vector<uint32_t>&mod);
    std::vector<uint32_t> mul(std::vector<uint32_t>&a,std::vector<uint32_t>&b,std::vector<uint32_t>&mod,std::vector<uint32_t>&mu);
    std::vector<uint32_t> mul_barrett(std::vector<uint32_t>&a,std::vector<uint32_t>&b);
    std::vector<uint32_t> karatsuba(std::vector<uint32_t>&a,std::vector<uint32_t>&b,std::vector<uint32_t>&mod,std::vector<uint32_t>&mu);
    std::vector<uint32_t> add(std::vector<uint32_t>&a,std::vector<uint32_t>&b);
    void karatsuba_add(std::vector<uint32_t>&rsu,std::vector<uint32_t>&a,int offset);
    std::vector<uint32_t> compute_mu( std::vector<uint32_t>& mod);
    void barrett_mod(std::vector<uint32_t>& a, std::vector<uint32_t>& mod, std::vector<uint32_t> & mu);
    void trim(std::vector<uint32_t> &a) { while(a.size()>1 && a.back()==0) a.pop_back(); }
    std::vector<uint32_t> lshift_bits(const std::vector<uint32_t>& v, unsigned s);
    void rshift_bits_inplace(std::vector<uint32_t>& v, unsigned s);
    // judge(int i);

};

#endif // PRIMENUMER_H
