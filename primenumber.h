#ifndef PRIMENUMER_H
#define PRIMENUMER_H
#include <vector>
#include<cstdint>
class PrimeNumer
{
private:
    std::vector<uint32_t> prime1;
    std::vector<uint32_t> prime2;
    std::vector<uint32_t>fixprime;
    std::vector<uint32_t>yu;
public:
    PrimeNumer();
    void generateOdd();
    bool miller_rabin();
    std::vector<uint32_t> findD(int &s);
    void div(std::vector<uint32_t>&d);
    std::vector<uint32_t> quickExp(std::vector<uint32_t> base,std::vector<uint32_t>&d,std::vector<uint32_t>&mod);
    int compare( std::vector<uint32_t>&a, std::vector<uint32_t>&mod);
    std::vector<uint32_t> sub( std::vector<uint32_t>&a, std::vector<uint32_t>&mod);
    void mod(std::vector<uint32_t>&a,std::vector<uint32_t>&mod);
    void mul_mod(std::vector<uint32_t>&rsu,std::vector<uint32_t>&base,std::vector<uint32_t>&mod);
    std::vector<uint32_t> mul(std::vector<uint32_t>&a,std::vector<uint32_t>&b);
    std::vector<uint32_t> karatsuba(std::vector<uint32_t>&a,std::vector<uint32_t>&b);
    std::vector<uint32_t> add(std::vector<uint32_t>&a,std::vector<uint32_t>&b);
    void karatsuba_add(std::vector<uint32_t>&rsu,std::vector<uint32_t>&a,int offset);

    // judge(int i);

};

#endif // PRIMENUMER_H
