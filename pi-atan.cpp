#include<iostream>
#include<gmp.h>
#include<math.h>
#include<cmath>
#include<string.h>
#include<chrono>

void arctan2(mpf_t ret, long long int n, unsigned int d, long long int numDigits){
    int bitCnt = numDigits * std::log2(10);

    mpf_init2(ret, bitCnt);

    mpf_t dsq;
    mpf_init2(dsq, bitCnt);
    mpf_set_ui(dsq, d);
    mpf_mul_ui(dsq, dsq, d);
    // std::cout << d << std::endl;
    // gmp_printf("%.50Ff", dsq);
    // std::cout << std::endl;


    mpf_t term;
    mpf_init2(term, bitCnt);
    mpf_init2(ret, bitCnt);
    mpf_set_ui(term, 1);
    mpf_div_ui(term, term, d);
    // std::cout << bitCnt << std::endl;
    // gmp_printf("%.50Ff", term);
    std::cout << std::endl;

    mpf_t total;
    mpf_init2(total, bitCnt);
    mpf_set(total, term);

    mpf_t totalTerm;
    mpf_init2(totalTerm, bitCnt);
    mpf_set(totalTerm, term);

    // gmp_printf("%.*Ff", 2*numDigits, total);
    // std::cout << std::endl << std::endl;

    mpf_t error;
    mpf_init2(error, bitCnt);
    mpf_set_ui(error, 1);
    mpf_div_2exp(error, error, bitCnt);

    // gmp_printf("%.50Ff", error);
    // std::cout << std::endl;

    long long int d2 = 0;
    mpf_t absTerm;
    mpf_init2(absTerm, bitCnt);
    mpf_set(absTerm, term);
    while(mpf_cmp(absTerm, error) >= 0){
        d2 += 1;
        mpf_mul_ui(term, term, n*n);
        mpf_div(term, term, dsq);
        mpf_neg(term, term);

        mpf_set(totalTerm, term);
        mpf_div_ui(totalTerm, totalTerm, 2*d2 + 1);
        // gmp_printf("%.*Ff", 2*numDigits, totalTerm);
        // std::cout << std::endl;
        mpf_add(total, total, totalTerm);
        // gmp_printf("%.*Ff", 2*numDigits, total);
        // std::cout << std::endl << std::endl;
        mpf_abs(absTerm, term);
        if(d2 % 10000 == 0){
            std::cout << d2 << " ";
        }
    }
    std::cout << "arctan of " << n << " / " << d << " took " << d2 << " iterations." << std::endl;

    // gmp_printf("%.*Ff", numDigits, total);
    // std::cout << std::endl;
    mpf_clear(dsq);
    mpf_clear(term);
    mpf_clear(totalTerm);
    mpf_set(ret, total);
    mpf_clear(total);
    mpf_clear(absTerm);
}

void arctan(mpz_t ret, long long int n, long long int d, long long int numDigits){

    int bitCnt = numDigits * std::log2(10);

    mpz_init2(ret, bitCnt);

    mpz_t dsq;
    mpz_init2(dsq, bitCnt);
    mpz_set_ui(dsq, d);
    mpz_mul_ui(dsq, dsq, d);
    // std::cout << d << std::endl;
    // gmp_printf("%.50Ff", dsq);
    // std::cout << std::endl;


    mpz_t term;
    mpz_init2(term, bitCnt);
    mpz_init2(ret, bitCnt);

    mpz_ui_pow_ui(term, 10, numDigits);//term = 10^numDigits
    mpz_mul_ui(term, term, n);
    // std::cout <<"d: " << d << std::endl;
    mpz_fdiv_q_ui(term, term, d);//term = 10^numDigits * n // d
    // gmp_printf("%Zd", term);
    // std::cout << std::endl;

    mpz_t total;
    mpz_init2(total, bitCnt);
    mpz_set(total, term);

    mpz_t totalTerm;
    mpz_init2(totalTerm, bitCnt);
    mpz_set(totalTerm, term);
 
    long long int d2 = 0;
    // gmp_printf("%Zd", term);
    // std::cout << std::endl;
    while(mpz_sgn(term) != 0){
        d2 += 1;
        // term = term * (-n*n)//(d*d)
        mpz_neg(term, term);
        mpz_fdiv_q(term, term, dsq);
        // total += term // (2*d2 + 1)
        mpz_fdiv_q_ui(totalTerm, term, 2*d2 + 1);
        mpz_add(total, total, totalTerm);
        //std::cout << d2 << std::endl;
        // gmp_printf("%Zd", term);
        // std::cout << std::endl;
        if(d2 % 10000 == 0){
            std::cout << d2 << " ";
        }
    }
    std::cout << std::endl << std::endl << "arctan of " << n << " / " << d << " took " << d2 << " iterations." << std::endl;
    mpz_clear(term);
    mpz_clear(totalTerm);
    mpz_init_set(ret, total);
    mpz_clear(total);
    mpz_clear(dsq);
}

void arctan3(mpf_t ret, int n, long long int d, long long int numDigits){

    int bitCnt = numDigits * std::log2(10);

    mpf_init2(ret, bitCnt);

    mpz_t dz;
    mpz_init_set_ui(dz, d);

    mpq_t dsq;
    mpq_init(dsq);
    mpq_set_ui(dsq, 1, 1);
    mpq_set_den(dsq, dz);
    mpq_mul(dsq, dsq, dsq);
    // std::cout << d << std::endl;
    // gmp_printf("%.50Ff", dsq);
    // std::cout << std::endl;

    mpq_t odd;
    mpq_init(odd);
    mpq_set_ui(odd, 1, 1);

    mpq_t term;
    mpq_init(term);
    mpq_set_ui(term, 1, 1);
    mpq_set_den(term, dz);

    mpq_t total;
    mpq_init(total);
    mpq_set(total, term);

    mpq_t totalTerm;
    mpq_init(totalTerm);
    mpq_set(totalTerm, term);

    mpq_t epsilon;
    mpq_init(epsilon);
    mpq_set_ui(epsilon, 1, 1);
    mpz_t epsilonDenom;
    mpz_init(epsilonDenom);
    mpz_ui_pow_ui(epsilonDenom, 10, numDigits);
    mpq_set_den(epsilon, epsilonDenom);
    mpz_clear(epsilonDenom);

    mpq_t two;
    mpq_init(two);
    mpq_set_ui(two, 2, 1);

    bool isNeg = false;
    long long int d2 = 0;
    // gmp_printf("%Zd", term);
    // std::cout << std::endl;
    while(true){
        d2 += 1;
        mpq_add(odd, odd, two);
        // term = term * (-n*n)//(d*d)
        mpq_neg(term, term);
        isNeg = !isNeg;
        mpq_mul(term, term, dsq);
        // total += term // (2*d2 + 1)
        mpq_div(totalTerm, term, odd);
        mpq_add(total, total, totalTerm);
        //std::cout << d2 << std::endl;
        // gmp_printf("%Zd", term);
        // std::cout << std::endl;
        if(isNeg){
            mpq_neg(term, term);
            if(mpq_cmp(term, epsilon) < 0) {
                break;
            }
            mpq_neg(term, term);
        } else if(mpq_cmp(term, epsilon) < 0) {
            break;
        }
        if(d2 % 1000 == 0){
            std::cout << d2 << " ";
            mpq_canonicalize(total);
        }
    }

    mpz_t num;
    mpz_init(num);
    mpq_get_num(num, total);
    mpz_t den;
    mpz_init(den);
    mpq_get_den(den, total);
    mpf_set_q(ret, total);
    mpz_clear(num);
    mpz_clear(den);
    std::cout << std::endl << std::endl << "arctan of " << n << " / " << d << " took " << d2 << " iterations." << std::endl;
    mpz_clear(dz);
    mpq_clear(term);
    mpq_clear(totalTerm);
    mpq_clear(total);
    mpq_clear(epsilon);
    mpq_clear(dsq);
    mpq_clear(odd);
    mpq_clear(two);
}


void calculatePi(mpf_t ret, long long int numDigits){
    int bitCnt = numDigits * std::log2(10);
    mpf_t pi;
    mpf_init2(pi, bitCnt);
    mpf_set_ui(pi, 1);

    mpf_init2(ret, bitCnt);
    mpf_set_ui(ret, 1);

    mpf_t a5;
    arctan2(a5, 1, 5, numDigits);

    mpf_t a239;
    arctan2(a239, 1, 239, numDigits);

        
    // gmp_printf("%.*Ff", numDigits, a5);
    // std::cout << std::endl;
        
    // gmp_printf("%.*Ff", numDigits, a239);
    // std::cout << std::endl;

    // pi = 4 * (4 * a5 - a239)
    mpf_set(pi, a5);
    mpf_mul_ui(pi, pi, 4);
    mpf_sub(pi, pi, a239);
    mpf_mul_ui(pi, pi, 4);

    // gmp_printf("%.*Ff", numDigits, pi);
    // std::cout << std::endl;

    mpf_set(ret, pi);

    mpf_clear(pi);
    mpf_clear(a5);
    mpf_clear(a239);
}

void calculatePiInt(mpz_t ret, long long int numDigits){
    mpz_t pi;
    mpz_init_set_ui(pi, 1);

    mpz_init_set_ui(ret, 1);

    mpz_t a5;
    arctan(a5, 1, 5, numDigits);

    mpz_t a239;
    arctan(a239, 1, 239, numDigits);

        
    // gmp_printf("%.*Ff", numDigits, a5);
    // std::cout << std::endl;
        
    // gmp_printf("%.*Ff", numDigits, a239);
    // std::cout << std::endl;

    // pi = 4 * (4 * a5 - a239)
    mpz_set(pi, a5);
    mpz_mul_ui(pi, pi, 4);
    mpz_sub(pi, pi, a239);
    mpz_mul_ui(pi, pi, 4);

    // gmp_printf("%.*Ff", numDigits, pi);
    // std::cout << std::endl;

    mpz_set(ret, pi);

    mpz_clear(pi);
    mpz_clear(a5);
    mpz_clear(a239);
}

void calculatePi2(mpf_t ret, long long int numDigits){
    int bitCnt = numDigits * std::log2(10);
    mpf_t pi;
    mpf_init2(pi, bitCnt);
    mpf_set_ui(pi, 1);

    mpf_init2(ret, bitCnt);
    mpf_set_ui(ret, 1);

    mpf_t a239;
    arctan2(a239, 1, 239, numDigits);
    mpf_mul_ui(a239, a239, 183);

    mpf_t a1023;
    arctan2(a1023, 1, 1023, numDigits);
    mpf_mul_ui(a1023, a1023, 32);

    mpf_t a5832;
    arctan2(a5832, 1, 5832, numDigits);
    mpf_mul_ui(a5832, a5832, 68);

    mpf_t a110443;
    arctan2(a110443, 1, 110443, numDigits);
    mpf_mul_ui(a110443, a110443, 12);

    mpf_t a4841182;
    arctan2(a4841182, 1, 4841182, numDigits);
    mpf_mul_ui(a4841182, a4841182, 12);

    mpf_t a6826318;
    arctan2(a6826318, 1, 6826318, numDigits);    
    mpf_mul_ui(a6826318, a6826318, 100);
        
    // gmp_printf("%.*Ff", numDigits, a5);
    // std::cout << std::endl;
        
    // gmp_printf("%.*Ff", numDigits, a239);
    // std::cout << std::endl;

    mpf_set(pi, a239);
    mpf_add(pi, pi, a1023);
    mpf_sub(pi, pi, a5832);
    mpf_add(pi, pi, a110443);
    mpf_sub(pi, pi, a4841182);
    mpf_sub(pi, pi, a6826318);
    mpf_mul_ui(pi, pi, 4);


    // gmp_printf("%.*Ff", numDigits, pi);
    // std::cout << std::endl;

    mpf_set(ret, pi);

    mpf_clear(pi);
    mpf_clear(a239);
    mpf_clear(a1023);
    mpf_clear(a5832);
    mpf_clear(a110443);
    mpf_clear(a4841182);
    mpf_clear(a6826318);
}


void calculatePi2Int(mpz_t ret, long long int numDigits){
    int extraDigits = 1000;
    numDigits += extraDigits;
    mpz_t pi;
    mpz_init_set_ui(pi, 1);

    mpz_init_set_ui(ret, 1);

    mpz_t a239;
    arctan(a239, 1, 239, numDigits);
    mpz_mul_ui(a239, a239, 183);

    mpz_t a1023;
    arctan(a1023, 1, 1023, numDigits);
    mpz_mul_ui(a1023, a1023, 32);

    mpz_t a5832;
    arctan(a5832, 1, 5832, numDigits);
    mpz_mul_ui(a5832, a5832, 68);

    mpz_t a110443;
    arctan(a110443, 1, 110443, numDigits);
    mpz_mul_ui(a110443, a110443, 12);

    mpz_t a4841182;
    arctan(a4841182, 1, 4841182, numDigits);
    mpz_mul_ui(a4841182, a4841182, 12);

    mpz_t a6826318;
    arctan(a6826318, 1, 6826318, numDigits);    
    mpz_mul_ui(a6826318, a6826318, 100);
        
    // gmp_printf("%.*Ff", numDigits, a5);
    // std::cout << std::endl;
        
    // gmp_printf("%.*Ff", numDigits, a239);
    // std::cout << std::endl;

    mpz_set(pi, a239);
    mpz_add(pi, pi, a1023);
    mpz_sub(pi, pi, a5832);
    mpz_add(pi, pi, a110443);
    mpz_sub(pi, pi, a4841182);
    mpz_sub(pi, pi, a6826318);
    mpz_mul_ui(pi, pi, 4);

    mpz_t extra;
    mpz_init_set_ui(extra, 1);
    mpz_ui_pow_ui(extra, 10, extraDigits);

    mpz_fdiv_q(pi, pi, extra);

    // gmp_printf("%.*Ff", numDigits, pi);
    // std::cout << std::endl;

    mpz_set(ret, pi);

    mpz_clear(extra);
    mpz_clear(pi);
    mpz_clear(a239);
    mpz_clear(a1023);
    mpz_clear(a5832);
    mpz_clear(a110443);
    mpz_clear(a4841182);
    mpz_clear(a6826318);
}


int main(){
    bool isFloat = true;
    long long int numDigits;
    std::cout << "Enter number of digits of pi to find:" << std::endl;
    std::cin >> numDigits;
    std::cout << std::endl << std::endl;
    if(isFloat){
        mpf_t pi;
        
        auto begin = std::chrono::steady_clock::now();

        calculatePi2(pi, numDigits + 1);

        auto end = std::chrono::steady_clock::now();
        int seconds = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() % 60;
        int minutes = std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() % 60;
        int hours = std::chrono::duration_cast<std::chrono::hours>(end - begin).count();

        // gmp_printf("%Zd", pi);
        // std::cout << std::endl;

        std::cout << std::endl << "Calculated " << numDigits << " of pi in " << hours << " hours, " << minutes << " minutes, " << seconds << "seconds" << std::endl;

        std::string s = std::to_string(numDigits).append("Floatcpp.txt");

        FILE *file;
        file = fopen(s.c_str(), "wt");
        gmp_fprintf(file, "%.*Ff", numDigits, pi);
        fclose(file);

        mpf_clear(pi);
    } else {
        mpz_t pi;
        
        auto begin = std::chrono::steady_clock::now();

        calculatePi2Int(pi, numDigits);

        auto end = std::chrono::steady_clock::now();
        int seconds = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() % 60;
        int minutes = std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() % 60;
        int hours = std::chrono::duration_cast<std::chrono::hours>(end - begin).count();

        // gmp_printf("%Zd", pi);
        // std::cout << std::endl;

        std::cout << std::endl << "Calculated " << numDigits << " of pi in " << hours << " hours, " << minutes << " minutes, " << seconds << "seconds" << std::endl;

        std::string s = std::to_string(numDigits).append("cpp.txt");

        FILE *file;
        file = fopen(s.c_str(), "wt");
        gmp_fprintf(file, "%Zd", pi);
        fclose(file);

        mpz_clear(pi);
    }
}
