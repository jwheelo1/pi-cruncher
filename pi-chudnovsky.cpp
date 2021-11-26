#include<iostream>
#include<gmp.h>
#include<math.h>
#include<cmath>
#include<string.h>
#include<chrono>

mpf_t ret;

mpz_t Pab;
mpz_t Qab;
mpz_t Tab;
mpz_t x3;
mpz_t t;
mpz_t gcd;

void bs(long int a, long int b, int level, bool showAll){
    if(b - a == 1){
        if( a == 0 ){
            mpz_set_ui(Pab, 1);
            mpz_set_ui(Qab, 1);
        } else  {
            mpz_set_ui(Pab, 6*a-5);
            mpz_mul_ui(Pab, Pab, 2*a-1);
            mpz_mul_ui(Pab, Pab, 6*a-1);
            mpz_set_ui(Qab, a);
            mpz_mul_ui(Qab, Qab, a);
            mpz_mul_ui(Qab, Qab, a);
            mpz_mul(Qab, Qab, x3);
        }

        mpz_set_ui(t, a);

        mpz_mul_ui(t, t, 545140134);
        mpz_add_ui(t, t, 13591409);
        mpz_mul(Tab, Pab, t);
        if( b % 2 ){
            mpz_neg(Tab, Tab);
        }

    } else {
        long int m = a+(b-a)*0.5224;

        if(level <= 6 || showAll){
            std::string spaces(level + 2, ' ');
            std::cout << spaces << "at level " << level << std::endl;
        }

        bs(a, m, level + 1, false);

        mpz_t Pam;
        mpz_t Qam;
        mpz_t Tam;
        mpz_init_set(Pam, Pab);
        mpz_init_set(Qam, Qab);
        mpz_init_set(Tam, Tab);

        
        if(level == 0){
            showAll = true;
        }

        bs(m, b, level + 1, showAll);

        if(level >= 4){
            mpz_gcd(gcd, Qab, Pam);
            mpz_divexact(Qab, Qab, gcd);
            mpz_divexact(Pam, Pam, gcd);
        }

        mpz_mul(Pab, Pam, Pab);
        mpz_mul(Tam, Tam, Qab);
        mpz_mul(Qab, Qam, Qab);
        mpz_mul(Pam, Pam, Tab);
        mpz_add(Tab, Tam, Pam);

        mpz_clear(Pam);
        mpz_clear(Qam);
        mpz_clear(Tam);
        if(level <= 6 || showAll){
            std::string spaces(level + 2, ' ');
            std::cout << spaces << "done" << std::endl;
        }
    }
}


void getSumBS(long int numDigits){
    long long int bitCnt = ceil(numDigits * std::log2(10));
    
    std::cout << " before inits " << std::endl;
    mpz_init(x3);
    mpz_set_ui(x3, 640320);
    mpz_mul_ui(x3, x3, 640320);
    mpz_mul_ui(x3, x3, 640320);
    mpz_fdiv_q_ui(x3, x3, 24);
    std::cout << " after first inits " << std::endl;
    mpz_init(Pab);
    mpz_init(Qab);
    mpz_init(Tab);
    mpz_init_set_ui(gcd, 1);
    
    mpz_init(t);

    std::cout << " after second inits " << std::endl;

    int N = floor(numDigits / 14.181647462725477655525521678181770863769125289828726959816854332 + 1);

    bs(0, N, 0, false);

    mpz_clear(Pab);
    mpz_clear(x3);
    mpz_clear(t);

    mpz_gcd(gcd, Qab, Tab);
    // gmp_printf("%Zd\n", Qab);
    mpz_divexact(Qab, Qab, gcd);
    mpz_divexact(Tab, Tab, gcd);
    // gmp_printf("%Zd\n", Qab);
    size_t sq = mpz_sizeinbase(Qab, 2);
    size_t st = mpz_sizeinbase(Tab, 2);
    std::cout << sq << "    " << st << std::endl;

    mpz_clear(gcd);

    std::cout << " after bs " << std::endl;
    mpf_t Tab_f;
    mpf_init2(Tab_f, st+1);
    mpf_set_z(Tab_f, Tab);
    std::cout << "init Tab_f" << std::endl;

    mpz_clear(Tab);

    mpf_t Qab_f;
    mpf_init2(Qab_f, sq+1);
    mpf_set_z(Qab_f, Qab);
    std::cout << "init Qab_f" << std::endl;

    mpz_clear(Qab);

    mpf_div(Qab_f, Qab_f, Tab_f);
    std::cout << "qab = qab / tab" << std::endl;
    mpf_clear(Tab_f);

    mpf_mul(ret, ret, Qab_f);
    std::cout << "ret = Tab_f / Qab_f" << std::endl;
    mpf_neg(ret, ret);
    std::cout << "neg ret" << std::endl;
    mpf_clear(Qab_f);
    std::cout << "clears" << std::endl;
}

void getSum(mpf_t ret, long long int numDigits){
    int bitCnt = numDigits * std::log2(10);
    //int bitCnt=1000000;

    long long int q = 0;
    // initialize values to recurse on
    const char l_init[] = {'1','3','5','9','1','4','0','9','@','0',0};
    mpf_t L_q;
    mpf_init2(L_q, bitCnt);
    mpf_set_str(L_q, l_init, 10);
    mpf_t X_q;
    mpf_init2(X_q, bitCnt);
    mpf_set_ui(X_q, 1);
    mpf_t M_q;
    mpf_init2(M_q, bitCnt);
    mpf_set_ui(M_q, 1);
    mpf_t K_q;
    mpf_init2(K_q, bitCnt);
    mpf_set_si(K_q, -6);
    
    const char l_step[] = {'5','4','5','1','4','0','1','3','4','@','0',-1};
    mpf_t L_step;
    mpf_init(L_step);
    mpf_set_str(L_step, l_step, 10);

    gmp_printf("L: %.20Ff", L_step);
    std::cout << std::endl;

    const char x_step[] = {'2','6','2','5','3','7','4','1','2','6','4','0','7','6','8','0','0','0','@','0',-1};
    mpf_t X_step;
    mpf_init(X_step);
    mpf_set_str(X_step, x_step, 10);
    mpf_neg(X_step, X_step);
    gmp_printf("X: %.20Ff", X_step);
    std::cout << std::endl << std::endl;

    mpf_t M_term_num;
    mpf_init2(M_term_num, bitCnt);

    mpf_t totalterm;
    mpf_init2(totalterm, bitCnt);
    mpf_set_ui(totalterm, 1);
    mpf_mul(totalterm, totalterm, M_q);
    mpf_mul(totalterm, totalterm, L_q);
    mpf_div(totalterm, totalterm, X_q);

    mpf_t total;
    mpf_init2(total, bitCnt);
    mpf_set(total, totalterm);

    mpf_set_ui(totalterm, 1);

    while(q < numDigits){
        q += 1;
        // L_q = L_q + L_step
        mpf_add(L_q, L_q, L_step);
        // X_q = X_q * X_step
        mpf_mul(X_q, X_q, X_step);
        // gmp_printf("L: %.20Ff, X: %.20Ff", L_q, X_q);
        // std::cout << std::endl;
        // M = M * (K^3 - 16K) / q^3
        mpf_add_ui(K_q, K_q, 12);

        mpf_set(M_term_num, K_q);
        mpf_mul(M_term_num, M_term_num, K_q);
        mpf_sub_ui(M_term_num, M_term_num, 16);
        mpf_mul(M_term_num, M_term_num, K_q);
        mpf_div_ui(M_term_num, M_term_num, q);
        mpf_div_ui(M_term_num, M_term_num, q);
        mpf_div_ui(M_term_num, M_term_num, q);
        mpf_mul(M_q, M_q, M_term_num);

        // term = (M_q * L_q) / X_q
        mpf_mul(totalterm, totalterm, M_q);
        // gmp_printf("t1: %.10Ff\n", totalterm);
        mpf_mul(totalterm, totalterm, L_q);
        // gmp_printf("t2: %.10Ff\n", totalterm);
        mpf_div(totalterm, totalterm, X_q);
        // gmp_printf("t3: %.100Ff\n", totalterm);

        //gmp_printf("%.10Ff", totalterm);

        mpf_add(total, total, totalterm);
        mpf_set_ui(totalterm, 1);

        // gmp_printf("%.*Ff", q, total);
        // std::cout << std::endl;
        if( q % 1000 == 0 ){
            std::cout << q << " ";
        }
    }

    std::cout << std::endl;

    mpf_set(ret, total);
    mpf_clear(L_q);
    mpf_clear(X_q);
    mpf_clear(M_q);
    mpf_clear(K_q);
    mpf_clear(L_step);
    mpf_clear(X_step);
    mpf_clear(M_term_num);
    mpf_clear(total);
    mpf_clear(totalterm);
}

int main() {
    long long int numDigits;
    std::cout << "Enter number of digits of pi to find:" << std::endl;
    std::cin >> numDigits;
    std::cout << std::endl << std::endl;

    long long int bitCnt = ceil(numDigits * std::log2(10));
    std::cout << numDigits << " " << bitCnt << std::endl;

    auto begin = std::chrono::steady_clock::now();

    mpf_init2(ret, bitCnt);
    mpf_sqrt_ui(ret, 10005);
    mpf_mul_ui(ret, ret, 426880);

    getSumBS(numDigits);
    std::cout << "after getSumBS" << std::endl;

    auto end = std::chrono::steady_clock::now();
    int seconds = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() % 60;
    int minutes = std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() % 60;
    int hours = std::chrono::duration_cast<std::chrono::hours>(end - begin).count();

    // gmp_printf("%Zd", pi);
    // std::cout << std::endl;

    std::cout << std::endl << "Calculated " << numDigits << " of pi in " << hours << " hours, " << minutes << " minutes, " << seconds << "seconds" << std::endl;

    std::string s = std::to_string(numDigits).append("Chudnovsky.txt");

    //gmp_printf("pi? %.*Ff", numDigits, pi);
    //std::cout << std::endl;

    FILE *file;
    file = fopen(s.c_str(), "wt");
    gmp_fprintf(file, "%.*Ff", numDigits, ret);
    fclose(file);

    mpf_clear(ret);
    return 0;
}
