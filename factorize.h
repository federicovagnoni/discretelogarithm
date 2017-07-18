//
// Created by federico on 18/07/17.
//

#ifndef DISCRETELOGARITHMS_FACTORIZE_H
#define DISCRETELOGARITHMS_FACTORIZE_H

void eratosthenes(mpz_t n, struct list *primelist);
struct list *isBsmooth(mpz_t n, mpz_t B, struct list *primelist);
void gcdext(mpz_t a_origin, mpz_t b_origin, mpz_t x, mpz_t y, mpz_t z);

#endif //DISCRETELOGARITHMS_FACTORIZE_H
