//
// Created by federico on 18/07/17.
//

#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include "utils.h"


int limit_num = 1;
int limit_denom = 1;

void g(mpz_t result, mpz_t x, mpz_t a, mpz_t modulus) {
    mpz_t pow;
    mpz_init(pow);
    mpz_pow_ui(pow, x, 2);
    mpz_add(pow, pow, a);
    mpz_mod(pow, pow, modulus);
    mpz_set(result, pow);
    mpz_clear(pow);
    return;
}

void choose_seed(mpz_t x, mpz_t y, mpz_t a, mpz_t n, gmp_randstate_t *randstate) {
    mpz_urandomm(x, *randstate, n);
    mpz_set(y, x);
    mpz_urandomm(a, *randstate, n);
}


void pollard_rho_opt(mpz_t factor, mpz_t n, mpz_t B, mpz_t x, mpz_t y) {
    mpz_t a, value, gcd, counter, limit;
    mpz_init_set_ui(a, 1);
    mpz_init(value);
    mpz_init(gcd);
    mpz_init_set_ui(counter, 0);
    mpz_init(limit);
    mpz_sqrt(limit, B);
    mpz_mul_ui(limit, limit, (unsigned long) limit_num);
    mpz_div_ui(limit, limit, (unsigned long) limit_denom);

    while (1) {

        g(x, x, a, n);
        g(y, y, a, n);
        g(y, y, a, n);
        mpz_sub(value, x, y);
        mpz_abs(value, value);
        mpz_gcd(gcd, value, n);

        if (mpz_cmp_ui(gcd, 1) == 0) {
            mpz_set(factor, n);
            mpz_clear(a);
            mpz_clear(value);
            mpz_clear(gcd);
            mpz_clear(counter);
            mpz_clear(limit);
            return;
        }
        if (mpz_cmp(gcd, n) == 0) {
            mpz_set(factor, n);
            mpz_clear(a);
            mpz_clear(value);
            mpz_clear(gcd);
            mpz_clear(counter);
            mpz_clear(limit);
            return;
        }

        mpz_set(factor, gcd);
        mpz_clear(a);
        mpz_clear(value);
        mpz_clear(gcd);
        mpz_clear(counter);
        mpz_clear(limit);
        return;
    }
}

/*
struct element *addorderelem(mpz_t value, struct list *list) {
    struct element *elem = init_elem(value);
    struct element *temp;
    for (temp = list->HEAD; temp != NULL; temp = temp->next) {
        if (mpz_cmp(value, temp->value) > 0) {
            if (temp->next != NULL) {
                if (mpz_cmp(value, temp->next->value) < 0) {
                    temp->next = elem;
                    elem->prev = temp;
                    elem->next = temp->next;
                    temp->next->prev = elem;
                    list->count++;
                    return elem;
                }
            } else {
                list->TAIL = elem;
                list->TAIL->prev = temp;
                temp->next = elem;
                list->count++;
                return elem;
            }
        }
    }
    if (list->HEAD == NULL) {
        list->HEAD = elem;
        list->TAIL = elem;
        list->count++;
        return elem;
    }
}
 */

void eratosthenes(mpz_t n, struct list *primelist) {

    mpz_t i, j, mul;

    mpz_init(i);
    mpz_init(j);
    mpz_init(mul);

    mpz_t *array;
    array = malloc(sizeof(mpz_t) * mpz_get_ui(n));
    if (array == NULL) {
        fprintf(stderr, "Error in malloc\n");
        exit(EXIT_FAILURE);
    }

    for (mpz_set_ui(i, 0); mpz_cmp(i, n) < 0; mpz_add_ui(i, i, 1)) {
        mpz_init_set(array[mpz_get_ui(i)], i);
    }

    for (mpz_set_ui(i, 2); mpz_cmp(i, n) < 0; mpz_add_ui(i, i, 1)) {
        if (mpz_cmp_ui(array[mpz_get_ui(i)], 0) == 0) {
            continue;
        }
        addelem(i, primelist);
        for (mpz_set_ui(j, 2); mpz_cmp(j, n) < 0; mpz_add_ui(j, j, 1)) {
            mpz_mul(mul, i, j);
            if (mpz_cmp(mul, n) >= 0) {
                break;
            } else {
                mpz_set_ui(array[mpz_get_ui(mul)], 0);
            }
        }
    }

    for (mpz_set_ui(i, 0); mpz_cmp(i, n) < 0; mpz_add_ui(i, i, 1)) {
        mpz_clear(array[mpz_get_ui(i)]);
    }
    mpz_clear(i);
    mpz_clear(j);
    mpz_clear(mul);
    free(array);
    return;
}

void trialdivision(mpz_t factor, mpz_t n, struct list *primelist) {

    struct element *temp = primelist->HEAD;

    mpz_t r, root;
    mpz_init(r);
    mpz_init(root);
    mpz_sqrt(root, n);

    while (temp != NULL && mpz_cmp(temp->value, root) <= 0) {
        mpz_mod(r, n, temp->value);
        if (mpz_cmp_ui(r, 0) == 0) {
            mpz_clear(r);
            mpz_clear(root);
            mpz_set(factor, temp->value);
            return;
        }
        temp = temp->next;
    }
    mpz_clear(r);
    mpz_clear(root);
    mpz_set(factor, n);
    return;
}

struct list *isBsmooth(mpz_t n, mpz_t B, struct list *primelist) {

    struct list *newlist = init_list();

    mpz_t sfactor, uno;
    mpz_t todivide;
    mpz_t factor;
    mpz_t x, y;
    mpz_init_set_ui(x, 2);
    mpz_init_set_ui(y, 2);

    mpz_init_set(todivide, n);
    mpz_init_set_ui(uno, 1);

    mpz_init(factor);
    mpz_init(sfactor);


    while (mpz_cmp(todivide, uno) != 0) {

        pollard_rho_opt(factor, todivide, B, x, y);

        if (mpz_cmp_si(factor, -1) == 0) {
            mpz_clear(factor);
            mpz_clear(todivide);
            mpz_clear(uno);
            mpz_clear(x);
            mpz_clear(y);
            freelist(newlist);
            return NULL;
        }

        while (mpz_cmp_ui(factor, 1) != 0) {
            trialdivision(sfactor, factor, primelist);
            if (mpz_cmp(sfactor, B) > 0) {
                mpz_clear(factor);
                mpz_clear(todivide);
                mpz_clear(uno);
                mpz_clear(x);
                mpz_clear(y);
                freelist(newlist);
                return NULL;
            }
            mpz_div(factor, factor, sfactor);
            mpz_div(todivide, todivide, sfactor);
            addnewfactor(sfactor, newlist);
        }
    }
    mpz_clear(factor);
    mpz_clear(todivide);
    mpz_clear(uno);
    mpz_clear(x);
    mpz_clear(y);
    return newlist;
}


void gcdext(mpz_t a_origin, mpz_t b_origin, mpz_t x, mpz_t y, mpz_t z) {

    // Base Case
    mpz_set_ui(x, 0);
    mpz_t lastx, lasty, temp, quotient, term, root, a, b, invers;
    mpz_init(root);
    mpz_init(invers);
    mpz_sqrt(root, a_origin);
    mpz_init_set(a, a_origin);
    mpz_init_set(b, b_origin);
    size_t rightsize = mpz_sizeinbase(root, 10);
    mpz_init(term);
    mpz_init(quotient);
    mpz_init(temp);
    mpz_init_set_ui(lastx, 1);
    mpz_set_ui(y, 1);
    mpz_init_set_ui(lasty, 0);

    while (mpz_cmp_ui(b, 0) != 0) {

        mpz_set(temp, b);
        mpz_div(quotient, a, b);
        mpz_mod(b, a, b);
        mpz_set(a, temp);

        mpz_set(temp, x);
        mpz_mul(term, quotient, x);
        mpz_sub(x, lastx, term);

        mpz_set(lastx, temp);

        mpz_set(temp, y);
        mpz_mul(term, quotient, y);
        mpz_sub(y, lasty, term);

        mpz_set(lasty, temp);

        if (mpz_cmp_ui(a, 1) == 0) {
            mpz_set_ui(x, 0);
            mpz_set_ui(y, 0);
            mpz_set_ui(z, 0);
            mpz_clear(lastx);
            mpz_clear(lasty);
            mpz_clear(temp);
            mpz_clear(quotient);
            mpz_clear(term);
            mpz_clear(root);
            mpz_clear(a);
            mpz_clear(b);
            return;
        }

        if (mpz_sizeinbase(lasty, 10) == rightsize || mpz_sizeinbase(a, 10) <= rightsize) {
            if (mpz_cmp_ui(lasty, 0) > 0) {
                mpz_set(x, lastx);
                mpz_set(y, lasty);
                mpz_set(z, a);
                mpz_clear(lastx);
                mpz_clear(lasty);
                mpz_clear(temp);
                mpz_clear(quotient);
                mpz_clear(term);
                mpz_clear(root);
                mpz_clear(a);
                mpz_clear(b);
                return;
            }
        }
    }
}
