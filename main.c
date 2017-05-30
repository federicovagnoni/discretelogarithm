#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
//#include <pari/pari.h>

void g(mpz_t result, mpz_t x, mpz_t n) {
    mpz_t pow;
    mpz_init(pow);
    mpz_pow_ui(pow, x, 2);
    mpz_mul_ui(pow, pow, 1);
    mpz_add_ui(pow, pow, 1);
    mpz_mod(pow, pow, n);
    mpz_set(result, pow);
    return;
}

void pollard_rho(mpz_t factor, mpz_t n) {
    mpz_t x, y, d, sub;
    mpz_init(sub);
    mpz_init_set_ui(x, 1);
    mpz_init_set_ui(y, 1);
    mpz_init_set_ui(d, 1);

    while (mpz_cmp_ui(d, 1) == 0) {
        g(x, x, n);
        //printf("x = %.0lf\n", x);
        g(y, y, n);
        g(y, y, n);
        //printf("y = %.0lf\n", y);
        mpz_sub(sub, x, y);
        mpz_abs(sub, sub);
        mpz_gcd(d, sub, n);
        //printf("d = %.0lf\n", d);
    }
    if (mpz_cmp(d, n) == 0) {
        mpz_set_ui(factor, (unsigned long) -1);
        return;
    } else {
        mpz_set(factor, d);
        return;
    }
}

//double isprime(double n) {
//    double x;
//    double x0 = 2;
//    int i = 0;
//    while (n != 1) {
//        x = pollard_rho(x0, n);
//        if (x == -1) {
////            i++;
////            x0 = x0 * -2;
////            //printf("Cambio x0: %.0lf (n = %.0lf)\n", x0, n);
////            if (i == 10) {
////                //printf("%.0lf  end\n", n);
//            break;
////            }
//        } else {
//            return 0;
//            x0 = 2;
//            //printf("%.0lf   %.0lf\n", n, x);
//            n = n / x;
//        }
//    }
//    return 1;
//}


struct element {
    mpz_t value;
    int exp;
    struct element *next;
    struct element *prev;
};

struct list {
    struct element *HEAD;
    struct element *TAIL;
    int count;
};

struct element *init_elem(mpz_t value) {
    struct element *elem;
    elem = malloc(sizeof(struct element));

    if (elem == NULL) {
        fprintf(stderr, "Error in malloc\n");
        exit(EXIT_FAILURE);
    }

    mpz_init_set(elem->value, value);
    elem->exp = 1;
    elem->next = NULL;
    elem->prev = NULL;

    return elem;
}

struct list *init_list() {
    struct list *newlist;
    newlist = malloc(sizeof(struct list));

    if (newlist == NULL) {
        fprintf(stderr, "Error in malloc\n");
        exit(EXIT_FAILURE);
    }

    newlist->HEAD = NULL;
    newlist->TAIL = NULL;
    newlist->count = 0;

    return newlist;
}


void addelem(mpz_t value, struct list *list) {
    struct element *elem = init_elem(value);

    if (list->HEAD == NULL) {
        list->HEAD = elem;
        list->TAIL = elem;
    } else {
        list->TAIL->next = elem;
        elem->prev = list->TAIL;
        list->TAIL = elem;
    }
    list->count++;
};

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

void freelist(struct list *list) {
    struct element *temp = NULL;
    if (list->count == 1) {
        mpz_clear(list->HEAD->value);
        free(list->HEAD);
        free(list);
        return;
    }
    while (temp != list->HEAD) {
        temp = list->TAIL->prev;
        mpz_clear(list->TAIL->value);
        free(list->TAIL);
        list->TAIL = temp;
    }
//    mpz_clear(list->HEAD->value);
    free(list->HEAD);
    free(list);
}

void printlist(struct list *list) {
    struct element *temp = list->HEAD;
    while (temp != NULL) {
        mpz_out_str(stdout, 10, temp->value);
        printf("^%d", temp->exp);
        temp = temp->next;
        if (temp != NULL) {
            printf(" * ");
        }
    }
    printf("\n");
    free(temp);
}

void returnlist(struct list *linkedlist, mpz_t *returnedlist) {
    struct element *temp = linkedlist->HEAD;
    int i = 0;
    while (temp != NULL) {
        mpz_init_set(returnedlist[i], temp->value);
        temp = temp->next;
        i++;
    }
}


struct list *primes_in_range(mpz_t n, mpz_t m) {
    struct list *newlist = init_list();
    mpz_t i;
    mpz_t r;
    mpz_t en;
    mpz_init_set(en, n);
    mpz_init(r);
    mpz_init(i);
    while (mpz_cmp(en, m) <= 0) {
        mpz_set_si(i, 2);
        while (mpz_cmp(i, en) < 0) {
            mpz_mod(r, en, i);
            if (mpz_cmp_si(r, 0) == 0) {
                break;
            }
            mpz_add_ui(i, i, 1);
        }
        if (mpz_cmp(i, en) == 0) {
            addelem(en, newlist);
        }
        mpz_add_ui(en, en, 1);
    }
    mpz_clear(en);
    mpz_clear(r);
    mpz_clear(i);
    return newlist;
}

void trialdivison(mpz_t factor, mpz_t n, struct list *list) {
    struct element *temp = list->HEAD;
    mpz_t r;
    mpz_init(r);
    while (temp != NULL) {
//        printf("n = ");
//        mpz_out_str(stdout, 10, n);
//        printf("\ntempvalue = ");
//        mpz_out_str(stdout, 10, temp->value);
        mpz_mod(r, n, temp->value);
//        printf("\nr = ");
//        mpz_out_str(stdout, 10, r);
//        printf("\n");
        if (mpz_cmp_ui(r, 0) == 0) {
            mpz_clear(r);
            mpz_set(factor, temp->value);
            return;
        }
        temp = temp->next;
    }
    mpz_clear(r);
    free(temp);
    mpz_set(factor, n);
    return;
}


struct element *addnewfactor(mpz_t value, struct list *list) {
    struct element *temp;
    for (temp = list->HEAD; temp != NULL; temp = temp->next) {
        if (mpz_cmp(temp->value, value) == 0) {
            temp->exp++;
            return temp;
        }
    }
    addelem(value, list);
    return list->TAIL;
}

struct list *factorsbytrialdivision(mpz_t n) {
    struct list *newlist = init_list();
    mpz_t factor;
    mpz_init(factor);
    mpz_t todivide;
    mpz_init_set(todivide, n);
    struct list *primelist;
    mpz_t root;
    mpz_init(root);
    mpz_t uno;
    mpz_init_set_ui(uno, 1);
    while (mpz_cmp_ui(todivide, 1) != 0) {
        mpz_sqrt(root, todivide);
        primelist = primes_in_range(uno, root);
        trialdivison(factor, todivide, primelist);
        freelist(primelist);
        if (mpz_cmp_ui(factor, -1) == 0) {
            break;
        }
        addnewfactor(factor, newlist);
        mpz_div(todivide, todivide, factor);
    }
    mpz_clear(factor);
    mpz_clear(todivide);
    mpz_clear(root);
    return newlist;
}

struct list *factorsbypollard(mpz_t n) {
    struct list *newlist = init_list();
    mpz_t factor;
    mpz_init(factor);
    mpz_t todivide;
    mpz_init_set(todivide, n);
    mpz_t root;
    mpz_init(root);
    mpz_t uno,due;
    mpz_init_set_ui(uno, 1);
    mpz_init_set_ui(due, 2);
    mpz_t mod;
    mpz_init(mod);

    mpz_mod(mod, todivide, due);
    while (mpz_cmp_ui(mod, 0) == 0) {
        mpz_div(todivide, todivide, due);
        addnewfactor(due, newlist);
        mpz_mod(mod, todivide, due);
    }

    while (mpz_cmp_ui(todivide, 1) != 0) {

        pollard_rho(factor, todivide);

        if (mpz_cmp_ui(factor, -1) == 0) {
            addnewfactor(todivide, newlist);
            break;
        }

        addnewfactor(factor, newlist);
        mpz_div(todivide, todivide, factor);

    }
    mpz_clear(factor);
    mpz_clear(todivide);
    mpz_clear(root);
    mpz_clear(mod);
    mpz_clear(uno);
    mpz_clear(due);
    return newlist;
}


void removeelem(mpz_t value, struct list *list) {
    struct element *temp;
    for (temp = list->HEAD; temp != NULL; temp = temp->next) {
        if (mpz_cmp(temp->value, value) == 0) {
            if (temp->prev != NULL) {
                temp->prev->next = temp->next;
            } else {
                list->HEAD = temp->next;
            }
            if (temp->next != NULL) {
                temp->next->prev = temp->prev;
            } else {
                list->TAIL = temp->prev;
            }
            list->count--;
            free(temp);
            return;
        }
    }
}

int isBsmooth(mpz_t m, mpz_t B) {
    mpz_t todivide;
    mpz_init_set(todivide, m);
    mpz_t factor;
    mpz_init(factor);
    mpz_t r;
    mpz_init(r);
    mpz_sqrt(r, m);
    mpz_t uno;
    mpz_init_set_ui(uno, 1);
    struct list *list = primes_in_range(uno, r);
    /*
    printf("primi tra ");
    mpz_out_str(stdout,10, uno);
    printf(" e ");
    mpz_out_str(stdout,10, r);
    printf("\n");
    printlist(list);
     */
    while (mpz_cmp_si(todivide, 1) != 0) {
        trialdivison(factor, todivide, list);
        /*
        printf("%d\n", mpz_cmp(factor,B));
        printf("B = ");
        mpz_out_str(stdout,10,B);
        printf("\nfactor = ");
        mpz_out_str(stdout,10, factor);
        printf("\n");
         */
        if (mpz_cmp(factor, B) > 0) {
            //printf("Il numero %.0lf ha %.0lf come fattore: non è %.0lf-smooth\n", m, factor, B);
            mpz_clear(todivide);
            mpz_clear(factor);
            mpz_clear(r);
            mpz_clear(uno);
            freelist(list);
            return 0;
        }
        //printf("%.0lf ha fattore %.0lf\n", m, factor);
        //printf("%.0lf è un fattore\n", factor);
        mpz_div(todivide, todivide, factor);
        //printf("il numero diventerà %.0lf\n", todivide);
    }
    mpz_clear(factor);
    mpz_clear(r);
    mpz_clear(todivide);
    mpz_clear(uno);
    freelist(list);
    return 1;

//    struct list *list = factorsbypollard(m);
//    printlist(list);
//    printf("\n");
//    fflush(stdout);
//    if (mpz_cmp(list->TAIL->value, B) < 0) {
//        return 1;
//    }
//    return 0;
}

double main(int argc, char *argv[]) {

    unsigned long long p, q, a, B, lp, exponent, t;

    printf("Inserisci il numero da fattorizzare: ");
    scanf("%llu", &q);
    printf("q = %llu\n", q);

    p = 1 + 2 * q;
    printf("p = %llu\n", p);

    exponent = (unsigned long long int) sqrtl(logl(p) * logl(logl(p)));
    printf("exponent = %llu\n", exponent);
    lp = (unsigned long long int) expl(exponent);
    printf("L(p) = %llu\n", lp);
    printf("L(p)^c = %.0Lf\n", sqrtl(lp));
    B = (unsigned long long int) ceill(sqrtl(lp));
    B = B + B * 3/2;

    printf("B = %llu\n", B);
    mpz_t am, pm, pm1, exp, try, qm, due;
    mpz_init(am);
    mpz_init(exp);
    mpz_init(try);
    mpz_init_set_ui(pm, p);
    mpz_init_set_ui(qm, q);
    mpz_init_set_ui(due, 2);
    mpz_init(pm1);
    mpz_sub_ui(pm1, pm, 1);

    int count = 0;
    struct element *divisor;
    struct list *divisors = factorsbypollard(pm1); //init_list();
//    addnewfactor(qm, divisors);
//    addnewfactor(due, divisors);

    printlist(divisors);

    for (mpz_set_ui(am, 2); mpz_cmp(am, pm) < 0; mpz_add_ui(am, am, 1)) {
        for (divisor = divisors->HEAD; divisor != NULL; divisor = divisor->next) {
            mpz_div(exp, pm1, divisor);
            mpz_powm(try, am, exp, pm);
            if (mpz_cmp_ui(try, 1) == 0) {
                break;
            } else {
                count++;
            }
        }
        if (count == divisors->count) {
            mpz_out_str(stdout, 10, am);
            printf(" è una radice primitiva\n");
            break;
        } else {
            count = 0;
        }

    }


    printf("Inserisci una radice primitiva (e base del logaritmo): ");
    scanf("%llu", &a);

    printf("Inserisci l'argomento del logaritmo: ");
    scanf("%llu", &t);
    mpz_t Bm, tm;
    mpz_init_set_ui(Bm, B);
    mpz_set_ui(am, a);
    mpz_init_set_ui(tm, t);

    mpz_t uno;
    mpz_init_set_ui(uno, 1);
    struct list *primelist = primes_in_range(uno, Bm);

    printf("La lista dei primi compresi tra 1 e ");
    mpz_out_str(stdout, 10, Bm);
    printf(" è: ");
    printlist(primelist);
    printf("\n");

    mpz_t listprime[primelist->count];
    returnlist(primelist, listprime);

    mpz_t pm2, logpm, p2m;
    mpz_init(pm2);
    mpz_init(logpm);
    mpz_sub_ui(pm2, pm, 2);
    mpz_init(p2m);
    mpz_mul_ui(p2m, pm, 2);
    mpz_set_ui(logpm, (unsigned long) ceill(log(p)));

    int i = 0;
    mpz_t randexp;
    srand(time(NULL));

    mpz_t candidate;
    mpz_t res;
    mpz_init(res);

    struct list *candidatefactors;
    gmp_randstate_t randstate;
    gmp_randinit_default(randstate);
    gmp_randseed_ui(randstate, (unsigned long) time(NULL));

//    struct list *newlist;
    struct element *temp1, *temp2;

    mpz_t *matrix = NULL;
    int l, k, z;
    mpz_t u;
    mpz_init(u);
    mpz_t base;
    mpz_init(base);
    mpz_init(randexp);
    mpz_init(candidate);
    int x, y, s;
    mpz_t randexp1, randexp2;
    mpz_t candidate1, candidate2;
    mpz_init(randexp1);
    mpz_init(randexp2);
    mpz_init(candidate1);
    mpz_init(candidate2);

    mpz_t temp;
    mpz_init(temp);

    mpz_t invers;
    mpz_init(invers);

    int SIZE = primelist->count;
    int RELATIONS = B + 4;
    int ROWSIZE = primelist->count + 1;
    mpz_t c;
    mpz_t sol[SIZE];
    mpz_init(c);

    int isredundant = 0, equals = 0;
    while (i < primelist->count) {


        for (mpz_set_ui(u, 1); mpz_cmp(u, pm1) <= 0; mpz_add_ui(u, u, 1)) {

            mpz_powm(candidate, am, u, pm);
//            mpz_pow_ui(res, am, mpz_get_ui(u));

//            mpz_out_str(stdout, 10, res);
//            printf("\n");
//            mpz_out_str(stdout, 10, u);
//            printf("\n");
//            fflush(stdout);
//
            mpz_out_str(stdout, 10, candidate);
            printf(" (pari a ");
            mpz_out_str(stdout, 10, am);
            printf(" alla ");
            mpz_out_str(stdout, 10, u);
            printf(") in prova\n");
            fflush(stdout);

            if (mpz_cmp(candidate, uno) != 0) {// && mpz_cmp(res, pm) > 0) {
                printf("Vedo se è Bsmooth\n");
                fflush(stdout);
                if (isBsmooth(candidate, Bm) == 1) { //&& mpz_cmp(res, pm) > 0) {// && mpz_cmp(res, p2m) < 0) {
                    mpz_out_str(stdout, 10, candidate);
                    printf(" (pari a ");
                    mpz_out_str(stdout, 10, am);
                    printf(" alla ");
                    mpz_out_str(stdout, 10, u);
                    printf(") è %llu-smooth\n", B);
                    fflush(stdout);

                    candidatefactors = factorsbytrialdivision(candidate);
                    printlist(candidatefactors);

                    matrix = realloc(matrix, sizeof(mpz_t) * ROWSIZE * (i + 1));

                    int j;
                    for (j = 0; j < ROWSIZE; j++) {
                        mpz_init_set_ui(*(matrix + i * ROWSIZE + j), 0);
                    }


                    j = 0;
                    for (temp1 = primelist->HEAD; temp1 != NULL; temp1 = temp1->next) {
                        for (temp2 = candidatefactors->HEAD; temp2 != NULL; temp2 = temp2->next) {
                            if (mpz_cmp(temp1->value, temp2->value) == 0) {
                                mpz_set_ui(*(matrix + i * ROWSIZE + j), (unsigned long) temp2->exp);
                            }
                        }
                        j++;
                    }

                    mpz_set(*(matrix + i * ROWSIZE + ROWSIZE - 1), u);

                    freelist(candidatefactors);

                    // Vedo se ci sono righe uguali
                    equals = 0;
                    isredundant = 0;

                    if (i > 0) {
                        for (l = i - 1; l >= 0; l--) {
                            for (j = 0; j < primelist->count; j++) {
                                if (mpz_cmp(*(matrix + i * ROWSIZE + j), *(matrix + l * ROWSIZE + j)) == 0) {
                                    equals++;
                                    if (equals == primelist->count) {
                                        isredundant = 1;
                                    }
                                }
                            }
                            equals = 0;
                        }
                    }

                    // Se ci sono esco e provo con un'altra combinazione
                    if (isredundant == 1) {
                        printf("Ci sono righe uguali!\n");
                        continue;
                    }


                    int zeros = 0;
                    isredundant = 0;

                    for (x = 0; x < ROWSIZE - 1; x++) {
                        if (mpz_cmp_ui(*(matrix + i * ROWSIZE + x), 0) != 0) {
                            if (mpz_invert(invers, *(matrix + i * ROWSIZE + x), pm1) == 0) {
                                isredundant = 1;
                            } else {
                                isredundant = 0;
                            }
                        } else {
                            zeros++;
                        }
                    }

                    if (zeros == ROWSIZE - 2 && isredundant == 0) {
                        for (x = 0; x < ROWSIZE; x++) {
                            mpz_mul(*(matrix + i * ROWSIZE + x), *(matrix + i * ROWSIZE + x), invers);
                            mpz_mod(*(matrix + i * ROWSIZE + x), *(matrix + i * ROWSIZE + x), pm1);
                            isredundant = 0;
                        }
                    } else {
                        isredundant = 1;
                    }

                    if (isredundant == 1) {
                        printf("E' una riga inutile!\n");
                        continue;
                    }

                    equals = 0;
                    isredundant = 0;

                    if (i > 0) {
                        for (l = i - 1; l >= 0; l--) {
                            for (j = 0; j < primelist->count; j++) {
                                if (mpz_cmp(*(matrix + i * ROWSIZE + j), *(matrix + l * ROWSIZE + j)) == 0) {
                                    equals++;
                                    if (equals == primelist->count) {
                                        isredundant = 1;
                                    }
                                }
                            }
                            equals = 0;
                        }
                    }

                    if (isredundant == 1) {
                        printf("Ci sono ancora righe uguali!\n");
                        continue;
                    }
                    for (k = 0; k < i + 1; k++) {
                        printf("[ ");
                        for (j = 0; j < ROWSIZE; j++) {
                            mpz_out_str(stdout, 10, *(matrix + k * ROWSIZE + j));
                            printf(" ");
                        }
                        printf("]\n");
                    }
                    printf("\n");

                    i++;
                    if (i == RELATIONS) {
                        break;
                    }
                }
            }
            if (i == RELATIONS) {
                break;
            }
        }
        if (i == RELATIONS) {
            break;
        }
    }

    mpz_clear(candidate1);
    mpz_clear(candidate2);
    mpz_clear(randexp1);
    mpz_clear(randexp2);
//    mpz_clear(randexp);
//    mpz_clear(candidate);
//    mpz_clear(res);
//    mpz_clear(base);

    int j;

    for (k = 0; k < primelist->count; k++) {
        printf("[ ");
        for (j = 0; j < ROWSIZE; j++) {
            mpz_out_str(stdout, 10, *(matrix + k * ROWSIZE + j));
            printf(" ");
        }
        printf("]\n");
    }
    printf("\n");


    for (mpz_set_ui(u, 1); mpz_cmp(u, pm1) <= 0; mpz_add_ui(u, u, 1)) {

        mpz_powm(candidate, am, u, pm);
        mpz_mul_ui(candidate, candidate, t);
        mpz_mod(candidate, candidate, pm);
        mpz_set(res, am);

//        for (s = 0; s < (int) mpz_get_ui(u) - 1; s++) {
//            mpz_mul(res, res, am);
//        }

        if (isBsmooth(candidate, Bm) == 1 &&
            mpz_cmp(candidate, uno) != 0) { //&& mpz_cmp(res, pm) > 0) {// && mpz_cmp(res, p2m) < 0) {
//            mpz_out_str(stdout, 10, candidate);
//            printf(" (pari a %llu *", t);
//            mpz_out_str(stdout, 10, am);
//            printf(" alla ");
//            mpz_out_str(stdout, 10, u);
//            printf(") è %llu-smooth\n", B);
//            fflush(stdout);

            candidatefactors = factorsbytrialdivision(candidate);
            printlist(candidatefactors);

            struct element *factor, *prime;
            mpz_t result;
            mpz_init_set_ui(result, 0);
            int pos = 0;
            int row;
            int times;

            for (factor = candidatefactors->HEAD; factor != NULL; factor = factor->next) {
                for (prime = primelist->HEAD; prime != NULL; prime = prime->next) {
                    if (mpz_cmp(prime->value, factor->value) == 0) {
                        break;
                    }
                    pos++;
                }
                for (row = 0; row < primelist->count; row++) {
//                    printf("Vedo posizione [%d,%d]\n", row, pos);
                    if (mpz_cmp_ui(*(matrix + row * ROWSIZE + pos), 1) == 0) {
//                        printf("Trovato un 1 in [%d,%d] per il fattore ", row, pos);
//                        mpz_out_str(stdout, 10, factor->value);
//                        printf("\n");
//                        printf("Il valore corrispondente è: ");
//                        mpz_out_str(stdout, 10, *(matrix + row * ROWSIZE + ROWSIZE - 1));
//                        printf("\n");
                        for (times = 0; times < factor->exp; times++) {
                            mpz_add(result, result, *(matrix + row * ROWSIZE + ROWSIZE - 1));
//                            mpz_out_str(stdout, 10, result);
//                            printf("\n");
                        }
                    }
                }
                pos = 0;
            }
            mpz_sub(result, result, u);
            mpz_mod(result, result, pm1);
            printf("Il risultato è: ");
            mpz_out_str(stdout, 10, result);
            printf("\n");
            break;
        }
    }


//    mpz_invert(invers, *(matrix + (SIZE - 1) * ROWSIZE + (SIZE- 1)), pm1);
//    mpz_mul(sol[SIZE], *(matrix + (SIZE - 1) * ROWSIZE + (SIZE)), invers);
//
//    mpz_t sum, num, denom;
//    mpz_init(num);
//    mpz_init(denom);
//    mpz_init_set_ui(sum ,0);
//    /* this loop is for backward substitution*/
//    for(i = SIZE; i > 0; i--)
//    {
//        mpz_set(sum, *(matrix + i * ROWSIZE + SIZE));
//        for(j = i+1; j < SIZE; j++)
//        {
//            mpz_submul(sum, sol[j],*(matrix + i * ROWSIZE + j));
//            // sum = sum - x[j] * *(matrix + i * ROWSIZE + j))
//        }
//        mpz_invert(denom, *(matrix + i * ROWSIZE + i), pm1);
//        mpz_sub(num, *(matrix + i * ROWSIZE + SIZE), sum);
//        mpz_mul(sol[i], num, denom);
//        printf("\n x%d =>", i+1);
//        mpz_out_str(stdout,10,sol[i]);
//    }


    mpz_clear(pm2);


    return EXIT_SUCCESS;
}




