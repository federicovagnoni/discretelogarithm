#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <time.h>

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



void g(mpz_t result, mpz_t x, mpz_t n) {
    mpz_t pow;
    mpz_init(pow);
    mpz_pow_ui(pow, x, 2);
    mpz_mul_ui(pow, pow, 1);
    mpz_add_ui(pow, pow, 1);
    mpz_mod(pow, pow, n);
    mpz_set(result, pow);
    mpz_clear(pow);
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
        mpz_clear(x);
        mpz_clear(y);
        mpz_clear(d);
        mpz_clear(sub);
        mpz_set_si(factor, -1);
        return;
    } else {
        mpz_set(factor, d);
        mpz_clear(x);
        mpz_clear(y);
        mpz_clear(d);
        mpz_clear(sub);
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
    return;
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
    mpz_t i, r, en;

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

        mpz_mod(r, n, temp->value);

        if (mpz_cmp_ui(r, 0) == 0) {
            mpz_clear(r);
            mpz_set(factor, temp->value);
            return;
        }
        temp = temp->next;
    }
    mpz_clear(r);
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

void factorsbytrialdivision(mpz_t n, struct list *list, struct list *primelist) {

    mpz_t factor, todivide, uno, root;

    mpz_init(factor);
    mpz_init(root);
    mpz_init_set(todivide, n);
    mpz_init_set_ui(uno, 1);
    mpz_sqrt(root, todivide);

    while (mpz_cmp_ui(todivide, 1) != 0) {
        trialdivison(factor, todivide, primelist);
        if (mpz_cmp_si(factor, -1) == 0) {
            break;
        }
        addnewfactor(factor, list);
        mpz_div(todivide, todivide, factor);
    }
    mpz_clear(factor);
    mpz_clear(todivide);
    mpz_clear(root);
    return;
}

struct list *factorsbypollard(mpz_t n) {

    struct list *newlist = init_list();
    struct list *primelist;
    mpz_t factor, todivide, uno, due, mod, root;

    mpz_init(factor);
    mpz_init(mod);
    mpz_init(root);
    mpz_init_set(todivide, n);
    mpz_init_set_ui(uno, 1);
    mpz_init_set_ui(due, 2);

    mpz_mod(mod, todivide, due);

    while (mpz_cmp_ui(mod, 0) == 0) {
        mpz_div(todivide, todivide, due);
        addnewfactor(due, newlist);
        mpz_mod(mod, todivide, due);
    }

    while (mpz_cmp_ui(todivide, 1) != 0) {

        pollard_rho(factor, todivide);

        if (mpz_cmp_si(factor, -1) == 0) {
            addnewfactor(todivide, newlist);
            break;
        }
        mpz_sqrt(root, factor);
        primelist = primes_in_range(uno, root);
        factorsbytrialdivision(factor, newlist, primelist);
        freelist(primelist);
        mpz_div(todivide, todivide, factor);

    }
    mpz_clear(factor);
    mpz_clear(todivide);
    mpz_clear(mod);
    mpz_clear(uno);
    mpz_clear(due);
    return newlist;
}


struct list* isBsmooth(mpz_t m, mpz_t B) {

    struct list *factors = factorsbypollard(m);
//    printlist(factors);
//    fflush(stdout);
    struct element *elem;
    for (elem = factors->HEAD; elem != NULL; elem = elem->next) {
        if (mpz_cmp(elem->value, B) > 0) {
            freelist(factors);
            return NULL;
        }
    }
    return factors;
}

double main(int argc, char *argv[]) {

    unsigned long long p, q, a, B, lp, exponent, t;
//    mpz_t b, e, primo;
//    mpz_init_set_ui(b, 6);
//    mpz_init_set_ui(e, 58756221822);
//    mpz_init_set_ui(primo, 98764321261);
//    mpz_powm(b, b, e, primo);
//    mpz_out_str(stdout, 10, b);

    printf("Inserisci il numero da fattorizzare: ");
    scanf("%llu", &q);
    printf("q = %llu\n", q);

    mpz_t am, pm, pm1, exp, try, qm, uno, due, Bm, tm, pm2, logpm, p2m;

    mpz_init(am);
    mpz_init(exp);
    mpz_init(try);
    mpz_init(pm);
    mpz_init(pm1);
    mpz_init(pm2);
    mpz_init(logpm);
    mpz_init(p2m);

    mpz_init_set_ui(qm, q);
    mpz_init_set_ui(uno, 1);
    mpz_init_set_ui(due, 2);

    p = 2*q + 1;

    mpz_mul(pm, qm, due);
    mpz_add(pm, pm, uno);
    printf("p = ");
    mpz_out_str(stdout, 10, pm);
    printf("\n");

    exponent = (unsigned long long int) sqrtl(logl(p) * logl(logl(p)));
    printf("exponent = %llu\n", exponent);
    lp = (unsigned long long int) expl(exponent);
    printf("L(p) = %llu\n", lp);
    printf("L(p)^c = %.0Lf\n", sqrtl(lp));
    B = (unsigned long long int) ceill(sqrtl(lp));

    B = B + B * 3/2;
    printf("B = %llu\n", B);

    mpz_init_set_ui(Bm, B);

    mpz_sub_ui(pm1, pm, 1);

//    printf("Calcolo i primi...");
//    fflush(stdout);
//    struct list *primes = primes_in_range(uno, root);
//    printf("done\n");
//    fflush(stdout);

    int count = 0;
    struct element *divisor;
    printf("Calcolo i divisori...");
    fflush(stdout);
    struct list *divisors = factorsbypollard(pm1); //init_list();
    printf("done\n\n");
    fflush(stdout);

    printf("p - 1 = ");
    printlist(divisors);
    printf("\n");
    fflush(stdout);

    for (mpz_set_ui(am, 2); mpz_cmp(am, pm) < 0; mpz_add_ui(am, am, 1)) {
        for (divisor = divisors->HEAD; divisor != NULL; divisor = divisor->next) {
            mpz_div(exp, pm1, divisor->value);
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

    freelist(divisors);
    free(divisor);


    printf("Inserisci una radice primitiva (e base del logaritmo): ");
    scanf("%llu", &a);

    mpz_set_ui(am, a);

    printf("Inserisci l'argomento del logaritmo: ");
    scanf("%llu", &t);

    mpz_init_set_ui(tm, t);

    struct list *primelist = primes_in_range(uno, Bm);

    printf("Tra 1 e ");
    mpz_out_str(stdout, 10, Bm);
    printf(" ci sono %d numeri primi\n", primelist->count);

    mpz_t listprime[primelist->count];
    returnlist(primelist, listprime);


    mpz_sub_ui(pm2, pm, 2);
    mpz_mul_ui(p2m, pm, 2);
    mpz_set_ui(logpm, (unsigned long) ceill(logl(p)));

    int i = 0;
    int j, l, k, z;
    int RELATIONS = primelist->count;
    int ROWSIZE = primelist->count + 1;

    int isredundant = 0, equals = 0;

    mpz_t candidate, res, u, invers, c, term;
    mpz_init(res);
    mpz_init(u);
    mpz_init(invers);
    mpz_init(candidate);
    mpz_init(c);
    mpz_init(term);

    struct list *candidatefactors;

    struct element *temp1, *temp2;

    mpz_t *matrix = NULL;

    int useless = 0;

    gmp_randstate_t rand;
    gmp_randinit_default(rand);

    gmp_randseed_ui(rand, (unsigned long) time(NULL));

    //
    printf("Calcolo delle %d relazioni...\n", RELATIONS);
    while (i < primelist->count) {

        for (mpz_set_ui(u, 1); mpz_cmp(u, pm) < 0; mpz_add_ui(u, u, 1)) {
//            mpz_urandomm(u, rand, pm);
//            mpz_add_ui(u, u, 1);

            mpz_powm(candidate, am, u, pm);
            if (useless == 0) {
                mpz_pow_ui(res, am, mpz_get_ui(u));
                if (mpz_cmp(res, pm) > 0)
                    useless = 1;
                    mpz_set_ui(res, 1);
            }

//            mpz_out_str(stdout, 10, res);
//            printf(" ovvero ");
//            mpz_out_str(stdout, 10, candidate);
//            printf(" (pari a ");
//            mpz_out_str(stdout, 10, am);
//            printf(" alla ");
//            mpz_out_str(stdout, 10, u);
//            printf(") in prova\n");
//            fflush(stdout);
//            mpz_out_str(stdout, 10, u);
//            printf("\n");

            if (mpz_cmp(res, pm) > 0 || useless == 1) {

                candidatefactors = isBsmooth(candidate, Bm);

                if (candidatefactors != NULL) {

//                    mpz_out_str(stdout, 10, candidate);
//                    printf(" (pari a ");
//                    mpz_out_str(stdout, 10, am);
//                    printf(" alla ");
//                    mpz_out_str(stdout, 10, u);
//                    printf(") è %llu-smooth\n", B);
//                    fflush(stdout);


                    matrix = realloc(matrix, sizeof(mpz_t) * ROWSIZE * (i + 1));

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

//                    for (k = 0; k < i + 1; k++) {
//                        printf("[ ");
//                        for (j = 0; j < ROWSIZE; j++) {
//                            mpz_out_str(stdout, 10, *(matrix + k * ROWSIZE + j));
//                            printf(" ");
//                        }
//                        printf("]\n");
//                    }
//                    printf("\n");

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
                        continue;
                    }


                    //altrimenti effettuo l'eliminazione di Gauss
                    for (k = 0; k < i + 1; k++) {
                        if (mpz_invert(invers, *(matrix + k * ROWSIZE + k), pm1) == 0) {
                            isredundant = 1;
                            break;
                        }
                        for (z = k + 1; z < i + 1; z++) {
                            mpz_mul(c, *(matrix + z * ROWSIZE + k), invers);
                            for (l = 0; l < ROWSIZE; l++) {
                                mpz_mul(term, c, *(matrix + k * ROWSIZE + l));
                                mpz_sub(*(matrix + z * ROWSIZE + l), *(matrix + z * ROWSIZE + l), term);
                                mpz_mod(*(matrix + z * ROWSIZE + l), *(matrix + z * ROWSIZE + l), pm1);
                            }
                        }
                    }

                    if (isredundant == 1) {
                        continue;
                    }

                    equals = 0;
                    isredundant = 0;

                    // Se la nuova riga è comunque uguale alle altre, la rimuovo
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
                        continue;
                    }

//                    for (k = 0; k < i + 1; k++) {
//                        printf("[ ");
//                        for (j = 0; j < ROWSIZE; j++) {
//                            mpz_out_str(stdout, 10, *(matrix + k * ROWSIZE + j));
//                            printf(" ");
//                        }
//                        printf("]\n");
//                    }
//                    printf("\n");

                    printf("%d di %d relazioni trovate\n", i + 1, RELATIONS);
                    fflush(stdout);

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
        mpz_set_ui(u, 0);
        printf("Ricomincio il calcolo...\n");
    }
    printf("done!\n\n");



//    for (k = 0; k < primelist->count; k++) {
//        printf("[ ");
//        for (j = 0; j < ROWSIZE; j++) {
//            mpz_out_str(stdout, 10, *(matrix + k * ROWSIZE + j));
//            printf(" ");
//        }
//        printf("]\n");
//    }
//    printf("\n");

    mpz_t sol[RELATIONS];
    for (j = 0; j < RELATIONS; j++) {
        mpz_init(sol[j]);
    }

    mpz_invert(invers, *(matrix + (RELATIONS - 1) * ROWSIZE + (RELATIONS - 1)), pm1);
    mpz_mul(sol[RELATIONS - 1], *(matrix + (RELATIONS - 1) * ROWSIZE + RELATIONS), invers);
    mpz_mod(sol[RELATIONS - 1], sol[RELATIONS - 1], pm1);

    mpz_t sum, num, denom;
    mpz_init(num);
    mpz_init(denom);
    mpz_init(sum);

    /* this loop is for backward substitution*/
    printf("Calcolo della matrice dei coefficienti...");
    for(i = RELATIONS - 2; i >= 0; i--)
    {
        mpz_set_ui(sum, 0);
        for(j = i + 1; j < RELATIONS; j++)
        {
            mpz_addmul(sum, sol[j], *(matrix + i * ROWSIZE + j));
            // sum = sum + x[j] * *(matrix + i * ROWSIZE + j))
        }
        mpz_invert(denom, *(matrix + i * ROWSIZE + i), pm1);
        mpz_sub(num, *(matrix + i * ROWSIZE + RELATIONS), sum);
        mpz_mul(sol[i], num, denom);
        mpz_mod(sol[i], sol[i], pm1);
    }

    for (k = 0; k < RELATIONS; k++) {
        for (j = 0; j < RELATIONS; j++) {
            if (k == j) {
                mpz_set_ui(*(matrix + k * ROWSIZE + j), 1);
                mpz_set(*(matrix + k * ROWSIZE + ROWSIZE - 1), sol[k]);
            } else {
                mpz_set_ui(*(matrix + k * ROWSIZE + j), 0);
            }
        }
    }

    printf("done\n\n");
    fflush(stdout);

//    for (k = 0; k < primelist->count; k++) {
//        printf("[ ");
//        for (j = 0; j < ROWSIZE; j++) {
//            mpz_out_str(stdout, 10, *(matrix + k * ROWSIZE + j));
//            printf(" ");
//        }
//        printf("]\n");
//    }
//    printf("\n");

    mpz_t result;
    mpz_init(result);

    printf("Calcolo della soluzione...");
    fflush(stdout);

    for (mpz_set_ui(u, 1); mpz_cmp(u, pm1) <= 0; mpz_add_ui(u, u, 1)) {

        mpz_powm(candidate, am, u, pm);
        mpz_mul_ui(candidate, candidate, t);
        mpz_mod(candidate, candidate, pm);
        mpz_set(res, am);

        candidatefactors = isBsmooth(candidate, Bm);

        if (candidatefactors != NULL && mpz_cmp(candidate, uno) != 0) {

            struct element *factor, *prime;
            mpz_set_ui(result, 0);
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
                    if (mpz_cmp_ui(*(matrix + row * ROWSIZE + pos), 1) == 0) {

                        for (times = 0; times < factor->exp; times++) {
                            mpz_add(result, result, *(matrix + row * ROWSIZE + ROWSIZE - 1));

                        }
                    }
                }
                pos = 0;
            }
            freelist(candidatefactors);
            mpz_sub(result, result, u);
            mpz_mod(result, result, pm1);
            printf("done\n\nIl risultato è: ");
            mpz_out_str(stdout, 10, result);
            break;
        }
    }

    return EXIT_SUCCESS;
}




