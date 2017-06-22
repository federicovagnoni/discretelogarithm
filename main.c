#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <zconf.h>
#include <mpfr.h>
#include <pthread.h>

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

void pollard_rho(mpz_t factor, mpz_t n) {
    mpz_t x, y, d, sub, a;
    mpz_init(sub);
    mpz_init(a);
    mpz_init_set_ui(x, 1);
    mpz_init_set_ui(y, 1);
    mpz_init_set_ui(d, 1);
    gmp_randstate_t randstate;
    gmp_randinit_mt(randstate);
    gmp_randseed_ui(randstate, (unsigned long) time(NULL));
    choose_seed(x, y, a, n, &randstate);

    while (mpz_cmp_ui(d, 1) == 0) {
        g(x, x, a, n);
        //printf("x = %.0lf\n", x);
        g(y, y, a, n);
        g(y, y, a, n);
        //printf("y = %.0lf\n", y);
        mpz_sub(sub, x, y);
        mpz_abs(sub, sub);
        mpz_gcd(d, sub, n);
        //printf("d = %.0lf\n", d);
        if (mpz_cmp(d, n) == 0) {
            mpz_clear(x);
            mpz_clear(y);
            mpz_clear(d);
            mpz_clear(sub);
            mpz_set_si(factor, -1);
            return;
        }
        if (mpz_cmp_ui(d, 1) == 0) {
            continue;
        }
    }
    mpz_set(factor, d);
    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(d);
    mpz_clear(sub);
    return;
}


void pollard_rho_opt(mpz_t factor, mpz_t n, mpz_t B) {
    mpz_t x, y, a, value, gcd, counter, root;
    mpz_init(x);
    mpz_init(y);
    mpz_init(a);
    mpz_init(value);
    mpz_init(gcd);
    mpz_init_set_ui(counter, 0);
    mpz_init_set(root, B);
//    mpz_sqrt(root, B);
//    mpz_add_ui(root, root, 15);
    gmp_randstate_t randstate;
    gmp_randinit_mt(randstate);
    gmp_randseed_ui(randstate, (unsigned long) time(NULL));
    choose_seed(x, y, a, n, &randstate);
    while (1 == 1) {
        g(x, x, a, n);
        g(y, y, a, n);
        g(y, y, a, n);
        mpz_sub(value, x, y);
        mpz_abs(value, value);
        mpz_gcd(gcd, value, n);
        if (mpz_cmp_ui(gcd, 1) == 0) {
            if (mpz_cmp(counter, root) == 0) {
//                printf("Guarda tanto non è Bsmooth...\n");
//                fflush(stdout);
                mpz_set_si(factor, -2);
                mpz_clear(x);
                mpz_clear(y);
                mpz_clear(a);
                mpz_clear(value);
                mpz_clear(gcd);
                mpz_clear(counter);
                mpz_clear(root);
                gmp_randclear(randstate);
                return;
            }
            mpz_add_ui(counter, counter, 1);
            choose_seed(x, y, a, n, &randstate);
            continue;
        }
        if (mpz_cmp(gcd, n) == 0) {
            mpz_set(factor, n);
            mpz_clear(x);
            mpz_clear(y);
            mpz_clear(a);
            mpz_clear(value);
            mpz_clear(gcd);
            mpz_clear(counter);
            mpz_clear(root);
            gmp_randclear(randstate);
            return;
        }
        mpz_set(factor, gcd);
        mpz_clear(x);
        mpz_clear(y);
        mpz_clear(a);
        mpz_clear(value);
        mpz_clear(gcd);
        mpz_clear(counter);
        mpz_clear(root);
        gmp_randclear(randstate);
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
    mpz_t r, root;
    mpz_init(r);
    mpz_init(root);
    mpz_sqrt(root, n);
    while (temp != NULL && mpz_cmp(temp->value, root) <= 0) {
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
        if (mpz_cmp(factor, todivide) == 0) {
            addnewfactor(factor, list);
            mpz_clear(factor);
            mpz_clear(todivide);
            mpz_clear(root);
            mpz_clear(uno);
            return;
        }
        addnewfactor(factor, list);
        mpz_div(todivide, todivide, factor);
    }
    mpz_clear(factor);
    mpz_clear(todivide);
    mpz_clear(root);
    mpz_clear(uno);
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
//        mpz_sqrt(root, factor);
//        primelist = primes_in_range(uno, root);
//        factorsbytrialdivision(factor, newlist, primelist);
//        freelist(primelist);
        mpz_div(todivide, todivide, factor);

    }
    mpz_clear(factor);
    mpz_clear(todivide);
    mpz_clear(mod);
    mpz_clear(uno);
    mpz_clear(due);
    return newlist;
}

struct list *isBsmooth2(mpz_t n, struct list *primelist, mpz_t B) {

    struct list *newlist = init_list();
    struct element *elem;
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

    mpz_clear(mod);
    mpz_clear(due);

    while (mpz_cmp(todivide, uno) != 0) {

        pollard_rho_opt(factor, todivide, B);

        if (mpz_cmp_si(factor, -2) == 0) {
            mpz_clear(factor);
            mpz_clear(todivide);
            mpz_clear(uno);
            freelist(newlist);
            return NULL;
        }
//        printf("Trovato fattore");
//        mpz_out_str(stdout, 10, factor);
//        printf("\n");
        factorsbytrialdivision(factor, newlist, primelist);

        for (elem = newlist->TAIL; elem != NULL; elem = elem->prev) {
            if (mpz_cmp(elem->value, B) > 0) {
                mpz_clear(factor);
                mpz_clear(todivide);
                mpz_clear(uno);
                freelist(newlist);
                return NULL;
            }
        }

        mpz_div(todivide, todivide, factor);

    }
    mpz_clear(factor);
    mpz_clear(todivide);
    mpz_clear(uno);
    return newlist;
}


struct list *isBsmooth(mpz_t m, mpz_t B) {

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
//        mpz_out_str(stdout, 10, lastx);
//        printf(" * ");
//        mpz_out_str(stdout, 10, a_origin);
//        printf(" + ");
//        mpz_out_str(stdout, 10, lasty);
//        printf(" * ");
//        mpz_out_str(stdout, 10, b_origin);
//        printf(" = ");
//        mpz_out_str(stdout, 10, a);
//        printf("\n");
        if (mpz_cmp_ui(a, 1) == 0) {
            mpz_set_ui(x, 0);
            mpz_set_ui(y, 0);
            mpz_set_ui(z, 0);
            return;
        }
//        mpz_invert(invers, lasty, a_origin);

        if (mpz_sizeinbase(lasty, 10) == rightsize || mpz_sizeinbase(a, 10) <= rightsize) {
            if (mpz_cmp_ui(lasty, 0) > 0) {
//            printf("x = ");
//            mpz_out_str(stdout, 10, lastx);
                mpz_set(x, lastx);
//            printf("\ny = ");
//            mpz_out_str(stdout, 10, lasty);
                mpz_set(y, lasty);
//            printf("\ngcd = ");
//            mpz_out_str(stdout, 10, a);
                mpz_set(z, a);
//            printf("\n");
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

pthread_rwlock_t primelistlock = PTHREAD_RWLOCK_INITIALIZER;
pthread_rwlock_t randlock = PTHREAD_RWLOCK_INITIALIZER;

pthread_mutex_t counter_mutex = PTHREAD_MUTEX_INITIALIZER;

gmp_randstate_t randstate;

int row_counter = 0;

struct tparams {
    mpz_t pm;
    mpz_t qm;
    mpz_t am;
    mpz_t Bm;
    mpz_t *matrix;
    int ROWSIZE;
    struct list *primelist;
};


void *findrelation(void *arg) {

    struct tparams *params = (struct tparams *) arg;
    mpz_t x, y, gcd, u, candidate, invers, c, term;
    int j, k, z, l, isredundant, equals;
    struct element *temp1, *temp2;

    int RELATIONS = params->ROWSIZE - 1;

    struct list *candidatefactorsofgcd;
    struct list *candidatefactorsofy;

    mpz_init(u);
    mpz_init(candidate);


    for (;;) {

        pthread_rwlock_wrlock(&randlock);
        mpz_urandomm(u, randstate, params->pm);
        pthread_rwlock_unlock(&randlock);

        mpz_add_ui(u, u, 1);
        mpz_powm(candidate, params->am, u, params->pm);

        if (mpz_sizeinbase(candidate, 10) < mpz_sizeinbase(params->pm, 10)  /*|| useless == 1*/) {

            mpz_init(x);
            mpz_init(y);
            mpz_init(gcd);

            gcdext(params->pm, candidate, x, y, gcd);

            if (mpz_cmp_ui(x, 0) == 0 && mpz_cmp_ui(y, 0) == 0 && mpz_cmp_ui(gcd, 0) == 0) {
                mpz_clear(x);
                mpz_clear(y);
                mpz_clear(gcd);
                continue;
            }

//            printf("gcd = ");
//            mpz_out_str(stdout, 10, gcd);
//            printf("\n");
//            printf("x = ");
//            mpz_out_str(stdout, 10, x);
//            printf("\n");
//            printf("y = ");
//            mpz_out_str(stdout, 10, y);
//            printf("\n");
//            printf("candidate = ");
//            mpz_out_str(stdout, 10, candidate);
//            printf(" pari a ");
//            mpz_out_str(stdout, 10, params->am);
//            printf(" ^ ");
//            mpz_out_str(stdout, 10, u);
//            printf("\n");
//            fflush(stdout);

            pthread_rwlock_rdlock(&primelistlock);

            candidatefactorsofgcd = isBsmooth2(gcd, params->primelist, params->Bm);
            candidatefactorsofy = isBsmooth2(y, params->primelist, params->Bm);

            pthread_rwlock_unlock(&primelistlock);

            // if it is
            if (candidatefactorsofgcd != NULL && candidatefactorsofy != NULL) {

                pthread_mutex_lock(&counter_mutex);

                if (row_counter == RELATIONS) {
                    pthread_mutex_unlock(&counter_mutex);
                    freelist(candidatefactorsofgcd);
                    freelist(candidatefactorsofy);
                    mpz_clear(x);
                    mpz_clear(y);
                    mpz_clear(gcd);
                    mpz_clear(u);
                    mpz_clear(candidate);
                    pthread_exit(NULL);
                }
//                printf("SONO BSMOOTH!!!");
//                fflush(stdout);
//                printlist(candidatefactorsofgcd);
//                printf("\n");
//                printlist(candidatefactorsofy);
//                printf("\n");
//                printf("gcd = ");
//                mpz_out_str(stdout, 10, gcd);
//                printf("\n");
//                printf("x = ");
//                mpz_out_str(stdout, 10, x);
//                printf("\n");
//                printf("y = ");
//                mpz_out_str(stdout, 10, y);
//                printf("\n");
//                printf("candidate = ");
//                mpz_out_str(stdout, 10, candidate);
//                printf(" pari a ");
//                mpz_out_str(stdout, 10, params->am);
//                printf(" ^ ");
//                mpz_out_str(stdout, 10, u);
//                printf("\n");


                // allocate the right amount of memory by adding a new row of mpz_t.
                // A row is formed by the coefficient related to a prime in the factor base
                // so the row is B/ln(B) numbers long + 1 for the exponent used

                // initialize all the elements to 0
                for (j = 0; j < params->ROWSIZE; j++) {
                    mpz_init_set_ui(*(params->matrix + row_counter * params->ROWSIZE + j), 0);
                }


                j = 0;
                // for every element in the factor base
                pthread_rwlock_rdlock(&primelistlock);
                for (temp1 = params->primelist->HEAD; temp1 != NULL; temp1 = temp1->next) {
                    // for every factor of the candidate (which is B-smooth)
                    for (temp2 = candidatefactorsofgcd->HEAD; temp2 != NULL; temp2 = temp2->next) {
                        // if the primes are equal
                        if (mpz_cmp(temp1->value, temp2->value) == 0) {
                            // set the related coefficient in the matrix equals to the exponent of that factor
                            mpz_set_si(*(params->matrix + row_counter * params->ROWSIZE + j), temp2->exp);
                        }
                    }
                    j++;
                }

                j = 0;
                for (temp1 = params->primelist->HEAD; temp1 != NULL; temp1 = temp1->next) {
                    // for every factor of the candidate (which is B-smooth)
                    for (temp2 = candidatefactorsofy->HEAD; temp2 != NULL; temp2 = temp2->next) {
                        // if the primes are equal
                        if (mpz_cmp(temp1->value, temp2->value) == 0) {
                            // set the related coefficient in the matrix equals to the exponent of that factor
                            mpz_sub_ui(*(params->matrix + row_counter * params->ROWSIZE + j), *(params->matrix + row_counter * params->ROWSIZE + j),
                                       (unsigned long) temp2->exp);
                            mpz_mod(*(params->matrix + row_counter * params->ROWSIZE + j), *(params->matrix + row_counter * params->ROWSIZE + j), params->qm);
                        }
                    }
                    j++;
                }
                // set the last element of the row equals to the exponent of the primitive root
                mpz_set(*(params->matrix + row_counter * params->ROWSIZE + params->ROWSIZE - 1), u);

                freelist(candidatefactorsofgcd);
                freelist(candidatefactorsofy);

                pthread_rwlock_unlock(&primelistlock);

                // If it's not redundat I can proceed with the Gauss elimination
                isredundant = 0;

                mpz_init(invers);
                mpz_init(c);
                mpz_init(term);
                for (k = 0; k < row_counter + 1; k++) {
                    // if the coefficient in the position [k, k] is not invertible then the row is useless
                    if (mpz_invert(invers, *(params->matrix + k * params->ROWSIZE + k), params->qm) == 0) {
                        isredundant = 1;
                        break;
                    }
                    for (z = k + 1; z < row_counter + 1; z++) {
                        mpz_mul(c, *(params->matrix + z * params->ROWSIZE + k), invers);
                        for (l = 0; l < params->ROWSIZE; l++) {
                            mpz_mul(term, c, *(params->matrix + k * params->ROWSIZE + l));
                            mpz_sub(*(params->matrix + z * params->ROWSIZE + l), *(params->matrix + z * params->ROWSIZE + l), term);
                            mpz_mod(*(params->matrix + z * params->ROWSIZE + l), *(params->matrix + z * params->ROWSIZE + l), params->qm);
                        }
                    }
                }
                mpz_clear(invers);
                mpz_clear(c);
                mpz_clear(term);


                if (isredundant == 1) {
                    for (j = 0; j < params->ROWSIZE; j++) {
                        mpz_clear(*(params->matrix + row_counter * params->ROWSIZE + j));
                    }
                    pthread_mutex_unlock(&counter_mutex);
                    continue;
                }



                equals = 0;

                // If I obtain a row equals to another, than it is useless
                if (row_counter > 0) {
                    for (l = row_counter - 1; l >= 0; l--) {
                        for (j = 0; j < RELATIONS; j++) {
                            if (mpz_cmp(*(params->matrix + row_counter * params->ROWSIZE + j), *(params->matrix + l * params->ROWSIZE + j)) == 0) {
                                equals++;
                                if (equals == RELATIONS) {
                                    isredundant = 1;
                                }
                            }
                        }
                        equals = 0;
                    }
                }

                if (isredundant == 1) {
                    for (j = 0; j < params->ROWSIZE; j++) {
                        mpz_clear(*(params->matrix + row_counter * params->ROWSIZE + j));
                    }
                    pthread_mutex_unlock(&counter_mutex);
                    continue;
                }

//                for (k = 0; k < row_counter + 1; k++) {
//                    printf("[ ");
//                    for (j = 0; j < params->ROWSIZE; j++) {
//                        mpz_out_str(stdout, 10, *(params->matrix + k * params->ROWSIZE + j));
//                        printf(" ");
//                    }
//                    printf("]\n");
//                }
//                printf("\n");

                printf("%d di %d relazioni trovate\n", row_counter + 1, RELATIONS);
                row_counter++;
                mpz_clear(x);
                mpz_clear(y);
                mpz_clear(gcd);

                if (row_counter == RELATIONS) {
                    pthread_mutex_unlock(&counter_mutex);
                    mpz_clear(u);
                    mpz_clear(candidate);
                    pthread_exit(NULL);
                }

                pthread_mutex_unlock(&counter_mutex);
            }
        }
    }
}


double main(int argc, char *argv[]) {

//    mpz_t b, e, primo;
//    mpz_init_set_ui(b, 6);
//    mpz_init_set_ui(e, 25257826019);
//    mpz_init_set_ui(primo, 98764321261);
//    mpz_powm(b, b, e, primo);
//    mpz_out_str(stdout, 10, b);


    mpz_t am, pm, pm1, exp, try, qm, uno, due, Bm, tm, pm2, p2m;
    mpfr_t logpm, mpfrpm, mpfrexponent;

    mpz_init(am);
    mpz_init(exp);
    mpz_init(try);
    mpz_init(pm);
    mpz_init(pm1);
    mpz_init(pm2);
    mpfr_init(logpm);
    mpz_init(p2m);
    mpz_init(qm);
    mpz_init(tm);

    printf("Inserisci il numero da fattorizzare: ");
    mpz_inp_str(qm, STDIN_FILENO, 10);
    printf("q = ");
    mpz_out_str(stdout, 10, qm);
    printf("\n");

//    mpz_init_set_ui(qm, q);
    mpz_init_set_ui(uno, 1);
    mpz_init_set_ui(due, 2);
    mpfr_init(mpfrpm);
    mpz_mul(pm, qm, due);
    mpz_add(pm, pm, uno);
    printf("p = ");
    mpz_out_str(stdout, 10, pm);
    printf("\n");

    mpz_t exponent, lp;
    mpz_init(exponent);
    mpz_init(lp);
    mpfr_t mpfrlp, tre;

    // Find B-best
    mpfr_init(mpfrlp);
    mpfr_rnd_t rnd;
    mpfr_set_default_rounding_mode(rnd);
    mpfr_init_set_ui(tre, 3, rnd);
    mpfr_set_z(mpfrpm, pm, rnd);
    mpfr_log(logpm, mpfrpm, rnd);
    mpfr_t loglogpm, product;
    mpfr_init(loglogpm);
    mpfr_init(product);

    mpfr_log(loglogpm, logpm, rnd);
    mpfr_mul(product, logpm, loglogpm, rnd);
    mpfr_div(product, product, tre, rnd);

    mpfr_init(mpfrexponent);

    mpfr_sqrt(mpfrexponent, product, rnd);
    mpfr_get_z(exponent, mpfrexponent, rnd);


    printf("exponent = ");
    mpz_out_str(stdout, 10, exponent);
    printf("\n");

    mpfr_exp(mpfrlp, mpfrexponent, rnd);
    mpfr_get_z(lp, mpfrlp, rnd);

    printf("L(p) = ");
    mpz_out_str(stdout, 10, lp);
    printf("\n");

    mpz_init_set(Bm, lp);
    printf("B = ");
    mpz_out_str(stdout, 10, Bm);
    printf("\n");


    mpz_sub_ui(pm1, pm, 1);

//    mpz_t factord, number;
//    mpz_init(factord);
//    mpz_init_set_ui(number, 160);
//
//    struct list *bau = isBsmooth2(number, Bm);
//    exit(EXIT_SUCCESS);

    int count = 0;
    int i = 0;
    int j, l, k, z;

    struct element *divisor;
    printf("Calcolo i divisori di p - 1...");
    fflush(stdout);
    struct list *divisors = init_list();/*= factorsbypollard(pm1)*/;
    printf("done\n\n");
    fflush(stdout);

    printf("p - 1 = ");
    addelem(qm, divisors);
    addelem(due, divisors);
    printlist(divisors);
    printf("\n");
    fflush(stdout);

    // Find a primitive root for Z_{p} (valid for p prime!)
    for (mpz_set_ui(am, 25); mpz_cmp(am, pm) < 0; mpz_add_ui(am, am, 1)) {
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
    mpz_inp_str(am, STDIN_FILENO, 10);

    printf("Inserisci l'argomento del logaritmo: ");
    mpz_inp_str(tm, STDIN_FILENO, 10);

    time_t start = time(NULL);
    // Calculate primes from 1 to B: it will be our factor base
    struct list *primelist = primes_in_range(uno, Bm);
    struct list *candidatefactorsofgcd;
    struct list *candidatefactorsofy;
//    struct list *factorsofam;
//    factorsofam = factorsbypollard(am);
    struct element *temp1, *temp2;




    // convert the dinamic list of mpz_t in a common array of mpz_t
    mpz_t listprime[primelist->count];
    returnlist(primelist, listprime);

    mpz_sub_ui(pm2, pm, 2);
    mpz_mul_ui(p2m, pm, 2);

    int RELATIONS = primelist->count;
    int ROWSIZE = primelist->count + 1;

    // Initialize matrix of coefficients
    mpz_t *matrix = malloc(RELATIONS * sizeof(mpz_t) * ROWSIZE);

    int isredundant = 0, equals = 0, useless = 0;

    mpz_t candidate, res, u, invers, c, term, x, y, gcd;
    mpz_init(res);
    mpz_init(candidate);

    mpz_init(u);


    // Initialize variable for random exponents
    gmp_randinit_mt(randstate);
    gmp_randseed_ui(randstate, (unsigned long) time(NULL));

//    mpz_t maxfactor;
//    mpz_init_set_ui(maxfactor, 1);
//
//
//    // Find the special relation
//    for (;;) {
//
//        mpz_urandomm(u, randstate, pm);
//        mpz_powm(candidate, am, u, pm);
//        mpz_mul(candidate, candidate, tm);
//        mpz_mod(candidate, candidate, pm);
//
//
//        if (mpz_sizeinbase(candidate, 10) < mpz_sizeinbase(pm, 10)  /*|| useless == 1*/) {
//
//            mpz_init(x);
//            mpz_init(y);
//            mpz_init(gcd);
//
//            gcdext(pm, candidate, x, y, gcd);
//
//            if (mpz_cmp_ui(x, 0) == 0 && mpz_cmp_ui(y, 0) == 0 && mpz_cmp_ui(gcd, 0) == 0) {
//                mpz_clear(x);
//                mpz_clear(y);
//                mpz_clear(gcd);
//            }
//
//
//            candidatefactorsofgcd = isBsmooth2(gcd, primelist, Bm);
//            candidatefactorsofy = isBsmooth2(y, primelist, Bm);
//
//
//            if (candidatefactorsofgcd != NULL && candidatefactorsofy != NULL) {
//                printf("Fattori di gcd: ");
//                printlist(candidatefactorsofgcd);
//                printf("\nFattori di y: ");
//                printlist(candidatefactorsofy);
//                printf("\n");
//                fflush(stdout);
//                for (factor = candidatefactorsofgcd->HEAD; factor != NULL; factor = factor->next) {
//                    addnewfactor(factor->value, candidatefactors);
//                    if (mpz_cmp(factor->value, maxfactor) > 0) {
//                        mpz_set(maxfactor, factor->value);
//                    }
//                    for (temp1 = candidatefactors->HEAD; temp1 != NULL; temp1 = temp1->next) {
//                        if (mpz_cmp(factor->value, temp1->value) == 0) {
//                            temp1->exp = factor->exp;
//                        }
//                    }
//                }
//
//                for (factor = candidatefactorsofy->HEAD; factor != NULL; factor = factor->next) {
//                    if (mpz_cmp(factor->value, maxfactor) > 0) {
//                        mpz_set(maxfactor, factor->value);
//                    }
//                    addnewfactor(factor->value, candidatefactors);
//                    for (temp1 = candidatefactors->HEAD; temp1 != NULL; temp1 = temp1->next) {
//                        if (mpz_cmp(factor->value, temp1->value) == 0) {
//                            temp1->exp = temp1->exp - factor->exp - 1;
//                        }
//                    }
//                }
//
//                printf("Fattori risultanti: ");
//                printlist(candidatefactors);
//                printf("\n");
//
//                mpz_set(Bm, maxfactor);
//                break;
//            }
//            mpz_clear(x);
//            mpz_clear(y);
//            mpz_clear(gcd);
//        }
//    }
//
//    freelist(primelist);
//    primelist = primes_in_range(uno, Bm);
//
//    RELATIONS = primelist->count;
//    ROWSIZE = primelist->count + 1;

    printf("Tra 1 e ");
    mpz_out_str(stdout, 10, Bm);
    printf(" ci sono %d numeri primi\n", RELATIONS);

    int MAX_THREAD = 4;
    pthread_t tids[MAX_THREAD];
    struct tparams *params[MAX_THREAD];

    for (;;) {
        for (int f = 0; f < MAX_THREAD; f++) {
            // put into a function
            params[f] = realloc(NULL, sizeof(struct tparams));
            mpz_init_set(params[f]->am, am);
            mpz_init_set(params[f]->Bm, Bm);
            mpz_init_set(params[f]->pm, pm);
            mpz_init_set(params[f]->qm, qm);
            params[f]->matrix = matrix;
            params[f]->ROWSIZE = ROWSIZE;
            params[f]->primelist = primelist;

            if (pthread_create(&tids[f], NULL, findrelation, (void *) params[f]) != 0) {
                fprintf(stderr, "Error in pthread_create()\n");
                exit(EXIT_FAILURE);
            }

        }
        for (int f = 0; f < MAX_THREAD; f++) {
            if (pthread_join(tids[f], NULL) != 0) {
                fprintf(stderr, "Error in pthread_join()\n");
                exit(EXIT_FAILURE);
            }
        }


//        for (k = 0; k < RELATIONS; k++) {
//            printf("[ ");
//            for (j = 0; j < ROWSIZE; j++) {
//                mpz_out_str(stdout, 10, *(matrix + k * ROWSIZE + j));
//                printf(" ");
//            }
//            printf("]\n");
//        }
//        printf("\n");

        break;
    }

    mpz_init(invers);
    mpz_init(c);
    mpz_init(term);
    // Initialize the vector of solutions
    mpz_t sol[RELATIONS];
    for (j = 0; j < RELATIONS; j++) {
        mpz_init(sol[j]);
    }

    // We will use the backward substitution
//    printf("Inverto ");
//    mpz_out_str(stdout, 10, *(matrix + (RELATIONS - 1) * ROWSIZE + (RELATIONS - 1)));
//    printf(" e viene: ");
    mpz_invert(invers, *(matrix + (RELATIONS - 1) * ROWSIZE + (RELATIONS - 1)), qm);
//    mpz_out_str(stdout, 10, invers);
//    printf("\n");
    mpz_mul(sol[RELATIONS - 1], *(matrix + (RELATIONS - 1) * ROWSIZE + RELATIONS), invers);
//    printf("e lo moltiplico per ");
//    mpz_out_str(stdout, 10, *(matrix + (RELATIONS - 1) * ROWSIZE + RELATIONS));
//    printf("\n");
    mpz_mod(sol[RELATIONS - 1], sol[RELATIONS - 1], qm);
//    printf("viene");
//    mpz_out_str(stdout, 10, sol[RELATIONS - 1]);
//    printf("\n");

    mpz_t sum, num, denom;
    mpz_init(num);
    mpz_init(denom);
    mpz_init(sum);

    /* this loop is for backward substitution*/
    printf("Calcolo della matrice dei coefficienti...");
    for (i = RELATIONS - 2; i >= 0; i--) {
        mpz_set_ui(sum, 0);
        for (j = i + 1; j < RELATIONS; j++) {
            mpz_addmul(sum, sol[j], *(matrix + i * ROWSIZE + j));
            // sum = sum + x[j] * *(matrix + i * ROWSIZE + j))
        }
        mpz_invert(denom, *(matrix + i * ROWSIZE + i), qm);
        mpz_sub(num, *(matrix + i * ROWSIZE + RELATIONS), sum);
        mpz_mul(sol[i], num, denom);
        mpz_mod(sol[i], sol[i], qm);
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

//    for (k = 0; k < RELATIONS; k++) {
//        printf("[ ");
//        for (j = 0; j < ROWSIZE; j++) {
//            mpz_out_str(stdout, 10, *(matrix + k * ROWSIZE + j));
//            printf(" ");
//        }
//        printf("]\n");
//    }
//    printf("\n");

    mpz_t value;

    for (k = 0; k < RELATIONS; k++) {
        mpz_init(value);
        mpz_powm(value, am, *(matrix + k * ROWSIZE + ROWSIZE - 1), pm);
//        mpz_out_str(stdout, 10, value);
//        printf("\n");
        if (mpz_cmp(value, listprime[k]) != 0) {
            mpz_add(*(matrix + k * ROWSIZE + ROWSIZE - 1), *(matrix + k * ROWSIZE + ROWSIZE - 1), qm);
        }
        mpz_clear(value);

    }

//    for (k = 0; k < RELATIONS; k++) {
//        printf("[ ");
//        for (j = 0; j < ROWSIZE; j++) {
//            mpz_out_str(stdout, 10, *(matrix + k * ROWSIZE + j));
//            printf(" ");
//        }
//        printf("]\n");
//    }
//    printf("\n");


    struct list *candidatefactors;
    struct element *factor, *prime;
    int pos;
    int row;
    int times;

    mpz_t result;
    mpz_init(result);

    // I have the complete matrix of coefficient
    // Now I have to multiplicate the argument of the logarithm with a random power of the primitive root
    // until I find a B-smooth number. The I can resolve the system
    printf("Calcolo della soluzione...");
    fflush(stdout);


    for (mpz_set_ui(u, 1); mpz_cmp(u, pm1) <= 0; mpz_add_ui(u, u, 1)) {

        mpz_powm(candidate, am, u, pm);
        mpz_mul(candidate, candidate, tm);
        mpz_mod(candidate, candidate, pm);
        mpz_set(res, am);

        candidatefactors = isBsmooth2(candidate, primelist, Bm);

        if (candidatefactors != NULL && mpz_cmp(candidate, uno) != 0) {

            mpz_set_ui(result, 0);

            pos = 0;
            // for every factor of the number obtained
            for (factor = candidatefactors->HEAD; factor != NULL; factor = factor->next) {
                // for every prime in the factor base
                for (prime = primelist->HEAD; prime != NULL; prime = prime->next) {
                    // if the numbers are equal stop the cycle
                    if (mpz_cmp(prime->value, factor->value) == 0) {
                        break;
                    }
                    // else increment the position
                    pos++;
                }
                // cycle on the row
                for (row = 0; row < primelist->count; row++) {
                    // if I find a 1 in position [row, pos]
                    if (mpz_cmp_ui(*(matrix + row * ROWSIZE + pos), 1) == 0) {
                        // I add factor->exp times the related exponent to solution
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
            printf("\nEseguito in %ld\n", time(NULL) - start);
            break;
        }
    }



    return EXIT_SUCCESS;
}




