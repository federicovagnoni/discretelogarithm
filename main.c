#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
//#include <pari/pari.h>

double g(double x, double n) {
    return fmod((pow(x, 2) + 1), n);
}

/* Standard C Function: Greatest Common Divisor */
double gcd(double a, double b) {
    double temp;
    while (b != 0) {
        temp = fmod(a, b);
        a = b;
        b = temp;
    }
    return a;
}

double pollard_rho(double x0, double n) {
    double x = x0;
    double y = 2;
    double d = 1;
    while (d == 1) {
        x = g(x, n);
        //printf("x = %.0lf\n", x);
        y = g(g(y, n), n);
        //printf("y = %.0lf\n", y);
        d = gcd(fabs(x - y), n);
        //printf("d = %.0lf\n", d);
    }
    if (d == n)
        return -1;
    else
        return d;
}

double isprime(double n) {
    double x;
    double x0 = 2;
    int i = 0;
    while (n != 1) {
        x = pollard_rho(x0, n);
        if (x == -1) {
//            i++;
//            x0 = x0 * -2;
//            //printf("Cambio x0: %.0lf (n = %.0lf)\n", x0, n);
//            if (i == 10) {
//                //printf("%.0lf  end\n", n);
            break;
//            }
        } else {
            return 0;
            x0 = 2;
            //printf("%.0lf   %.0lf\n", n, x);
            n = n / x;
        }
    }
    return 1;
}


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
        addnewfactor(factor, newlist);
        mpz_div(todivide, todivide, factor);
    }
    mpz_clear(factor);
    mpz_clear(todivide);
    mpz_clear(root);
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
    B = B + B * 3 / 2;

    printf("B = %llu\n", B);

    printf("Inserisci una radice primitiva (e base del logaritmo): ");
    scanf("%llu", &a);

    printf("Inserisci l'argomento del logaritmo: ");
    scanf("%llu", &t);
    mpz_t Bm, am, pm, tm;
    mpz_init_set_ui(Bm, B);
    mpz_init_set_ui(am, a);
    mpz_init_set_ui(pm, p);
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
//    mpz_out_str(stdout,10,logpm);
//    fflush(stdout);

    int i = 0;
    mpz_t randexp;
    srand(time(NULL));

    mpz_t candidate;
    mpz_t res;
    mpz_init(res);

    struct list *candidatefactors;
    gmp_randstate_t randstate;
    gmp_randinit_default(randstate);

//    struct list *newlist;
    struct element *temp1, *temp2;

    mpz_t *matrix = NULL;
    int k, u, z;
    mpz_t base;
    mpz_init(base);
    mpz_init(randexp);
    mpz_init(candidate);
    int y, s;
    mpz_t randexp1, randexp2;
    mpz_t candidate1, candidate2;
    mpz_init(randexp1);
    mpz_init(randexp2);
    mpz_init(candidate1);
    mpz_init(candidate2);

    mpz_t temp, pm1;
    mpz_init(temp);
    mpz_init(pm1);
    mpz_sub_ui(pm1, pm, 1);
    mpz_t invers;
    mpz_init(invers);

    int SIZE = primelist->count;
    int ROWSIZE = primelist->count * 2;
    mpz_t c;
    mpz_t x[SIZE];
    mpz_init(c);


    while (i < primelist->count) {
        for (u = 0; u < primelist->count; u++) {
            mpz_set(base, listprime[u]);
            mpz_set(res, base);
//            printf("base = ");
//            mpz_out_str(stdout,10, base);
//            printf("\n");
            //mpz_set_ui(randexp, rand() % (p - 2));
            for (z = 1; mpz_cmp_ui(logpm, z) >= 0; z++) {
//                printf("Exp = %d\n", z);
                //mpz_urandomm(randexp, randstate, logpm);
//                printf("exp = %d\n", z);
                mpz_set_ui(randexp, z);
//                mpz_add_ui(randexp, randexp, 1);
                //        mpz_set(base, listprime[rand() % (primelist->count)]);
                mpz_powm(candidate, base, randexp, pm);

//                printf("candidate = ");
//                mpz_out_str(stdout,10,candidate);
//                printf("\n");

                for (k = 0; k < (int) mpz_get_ui(randexp) - 1; k++) {
                    mpz_mul(res, res, base);
                }
//                printf("Res = ");
//                mpz_out_str(stdout,10,res);
//                printf("\n");



                if (isBsmooth(candidate, Bm) == 1 && mpz_cmp(candidate, uno) != 0 && mpz_cmp(res, pm) > 0) {
                    mpz_out_str(stdout, 10, candidate);
                    printf(" (pari a ");
                    mpz_out_str(stdout, 10, base);
                    printf(" alla ");
                    mpz_out_str(stdout, 10, randexp);
                    printf(") è %llu-smooth\n", B);
                    fflush(stdout);
                    candidatefactors = factorsbytrialdivision(candidate);


                    matrix = realloc(matrix, sizeof(mpz_t) * (primelist->count) * (i + 1) * 2);

                    int j;
                    for (j = 0; j < primelist->count * 2; j++) {
                        mpz_init_set_ui(*(matrix + i * (primelist->count * 2) + j), 0);
                    }


                    j = 0;
                    for (temp1 = primelist->HEAD; temp1 != NULL; temp1 = temp1->next) {
                        for (temp2 = candidatefactors->HEAD; temp2 != NULL; temp2 = temp2->next) {
                            if (mpz_cmp(temp1->value, temp2->value) == 0) {
                                mpz_set_ui(*(matrix + i * primelist->count * 2 + j), (unsigned long) temp2->exp);
                            }
                        }
                        j++;
                    }
                    j = primelist->count;
                    for (temp1 = primelist->HEAD; temp1 != NULL; temp1 = temp1->next) {
                        if (mpz_cmp(temp1->value, base) == 0) {
//                                printf("Aggiungo %d in posizione %d\n", (int) mpz_get_ui(randexp), i *primelist->count * 2 + j);
                            mpz_set(*(matrix + i * primelist->count * 2 + j), randexp);
                        }
                        j++;
                    }

                    for (j = primelist->count; j < primelist->count * 2; j++) {
                        if (j >= primelist->count) {
                            mpz_sub(*(matrix + i * primelist->count * 2 + (j - primelist->count)),
                                    *(matrix + i * primelist->count * 2 + (j - primelist->count)),
                                    *(matrix + i * primelist->count * 2 + j));
                            mpz_mod(*(matrix + i * primelist->count * 2 + (j - primelist->count)),
                                    *(matrix + i * primelist->count * 2 + (j - primelist->count)), pm1);
                            mpz_set_ui(*(matrix + i * primelist->count * 2 + j), 0);
                        }
                    }

                    mpz_invert(invers, *(matrix + i * primelist->count * 2 + i), pm1);
                    if (mpz_cmp_ui(invers, 0) == 0) {
                        continue;
                    }

                    //*(matrix + i * (primelist->count + 1) + primelist->count * 2 - 1) = (int) mpz_get_ui(randexp);
                    freelist(candidatefactors);
                    i++;

                }
                if (i == primelist->count) {
                    break;
                }
                mpz_set(res, base);
            }
        }



        // search for non trivial relations
        for (u = 0; u < primelist->count; u++) { //first factor
            for (k = u + 1; k < primelist->count; k++) { //second factor
                for (y = 1; mpz_cmp_ui(logpm, y) >= 0; y++) { //first exponent
                    for (z = 1; mpz_cmp_ui(logpm, z) >= 0; z++) { //second exponent
                        mpz_set_ui(randexp1, y);
                        mpz_set_ui(randexp2, z);
                        mpz_powm(candidate1, listprime[u], randexp1, pm);
                        mpz_set(res, listprime[u]);
                        mpz_powm(candidate2, listprime[k], randexp2, pm);
                        mpz_mul(candidate, candidate1, candidate2);
                        mpz_mod(candidate, candidate, pm);
                        for (s = 0; s < (int) mpz_get_ui(randexp1) - 1; s++) {
                            mpz_mul(res, res, listprime[u]);
                        }
                        for (s = 0; s < (int) mpz_get_ui(randexp2) - 1; s++) {
                            mpz_mul(res, res, listprime[k]);
                        }
                        if (isBsmooth(candidate, Bm) == 1 && mpz_cmp(candidate, uno) != 0 && mpz_cmp(res, pm) > 0 &&
                            mpz_cmp(res, p2m) < 0) {
                            mpz_out_str(stdout, 10, candidate);
                            printf(" (pari a ");
                            mpz_out_str(stdout, 10, listprime[u]);
                            printf(" alla ");
                            mpz_out_str(stdout, 10, randexp1);
                            printf(" * ");
                            mpz_out_str(stdout, 10, listprime[k]);
                            printf(" alla ");
                            mpz_out_str(stdout, 10, randexp2);
                            printf(") è %llu-smooth\n", B);
                            fflush(stdout);
                            candidatefactors = factorsbytrialdivision(candidate);


                            matrix = realloc(matrix, sizeof(mpz_t) * (primelist->count) * (i + 1) * 2);

                            int j;
                            for (j = 0; j < primelist->count * 2; j++) {
                                mpz_init_set_ui(*(matrix + i * (primelist->count * 2) + j), 0);
                            }


                            j = 0;
                            for (temp1 = primelist->HEAD; temp1 != NULL; temp1 = temp1->next) {
                                for (temp2 = candidatefactors->HEAD; temp2 != NULL; temp2 = temp2->next) {
                                    if (mpz_cmp(temp1->value, temp2->value) == 0) {
                                        mpz_set_ui(*(matrix + i * primelist->count * 2 + j),
                                                   (unsigned long) temp2->exp);
                                    }
                                }
                                j++;
                            }
                            j = primelist->count;
                            for (temp1 = primelist->HEAD; temp1 != NULL; temp1 = temp1->next) {
                                if (mpz_cmp(temp1->value, listprime[u]) == 0) {
//                                printf("Aggiungo %d in posizione %d\n", (int) mpz_get_ui(randexp), i *primelist->count * 2 + j);
                                    mpz_set(*(matrix + i * primelist->count * 2 + j), randexp1);
                                }
                                if (mpz_cmp(temp1->value, listprime[k]) == 0) {
                                    mpz_set(*(matrix + i * primelist->count * 2 + j), randexp2);
                                }
                                j++;
                            }

                            for (j = primelist->count; j < primelist->count * 2; j++) {
                                if (j >= primelist->count) {
                                    mpz_sub(*(matrix + i * primelist->count * 2 + (j - primelist->count)),
                                            *(matrix + i * primelist->count * 2 + (j - primelist->count)),
                                            *(matrix + i * primelist->count * 2 + j));
                                    mpz_mod(*(matrix + i * primelist->count * 2 + (j - primelist->count)),
                                            *(matrix + i * primelist->count * 2 + (j - primelist->count)), pm1);
                                    mpz_set_ui(*(matrix + i * primelist->count * 2 + j), 0);
                                }
                            }


                            mpz_invert(invers, *(matrix + i * primelist->count * 2 + i), pm1);
                            if (mpz_cmp_ui(invers, 0) == 0) {
                                continue;
                            }


                            //*(matrix + i * (primelist->count + 1) + primelist->count * 2 - 1) = (int) mpz_get_ui(randexp);
                            freelist(candidatefactors);
                            i++;

                        }
                        if (i == primelist->count) {
                            break;
                        }
                    }
                }
            }
        }
    }
    mpz_clear(candidate1);
    mpz_clear(candidate2);
    mpz_clear(randexp1);
    mpz_clear(randexp2);
    mpz_clear(randexp);
    mpz_clear(candidate);
    mpz_clear(res);
    mpz_clear(base);

    int j;

    for (k = 0; k < primelist->count; k++) {
        printf("[ ");
        for (j = 0; j < primelist->count * 2; j++) {
            if (j == primelist->count) {
                printf("] = [ ");
            }
            mpz_out_str(stdout, 10, *(matrix + k * primelist->count * 2 + j));
            printf(" ");
        }
        printf("]\n");
    }
    printf("\n");



    for (i = 0; i < SIZE; i++) {
        for (j = i + 1; j < SIZE; j++) {
            mpz_invert(temp, *(matrix + i * ROWSIZE + i), pm1);
//            if (mpz_cmp_ui(*(matrix + i * ROWSIZE + i), 0) == 0) {
//                i--;
//                break;
//            }
            printf("L'inverso di ");
            mpz_out_str(stdout, 10, *(matrix + i * ROWSIZE + i));
            printf(" modulo ");
            mpz_out_str(stdout, 10, pm1);
            printf(" è :");
            mpz_out_str(stdout, 10, temp);
            printf("\n");
            mpz_mul(c, *(matrix + j * ROWSIZE + i), temp);
            if (mpz_cmp(*(matrix + i * ROWSIZE + i), *(matrix + j * ROWSIZE + i)) == 0) {
                mpz_set_ui(c, 1);
            }

            printf("c = ");
            mpz_out_str(stdout, 10, c);
            printf("\n");
            for (k = 0; k < SIZE; k++) {
                mpz_mul(temp, c, *(matrix + i * ROWSIZE + k));
                mpz_sub(*(matrix + j * ROWSIZE + k), *(matrix + j * ROWSIZE + k), temp);
                mpz_mod(*(matrix + j * ROWSIZE + k), *(matrix + j * ROWSIZE + k), pm1);
            }
        }
    }


    for (k = 0; k < primelist->count; k++) {
        printf("[ ");
        for (j = 0; j < primelist->count * 2; j++) {
            if (j == primelist->count) {
                printf("] = [ ");
            }
            mpz_out_str(stdout, 10, *(matrix + k * primelist->count * 2 + j));
            printf(" ");
        }
        printf("]\n");
    }
    printf("\n");

//    for (k = 0; k < primelist->count; k++) {
//        printf("[ ");
//        for (j = 0; j < primelist->count * 2; j++) {
//            if (j == primelist->count) {
//                printf("] = [ ");
//            }
//            printf("%d ", *(matrix + k * primelist->count * 2 + j));
//
//        }
//        printf("]\n");
//    }
//    printf("\n");
//
//    float sum;
//    /* this loop is for backward substitution*/
//    for(i = 0; i < SIZE; i++)
//    {
//        sum = *(matrix + i * ROWSIZE + SIZE);
//        for(j = 0; j < i; j++)
//        {
//            sum = sum - x[j] * (float) *(matrix + i * ROWSIZE + j);
//        }
//        x[i] = sum / *(matrix + i * ROWSIZE + i);
//        printf("\n x%d => %f", i+1, x[i]);
//    }


    mpz_clear(pm2);


    return EXIT_SUCCESS;
}




