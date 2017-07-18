#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <zconf.h>
#include <mpfr.h>
#include <pthread.h>
#include "utils.h"
#include "factorize.h"

float denom = 4;

pthread_rwlock_t primelistlock = PTHREAD_RWLOCK_INITIALIZER;
pthread_mutex_t randlock = PTHREAD_MUTEX_INITIALIZER;

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
    struct list *primelist_tob;
};


void *findrelation(void *arg) {

    struct tparams *params = (struct tparams *) arg;
    mpz_t x, y, gcd, u, candidate, invers, c, term;
    int j, k, l, isredundant;
    struct element *temp1, *temp2;

    int RELATIONS = params->ROWSIZE - 1;

    struct list *candidatefactorsofgcd;
    struct list *candidatefactorsofy;


    for (;;) {
        mpz_init(u);
        mpz_init(candidate);

        pthread_mutex_lock(&randlock);
        mpz_urandomm(u, randstate, params->pm);
        pthread_mutex_unlock(&randlock);

        mpz_add_ui(u, u, 1);
        mpz_powm(candidate, params->am, u, params->pm);

        if (mpz_sizeinbase(candidate, 10) < mpz_sizeinbase(params->pm, 10) /*|| useless == 1*/) {

            mpz_init(x);
            mpz_init(y);
            mpz_init(gcd);

            gcdext(params->pm, candidate, x, y, gcd);

            if (mpz_cmp_ui(x, 0) == 0 && mpz_cmp_ui(y, 0) == 0 && mpz_cmp_ui(gcd, 0) == 0) {
                mpz_clear(x);
                mpz_clear(y);
                mpz_clear(gcd);
                mpz_clear(candidate);
                mpz_clear(u);
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

            pthread_rwlock_rdlock(&primelistlock);

            candidatefactorsofgcd = isBsmooth(gcd, params->Bm, params->primelist_tob);
            if (candidatefactorsofgcd == NULL) {
                mpz_clear(x);
                mpz_clear(y);
                mpz_clear(gcd);
                mpz_clear(candidate);
                mpz_clear(u);
                continue;
            }

            candidatefactorsofy = isBsmooth(y, params->Bm, params->primelist_tob);
            if (candidatefactorsofy == NULL) {
                freelist(candidatefactorsofgcd);
                mpz_clear(x);
                mpz_clear(y);
                mpz_clear(gcd);
                mpz_clear(candidate);
                mpz_clear(u);
                continue;
            }

            pthread_rwlock_unlock(&primelistlock);


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
            for (temp1 = params->primelist_tob->HEAD; temp1 != NULL; temp1 = temp1->next) {
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
            for (temp1 = params->primelist_tob->HEAD; temp1 != NULL; temp1 = temp1->next) {
                // for every factor of the candidate (which is B-smooth)
                for (temp2 = candidatefactorsofy->HEAD; temp2 != NULL; temp2 = temp2->next) {
                    // if the primes are equal
                    if (mpz_cmp(temp1->value, temp2->value) == 0) {
                        // set the related coefficient in the matrix equals to the exponent of that factor
                        mpz_sub_ui(*(params->matrix + row_counter * params->ROWSIZE + j),
                                   *(params->matrix + row_counter * params->ROWSIZE + j),
                                   (unsigned long) temp2->exp);
                        mpz_mod(*(params->matrix + row_counter * params->ROWSIZE + j),
                                *(params->matrix + row_counter * params->ROWSIZE + j), params->qm);
                    }
                }
                j++;
            }
            // set the last element of the row equals to the exponent of the primitive root
            mpz_set(*(params->matrix + row_counter * params->ROWSIZE + params->ROWSIZE - 1), u);

//            for (int n = 0; n < params->ROWSIZE; n++) {
//                mpz_out_str(stdout, 10, *(params->matrix + row_counter * params->ROWSIZE + n));
//                printf(" ");
//            }
//            printf("\n");

//            printf("gcd = ");
//            mpz_out_str(stdout, 10, gcd);
//            printf("\ny = ");
//            mpz_out_str(stdout, 10, y);
//            printf("\nFattori di gcd = ");
//            printlist(candidatefactorsofgcd);
//            printf("Fattori di y = ");
//            printlist(candidatefactorsofy);

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
                if (row_counter >= k + 1) {
                    mpz_mul(c, *(params->matrix + row_counter * params->ROWSIZE + k), invers);
                    for (l = 0; l < params->ROWSIZE; l++) {
                        mpz_mul(term, c, *(params->matrix + k * params->ROWSIZE + l));
                        mpz_sub(*(params->matrix + row_counter * params->ROWSIZE + l),
                                *(params->matrix + row_counter * params->ROWSIZE + l), term);
                        mpz_mod(*(params->matrix + row_counter * params->ROWSIZE + l),
                                *(params->matrix + row_counter * params->ROWSIZE + l), params->qm);
                    }
                }
            }

            mpz_clear(invers);
            mpz_clear(c);
            mpz_clear(term);
//                for (k = 0; k < row_counter + 1; k++) {
//                    printf("[ ");
//                    for (j = 0; j < params->ROWSIZE; j++) {
//                        mpz_out_str(stdout, 10, *(params->matrix + k * params->ROWSIZE + j));
//                        printf(" ");
//                    }
//                    printf("]\n");
//                }
//                printf("\n");

            if (isredundant == 1) {
                for (j = 0; j < params->ROWSIZE; j++) {
                    mpz_clear(*(params->matrix + row_counter * params->ROWSIZE + j));
                }
                pthread_mutex_unlock(&counter_mutex);
                mpz_clear(x);
                mpz_clear(y);
                mpz_clear(gcd);
                mpz_clear(candidate);
                mpz_clear(u);
                continue;
            }

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

        mpz_clear(u);
        mpz_clear(candidate);

    }
}


double main(int argc, char *argv[]) {

    mpz_t b, e, primo;
    mpz_init_set_ui(b, 5);
    char *str = "173967506280653770929324";
    mpz_init_set_str(e, str, 10);
    mpz_out_str(stdout, 10, e);
    printf("\n");
    char *ostr = "200000000000000000002967";
    mpz_init_set_str(primo, ostr, 10);
    mpz_out_str(stdout, 10, primo);
    printf("\n");
    mpz_powm(b, b, e, primo);
    mpz_out_str(stdout, 10, b);
    printf("\n");

    mpz_t am, pm, pm1, exp, try, qm, uno, due, Bm, tm, pm2, p2m, exponent, lp;
    mpfr_t logpm, mpfrpm, mpfrexponent, mpfrlp, coeff, loglogpm, product;
    mpfr_rnd_t rnd;

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
    mpz_init_set_ui(uno, 1);
    mpz_init_set_ui(due, 2);
    mpfr_init(mpfrpm);
    mpz_init(exponent);
    mpz_init(lp);
    mpfr_init(mpfrlp);
    mpfr_init(loglogpm);
    mpfr_init(product);
    mpfr_init(mpfrexponent);



    printf("Inserisci il numero da fattorizzare: ");
    mpz_inp_str(qm, STDIN_FILENO, 10);
    printf("q = ");
    mpz_out_str(stdout, 10, qm);
    printf("\n");

    mpz_mul(pm, qm, due);
    mpz_add(pm, pm, uno);
    printf("p = ");
    mpz_out_str(stdout, 10, pm);
    printf("\n");

    // Find B-best
    mpfr_set_default_rounding_mode(rnd);
    mpfr_init_set_d(coeff, denom, rnd);
    mpfr_set_z(mpfrpm, pm, rnd);
    mpfr_log(logpm, mpfrpm, rnd);


    mpfr_log(loglogpm, logpm, rnd);
    mpfr_mul(product, logpm, loglogpm, rnd);
    mpfr_div(product, product, coeff, rnd);


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

    int count = 0;
    int i = 0;
    int j, k;

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
    for (mpz_set_ui(am, 5); mpz_cmp(am, pm) < 0; mpz_add_ui(am, am, 1)) {
        for (divisor = divisors->HEAD; divisor != NULL; divisor = divisor->next) {
//            mpz_div(exp, pm1, divisor->value);
            mpz_powm(try, am, divisor->value, pm);
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
    struct list *primelist_tob = init_list(); //= primes_in_range(uno, Bm);
    eratosthenes(Bm, primelist_tob);

    struct list *candidatefactorsofgcd;
    struct list *candidatefactorsofy;
    struct element *temp1;


    // convert the dinamic list of mpz_t in a common array of mpz_t
    mpz_t listprime[primelist_tob->count];
    returnlist(primelist_tob, listprime);

    mpz_sub_ui(pm2, pm, 2);
    mpz_mul_ui(p2m, pm, 2);

    int RELATIONS = primelist_tob->count;
    int ROWSIZE = primelist_tob->count + 1;

    // Initialize matrix of coefficients
    mpz_t *matrix = malloc(RELATIONS * sizeof(mpz_t) * ROWSIZE);

    mpz_t candidate, res, u, invers, c, term, x, y, gcd;
    mpz_init(res);
    mpz_init(candidate);

    mpz_init(u);


    // Initialize variable for random exponents
    gmp_randinit_mt(randstate);
    gmp_randseed_ui(randstate, (unsigned long) time(NULL));

    printf("Tra 1 e ");
    mpz_out_str(stdout, 10, Bm);
    printf(" ci sono %d numeri primi\n", RELATIONS);

    struct list *candidatefactors = init_list();
    struct element *factor, *prime;

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
            params[f]->primelist_tob = primelist_tob;

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
    mpz_invert(invers, *(matrix + (RELATIONS - 1) * ROWSIZE + (RELATIONS - 1)), qm);
    mpz_mul(sol[RELATIONS - 1], *(matrix + (RELATIONS - 1) * ROWSIZE + RELATIONS), invers);
    mpz_mod(sol[RELATIONS - 1], sol[RELATIONS - 1], qm);

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
        mpz_out_str(stdout, 10, value);
        printf("\n");
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

        if (mpz_sizeinbase(candidate, 10) < mpz_sizeinbase(pm, 10)  /*|| useless == 1*/) {

            mpz_init(x);
            mpz_init(y);
            mpz_init(gcd);

            gcdext(pm, candidate, x, y, gcd);

            if (mpz_cmp_ui(x, 0) == 0 && mpz_cmp_ui(y, 0) == 0 && mpz_cmp_ui(gcd, 0) == 0) {
                mpz_clear(x);
                mpz_clear(y);
                mpz_clear(gcd);
                continue;
            }


            candidatefactorsofgcd = isBsmooth(gcd, Bm, primelist_tob);
            candidatefactorsofy = isBsmooth(y, Bm, primelist_tob);


            // if it is
            if (candidatefactorsofgcd != NULL && candidatefactorsofy != NULL) {

                printf("candidate = ");
                mpz_out_str(stdout, 10, candidate);
                printf("\ngcd = ");
                mpz_out_str(stdout, 10, gcd);
                printf("\ny = ");
                mpz_out_str(stdout, 10, y);
                printf("\n");
                printf("fattori di gcd = ");
                printlist(candidatefactorsofgcd);
                printf("\nfattori di y = ");
                printlist(candidatefactorsofy);
                printf("\n");


                mpz_set_ui(result, 0);
                for (factor = candidatefactorsofgcd->HEAD; factor != NULL; factor = factor->next) {
                    addnewfactor(factor->value, candidatefactors);
                    for (temp1 = candidatefactors->HEAD; temp1 != NULL; temp1 = temp1->next) {
                        if (mpz_cmp(factor->value, temp1->value) == 0) {
                            temp1->exp = factor->exp;
                        }
                    }
                }

                for (factor = candidatefactorsofy->HEAD; factor != NULL; factor = factor->next) {
                    addnewfactor(factor->value, candidatefactors);
                    for (temp1 = candidatefactors->HEAD; temp1 != NULL; temp1 = temp1->next) {
                        if (mpz_cmp(factor->value, temp1->value) == 0) {
                            temp1->exp = temp1->exp - factor->exp - 1;
                        }
                    }
                }

//                printf("I fattori risultanti sono = ");
//                printlist(candidatefactors);
//                printf("\n");


                pos = 0;
                // for every factor of the number obtained
                for (factor = candidatefactors->HEAD; factor != NULL; factor = factor->next) {
                    // for every prime in the factor base
                    for (prime = primelist_tob->HEAD; prime != NULL; prime = prime->next) {
                        // if the numbers are equal stop the cycle
                        if (mpz_cmp(prime->value, factor->value) == 0) {
                            break;
                        }
                        // else increment the position
                        pos++;
                    }
                    // cycle on the row
                    for (row = 0; row < primelist_tob->count; row++) {
                        // if I find a 1 in position [row, pos]
                        if (mpz_cmp_ui(*(matrix + row * ROWSIZE + pos), 1) == 0) {
                            // I add factor->exp times the related exponent to solution
                            if (factor->exp > 0) {
//                                printf("%d volte ", factor->exp);
//                                mpz_out_str(stdout, 10, *(matrix + row * ROWSIZE + ROWSIZE - 1));
//                                printf("\n");
                                for (times = 0; times < factor->exp; times++) {
                                    mpz_add(result, result, *(matrix + row * ROWSIZE + ROWSIZE - 1));

                                }
                            } else {
//                                printf("%d volte ", factor->exp);
//                                mpz_out_str(stdout, 10, *(matrix + row * ROWSIZE + ROWSIZE - 1));
//                                printf("\n");
                                for (times = 0; times < -factor->exp; times++) {
                                    mpz_sub(result, result, *(matrix + row * ROWSIZE + ROWSIZE - 1));
                                }
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
                printf("\nEseguito in %ld secondi\n", time(NULL) - start);
                break;
            }
        }
    }


    return EXIT_SUCCESS;
}




