#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
//#include <pari/pari.h>

double g(double x, double n) {
	return fmod((pow(x,2) + 1), n);
}

/* Standard C Function: Greatest Common Divisor */
double gcd (double a, double b ) {
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
        }
        else {
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

void freelist(struct list *list) {
    struct element *temp = NULL;
    if (list->count == 1) {
        free(list->HEAD);
        free(list);
        return;
    }
    while (temp != list->HEAD) {
        temp = list->TAIL->prev;
        free(list->TAIL);
        list->TAIL = temp;
    }
    free(list->HEAD);
    free(list);
}

void printlist(struct list *list) {
    struct element *temp = list->HEAD;
    while (temp != NULL) {
        mpz_out_str(stdout,10, temp->value);
        printf("(alla %d)\n", temp->exp);
        temp = temp->next;
    }
    free(temp);
}

//double *returnlist(struct list *linkedlist) {
//    double *list;
//    list = malloc(sizeof(double)*(linkedlist->count));
//    if (list == NULL) {
//        fprintf(stderr, "Error in malloc\n");
//        exit(EXIT_FAILURE);
//    }
//    struct element *temp = linkedlist->HEAD;
//    int i = 0;
//    while (temp != NULL) {
//        list[i] = temp->value;
//        temp = temp->next;
//        i++;
//    }
//    return list;
//}


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
        if (mpz_cmp(i,en) == 0) {
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


void addnewfactor(mpz_t value, struct list *list) {
    struct element *temp;
    for (temp = list->HEAD; temp != NULL; temp = temp->next) {
        if (mpz_cmp(temp->value, value)) {
            temp->exp++;
            return;
        }
    }
    free(temp);
    addelem(value, list);
    return;
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
    while (!mpz_cmp_si(todivide, 1)) {
        mpz_sqrt(root,todivide);
        primelist = primes_in_range((__mpz_struct *) 1, root);
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
//    printf("primi tra ");
//    mpz_out_str(stdout,10, uno);
//    printf(" e ");
//    mpz_out_str(stdout,10, r);
//    printf("\n");
//    printlist(list);
    while (mpz_cmp_si(todivide, 1) != 0) {
        trialdivison(factor, todivide, list);
        //printf("%d\n", mpz_cmp(factor,B));
//        printf("B = ");
//        mpz_out_str(stdout,10,B);
//        printf("\nfactor = ");
//        mpz_out_str(stdout,10, factor);
//        printf("\n");
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

double main(int argc, char* argv[]) {

    unsigned long long p, q, a, B, lp, exponent;

    printf("Inserisci il numero da fattorizzare: ");
	scanf("%llu",&q);
	printf("q = %llu\n", q);

    p = 1 + 2*q;
    printf("p = %llu\n", p);

    exponent = (unsigned long long int) sqrtl(logl(p) * logl(logl(p)));
    printf("exponent = %llu\n", exponent);
    lp = (unsigned long long int) expl(exponent);
    printf("L(p) = %llu\n", lp);
    printf("L(p)^c = %.0Lf\n", sqrtl(lp));
    B = (unsigned long long int) roundl(sqrtl(lp));

    printf("B = %llu\n", B);

    printf("Inserisci una radice primitiva: ");
    scanf("%llu",&a);
    mpz_t Bm, am, pm;
    mpz_init_set_ui(Bm, B);
    mpz_init_set_ui(am, a);
    mpz_init_set_ui(pm, p);

    mpz_t uno;
    mpz_init_set_ui(uno, 1);
    struct list *primelist = primes_in_range(uno, Bm);
    //printlist(primelist);
    struct element *temp;

//    struct list *factorsofp_1 = factorsbytrialdivision(p-1);
//    printlist(factorsofp_1);
//
//    for (temp = factorsofp_1->HEAD; temp != NULL; temp = temp->next) {
//        printf("%.0lf\n",pow(a,(p-1)/temp->value));
//        printf("%.0lf (pari a %.0lf alla %.0lf) è %.0lf-smooth\n", fmod(pow(a,(p-1)/temp->value), p), a, (p-1)/temp->value, B);
//    }

//    mpz_t i;
//    mpz_init_set_ui(i, 1);
    mpz_t pm2;
    mpz_init(pm2);
    mpz_sub_ui(pm2, pm, 2);
    int i;
    mpz_t randexp;
    mpz_init(randexp);
    srand(time(NULL));

    for (i = 0; i < 50; i++) {
        mpz_set_ui(randexp, rand() % (p - 2));
        mpz_t candidate;
        mpz_init(candidate);
        mpz_powm(candidate, am, randexp, pm);
        if (isBsmooth(candidate, Bm) == 1) {
            mpz_out_str(stdout, 10, candidate);
            printf(" (pari a ");
            mpz_out_str(stdout, 10, am);
            printf(" alla ");
            mpz_out_str(stdout, 10, randexp);
            printf(") è %llu-smooth\n", B);
        } else {
            mpz_out_str(stdout, 10, candidate);
            printf(" (pari a ");
            mpz_out_str(stdout, 10, am);
            printf(" alla ");
            mpz_out_str(stdout, 10, randexp);
            printf(") NON è %llu-smooth\n", B);
        }

    }

//    for (temp = primelist->HEAD; temp != NULL; temp = temp->next) {
//        //printf("%.0lf\n", temp->value);
//        double candidate = fmod(pow(temp->value,3), p);
//        if (isBsmooth(candidate, B) == 1) {
//            printf("%.0lf (pari a %.0lf alla %.0lf) è %.0lf-smooth\n", candidate, temp->value, (double) 3, B);
//        }
//    }

    return EXIT_SUCCESS;
}



