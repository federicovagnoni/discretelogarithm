#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
    double value;
    struct element *next;
    struct element *prev;
};

struct list {
    struct element *HEAD;
    struct element *TAIL;
    int count;
};

struct element *init_elem(double value) {
    struct element *elem;
    elem = malloc(sizeof(struct element));

    if (elem == NULL) {
        fprintf(stderr, "Error in malloc\n");
        exit(EXIT_FAILURE);
    }

    elem->value = value;
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


void addelem(double value, struct list *list) {
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
    //TODO da finire
    temp = list->TAIL->prev;

    free(list->TAIL);
}

void printlist(struct list *list) {
    struct element *temp = list->HEAD;
    while (temp != NULL) {
        printf("%.0lf\n", temp->value);
        temp = temp->next;
    }
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


struct list *primes_in_range(double n, double m) {
    struct list *newlist = init_list();
    double i;
    while (n <= m) {
        i = 2;
        while (i < n) {
            if (fmod(n,i) == 0) {
                break;
            }
            i++;
        }
        if (i == n) {
            addelem(n, newlist);
        }
        n++;
    }
    return newlist;
}

double trialdivison(double n, struct list *list) {
    struct element *temp = list->HEAD;
    while (temp != NULL) {
        if (fmod(n, temp->value) == 0) {
            return temp->value;
        }
        temp = temp->next;
    }
}


int isBsmooth(double m, double B) {
    double todivide = m;
    double factor;
    struct list *list = primes_in_range(1, round(sqrt(m)));
    while (todivide != 1) {
        factor = trialdivison(todivide, list);
        if (factor == -1) {
            factor = todivide;
        }
        if (factor > B) {
            printf("Il numero %.0lf ha %.0lf come fattore: non è %.0lf-smooth\n", m, factor, B);
            return 0;
        }
        //printf("%.0lf è un fattore\n", factor);
        todivide = todivide / factor;
        //printf("il numero diventerà %.0lf\n", todivide);
    }
    return 1;
}

double main(int argc, char* argv[]) {

    double p, q, a, B, lp, exponent;

    printf("Inserisci il numero da fattorizzare: ");
	scanf("%lf",&q);
	printf("q = %.0lf\n", q);

    p = 1 + 2*q;
    printf("p = %.0lf\n", p);

    exponent = sqrt(log(p)*log(log(p)));
    lp = exp(exponent);
    B = round(pow(lp, 1/sqrt(2) - 1/2));

    printf("B = %.0lf\n", B);

    printf("Inserisci una radice primitiva: ");
    scanf("%lf",&a);

    struct list *primelist = primes_in_range(1,B);


    return EXIT_SUCCESS;
}



