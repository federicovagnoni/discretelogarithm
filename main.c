#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double g(double x, double n) {
	return fmod((pow(x,2) + 1), n);
}

/* Standard C Function: Greatest Common Divisor */
double gcd (double a, double b )
{
  double temp;
  while (b != 0) {
      temp = fmod(a, b);
      a = b;
      b = temp;
  }
  return a;
}

double pollard_rho(double x0, double n)
{
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


double primes_in_range(double n, double m) {
    double i;
    while (n <= m) {
        i = 2;
        while (i <= n) {
            if (fmod(n,i) == 0) {
                break;
            }
            i++;
        }
        if (i == n) {
            //TODO aggiungi alla lista dinamica
        }
        n++;
    }
}

double main(int argc, char* argv[]) {

    double p;
    double q;
    printf("Inserisci il numero da fattorizzare: ");
	scanf("%lf",&q);
	printf("q = %.0lf\n", q);

    p = 1 + 2*q;
    printf("p = %.0lf\n", p);

    double B;
    double lp;
    double exponent = sqrt(log(p)*log(log(p)));
    lp = exp(exponent);
    B = round(pow(lp, 1/sqrt(2) - 1/2));

    printf("B = %lf\n", B);
    primes_in_range(1,B);

	return EXIT_SUCCESS;
}



