#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double g(double x, double n) {
	return fmod((pow(x,2) + 2), n);
}

/* Standard C Function: Greatest Common Divisor */
double gcd (double a, double b )
{
  double c;
  while (a != 0) {
     c = a;
     a = fmod(b, a);
     b = c;
  }
  return b;
}

double pollard_ro(double x0, double n)
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

double main() {
	double n;
	scanf("%lf",&n);
	printf("n = %.0lf\n", n);

	//double[20];
	double x = 1;
	double x0 = 2;
	int i = 0;
	while (n != 1) {
		x = pollard_ro(x0, n);
		if (x == -1) {
			i++;
			x0++;
			printf("i = %d   x0 = %.0lf n = %.0lf\n",i,x0,n);
			if (i == 10) {
				printf("%.0lf  end\n", n);
				break;
			}
		}
		else {
			x0 = 2;
			printf("%.0lf   %.0lf\n", n, x);
			n = n / x;
		}
	}
	return EXIT_SUCCESS;
}



