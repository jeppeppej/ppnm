#include <stdio.h>
#include <math.h>


#define forAll(i) for (int i = 0; i < n; i++)

void stepper_rk12 (
		int n,
		void f (double t, double y[], double dydt[]),
		double t, double y[], double h,
		double step[], double dy[]
		)
{
	double k0[n];
	f(t, y, k0);

	double y1[n];
	forAll(i) y1[i] = y[i] + (0.5*h)*k0[i];

	double k1[n];
	f(t + 0.5*h, y1, k1);

	forAll(i) step[i] = y[i] + h*k1[i];
	forAll(i) dy[i] = (k1[i] - k0[i])*h;
}

void driver_rk12 (
		int n,
		void f (double t, double y[], double dydt[]),
		double a, double b,
		double y[], double h,
		double acc, double eps
		)
{
		double t = a;

		while (t < b)
		{
			if (t + h > b)
			{
				h = b - t;
			}

			double steps[n], dy[n];
		       	stepper_rk12 (n, f, t, y, h, steps, dy);

			double sum = 0;
			forAll(i) sum += y[i]*y[i];
			double norm_y = sqrt(sum);

			sum = 0;
			forAll(i) sum += dy[i]*dy[i];
			double error = sqrt(sum);

			double tolerance = (acc + eps*norm_y)*sqrt(h/(b - a));

			if (error < tolerance)
			{
				t = t + h;
				forAll(i) y[i] = steps[i];
			}

			if (error > 0)
			{
				h *= 0.95*pow(tolerance/error, 0.25);
			}
			else
			{
				h *= 2;
			}
		}
}

void sir_model (
		double time_between_contacts,
		double t, double y[], double dydt[]
		)
{
	double population_size = 16*1e6;
	double recovery_time = 14;

	dydt[0] = -y[0]*y[1]/population_size/time_between_contacts;
	dydt[1] = y[0]*y[1]/population_size/time_between_contacts - y[1]/recovery_time;
	dydt[2] = y[1]/recovery_time;
}

int main ()
{
	double time_between_contacts[] = {1, 2, 4};
	int n = sizeof (time_between_contacts)/sizeof (time_between_contacts[0]);

	int a = 0;
	int b = 100;
	double y[n];
	double h = 0.1;
	double acc = 1e-2;
	double eps = 1e-2;

	FILE * susceptible = fopen ("s.txt", "w"); 
	FILE * infectious = fopen ("i.txt", "w");
	FILE * removed = fopen ("r.txt", "w"); 

	for (int t = 1; t < b; t++)
	{
		y[0] = 6*1e7;
		y[1] = 10;
		y[2] = 0;

		fprintf (susceptible, "%i ", t);
		fprintf (infectious, "%i ", t);
		fprintf (removed, "%i ", t);

		for (int i = 0; i < n; i++)
		{
			void sir_model_given_time_between_contacts (double t, double y[], double dydt[])
			{
				sir_model (time_between_contacts[i], t, y, dydt);
			}

			driver_rk12 (3, sir_model_given_time_between_contacts, (double) a, (double) t, y, h, acc, eps);

			fprintf (susceptible, "%g ", y[0]);
			fprintf (infectious, "%g ", y[1]);
			fprintf (removed, "%g ", y[2]);
		}

		fprintf (susceptible, "\n");
		fprintf (infectious, "\n");
		fprintf (removed, "\n");
	}

	fclose (susceptible);
	fclose (infectious);
	fclose (removed);

	return 0;
}
