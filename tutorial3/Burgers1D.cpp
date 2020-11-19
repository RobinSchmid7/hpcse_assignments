#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include <vector>

struct Simulation
{
    int N;
    double dx;
    double t;
    double * u;
    double * rhs;

    Simulation(int n)
    {
        N  = n;
        t  = 0.0;
        dx = 1.0 / N;

        for (int i = 0 ; i < N ; i++)
        {
            double x = (i + 0.5) * dx;
            u[i] = sin(2*M_PI*x);
        }
    }
    double compute_time_step()
    {
        double max_u = 0.0;
        for (int i = 0 ; i < N ; i++)
        {
            max_u = std::max(max_u, std::fabs(u[i]) );
        }
        
        assert (max_u > 1e-16);

        double dt = 0.5 * dx/max_u;
        return dt;
    }

    double upwind_derivative(double u1, double u2, double u3)
    {
        if (u2 > 0) return (u2-u1)/dx;
        else        return (u3-u2)/dx;
    }

    void advance(double dt)
    {
        for (int i = 0 ; i < N ; i++)
        {
            int ip1 = i + 1;
            int im1 = i - 1;

            double dudx = upwind_derivative(u[im1],u[i],u[ip1]);

            rhs[i] = u[i]*dudx;
        }

        for (int i = 0 ; i < N ; i++)
            u[i] -= dt * rhs[i];

        t += dt;
    }

    void save_to_file()
    {
        std::ofstream myfile;
        myfile.open("example.txt", std::fstream::app);
        myfile << t << " " ;
        for (int i = 0; i < N ; i++)
        {
            myfile << u[i] << " ";
        }
        myfile << "\n";
    }
};

int main()
{
  Simulation test(100);

  double t_max = 1.0;
  double t = 0.0;
  int j = 0;
  while (t < t_max)
  {
    double dt = test.compute_time_step();

    test.advance(dt);

    t += dt;

    if (j%10 == 0) test.save_to_file();
    j ++ ;
  }
  return 0;
}
