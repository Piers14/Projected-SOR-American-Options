#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cassert>
#include <fstream>

std::vector<std::vector<double>> AmericanOption(double (*IC)(double), double (*BC1)(double),
                                                double (*BC2)(double), double (*psi)(double),
                                                double sigma, double r, double R,
                                                double T, int timeSteps, int N);
void OutputMatrix(std::vector<std::vector<double>> vect);
std::vector<double> ProjectedSOR(std::vector<std::vector<double>> B,
                                 std::vector<double> f_x,std::vector<double> psi_x,
                                 std::vector<double> u_0, double w, bool Output);
double TwoNormError(std::vector<double> v1, std::vector<double> v2);

std::vector<std::vector<double>> AmericanOption(double (*IC)(double), double (*BC1)(double),
                                                double (*BC2)(double), double (*psi)(double),
                                                double sigma, double r, double R,
                                                double T, int timeSteps, int N)
{
    std::vector<double> x(N+1, 0); // spacial grid points
    std::vector<double> initValues(N+1, 0); // initial values
    double h = R / N; // mesh width
    for(int i = 0; i < N+1; i ++)
    {
        x[i] = i*h;
        initValues[i] = IC(x[i]);
    }
    double delta_t = T / (double)timeSteps; // time step width
    double lambda_1 = (sigma*sigma*delta_t)/(2*h*h);
    double lambda_2 = (r*delta_t)/(h);
    std::vector<std::vector<double>> solution(timeSteps+1, std::vector<double>(N+1, 0.0));
    for(int j = 0; j < N+1; j ++)
    {
        solution[0][j] = initValues[j]; // adding the initial values to solution matrix
    }
    // Setting up f and psi
    std::vector<double> f_x(N-1);
    std::vector<double> psi_x(N-1);
    for(int i = 0; i < N-1; i ++)
    {
        psi_x[i] = psi(x[i+1]);
    }

    std::vector<double> interiorSol(N-1);

    // Setting up matrix D
    std::vector<std::vector<double>> D(N-1, std::vector<double>(N-1));
    for(int i = 1; i < N-1; i ++)
    {
        D[i-1][i] = (-lambda_1 * x[i]*x[i] - lambda_2*x[i]);
        D[i][i-1] = (-lambda_1 * x[i+1]*x[i+1]);
    }
    for(int i = 0; i < N-1; i ++)
    {
        D[i][i] = (1 + r*delta_t + 2*lambda_1*x[i+1]*x[i+1] + lambda_2*x[i+1]);
    }
    std::vector<double> u_0(N-1, 0.0);
    //OutputMatrix(D);
    // Iterative Step
    for(int i = 1; i < timeSteps+1; i ++)
    {
        // First set up f_x = -v^n
        // end points
        solution[i][0] = BC1(delta_t * i);
        solution[i][N] = BC2(delta_t * i);
        // setting up RHS for SOR
        for(int j = 0; j < N-1; j ++)
        {
            f_x[j] = -solution[i-1][j+1];
        }
        f_x[0] -= lambda_1 * x[1] * x[1] * solution[i][0];
        f_x[N-2] -= solution[i][N] * (lambda_1 * x[N-1]*x[N-1] + lambda_2*x[N-1]);

        // Solves for -v^{n+1}
        interiorSol = ProjectedSOR(D, f_x, psi_x, u_0, 1.8, false);
        // Add interior solution into the solution matrix
        for(int j = 1; j < N; j ++)
        {
            solution[i][j] = -interiorSol[j-1];
            u_0[j-1] = interiorSol[j-1];
        }
    }
    return solution;
}

std::vector<double> ProjectedSOR(std::vector<std::vector<double>> B,
                                 std::vector<double> f_x,std::vector<double> psi_x,
                                 std::vector<double> u_0, double w, bool Output)
{
    assert(f_x.size() == B.size());
    int m = f_x.size(); // number of interior points
    std::vector<double> u_old(m);
    std::vector<double> u_new(m);
    std::vector<double> u_bar(m);
    for(int i = 0; i < m; i++)
    {
        u_old[i] = u_0[i];
    }

    double LSum = 0; // the left sum (j<i)
    double RSum = 0; // the right sum (j>i)
    bool repeat = true;
    int ctr = 0; // counter of iterations
    while(repeat)
    {
        for(int i = 0; i < m; i++)
        {
            LSum = RSum = 0;
            for(int j = 0; j < i; j ++)
            {
                LSum += B[i][j] * u_new[j];
            }
            for(int j = (m-1); j > i; j --)
            {
                RSum += B[i][j] * u_old[j];
            }
            u_bar[i] = (f_x[i] - LSum - RSum) / B[i][i];
            u_new[i] = std::min(psi_x[i], (w*u_bar[i] + (1.0-w)*u_old[i]));
        }
        ctr++;
        // Testing for convergence
        if(TwoNormError(u_old, u_new) < 10e-9)
        {
            repeat = false;
        }
        // Testing for max iterations
        if(ctr > 10e5)
        {
            repeat = false;
            std::cout << "Max iterations reached\n";
        }
        for(int i = 0; i < m; i++)
        {
            u_old[i] = u_new[i];
        }
        // Outputs each iteration if output = true (used for 4c)
        if(Output)
        {
            for(int i = 0; i < m; i++)
            {
                std::cout << std::setw(10) << std::setprecision(5)
                      << std::scientific << u_new[i] << " ";
            }
            std::cout << std::endl;
        }
    }
    return u_new;
}

double TwoNormError(std::vector<double> v1, std::vector<double> v2)
{
    assert(v1.size() == v2.size());
    double sum = 0;
    for(int i = 0; i < v1.size(); i++)
    {
        sum += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    }
    return sqrt(sum);
}

void OutputMatrix(std::vector<std::vector<double>> vect)
{
    for (int i = 0; i < vect.size(); i++) {
        for (int j = 0; j < vect[i].size(); j++)
            std::cout << std::setw(10) << std::setprecision(5)
                      << std::scientific << vect[i][j] << " ";
        std::cout << std::endl;
    }
}
