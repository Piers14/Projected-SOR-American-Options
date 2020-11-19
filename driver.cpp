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

// Example Boundary and Initial Conditions
double TestBC1(double t){return 100;} // boundary condition 1
double TestBC2(double t){return 0;} // boundary condition 2
double TestIC(double x){return std::max(100-x, 0.0);} // intitial condition
double psi(double x){return -TestIC(x);}


int main()
{
    // Example
    int N = 5;
    int m = 10;

    std::vector<std::vector<double>> test;
    test = AmericanOption(TestIC, TestBC1, TestBC2, psi, 0.5, 0.05, 300, 5, m, N);
    std::cout << "Example Approximation: \n";
    OutputMatrix(test);

    return 0;
}




