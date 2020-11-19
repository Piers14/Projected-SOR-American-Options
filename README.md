# Projected-SOR-American-Options
 
A routine which computes the price of an American Option under Black-Scholes assumptions. Given an intitial condition (the payoff) and boundary conditions, it employs the projected SOR method at each time-step to compute the value of the option at the next time step (across a range of spot prices).

The disretisation can be found in MATH_README.pdf

An example can be found in driver.cpp.
