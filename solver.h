//
// Created by rfiba on 19/12/18.
//

#ifndef FEM_SOLVER_H
#define FEM_SOLVER_H

#include "Grid.h"

const double eps = 1e-12;
bool gauss(int n, double ** AB, double * X);
void solve(double timeStep, int numbeOfIteration, Grid *grid);
#endif //FEM_SOLVER_H
