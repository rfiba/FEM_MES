#include <iostream>
#include <cmath>
#include "Node.h"
#include "Element.h"
#include "Grid.h"
#include "solver.h"

using namespace std;

int main() {
    unsigned short h,l;
    int simulationStep, simulationTime;
    double height, length, initialTemperature, ambientTemperature, alpha, specificHeat, conductivity, density;
    cin >> initialTemperature >> simulationTime >> simulationStep
        >> ambientTemperature >> alpha >> height >> length >> h >> l >> specificHeat >> conductivity >> density;

    Grid test = Grid(h,l, height, length,conductivity, initialTemperature, alpha,specificHeat, density,
                     ambientTemperature, 1 / sqrt(3));

    solve(simulationStep,simulationTime/simulationStep, &test);

    return 0;
}