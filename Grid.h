//
// Created by rfiba on 17/10/18.
//

#ifndef FEM_GRID_H
#define FEM_GRID_H


#include "Node.h"
#include "Element.h"

class Grid {
private:
    Node** nodes;
    Element** elements;
    unsigned short hOfGrid;
    unsigned short lOfGrid;
    double H;
    double L;
    double K;
    double initialTemperature;
    double PC;
    double alpha;
    double specificHeat;
    double density;
    double** matrixH;
    double** matrixC;
    double ambientTemperature;
    double *vectorP;
    double *vectorT0;
    void prepareNodes();
    void prepareElements();
public:
    Grid(int nH, int nL, double H, double L, double K, double initialTemperature, double alphaToAdd,
         double specificHeatToAdd, double densityToADD, double ambientTemperatureToAdd, double PC);
    ~Grid();
    void showGridByNodes();
    void showGridByElements();
    void agregateMatrixH();
    void showMatrixH();
    void prepareLocalMatricesH();
    void prepareLocalMatricesC();
    void agregateMatrixC();
    void prepareLocalVectorsP();
    void agregateVectorP();
    void showMatrixC();
    void showVectorP();
    void divideMatrixCbyTimeStep(double timeStep);
    void sumMatrixHandMatrixCbyTimeStep();
    void addBoundaryConditionOnElements();
    void sumVectorPandMatrixCbyTimeSteptimesTemperatures();
    void prepareVectorT0();
    void minusVectorP();
    int getHOfGrid();
    int getLOfGrid();
    double* getVectorP();
    double* getVectorT0();
    double** getMatrixH();
};


#endif //FEM_GRID_H
