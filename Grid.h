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
    double enviromentTempreture;
    double PC;
    double alpha;
    double** matrixH;
    void prepareNodes();
    void prepareElements();
public:
    Grid(int nH, int nL, double H, double L, double K, double enviromentTemperature, double alpha, double PC);
    ~Grid();
    void showGridByNodes();
    void showGridByElements();
    void agregateMatrixH();
    void showMatrixH();
    void prepareLocalMatricesH();

};


#endif //FEM_GRID_H
