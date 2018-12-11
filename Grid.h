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
    unsigned short H;
    unsigned short L;
    double K;
    double enviromentTempreture;
    double PC;

    void prepareNodes();
    void prepareElements();
public:
    Grid(int nH, int nL, int H, int L, double K, double enviromentTemperature, double PC);
    ~Grid();
    void showGridByNodes();
    void showGridByElements();
    void prepareLocalHMatrix();

};


#endif //FEM_GRID_H
