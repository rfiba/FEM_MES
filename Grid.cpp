//
// Created by rfiba on 17/10/18.
//

#include "Grid.h"

Grid::Grid(int nH, int nL, int HtoAdd, int LtoAdd, double KtoAdd, double enviromentTemperatureToAdd, double PCToAdd) {
    hOfGrid = nH;
    lOfGrid = nL;
    K = KtoAdd;
    H = HtoAdd;
    L = LtoAdd;
    PC = PCToAdd;
    enviromentTempreture = enviromentTemperatureToAdd;
    nodes = new Node*[hOfGrid];
    for(int i = 0; i < hOfGrid; i++)
        nodes[i] = new Node[lOfGrid];
    elements = new Element*[(hOfGrid-1)];
    for(int i = 0; i < hOfGrid-1; i++)
        elements[i] = new Element[lOfGrid-1];
    prepareNodes();
    prepareElements();
}

Grid::~Grid() {
    for(int i = 0; i < hOfGrid-1; ++i) {
        delete[] elements[i];
    }
    delete elements;
    
    for(int i = 0; i < hOfGrid; ++i) {
        delete[] nodes[i];
    }
    delete []nodes;
}



void Grid::prepareElements() {
    int iid = 0;

    for(int i = 0; i < lOfGrid-1; i++ )
    {
        for(int j = 0; j < hOfGrid-1; j++)
        {
            elements[j][i].addNodes(iid, &nodes[j][i], &nodes[j][i+1],&nodes[j+1][i+1] ,&nodes[j+1][i] );
            iid++;
        }
    }
}

void Grid::prepareNodes() {
    for(int i = 0; i < hOfGrid; i++)
    {
        for(int j = 0; j < lOfGrid; j++)
        {
            nodes[i][j].setXYT(((double)L/(lOfGrid-1))*j ,((double)H/(hOfGrid-1))*i , enviromentTempreture);
        }
    }
}

void Grid::showGridByNodes() {

    for(int i = hOfGrid-1; i >= 0; i--)
    {
        for(int j = 0; j < lOfGrid; j++)

            cout <<  nodes[i][j] << "\t\t";
        cout << endl;
    }
}

void Grid::showGridByElements() {
    for(int i = 0; i < lOfGrid-1; i++)
    {
        for(int j = 0; j < hOfGrid-1; j++)

            cout << elements[j][i] << endl;

    }
}

void Grid::prepareLocalHMatrix() {

    for(int i = 0; i < lOfGrid-1; i++ )
    {
        for(int j = 0; j < hOfGrid-1; j++)
        {
            elements[i][j].prepareMatrixH();
        }
    }
}
