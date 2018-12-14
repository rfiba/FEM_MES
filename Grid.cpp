//
// Created by rfiba on 17/10/18.
//

#include "Grid.h"

Grid::Grid(int nH, int nL, double HtoAdd, double LtoAdd, double KtoAdd, double enviromentTemperatureToAdd,
           double alphaToAdd, double specificHeatToAdd, double densityToADD, double PCToAdd) {
    hOfGrid = nH;
    lOfGrid = nL;
    K = KtoAdd;
    H = HtoAdd;
    L = LtoAdd;
    PC = PCToAdd;
    enviromentTempreture = enviromentTemperatureToAdd;
    alpha = alphaToAdd;
    specificHeat =specificHeatToAdd;
    density = densityToADD;

    nodes = new Node*[hOfGrid];
    for(int i = 0; i < hOfGrid; i++)
        nodes[i] = new Node[lOfGrid];

    matrixH = new double*[hOfGrid*lOfGrid];
    for(int i = 0; i < hOfGrid*lOfGrid; i++)
        matrixH[i] = new double[lOfGrid*hOfGrid];

    matrixC = new double*[hOfGrid*lOfGrid];
    for(int i = 0; i < hOfGrid*lOfGrid; i++)
        matrixC[i] = new double[lOfGrid*hOfGrid];

    elements = new Element*[(hOfGrid-1)];
    for(int i = 0; i < hOfGrid-1; i++)
        elements[i] = new Element[lOfGrid-1];

    prepareNodes();
    prepareElements();

    for(int i = 0; i < hOfGrid*lOfGrid; i++)
    {
        for(int j = 0; j < hOfGrid*lOfGrid; j++)
        {
            matrixH[i][j] = 0;
            matrixC[i][j] = 0;
        }
    }
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
            elements[j][i].addNodes(iid, &nodes[j][i], &nodes[j][i+1],&nodes[j+1][i+1] ,&nodes[j+1][i], i, hOfGrid, PC);
            if(i == 0)
                elements[j][i].setOutsideFlag(1, 3);
            if(j == 0)
                elements[j][i].setOutsideFlag(1, 0);
            if( i == lOfGrid-2)
                elements[j][i].setOutsideFlag(1, 1);
            if ( j==hOfGrid-2)
                elements[j][i].setOutsideFlag(1, 2);
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

void Grid::agregateMatrixH() {
    for(int i = 0; i < hOfGrid-1; i ++)
    {
        for(int j  = 0; j < lOfGrid-1; j++)
            elements[i][j].agregateMatrixH(matrixH);
    }
}

void Grid::showMatrixH() {
    for(int i = 0; i < hOfGrid*lOfGrid; i++)
    {
        for(int j = 0; j < hOfGrid*lOfGrid; j++)
            cout << matrixH[i][j] << " ";
        cout << endl;
    }
}

void Grid::prepareLocalMatricesH() {
    for(int i = 0; i < hOfGrid-1; i ++)
    {
        for(int j  = 0; j < lOfGrid-1; j++) {
            elements[i][j].prepareMatrixH(K);
            //elements[i][j].showJacobian();
            //elements[i][j].showReversedJacobiMatrix();
            //elements[i][j].showDetJacobian();
            //elements[i][j].showNdYPCmatrix();
            //elements[i][j].showNdYPCmatrix();
            //cout << endl << endl;
            //elements[i][j].showNdXPCmatrix();
            elements[i][j].addBoundaryCondition(alpha);

            //elements[i][j].showMatrixH();
        }
    }
}

void Grid::prepareLocalMatricesC() {
    for(int i = 0; i < hOfGrid-1; i ++)
    {
        for(int j  = 0; j < lOfGrid-1; j++) {
            elements[i][j].prepareMatrixC(specificHeat, density);
        }}
}

void Grid::agregateMatrixC() {

    for(int i = 0; i < hOfGrid-1; i ++)
    {
        for(int j  = 0; j < lOfGrid-1; j++)
            elements[i][j].agregateMatrixC(matrixC);
    }
}

void Grid::showMatrixC() {
    for(int i = 0; i < hOfGrid*lOfGrid; i++)
    {
        for(int j = 0; j < hOfGrid*lOfGrid; j++)
            cout << matrixC[i][j] << " ";
        cout << endl;
    }
}