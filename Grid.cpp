//
// Created by rfiba on 17/10/18.
//

#include "Grid.h"

Grid::Grid(int nH, int nL, double HtoAdd, double LtoAdd, double KtoAdd, double initialTemperatureToAdd,
           double alphaToAdd, double specificHeatToAdd, double densityToADD,double ambientTemperatureToAdd, double PCToAdd) {
    hOfGrid = nH;
    lOfGrid = nL;
    K = KtoAdd;
    H = HtoAdd;
    L = LtoAdd;
    PC = PCToAdd;
    initialTemperature = initialTemperatureToAdd;
    ambientTemperature = ambientTemperatureToAdd;
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
    vectorP = new double[lOfGrid*hOfGrid];
    vectorT0 = new double[lOfGrid*hOfGrid];
    for(int i = 0; i < hOfGrid*lOfGrid; i++)
    {
        vectorP[i] = 0;
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

    delete vectorP;
    delete vectorT0;
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
            nodes[i][j].setXYT(((double)L/(lOfGrid-1))*j ,((double)H/(hOfGrid-1))*i , initialTemperature);
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
            //elements[i][j].addBoundaryCondition(alpha);

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

void Grid::prepareLocalVectorsP() {
    for(int i = 0; i < hOfGrid-1; i ++)
    {
        for(int j  = 0; j < lOfGrid-1; j++)
            elements[i][j].prepareVectorP(PC,alpha, ambientTemperature);
    }
    elements[0][0].prepareVectorP(PC,alpha, ambientTemperature);
}

void Grid::agregateVectorP() {
    for(int i = 0; i < hOfGrid-1; i ++) {
        for (int j = 0; j < lOfGrid - 1; j++)
            elements[i][j].agregateVectorP(vectorP, hOfGrid*lOfGrid);
    }
}

void Grid::showVectorP() {
    cout << "{";
    for(int i = 0; i < hOfGrid*lOfGrid; i++)
        cout << vectorP[i] << " ";
    cout << "}\n";
}

void Grid::divideMatrixCbyTimeStep(double timeStep) {
    for(int i = 0; i < hOfGrid*lOfGrid; i++)
    {
        for(int j = 0; j < hOfGrid*lOfGrid; j++)
            matrixC[i][j] /= timeStep;
    }
}

void Grid::sumMatrixHandMatrixCbyTimeStep() {
    for(int i = 0; i < hOfGrid*lOfGrid; i++) {
        for (int j = 0; j < hOfGrid * lOfGrid; j++)
            matrixH[i][j] += matrixC[i][j];
    }
}

void Grid::addBoundaryConditionOnElements() {
    for(int i = 0; i < hOfGrid-1; i ++) {
        for (int j = 0; j < lOfGrid - 1; j++)
            elements[i][j].addBoundaryCondition(alpha);
    }
}

void Grid::sumVectorPandMatrixCbyTimeSteptimesTemperatures() {

    for(int i = 0; i < hOfGrid*lOfGrid; i++)
    {
        for(int j = 0; j < hOfGrid*lOfGrid;j++)
            vectorP[i]+=matrixC[i][j]*vectorT0[j];
    }

}

void Grid::prepareVectorT0() {
    for(int i = 0; i < hOfGrid-1; i ++) {
        for (int j = 0; j < lOfGrid - 1; j++)
            elements[i][j].getTemperaturesInVector(vectorT0);
    }


}

void Grid::minusVectorP() {
    for(int i = 0; i < hOfGrid*lOfGrid; i++)
        vectorP[i] *= -1;
}

int Grid::getHOfGrid() {
    return hOfGrid;
}

int Grid::getLOfGrid() {
    return lOfGrid;
}

double *Grid::getVectorP() {
    return vectorP;
}

double *Grid::getVectorT0() {
    return vectorT0;
}

double **Grid::getMatrixH() {
    return matrixH;
}

void Grid::setTemperatures(double *vector) {
    for(int i = 0; i < hOfGrid-1; i ++) {
        for (int j = 0; j < lOfGrid - 1; j++)
            elements[i][j].setTemperature(vector);
    }
}

void Grid::clearLocalandGlobalVectorMatrices() {
    for(int i = 0; i < hOfGrid-1; i ++) {
        for (int j = 0; j < lOfGrid - 1; j++)
            elements[i][j].clearVectorAndMatrices();
    }

    for(int i = 0; i < hOfGrid*lOfGrid; i++) {
        for (int j = 0; j < hOfGrid * lOfGrid; j++)
        {
            matrixC[i][j] = 0;
            matrixH[i][j] = 0;
        }
        vectorP[i] =0;
    }
}

