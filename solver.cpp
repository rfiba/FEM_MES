//
// Created by rfiba on 19/12/18.
//

#include "solver.h"
void solve(double timeStep, int numberOfIteration, Grid *grid)
{
    int h = grid->getHOfGrid();
    int l = grid->getLOfGrid();
    double **AB;
    double **matrixH;
    double *vectorP;
    double* result = new double [h*l];
    AB = new double*[h*l];
    for(int i = 0; i< h*l; i++)
        AB[i] = new double[h*l+1];
    for(int i = 0; i < numberOfIteration; i++)
    {
        if(1){
            grid->prepareLocalMatricesH();
            grid->addBoundaryConditionOnElements();
            grid->agregateMatrixH();
            grid->prepareLocalMatricesC();
            grid->agregateMatrixC();
            grid->prepareLocalVectorsP();
            grid->agregateVectorP();
        }
        grid->divideMatrixCbyTimeStep(timeStep);
        grid->sumMatrixHandMatrixCbyTimeStep();
        grid->prepareVectorT0();
        grid->sumVectorPandMatrixCbyTimeSteptimesTemperatures();
        //grid->minusVectorP();

        vectorP = grid->getVectorP();

        matrixH = grid->getMatrixH();
        for(int j = 0; j < h*l; j++)
        {
            for(int k = 0; k < h*l; k++) {
                AB[j][k] = matrixH[j][k];

            }
        }


        grid->showVectorP();
        for(int j = 0; j < h*l; j++)
            AB[j][h*l] = vectorP[j];

        for(int j = 0; j < h*l; j++) {
            for (int k = 0; k <= h * l; k++) {
                cout << AB[j][k] << " ";
            }
            cout << endl;
        }

        cout <<  gauss(h*l,AB , result) << endl;
        vectorP = grid->getVectorT0();
        for(int j = 0; j < h*l; j++)
            cout << result[j] << " ";
        cout << endl;
    }
}

bool gauss(int n, double ** AB, double * X)
{

    int i,j,k;
    double m,s;

    for(i = 0; i < n - 1; i++)
    {
        for(j = i + 1; j < n; j++)
        {
            if(fabs(AB[i][i]) < eps) return false;
            m = -AB[j][i] / AB[i][i];
            for(k = i + 1; k <= n; k++)
                AB[j][k] += m * AB[i][k];
        }
    }

    for(i = n - 1; i >= 0; i--)
    {
        s = AB[i][n];
        for(j = n - 1; j >= i + 1; j--)
            s -= AB[i][j] * X[j];
        if(fabs(AB[i][i]) < eps) return false;
        X[i] = s / AB[i][i];
    }
    return true;
}