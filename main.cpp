#include <iostream>
#include <cmath>
#include "Node.h"
#include "Element.h"
#include "Grid.h"
#include "solver.h"

using namespace std;

int main() {
    Grid test = Grid(4,4, 0.1, 0.1,25, 100, 300,700, 7800, 1200, 1 / sqrt(3));
    cout << sizeof(Grid) << endl;
    cout << sizeof(Element)*16 << endl;
    solve(50,1, &test);

    test.showGridByNodes();
    test.showGridByElements();
    test.prepareLocalMatricesH();
    test.addBoundaryConditionOnElements();
    test.agregateMatrixH();
    test.showMatrixH();
    cout << "_________" << endl;
    test.prepareLocalMatricesC();
    test.agregateMatrixC();
    test.showMatrixC();
    cout << "_________" << endl;
    test.prepareLocalVectorsP();
    test.agregateVectorP();
    test.showVectorP();
    test.divideMatrixCbyTimeStep(50);

    test.sumMatrixHandMatrixCbyTimeStep();
    test.showMatrixH();
    test.prepareVectorT0();
    test.sumVectorPandMatrixCbyTimeSteptimesTemperatures();
    test.showVectorP();

    Node nodes[4];
    nodes[0].setXYT(0,	0,	25);
    nodes[1].setXYT(0.0333333,	0,	25);
    nodes[2].setXYT(0.0333333,	0.0333333,	25);
    nodes[3].setXYT(0,	0.0333333,	25);
    Element testC = Element(0,nodes, 1 / sqrt(3));

    for(int i = 0; i < 4; i++)
        testC.setOutsideFlag(1,i);
    testC.calculateLengths();
    testC.prepareVectorP(1 / sqrt(3),300,1200);
    /*testC.prepareMatrixH(30);
    testC.addBoundaryCondition(25);
    //testC.prepareMatrixC();
    cout << "----\n";
    //testC.showMatrixC();
    cout << "----\n";
    //testC.prepareJacobian();
    //testC.showJacobian();
    //testC.reverseJacobiMatrix();
    cout << endl;
    //testC.showReversedJacobiMatrix();
    //testC.prepareDNdMatrix();
    //testC.showdNdMatrix();
    //testC.prepareNdXYPCmatrices();
    //testC.showNdYPCmatrix();
    //testC.prepareMatrixH();
    */cout << "-----------" << endl;
    //Grid test1 = Grid(1,1, 0.1, 0.1,25, 25, 300,700, 7800, 1 / sqrt(3));
    //test1.prepareLocalVectorP();
    //cout << sizeof(Element) << endl;
    return 0;
}