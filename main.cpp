#include <iostream>
#include <cmath>
#include "Node.h"
#include "Element.h"
#include "Grid.h"

using namespace std;

int main() {
    Grid test = Grid(4,4, 0.1, 0.1,25, 25, 300,700, 7800, 1 / sqrt(3));
    cout << sizeof(Grid) << endl;
    cout << sizeof(Element)*16 << endl;

    test.showGridByNodes();
    test.showGridByElements();
    test.prepareLocalMatricesH();
    test.agregateMatrixH();
    test.showMatrixH();
    cout << "_________" << endl;
    test.prepareLocalMatricesC();
    test.agregateMatrixC();
    test.showMatrixC();
    cout << "_________" << endl;
    Node nodes[4];
    nodes[0].setXYT(0,0, 20);
    nodes[1].setXYT(0.025,0, 20);
    nodes[2].setXYT(0.025,0.025, 20);
    nodes[3].setXYT(0,0.025, 20);
    Element testC = Element(0,nodes, 1 / sqrt(3));

    for(int i = 0; i < 4; i++)
        testC.setOutsideFlag(1,i);
    testC.calculateLengths();

    testC.prepareMatrixH(30);
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
    testC.showMatrixH();
    cout << sizeof(Element) << endl;
    return 0;
}