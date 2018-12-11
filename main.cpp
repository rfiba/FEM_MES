#include <iostream>
#include <cmath>
#include "Node.h"
#include "Element.h"
#include "Grid.h"

using namespace std;

int main() {
    Grid test = Grid(6,4, 40, 30,30, 20,1 / sqrt(3));
    test.showGridByNodes();
    cout << endl;
    test.showGridByElements();

    Node nodes[4];
    nodes[0].setXYT(0,0, 20);
    nodes[1].setXYT(0.025,0, 20);
    nodes[2].setXYT(0.025,0.025, 20);
    nodes[3].setXYT(0,0.025, 20);
    Element testC = Element(0,nodes, 1 / sqrt(3));
    testC.prepareMatrixH();
    testC.prepareMatrixC();
    cout << "----\n";
    testC.showMatrixC();
    cout << "----\n";
    /*testC.prepareJacobian();
    testC.showJacobian();
    testC.reverseJacobiMatrix();
    cout << endl;
    testC.showReversedJacobiMatrix();
    testC.prepareDNdMatrix();
    testC.showdNdMatrix();
    testC.prepareNdXYPCmatrices();
    testC.showNdYPCmatrix();
    testC.prepareMatrixH();*/
    testC.showMatrixH();
    cout << sizeof(Element) << endl;
    return 0;
}