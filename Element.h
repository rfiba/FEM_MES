//
// Created by rfiba on 12/10/18.
//

#ifndef FEM_ELEMENT_H
#define FEM_ELEMENT_H

#include "Node.h"
#include "functions.h"


class Element {
private:
    int globalNodesID[4];
    unsigned short id;
    Node* nodes[4];
    double jacobiMatrix[4][4];
    double detJacobian[4];
    double reversedJacobiMatrix[4][4];
    double dNdXmatrix[4][4];
    double dNdYmatrix[4][4];
    double dNdXPCmatrix[4][4][4];
    double dNdYPCmatrix[4][4][4];
    double matrixH[4][4];
    double matrixC[4][4];
    double shapeFunctionMatrix[4][4];
    double shapeDKsiFunctionMatrix[4][4];
    double shapeDEtaFunctionMatrix[4][4];
    double NNPCMatrix[4][4][4];
    double PC;

public:
    Element(){};
    Element(unsigned short idToAdd, Node nodesToAdd[4],double PC);
    void addNodes(unsigned short idToAdd, Node* aToAdd, Node* bToAdd, Node* cToAdd, Node* dToAdd);
    friend ostream &operator<<(ostream &output, Element &toShow);
    void setID(unsigned short idToAdd);
    void showJacobian();
    void prepareJacobian();
    void reverseJacobiMatrix();
    void showReversedJacobiMatrix();
    void prepareDNdMatrix();
    void showdNdMatrix();
    void showNdXPCmatrix();
    void prepareNdXYPCmatrices();
    void showNdYPCmatrix();
    void prepareMatrixH();
    void prepareMatrixC();

    void showMatrixH();
    void showMatrixC();
};


#endif //FEM_ELEMENT_H