//
// Created by rfiba on 12/10/18.
//

#include "Element.h"

Element::Element(unsigned short idToAdd, Node *nodesToAdd, double PC) {
    id = idToAdd;

    for(int i = 0; i < 4; i++)
        nodes[i] = &nodesToAdd[i];

    double PCArray[] = {-PC,-PC, PC, -PC, PC, PC, -PC, PC};
    double (*pointerShapeF[])(double, double)={N1,N2,N3,N4};
    double (*pointerKsiDerivateF[])(double)={N1_ksi_derivative,N2_ksi_derivative,N3_ksi_derivative,N4_ksi_derivative};
    double (*pointerEtaDerivateF[])(double)={N1_eta_derivative,N2_eta_derivative,N3_eta_derivative,N4_eta_derivative};

    for(int i = 0, tmp =0; i<4;i++, tmp = 0)
    {
        for(int j = 0; j < 4; j++, tmp+=2)
        {
            shapeFunctionMatrix[j][i] = pointerShapeF[i](PCArray[tmp], PCArray[tmp+1]);
            shapeDKsiFunctionMatrix[i][j]  = pointerKsiDerivateF[i](PCArray[tmp+1]);
            shapeDEtaFunctionMatrix[i][j]  = pointerEtaDerivateF[i](PCArray[tmp]);
        }
    }

    for(int i = 0; i < 4; i++)
        outsideFlags[i] = false;
}

ostream &operator<<(ostream &output, Element &toShow) {
    output << toShow.id << endl << "{";
    for(Node* i: toShow.nodes)
        output << *i << ",\t";
    output << "}";
    return output;
}

void Element::addNodes(unsigned short idToAdd, Node *aToAdd, Node *bToAdd, Node *cToAdd, Node *dToAdd) {
    id = idToAdd;
    nodes[0] = aToAdd;
    nodes[1] = bToAdd;
    nodes[2] = cToAdd;
    nodes[3] = dToAdd;
}

void Element::setID(unsigned short idToAdd) {
    id = idToAdd;
}

void Element::prepareJacobian() {

    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++) {
            jacobiMatrix[i][j] = 0;
            switch(i){
                case 0:
                    for(int k = 0; k < 4; k++)
                        jacobiMatrix[i][j] += shapeDKsiFunctionMatrix[k][i] * nodes[k]->getX();
                    break;
                case 1:
                    for(int k = 0; k < 4; k++)
                        jacobiMatrix[i][j] += shapeDKsiFunctionMatrix[k][i] * nodes[k]->getY();
                    break;
                case 2:
                    for(int k = 0; k < 4; k++)
                        jacobiMatrix[i][j] += shapeDEtaFunctionMatrix[k][i] * nodes[k]->getX();
                    break;
                case 3:
                    for(int k = 0; k < 4; k++)
                        jacobiMatrix[i][j] += shapeDEtaFunctionMatrix[k][i] * nodes[k]->getY();
                    break;
            }
        }
    }
    for(int i = 0; i < 4; i++)
    {
        detJacobian[i] = jacobiMatrix[0][i]*jacobiMatrix[3][i] - jacobiMatrix[1][i]*jacobiMatrix[2][i];
    }
}

void Element::showJacobian() {
    for(int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++)
            cout << jacobiMatrix[i][j] << " ";
        cout << endl;
    }

    for(int i = 0; i < 4; i++)
        cout << detJacobian[i] << " ";
}

void Element::reverseJacobiMatrix() {
    for(int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++)
        {
            switch(j)
            {
                case 0:
                    reversedJacobiMatrix[j][i] = jacobiMatrix[3][i]/detJacobian[i];
                    break;
                case 1:
                    reversedJacobiMatrix[j][i] = -jacobiMatrix[1][i]/detJacobian[i];
                    break;
                case 2:
                    reversedJacobiMatrix[j][i] = -jacobiMatrix[2][i]/detJacobian[i];
                    break;
                case 3:
                    reversedJacobiMatrix[j][i] = jacobiMatrix[0][i]/detJacobian[i];
                    break;
            }
        }
    }
}

void Element::showReversedJacobiMatrix() {
    for(int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++)
            cout << reversedJacobiMatrix[i][j] << " ";
        cout << endl;
    }
}

void Element::prepareDNdMatrix() {
    for(int i = 0; i< 4; i++)
    {
        for(int j = 0; j <4; j++){
            dNdXmatrix[i][j] = reversedJacobiMatrix[0][j]*shapeDKsiFunctionMatrix[i][j]
                               + reversedJacobiMatrix[1][j] * shapeDEtaFunctionMatrix[i][j];
            dNdYmatrix[i][j] = reversedJacobiMatrix[2][j]*shapeDKsiFunctionMatrix[i][j]
                               + reversedJacobiMatrix[3][j] * shapeDEtaFunctionMatrix[i][j];
        }
    }
}

void Element::showdNdMatrix() {
    for(int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++)
            cout << dNdXmatrix[i][j] << " ";
        cout << endl;
    }
    for(int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++)
            cout << dNdYmatrix[i][j] << " ";
        cout << endl;
    }
}

void Element::prepareNdXYPCmatrices() {
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            for(int k = 0; k < 4; k++){
                dNdXPCmatrix[i][k][j] = dNdXmatrix[j][i]*dNdXmatrix[k][i];
                dNdYPCmatrix[i][k][j] = dNdYmatrix[j][i]*dNdYmatrix[k][i];
            }
        }
    }
}

void Element::showNdXPCmatrix() {
    for(int i = 0 ; i < 4; i ++){
        for(int j = 0; j < 4; j++)
        {
            for(int k = 0; k < 4; k++)
                cout << dNdXPCmatrix[i][j][k] << " ";
            cout << endl;
        }
        cout << "---" << endl;
    }
}

void Element::showNdYPCmatrix() {
    for(int i = 0 ; i < 4; i ++){
        for(int j = 0; j < 4; j++)
        {
            for(int k = 0; k < 4; k++)
                cout << dNdYPCmatrix[i][j][k] << " ";
            cout << endl;
        }
        cout << "---" << endl;
    }
}

void Element::prepareMatrixH() {
    prepareJacobian();
    reverseJacobiMatrix();
    prepareDNdMatrix();
    prepareNdXYPCmatrices();
    for(int i = 0; i <4; i++)
    {
        for(int j = 0; j <4; j ++)
            matrixH[i][j] = 0;
    }
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            for(int k = 0; k < 4; k++)
                matrixH[j][k] += (dNdXPCmatrix[i][j][k] + dNdYPCmatrix[i][j][k])*detJacobian[i];
        }
    }
}

void Element::showMatrixH(){
    for(int i = 0; i <4; i++)
    {
        for(int j = 0; j <4; j ++)
            cout << matrixH[i][j]*30 << " ";
        cout << endl;
    }
}

void Element::prepareMatrixC() {
    for(int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                NNPCMatrix[i][k][j] = detJacobian[i]*shapeFunctionMatrix[k][i]*shapeFunctionMatrix[j][i]*700*7800;
            }
        }
    }

    for(int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                matrixC[i][j] += NNPCMatrix[k][i][j];
            }
        }
    }
}

void Element::showMatrixC() {
    for(int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << matrixC[i][j] << " ";
        }
        cout << endl;
    }

}

void Element::setOutsideFlag(bool flagToAdd, int i ) {
    outsideFlags[i] = flagToAdd;
}

void Element::calculateLengths() {
    for(int i = 0; i < 4; i++)
    {
        sideLengths[i] =
                sqrt(pow(nodes[i]->getX()-nodes[(i+1)%4]->getX(),2) + pow(nodes[i]->getY()-nodes[(i+1)%4]->getY(),2));
    }
}



