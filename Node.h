//
// Created by rfiba on 12/10/18.
//

#ifndef FEM_NODE_H
#define FEM_NODE_H

#include <iostream>

using namespace std;

class Node {
private:
    double x;
    double y;
    double t;
public:
    Node(){};
    Node(double xToAdd, double yToAdd, double tToAdd): x(xToAdd), y(yToAdd), t(tToAdd){};

    friend ostream &operator<<(ostream &output, Node &toShow);
    friend istream &operator>>(istream &input, Node &toFill);
    void setXYT(double xToAdd,double yToAdd, double tToADD);
    double getX();
    double getY();
};


#endif //FEM_NODE_H
