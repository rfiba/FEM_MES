//
// Created by rfiba on 12/10/18.
//

#include "Node.h"

ostream &operator<<(ostream &output, Node &toShow) {
    return output << toShow.x << "\t" << toShow.y << "\t" << toShow.t;
}

istream &operator>>(istream &input, Node &toFill) {
    return input >> toFill.x >> toFill.y >> toFill.t;
}

void Node::setXYT(double xToAdd, double yToAdd, double tToADD) {
    x = xToAdd;
    y = yToAdd;
    t = tToADD;
}

double Node::getX() {
    return x;
}

double Node::getY() {
    return y;
}

double Node::getT() {
    return t;
}

void Node::setTemperature(double toSet) {
    t = toSet;
}
