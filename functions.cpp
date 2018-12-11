//
// Created by rfiba on 07/11/18.
//

#include "functions.h"




double N1(double ksi, double eta) {
    return 0.25*(1-ksi)*(1-eta);
}

double N2(double ksi, double eta) {
    return 0.25*(1+ksi)*(1-eta);
}

double N3(double ksi, double eta) {
    return 0.25*(1+ksi)*(1+eta);
}

double N4(double ksi, double eta) {
    return 0.25*(1-ksi)*(1+eta);
}

double N1_ksi_derivative(double eta) {
    return -0.25*(1-eta);
}

double N2_ksi_derivative(double eta) {
    return 0.25*(1-eta);
}

double N3_ksi_derivative(double eta) {
    return 0.25*(1+eta);
}

double N4_ksi_derivative(double eta) {
    return -0.25*(1+eta);
}

double N1_eta_derivative(double ksi) {
    return -0.25*(1-ksi);
}

double N2_eta_derivative(double ksi) {
    return -0.25*(1+ksi);
}

double N3_eta_derivative(double ksi) {
    return 0.25*(1+ksi);
}

double N4_eta_derivative(double ksi) {
    return 0.25*(1-ksi);
}


