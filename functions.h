//
// Created by rfiba on 05/12/18.
//

#ifndef FEM_FUNCTIONS_H
#define FEM_FUNCTIONS_H

#include <cmath>



double N1(double ksi, double eta) ;
double N2(double ksi, double eta) ;
double N3(double ksi, double eta) ;
double N4(double ksi, double eta) ;
double N1_ksi_derivative(double eta) ;
double N2_ksi_derivative(double eta) ;
double N3_ksi_derivative(double eta) ;
double N4_ksi_derivative(double eta) ;
double N1_eta_derivative(double ksi) ;
double N2_eta_derivative(double ksi) ;
double N3_eta_derivative(double ksi) ;
double N4_eta_derivative(double ksi) ;

#endif //FEM_FUNCTIONS_H
