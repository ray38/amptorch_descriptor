#include <math.h>
#include <stdio.h>
#include "atomistic_mcsh.h"



double calc_C1(double A, double B, double alpha, double beta){
    double temp = sqrt(M_PI / (alpha + beta));
    return A * B * temp * temp * temp;
}

double calc_C2(double alpha, double beta){
    return -1.0 * (alpha * beta / (alpha + beta));
}


double calc_lambda(double alpha, double beta){
    return  alpha / (alpha + beta);
}

double calc_gamma(double alpha, double beta){
    return alpha + beta;
}

double dx0dx(){
    return 1;
}

double dy0dy(){
    return 1;
}

double dz0dz(){
    return 1;
}

void calc_MCSH_0_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{   
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double m_0_1 = C1 * exp( C2 * r0_sqr);

    deriv[0] = dx0dx() * m_0_1 * (2.0 * C2 * x0);
    deriv[1] = dy0dy() * m_0_1 * (2.0 * C2 * y0);
    deriv[2] = dz0dz() * m_0_1 * (2.0 * C2 * z0);

    value[0] = m_0_1;
}

void calc_MCSH_1_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{   
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double lambda = calc_lambda(alpha, beta);

    double temp = C1 * lambda * exp( C2 * r0_sqr);

    double miu_1_1_1 = temp * x0;
    double miu_1_1_2 = temp * y0;
    double miu_1_1_3 = temp * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = dx0dx() * temp * (1.0 + 2.0 * C2 * x0 * x0);
    deriv[1] = dy0dy() * temp * (2 * C2 * x0 * y0);
    deriv[2] = dz0dz() * temp * (2 * C2 * x0 * z0);

    // dmiu2 dx/dy/dz
    deriv[3] = dx0dx() * temp * (2 * C2 * x0 * y0);
    deriv[4] = dy0dy() * temp * (1.0 + 2.0 * C2 * y0 * y0);
    deriv[5] = dz0dz() * temp * (2 * C2 * y0 * z0);

    // dmiu3 dx/dy/dz
    deriv[6] = dx0dx() * temp * (2 * C2 * x0 * z0);
    deriv[7] = dy0dy() * temp * (2 * C2 * y0 * z0);
    deriv[8] = dz0dz() * temp * (1.0 + 2.0 * C2 * z0 * z0);

    value[0] = miu_1_1_1;
    value[1] = miu_1_1_2;
    value[2] = miu_1_1_3;
}

void calc_MCSH_2_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{   
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    
    double lambda = calc_lambda(alpha, beta);

    double gamma = calc_gamma(alpha, beta);
    double C4 = (3/(2*gamma)) - 1;

    double temp = C1 * exp( C2 * r0_sqr);

    double miu_2_1_1 = temp * (C4 + 3*lambda*lambda*x0*x0);
    double miu_2_1_2 = temp * (C4 + 3*lambda*lambda*y0*y0);
    double miu_2_1_3 = temp * (C4 + 3*lambda*lambda*z0*z0);

    // dmiu1 dx/dy/dz
    deriv[0] = dx0dx() * temp * x0 * (6*lambda*lambda + 2*C2*C4 + 6*C2*lambda*lambda*x0*x0); 
    deriv[1] = dy0dy() * miu_2_1_1 * (2*C2*y0);
    deriv[2] = dz0dz() * miu_2_1_1 * (2*C2*z0);

    // dmiu2 dx/dy/dz
    deriv[3] = dx0dx() * miu_2_1_2 * (2*C2*x0);
    deriv[4] = dy0dy() * temp * x0 * (6*lambda*lambda + 2*C2*C4 + 6*C2*lambda*lambda*y0*y0); 
    deriv[5] = dz0dz() * miu_2_1_2 * (2*C2*z0);

    // dmiu3 dx/dy/dz
    deriv[6] = dx0dx() * miu_2_1_3 * (2*C2*x0);
    deriv[7] = dy0dy() * miu_2_1_3 * (2*C2*y0);
    deriv[8] = dz0dz() * temp * x0 * (6*lambda*lambda + 2*C2*C4 + 6*C2*lambda*lambda*z0*z0); 

    value[0] = miu_2_1_1;
    value[1] = miu_2_1_2;
    value[2] = miu_2_1_3;
}

void calc_MCSH_2_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{   
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    
    double lambda = calc_lambda(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr) * lambda * lambda * 3;

    double miu_2_2_1 = temp * x0 * y0;
    double miu_2_2_2 = temp * x0 * z0;
    double miu_2_2_3 = temp * y0 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = dx0dx() * temp * y0 * (1 + 2*C2*x0*x0); 
    deriv[1] = dy0dy() * temp * x0 * (1 + 2*C2*y0*y0); 
    deriv[2] = dz0dz() * miu_2_2_1 * (2*C2*z0);

    // dmiu2 dx/dy/dz
    deriv[3] = dx0dx() * temp * z0 * (1 + 2*C2*x0*x0);
    deriv[4] = dy0dy() * miu_2_2_2 * (2*C2*y0); 
    deriv[5] = dz0dz() * temp * x0 * (1 + 2*C2*z0*z0);

    // dmiu3 dx/dy/dz
    deriv[6] = dx0dx() * miu_2_2_3 * (2*C2*x0);
    deriv[7] = dy0dy() * temp * z0 * (1 + 2*C2*y0*y0);
    deriv[8] = dz0dz() * temp * y0 * (1 + 2*C2*z0*z0);

    value[0] = miu_2_2_1;
    value[1] = miu_2_2_2;
    value[2] = miu_2_2_3;
}

void calc_MCSH_3_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{   
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    
    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
    double x0_sqr = x0*x0;
    double y0_sqr = y0*y0;
    double z0_sqr = z0*z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (45/(2*gamma)) - 9;

    double temp = C1 * exp( C2 * r0_sqr) * lambda;

    double miu_3_1_1 = temp * x0 * (C3 + 15*lambda_sqr*x0_sqr);
    double miu_3_1_2 = temp * y0 * (C3 + 15*lambda_sqr*y0_sqr);
    double miu_3_1_3 = temp * z0 * (C3 + 15*lambda_sqr*z0_sqr);

    // dmiu1 dx/dy/dz
    deriv[0] = dx0dx() * temp * ((C3 + 45*lambda_sqr*x0_sqr) + 2*C2*x0_sqr*(C3+15*lambda_sqr*x0_sqr)); 
    deriv[1] = dy0dy() * miu_3_1_1 * (2*C2*y0); 
    deriv[2] = dz0dz() * miu_3_1_1 * (2*C2*z0);

    // dmiu2 dx/dy/dz
    deriv[3] = dx0dx() * miu_3_1_2 * (2*C2*x0);
    deriv[4] = dy0dy() * temp * ((C3 + 45*lambda_sqr*y0_sqr) + 2*C2*y0_sqr*(C3+15*lambda_sqr*y0_sqr)); 
    deriv[5] = dz0dz() * miu_3_1_2 * (2*C2*z0);

    // dmiu3 dx/dy/dz
    deriv[6] = dx0dx() * miu_3_1_3 * (2*C2*x0);
    deriv[7] = dy0dy() * miu_3_1_3 * (2*C2*y0);
    deriv[8] = dz0dz() * temp * ((C3 + 45*lambda_sqr*z0_sqr) + 2*C2*z0_sqr*(C3+15*lambda_sqr*z0_sqr));

    value[0] = miu_3_1_1;
    value[1] = miu_3_1_2;
    value[2] = miu_3_1_3;
}

void calc_MCSH_3_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{   
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    
    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
    double x0_sqr = x0*x0;
    double y0_sqr = y0*y0;
    double z0_sqr = z0*z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (15/(2*gamma)) - 3;

    double temp = C1 * exp( C2 * r0_sqr) * lambda;

    double miu_3_2_1 = temp * y0 * (C3 + 15*lambda_sqr*x0_sqr);
    double miu_3_2_2 = temp * x0 * (C3 + 15*lambda_sqr*y0_sqr);
    double miu_3_2_3 = temp * z0 * (C3 + 15*lambda_sqr*x0_sqr);
    double miu_3_2_4 = temp * x0 * (C3 + 15*lambda_sqr*z0_sqr);
    double miu_3_2_5 = temp * z0 * (C3 + 15*lambda_sqr*y0_sqr);
    double miu_3_2_6 = temp * y0 * (C3 + 15*lambda_sqr*z0_sqr);

    // dmiu1 dx/dy/dz
    deriv[0] = dx0dx() * temp * x0 * y0 * (30*lambda_sqr + 2*C2*(C3+15*lambda_sqr*x0_sqr));  
    deriv[1] = dy0dy() * (C3+15*lambda_sqr*x0_sqr) * (1+2*C2*y0_sqr);
    deriv[2] = dz0dz() * miu_3_2_1 * (2*C2*z0);

    // dmiu2 dx/dy/dz
    deriv[3] = dx0dx() * (C3+15*lambda_sqr*y0_sqr) * (1+2*C2*x0_sqr);
    deriv[4] = dy0dy() * temp * x0 * y0 * (30*lambda_sqr + 2*C2*(C3+15*lambda_sqr*y0_sqr));  
    deriv[5] = dz0dz() * miu_3_2_2 * (2*C2*z0);

    // dmiu3 dx/dy/dz
    deriv[6] = dx0dx() * temp * x0 * z0 * (30*lambda_sqr + 2*C2*(C3+15*lambda_sqr*x0_sqr));  
    deriv[7] = dy0dy() * miu_3_2_3 * (2*C2*y0);
    deriv[8] = dz0dz() * (C3+15*lambda_sqr*x0_sqr) * (1+2*C2*z0_sqr);

    // dmiu4 dx/dy/dz
    deriv[9] = dx0dx() * (C3+15*lambda_sqr*z0_sqr) * (1+2*C2*x0_sqr);
    deriv[10] = dy0dy() * miu_3_2_4 * (2*C2*y0);
    deriv[11] = dz0dz() * temp * x0 * z0 * (30*lambda_sqr + 2*C2*(C3+15*lambda_sqr*z0_sqr));  

    // dmiu5 dx/dy/dz
    deriv[12] = dx0dx() * miu_3_2_5 * (2*C2*x0);
    deriv[13] = dy0dy() * temp * y0 * z0 * (30*lambda_sqr + 2*C2*(C3+15*lambda_sqr*y0_sqr));  
    deriv[14] = dz0dz() * (C3+15*lambda_sqr*y0_sqr) * (1+2*C2*z0_sqr);

    // dmiu6 dx/dy/dz
    deriv[15] = dx0dx() * miu_3_2_6 * (2*C2*x0);
    deriv[16] = dy0dy() * (C3+15*lambda_sqr*z0_sqr) * (1+2*C2*y0_sqr);
    deriv[17] = dz0dz() * temp * y0 * z0 * (30*lambda_sqr + 2*C2*(C3+15*lambda_sqr*z0_sqr)); 

    value[0] = miu_3_2_1;
    value[1] = miu_3_2_2;
    value[2] = miu_3_2_3;
    value[3] = miu_3_2_4;
    value[4] = miu_3_2_5;
    value[5] = miu_3_2_6;
}

void calc_MCSH_3_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{   
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double lambda = calc_lambda(alpha, beta);

    double temp =  C1 * exp( C2 * r0_sqr) * lambda * lambda * lambda;
    double m_3_3 = temp * x0 * y0 * z0;

    deriv[0] = dx0dx() * temp * y0 * z0 * (1 + 2*C2*x0*x0);
    deriv[1] = dy0dy() * temp * x0 * z0 * (1 + 2*C2*y0*y0);
    deriv[2] = dz0dz() * temp * x0 * y0 * (1 + 2*C2*z0*z0);

    value[0] = m_3_3;
}

void calc_MCSH_4_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{   
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    
    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
    double x0_sqr = x0*x0;
    double y0_sqr = y0*y0;
    double z0_sqr = z0*z0;

    double lambda_4 = lambda_sqr * lambda_sqr;
    double x0_4 = x0_sqr*x0_sqr;
    double y0_4 = y0_sqr*y0_sqr;
    double z0_4 = z0_sqr*z0_sqr;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (315/gamma) - 90;
    double C4 = (315/(4*gamma*gamma)) - (45/gamma) + 9;

    double temp = C1 * exp( C2 * r0_sqr);

    double miu_4_1_1 = temp * (105*lambda_4*x0_4 + C3*lambda_sqr*x0_sqr + C4);
    double miu_4_1_2 = temp * (105*lambda_4*y0_4 + C3*lambda_sqr*y0_sqr + C4);
    double miu_4_1_3 = temp * (105*lambda_4*z0_4 + C3*lambda_sqr*z0_sqr + C4);

    // dmiu1 dx/dy/dz
    deriv[0] = dx0dx() * temp * (210*lambda_4*x0_4*x0 + (420*lambda_4 + 2*C2*C3*lambda_sqr)*x0_sqr*x0 + (2*C3*lambda_sqr + C4)*x0); 
    deriv[1] = dy0dy() * miu_4_1_1 * (2*C2*y0); 
    deriv[2] = dz0dz() * miu_4_1_1 * (2*C2*z0);

    // dmiu2 dx/dy/dz
    deriv[3] = dx0dx() * miu_4_1_2 * (2*C2*x0);
    deriv[4] = dy0dy() * temp * (210*lambda_4*y0_4*y0 + (420*lambda_4 + 2*C2*C3*lambda_sqr)*y0_sqr*y0 + (2*C3*lambda_sqr + C4)*y0); 
    deriv[5] = dz0dz() * miu_4_1_2 * (2*C2*z0);

    // dmiu3 dx/dy/dz
    deriv[6] = dx0dx() * miu_4_1_3 * (2*C2*x0);
    deriv[7] = dy0dy() * miu_4_1_3 * (2*C2*y0);
    deriv[8] = dz0dz() * temp * (210*lambda_4*z0_4*z0 + (420*lambda_4 + 2*C2*C3*lambda_sqr)*z0_sqr*z0 + (2*C3*lambda_sqr + C4)*z0);

    value[0] = miu_4_1_1;
    value[1] = miu_4_1_2;
    value[2] = miu_4_1_3;
}

void calc_MCSH_4_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{   
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    
    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
    double lambda_3 = lambda_sqr * lambda;
    double x0_sqr = x0*x0;
    double y0_sqr = y0*y0;
    double z0_sqr = z0*z0;
    double x0_3 = x0_sqr * x0;
    double y0_3 = y0_sqr * y0;
    double z0_3 = z0_sqr * z0;

    double x0_lambda_3 = x0_3 * lambda_3;
    double y0_lambda_3 = y0_3 * lambda_3;
    double z0_lambda_3 = z0_3 * lambda_3;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (315/(2*gamma)) - 45;

    double temp = C1 * exp( C2 * r0_sqr) * lambda;

    double tempx = C3 * lambda * x0 + 105 * x0_lambda_3;
    double tempy = C3 * lambda * y0 + 105 * y0_lambda_3;
    double tempz = C3 * lambda * z0 + 105 * z0_lambda_3;

    double miu_4_2_1 = temp * y0 * tempx;
    double miu_4_2_2 = temp * x0 * tempy;
    double miu_4_2_3 = temp * z0 * tempx;
    double miu_4_2_4 = temp * x0 * tempz;
    double miu_4_2_5 = temp * z0 * tempy;
    double miu_4_2_6 = temp * y0 * tempz;

    // dmiu1 dx/dy/dz
    deriv[0] = dx0dx() * temp * y0 * ((C3*lambda + 315*lambda_3*x0_sqr) + 2*C2*x0*tempx);  
    deriv[1] = dy0dy() * temp * tempx * (1 + 2*C2*y0_sqr);
    deriv[2] = dz0dz() * miu_4_2_1 * (2*C2*z0);

    // dmiu2 dx/dy/dz
    deriv[3] = dx0dx() * temp * tempy * (1 + 2*C2*x0_sqr);
    deriv[4] = dy0dy() * temp * x0 * ((C3*lambda + 315*lambda_3*y0_sqr) + 2*C2*y0*tempy);  
    deriv[5] = dz0dz() * miu_4_2_2 * (2*C2*z0);

    // dmiu3 dx/dy/dz
    deriv[6] = dx0dx() * temp * z0 * ((C3*lambda + 315*lambda_3*x0_sqr) + 2*C2*x0*tempx);
    deriv[7] = dy0dy() * miu_4_2_3 * (2*C2*y0);
    deriv[8] = dz0dz() * temp * tempx * (1 + 2*C2*z0_sqr);

    // dmiu4 dx/dy/dz
    deriv[9] = dx0dx() * temp * tempz * (1 + 2*C2*x0_sqr);
    deriv[10] = dy0dy() * miu_4_2_4 * (2*C2*y0);
    deriv[11] = dz0dz() * temp * x0 * ((C3*lambda + 315*lambda_3*z0_sqr) + 2*C2*z0*tempz); 

    // dmiu5 dx/dy/dz
    deriv[12] = dx0dx() * miu_4_2_5 * (2*C2*x0);
    deriv[13] = dy0dy() * temp * z0 * ((C3*lambda + 315*lambda_3*y0_sqr) + 2*C2*y0*tempy);  
    deriv[14] = dz0dz() * temp * tempy * (1 + 2*C2*z0_sqr);

    // dmiu6 dx/dy/dz
    deriv[15] = dx0dx() * miu_4_2_6 * (2*C2*x0);
    deriv[16] = dy0dy() * temp * tempz * (1 + 2*C2*y0_sqr);
    deriv[17] = dz0dz() * temp * y0 * ((C3*lambda + 315*lambda_3*z0_sqr) + 2*C2*z0*tempz);

    value[0] = miu_4_2_1;
    value[1] = miu_4_2_2;
    value[2] = miu_4_2_3;
    value[3] = miu_4_2_4;
    value[4] = miu_4_2_5;
    value[5] = miu_4_2_6;
}


void calc_MCSH_4_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{   
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    
    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
    double x0_sqr = x0*x0;
    double y0_sqr = y0*y0;
    double z0_sqr = z0*z0;

    double lambda_x0_sqr = lambda_sqr * x0_sqr;
    double lambda_y0_sqr = lambda_sqr * y0_sqr;
    double lambda_z0_sqr = lambda_sqr * z0_sqr;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (105/(2*gamma)) - 15;
    double C4 = (105/(4*gamma*gamma)) - (15/gamma) + 3;

    double temp = C1 * exp( C2 * r0_sqr);
    double temp1 = 105 * lambda_x0_sqr * lambda_y0_sqr + C3 * lambda_x0_sqr + C3 * lambda_y0_sqr + C4;
    double temp2 = 105 * lambda_x0_sqr * lambda_z0_sqr + C3 * lambda_x0_sqr + C3 * lambda_z0_sqr + C4;
    double temp3 = 105 * lambda_y0_sqr * lambda_z0_sqr + C3 * lambda_y0_sqr + C3 * lambda_z0_sqr + C4;

    double miu_4_3_1 = temp * temp1;
    double miu_4_3_2 = temp * temp2;
    double miu_4_3_3 = temp * temp3;

    // dmiu1 dx/dy/dz
    deriv[0] = dx0dx() * temp * (2*lambda_sqr*x0*(105*lambda_y0_sqr+C3) + 2*C2*x0*temp1); 
    deriv[1] = dy0dy() * temp * (2*lambda_sqr*y0*(105*lambda_x0_sqr+C3) + 2*C2*y0*temp1); 
    deriv[2] = dz0dz() * miu_4_3_1 * (2*C2*z0);

    // dmiu2 dx/dy/dz
    deriv[3] = dx0dx() * temp * (2*lambda_sqr*x0*(105*lambda_z0_sqr+C3) + 2*C2*x0*temp2); 
    deriv[4] = dy0dy() * miu_4_3_2 * (2*C2*y0);
    deriv[5] = dz0dz() * temp * (2*lambda_sqr*z0*(105*lambda_x0_sqr+C3) + 2*C2*z0*temp2); 

    // dmiu3 dx/dy/dz
    deriv[6] = dx0dx() * miu_4_3_3 * (2*C2*x0);
    deriv[7] = dy0dy() * temp * (2*lambda_sqr*y0*(105*lambda_z0_sqr+C3) + 2*C2*y0*temp3); 
    deriv[8] = dz0dz() * temp * (2*lambda_sqr*z0*(105*lambda_y0_sqr+C3) + 2*C2*z0*temp3);

    value[0] = miu_4_3_1;
    value[1] = miu_4_3_2;
    value[2] = miu_4_3_3;
}

void calc_MCSH_4_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{   
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    
    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
    double x0_sqr = x0*x0;
    double y0_sqr = y0*y0;
    double z0_sqr = z0*z0;

    double lambda_x0_sqr = lambda_sqr * x0_sqr;
    double lambda_y0_sqr = lambda_sqr * y0_sqr;
    double lambda_z0_sqr = lambda_sqr * z0_sqr;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (105/(2*gamma)) - 15;

    double temp = C1 * exp( C2 * r0_sqr) * lambda_sqr;

    double miu_4_4_1 = temp * y0 * z0 * (105*lambda_x0_sqr + C3);
    double miu_4_4_2 = temp * x0 * z0 * (105*lambda_y0_sqr + C3);
    double miu_4_4_3 = temp * x0 * y0 * (105*lambda_z0_sqr + C3);

    // dmiu1 dx/dy/dz
    deriv[0] = dx0dx() * temp * x0 * y0 * z0 * (210*C2*lambda_x0_sqr + 210*lambda_sqr + 2*C2*C3); 
    deriv[1] = dy0dy() * temp * z0 * (105*lambda_x0_sqr + C3) * (1 + 2*C2*y0_sqr);
    deriv[2] = dz0dz() * temp * y0 * (105*lambda_x0_sqr + C3) * (1 + 2*C2*z0_sqr);

    // dmiu2 dx/dy/dz
    deriv[3] = dx0dx() * temp * z0 * (105*lambda_y0_sqr + C3) * (1 + 2*C2*x0_sqr);
    deriv[4] = dy0dy() * temp * x0 * y0 * z0 * (210*C2*lambda_y0_sqr + 210*lambda_sqr + 2*C2*C3); 
    deriv[5] = dz0dz() * temp * x0 * (105*lambda_y0_sqr + C3) * (1 + 2*C2*z0_sqr);

    // dmiu3 dx/dy/dz
    deriv[6] = dx0dx() * temp * y0 * (105*lambda_z0_sqr + C3) * (1 + 2*C2*x0_sqr);
    deriv[7] = dy0dy() * temp * x0 * (105*lambda_z0_sqr + C3) * (1 + 2*C2*y0_sqr);
    deriv[8] = dz0dz() * temp * x0 * y0 * z0 * (210*C2*lambda_z0_sqr + 210*lambda_sqr + 2*C2*C3);

    value[0] = miu_4_4_1;
    value[1] = miu_4_4_2;
    value[2] = miu_4_4_3;
}

int get_mcsh_type(int mcsh_order, int group_num)
{
    if (mcsh_order == 0) {
        if (group_num == 1) {
            return 1;
        } else {
            return 0;
        }
    } else if (mcsh_order == 1) {
        if (group_num == 1) {
            return 2;
        } else {
            return 0;
        }
    } else if (mcsh_order == 2) {
        if (group_num == 1) {
            return 2;
        } else if (group_num == 2){
            return 2;
        } else {
            return 0;
        }
    } else if (mcsh_order == 3) {
        if (group_num == 1) {
            return 2;
        } else if (group_num == 2){
            return 3;
        } else if (group_num == 3){
            return 1;
        } else {
            return 0;
        }
    } else if (mcsh_order == 4) {
        if (group_num == 1) {
            return 2;
        } else if (group_num == 2){
            return 3;
        } else if (group_num == 3){
            return 2;
        } else if (group_num == 4){
            return 2;
        } else {
            return 0;
        }
    } else {
        return 0;
    }
}


AtomisticMCSHFunction get_mcsh_function(int mcsh_order, int group_num)
{
    if (mcsh_order == 0) {
        if (group_num == 1) {
            return calc_MCSH_0_1;
        } 
    } else if (mcsh_order == 1) {
        if (group_num == 1) {
            return calc_MCSH_1_1;
        } 
    } else if (mcsh_order == 2) {
        if (group_num == 1) {
            return calc_MCSH_2_1;
        } else if (group_num == 2){
            return calc_MCSH_2_2;    
        }
    } else if (mcsh_order == 3) {
        if (group_num == 1) {
            return calc_MCSH_3_1;
        } else if (group_num == 2){
            return calc_MCSH_3_2;
        } else if (group_num == 3){
            return calc_MCSH_3_3;
        }
    } else if (mcsh_order == 4) {
        if (group_num == 1) {
            return calc_MCSH_4_1;
        } else if (group_num == 2){
            return calc_MCSH_4_2;
        } else if (group_num == 3){
            return calc_MCSH_4_3;
        } else if (group_num == 4){
            return calc_MCSH_4_4;
        }
    }
}
