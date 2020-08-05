#include <math.h>

const int IMPLEMENTED_MCSH_TYPE[11][2] = {
    {0, 1},
    {1, 1},
    {2, 1},
    {2, 2},
    {3, 1},
    {3, 2},
    {3, 3},
    {4, 1},
    {4, 2},
    {4, 3},
    {4, 4}
}; // Change this when you implement new type!

const int NUM_IMPLEMENTED_TYPE = 11;

typedef void (*AtomisticMCSHFunction) ( double, double, double, double, double, double, double, double, double *, double *);
// typedef void (*MCSHFunctionType1) ( double, double, double, double, double, double, double, double *, double *);
// typedef void (*MCSHFunctionType2) ( double, double, double, double, double, double, double, double *, double *);
// typedef void (*MCSHFunctionType3) ( double, double, double, double, double, double, double, double *, double *);

double calc_C1(double A, double B, double alpha, double beta);

double calc_C2(double alpha, double beta);

double calc_lambda(double alpha, double beta);

double calc_gamma(double alpha, double beta);

double dx0dx();

double dy0dy();

double dz0dz();

void calc_MCSH_0_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_MCSH_1_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_MCSH_2_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_MCSH_2_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_MCSH_3_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_MCSH_3_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_MCSH_3_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_MCSH_4_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_MCSH_4_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_MCSH_4_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);
void calc_MCSH_4_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv);

int get_mcsh_type(int mcsh_order, int group_num);
AtomisticMCSHFunction get_mcsh_function(int mcsh_order, int group_num);