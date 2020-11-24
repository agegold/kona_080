/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6745838954115908166);
void inv_err_fun(double *nom_x, double *true_x, double *out_7460209616644160120);
void H_mod_fun(double *state, double *out_5659993508820358901);
void f_fun(double *state, double dt, double *out_4418064154231686945);
void F_fun(double *state, double dt, double *out_1551122113618590143);
void h_3(double *state, double *unused, double *out_4371216109079759136);
void H_3(double *state, double *unused, double *out_8475332632645043074);
void h_4(double *state, double *unused, double *out_3754451408422766119);
void H_4(double *state, double *unused, double *out_5260716124605144498);
void h_9(double *state, double *unused, double *out_1230359720489346662);
void H_9(double *state, double *unused, double *out_3817493116598694323);
void h_10(double *state, double *unused, double *out_4208734631934960756);
void H_10(double *state, double *unused, double *out_2704434025042125754);
void h_12(double *state, double *unused, double *out_5500115053748913198);
void H_12(double *state, double *unused, double *out_7367271636344354067);
void h_31(double *state, double *unused, double *out_2737052013050023636);
void H_31(double *state, double *unused, double *out_4692181861200934246);
void h_32(double *state, double *unused, double *out_1681446629171863588);
void H_32(double *state, double *unused, double *out_7278323470326339029);
void h_13(double *state, double *unused, double *out_5594889560193853760);
void H_13(double *state, double *unused, double *out_7025223638588436421);
void h_14(double *state, double *unused, double *out_1230359720489346662);
void H_14(double *state, double *unused, double *out_3817493116598694323);
void h_19(double *state, double *unused, double *out_6989142718390137094);
void H_19(double *state, double *unused, double *out_615186734503736858);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);