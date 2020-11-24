/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_3322668207476132172);
void inv_err_fun(double *nom_x, double *true_x, double *out_3454527342346442355);
void H_mod_fun(double *state, double *out_3918898060006261843);
void f_fun(double *state, double dt, double *out_2000204643949704607);
void F_fun(double *state, double dt, double *out_7131835208994045609);
void h_25(double *state, double *unused, double *out_2709173137952021138);
void H_25(double *state, double *unused, double *out_3349371623221784425);
void h_24(double *state, double *unused, double *out_7101462010527390649);
void H_24(double *state, double *unused, double *out_930663648763150672);
void h_30(double *state, double *unused, double *out_5930704338766520746);
void H_30(double *state, double *unused, double *out_5483473539538583911);
void h_26(double *state, double *unused, double *out_3927877053138189462);
void H_26(double *state, double *unused, double *out_2355221531717424848);
void h_27(double *state, double *unused, double *out_7710846405612866783);
void H_27(double *state, double *unused, double *out_4195891551701958599);
void h_29(double *state, double *unused, double *out_52418851207805658);
void H_29(double *state, double *unused, double *out_6058032153394120086);
void h_28(double *state, double *unused, double *out_3759231009639563890);
void H_28(double *state, double *unused, double *out_7948496126999645538);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
