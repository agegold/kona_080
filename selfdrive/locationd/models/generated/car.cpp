
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_3322668207476132172) {
   out_3322668207476132172[0] = delta_x[0] + nom_x[0];
   out_3322668207476132172[1] = delta_x[1] + nom_x[1];
   out_3322668207476132172[2] = delta_x[2] + nom_x[2];
   out_3322668207476132172[3] = delta_x[3] + nom_x[3];
   out_3322668207476132172[4] = delta_x[4] + nom_x[4];
   out_3322668207476132172[5] = delta_x[5] + nom_x[5];
   out_3322668207476132172[6] = delta_x[6] + nom_x[6];
   out_3322668207476132172[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_3454527342346442355) {
   out_3454527342346442355[0] = -nom_x[0] + true_x[0];
   out_3454527342346442355[1] = -nom_x[1] + true_x[1];
   out_3454527342346442355[2] = -nom_x[2] + true_x[2];
   out_3454527342346442355[3] = -nom_x[3] + true_x[3];
   out_3454527342346442355[4] = -nom_x[4] + true_x[4];
   out_3454527342346442355[5] = -nom_x[5] + true_x[5];
   out_3454527342346442355[6] = -nom_x[6] + true_x[6];
   out_3454527342346442355[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_3918898060006261843) {
   out_3918898060006261843[0] = 1.0;
   out_3918898060006261843[1] = 0.0;
   out_3918898060006261843[2] = 0.0;
   out_3918898060006261843[3] = 0.0;
   out_3918898060006261843[4] = 0.0;
   out_3918898060006261843[5] = 0.0;
   out_3918898060006261843[6] = 0.0;
   out_3918898060006261843[7] = 0.0;
   out_3918898060006261843[8] = 0.0;
   out_3918898060006261843[9] = 1.0;
   out_3918898060006261843[10] = 0.0;
   out_3918898060006261843[11] = 0.0;
   out_3918898060006261843[12] = 0.0;
   out_3918898060006261843[13] = 0.0;
   out_3918898060006261843[14] = 0.0;
   out_3918898060006261843[15] = 0.0;
   out_3918898060006261843[16] = 0.0;
   out_3918898060006261843[17] = 0.0;
   out_3918898060006261843[18] = 1.0;
   out_3918898060006261843[19] = 0.0;
   out_3918898060006261843[20] = 0.0;
   out_3918898060006261843[21] = 0.0;
   out_3918898060006261843[22] = 0.0;
   out_3918898060006261843[23] = 0.0;
   out_3918898060006261843[24] = 0.0;
   out_3918898060006261843[25] = 0.0;
   out_3918898060006261843[26] = 0.0;
   out_3918898060006261843[27] = 1.0;
   out_3918898060006261843[28] = 0.0;
   out_3918898060006261843[29] = 0.0;
   out_3918898060006261843[30] = 0.0;
   out_3918898060006261843[31] = 0.0;
   out_3918898060006261843[32] = 0.0;
   out_3918898060006261843[33] = 0.0;
   out_3918898060006261843[34] = 0.0;
   out_3918898060006261843[35] = 0.0;
   out_3918898060006261843[36] = 1.0;
   out_3918898060006261843[37] = 0.0;
   out_3918898060006261843[38] = 0.0;
   out_3918898060006261843[39] = 0.0;
   out_3918898060006261843[40] = 0.0;
   out_3918898060006261843[41] = 0.0;
   out_3918898060006261843[42] = 0.0;
   out_3918898060006261843[43] = 0.0;
   out_3918898060006261843[44] = 0.0;
   out_3918898060006261843[45] = 1.0;
   out_3918898060006261843[46] = 0.0;
   out_3918898060006261843[47] = 0.0;
   out_3918898060006261843[48] = 0.0;
   out_3918898060006261843[49] = 0.0;
   out_3918898060006261843[50] = 0.0;
   out_3918898060006261843[51] = 0.0;
   out_3918898060006261843[52] = 0.0;
   out_3918898060006261843[53] = 0.0;
   out_3918898060006261843[54] = 1.0;
   out_3918898060006261843[55] = 0.0;
   out_3918898060006261843[56] = 0.0;
   out_3918898060006261843[57] = 0.0;
   out_3918898060006261843[58] = 0.0;
   out_3918898060006261843[59] = 0.0;
   out_3918898060006261843[60] = 0.0;
   out_3918898060006261843[61] = 0.0;
   out_3918898060006261843[62] = 0.0;
   out_3918898060006261843[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_2000204643949704607) {
   out_2000204643949704607[0] = state[0];
   out_2000204643949704607[1] = state[1];
   out_2000204643949704607[2] = state[2];
   out_2000204643949704607[3] = state[3];
   out_2000204643949704607[4] = state[4];
   out_2000204643949704607[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_2000204643949704607[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_2000204643949704607[7] = state[7];
}
void F_fun(double *state, double dt, double *out_7131835208994045609) {
   out_7131835208994045609[0] = 1;
   out_7131835208994045609[1] = 0;
   out_7131835208994045609[2] = 0;
   out_7131835208994045609[3] = 0;
   out_7131835208994045609[4] = 0;
   out_7131835208994045609[5] = 0;
   out_7131835208994045609[6] = 0;
   out_7131835208994045609[7] = 0;
   out_7131835208994045609[8] = 0;
   out_7131835208994045609[9] = 1;
   out_7131835208994045609[10] = 0;
   out_7131835208994045609[11] = 0;
   out_7131835208994045609[12] = 0;
   out_7131835208994045609[13] = 0;
   out_7131835208994045609[14] = 0;
   out_7131835208994045609[15] = 0;
   out_7131835208994045609[16] = 0;
   out_7131835208994045609[17] = 0;
   out_7131835208994045609[18] = 1;
   out_7131835208994045609[19] = 0;
   out_7131835208994045609[20] = 0;
   out_7131835208994045609[21] = 0;
   out_7131835208994045609[22] = 0;
   out_7131835208994045609[23] = 0;
   out_7131835208994045609[24] = 0;
   out_7131835208994045609[25] = 0;
   out_7131835208994045609[26] = 0;
   out_7131835208994045609[27] = 1;
   out_7131835208994045609[28] = 0;
   out_7131835208994045609[29] = 0;
   out_7131835208994045609[30] = 0;
   out_7131835208994045609[31] = 0;
   out_7131835208994045609[32] = 0;
   out_7131835208994045609[33] = 0;
   out_7131835208994045609[34] = 0;
   out_7131835208994045609[35] = 0;
   out_7131835208994045609[36] = 1;
   out_7131835208994045609[37] = 0;
   out_7131835208994045609[38] = 0;
   out_7131835208994045609[39] = 0;
   out_7131835208994045609[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_7131835208994045609[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_7131835208994045609[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7131835208994045609[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7131835208994045609[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_7131835208994045609[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_7131835208994045609[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_7131835208994045609[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_7131835208994045609[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_7131835208994045609[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_7131835208994045609[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7131835208994045609[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7131835208994045609[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_7131835208994045609[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_7131835208994045609[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_7131835208994045609[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7131835208994045609[56] = 0;
   out_7131835208994045609[57] = 0;
   out_7131835208994045609[58] = 0;
   out_7131835208994045609[59] = 0;
   out_7131835208994045609[60] = 0;
   out_7131835208994045609[61] = 0;
   out_7131835208994045609[62] = 0;
   out_7131835208994045609[63] = 1;
}
void h_25(double *state, double *unused, double *out_2709173137952021138) {
   out_2709173137952021138[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3349371623221784425) {
   out_3349371623221784425[0] = 0;
   out_3349371623221784425[1] = 0;
   out_3349371623221784425[2] = 0;
   out_3349371623221784425[3] = 0;
   out_3349371623221784425[4] = 0;
   out_3349371623221784425[5] = 0;
   out_3349371623221784425[6] = 1;
   out_3349371623221784425[7] = 0;
}
void h_24(double *state, double *unused, double *out_7101462010527390649) {
   out_7101462010527390649[0] = state[4];
   out_7101462010527390649[1] = state[5];
}
void H_24(double *state, double *unused, double *out_930663648763150672) {
   out_930663648763150672[0] = 0;
   out_930663648763150672[1] = 0;
   out_930663648763150672[2] = 0;
   out_930663648763150672[3] = 0;
   out_930663648763150672[4] = 1;
   out_930663648763150672[5] = 0;
   out_930663648763150672[6] = 0;
   out_930663648763150672[7] = 0;
   out_930663648763150672[8] = 0;
   out_930663648763150672[9] = 0;
   out_930663648763150672[10] = 0;
   out_930663648763150672[11] = 0;
   out_930663648763150672[12] = 0;
   out_930663648763150672[13] = 1;
   out_930663648763150672[14] = 0;
   out_930663648763150672[15] = 0;
}
void h_30(double *state, double *unused, double *out_5930704338766520746) {
   out_5930704338766520746[0] = state[4];
}
void H_30(double *state, double *unused, double *out_5483473539538583911) {
   out_5483473539538583911[0] = 0;
   out_5483473539538583911[1] = 0;
   out_5483473539538583911[2] = 0;
   out_5483473539538583911[3] = 0;
   out_5483473539538583911[4] = 1;
   out_5483473539538583911[5] = 0;
   out_5483473539538583911[6] = 0;
   out_5483473539538583911[7] = 0;
}
void h_26(double *state, double *unused, double *out_3927877053138189462) {
   out_3927877053138189462[0] = state[7];
}
void H_26(double *state, double *unused, double *out_2355221531717424848) {
   out_2355221531717424848[0] = 0;
   out_2355221531717424848[1] = 0;
   out_2355221531717424848[2] = 0;
   out_2355221531717424848[3] = 0;
   out_2355221531717424848[4] = 0;
   out_2355221531717424848[5] = 0;
   out_2355221531717424848[6] = 0;
   out_2355221531717424848[7] = 1;
}
void h_27(double *state, double *unused, double *out_7710846405612866783) {
   out_7710846405612866783[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4195891551701958599) {
   out_4195891551701958599[0] = 0;
   out_4195891551701958599[1] = 0;
   out_4195891551701958599[2] = 0;
   out_4195891551701958599[3] = 1;
   out_4195891551701958599[4] = 0;
   out_4195891551701958599[5] = 0;
   out_4195891551701958599[6] = 0;
   out_4195891551701958599[7] = 0;
}
void h_29(double *state, double *unused, double *out_52418851207805658) {
   out_52418851207805658[0] = state[1];
}
void H_29(double *state, double *unused, double *out_6058032153394120086) {
   out_6058032153394120086[0] = 0;
   out_6058032153394120086[1] = 1;
   out_6058032153394120086[2] = 0;
   out_6058032153394120086[3] = 0;
   out_6058032153394120086[4] = 0;
   out_6058032153394120086[5] = 0;
   out_6058032153394120086[6] = 0;
   out_6058032153394120086[7] = 0;
}
void h_28(double *state, double *unused, double *out_3759231009639563890) {
   out_3759231009639563890[0] = state[5];
   out_3759231009639563890[1] = state[6];
}
void H_28(double *state, double *unused, double *out_7948496126999645538) {
   out_7948496126999645538[0] = 0;
   out_7948496126999645538[1] = 0;
   out_7948496126999645538[2] = 0;
   out_7948496126999645538[3] = 0;
   out_7948496126999645538[4] = 0;
   out_7948496126999645538[5] = 1;
   out_7948496126999645538[6] = 0;
   out_7948496126999645538[7] = 0;
   out_7948496126999645538[8] = 0;
   out_7948496126999645538[9] = 0;
   out_7948496126999645538[10] = 0;
   out_7948496126999645538[11] = 0;
   out_7948496126999645538[12] = 0;
   out_7948496126999645538[13] = 0;
   out_7948496126999645538[14] = 1;
   out_7948496126999645538[15] = 0;
}
}

extern "C"{
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
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
