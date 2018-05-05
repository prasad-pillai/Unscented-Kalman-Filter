#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


/**
* Initializes Unscented Kalman filter
* This is scaffolding, do not modify
*/
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.6;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // dimentaion of the state vector
  n_x_ = 5;

  lambda_ = 3- n_x_;

  //augmented state diamension
  n_aug_ = n_x_ + 2;

  //number of sigma points
  n_sig_ = 2*n_aug_ + 1;

  //declare weights
  weights_ = VectorXd(n_sig_);

  //declare X_sig_pred_
  X_sig_pred_ = MatrixXd(n_x_, n_sig_);

  //initial state vector
  x_ = VectorXd(n_x_);

  //initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  //set NIS
  nis_laser_ = 0;

  nis_radar_ = 0;

  is_initialized_ = false;

}

UKF::~UKF() {}

/**
* @param {MeasurementPackage} meas_package The latest measurement data of
* either radar or laser.
*/
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  //Steps
  //1. initialization
  //2. find delta_t
  //3. call Prediction()
  //4. call update depanding on the sensor

  // initialization
  if(!is_initialized_){
    x_ << 1,1,0,0,0;

    //first measurement saving
    float m1p1 = meas_package.raw_measurements_[0];
    float m1p2 = meas_package.raw_measurements_[1];

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
      x_[0] = m1p1*cos(m1p2);
      x_[1] = m1p1*sin(m1p2);
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
      x_[0] = m1p1;
      x_[1] = m1p2;
    }

    //initialize state uncertainity matrix
    P_ << 1, 0, 0, 0, 0,
    0, 1, 0, 0, 0,
    0, 0, 10, 0, 0,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 10;

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;

  }
  // 2. FInd delta time

  float del_t = (meas_package.timestamp_ - time_us_)/1000000.0;

  //set previous timestamp_
  time_us_ = meas_package.timestamp_;

  // Step 3 call Prediction()
  Prediction(del_t);

  // Step 4 update
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else {
    UpdateLidar(meas_package);
  }

  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}

/**
* Predicts sigma points, the state, and the state covariance matrix.
* @param {double} delta_t the change in time (in seconds) between the last
* measurement and this one.
*/
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // Steps
  // 1. Generate sigma points at time t using augmented sigma points to consider process noise
  // 2. calculate sigma points at time t+delta_t by passing them through the non linear function

  // initialization of variables

  // augmented state vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;
  x_aug[n_x_] = 0;
  x_aug[n_x_ + 1] = 0;

  // augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(P_.rows(), P_.cols()) = P_;

  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  double sqrt_lambda = sqrt(lambda_ + n_aug_);

  //create sigma point matrix at timestep Xt
  MatrixXd X_sig_aug = MatrixXd(n_aug_, n_sig_);
  X_sig_aug.fill(0.0);
  X_sig_aug.col(0) = x_aug;
  for (int i = 0; i < L.cols(); i++) {
    X_sig_aug.col(i + 1) = x_aug + (sqrt_lambda * L.col(i));
    X_sig_aug.col(n_aug_ + i + 1) = x_aug - (sqrt_lambda * L.col(i));
  }

  //create sigma point matrix at timestep Xt + delta_t
  X_sig_pred_.fill(0.0);

  //predict sigma points
  for (int i = 0; i< n_sig_; i++)
  {
    //extract values for better readability
    double p_x = X_sig_aug(0,i);
    double p_y = X_sig_aug(1,i);
    double v = X_sig_aug(2,i);
    double yaw = X_sig_aug(3,i);
    double yaw_dot = X_sig_aug(4,i);
    double nu_a = X_sig_aug(5,i);
    double nu_yaw_dd = X_sig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yaw_dot) > 0.001) {
      px_p = p_x + v/yaw_dot * ( sin (yaw + yaw_dot*delta_t) - sin(yaw));
      py_p = p_y + v/yaw_dot * ( cos(yaw) - cos(yaw+yaw_dot*delta_t) );
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yaw_dot*delta_t;
    double yawd_p = yaw_dot;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yaw_dd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yaw_dd*delta_t;

    //put predicted sigma point into right column of the matrix
    X_sig_pred_(0,i) = px_p;
    X_sig_pred_(1,i) = py_p;
    X_sig_pred_(2,i) = v_p;
    X_sig_pred_(3,i) = yaw_p;
    X_sig_pred_(4,i) = yawd_p;

  }

  // calculate weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i < n_sig_; i++) {
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  // calculate predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    x_ += weights_(i) * X_sig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {

    // state difference
    VectorXd x_diff = X_sig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ += weights_(i) * x_diff * x_diff.transpose() ;
  }
}

/**
* Updates the state and the state covariance matrix using a laser measurement.
* @param {MeasurementPackage} meas_package
*/
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  if (use_laser_ == true){

    std::cout << "LASER" << std::endl;

    int dim = 2;

    // measurement matrix for laser
    MatrixXd H_ = MatrixXd(dim, n_x_);
    H_.fill(0.0);
    H_.row(0)[0] = 1;
    H_.row(1)[1] = 1;

    VectorXd z = meas_package.raw_measurements_;

    //measurement covariance matrix for laser
    MatrixXd R_laser = MatrixXd(dim, dim);
    R_laser.fill(0.0);
    R_laser << std_laspx_*std_laspx_, 0,
    0, std_laspy_*std_laspy_;

    VectorXd y = z - H_ * x_;
    MatrixXd H_t = H_.transpose();
    MatrixXd S = (H_ * P_ * H_t) + R_laser;
    MatrixXd Si = S.inverse();
    MatrixXd K = P_ * H_t * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = Eigen::MatrixXd::Identity(x_size, x_size);
    P_ = (I - (K * H_)) * P_;

    // calculate NIS for Laser
    nis_laser_ = y.transpose() * Si * y;
  }
}

/**
* Updates the state and the state covariance matrix using a radar measurement.
* @param {MeasurementPackage} meas_package
*/
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  // measurement diamension

  if (use_radar_){
    std::cout << "RADAR" << std::endl;


    int dim = 3;

    double weight_0 = lambda_/(lambda_+n_aug_);
    weights_(0) = weight_0;
    for (int i=1; i<n_sig_; i++) {
      double weight = 0.5/(n_aug_+lambda_);
      weights_(i) = weight;
    }

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(dim, n_sig_);

    //transform sigma points into measurement space
    for (int i = 0; i < n_sig_; i++) {

      // extract values for better readibility
      double p_x = X_sig_pred_(0,i);
      double p_y = X_sig_pred_(1,i);
      double v  = X_sig_pred_(2,i);
      double yaw = X_sig_pred_(3,i);

      double v1 = cos(yaw)*v;
      double v2 = sin(yaw)*v;

      // measurement model
      Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);
      Zsig(1,i) = atan2(p_y,p_x);
      Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);
    }

    //mean predicted measurement
    VectorXd z_pred = VectorXd(dim);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred += weights_(i) * Zsig.col(i);
    }

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(dim,dim);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
      //difference
      VectorXd z_diff = Zsig.col(i) - z_pred;

      //angle normalization
      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

      S += weights_(i) * z_diff * z_diff.transpose();
    }


    MatrixXd R_radar = MatrixXd(dim, dim);
    R_radar.fill(0.0);
    R_radar << std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0, std_radrd_*std_radrd_;


    S += R_radar;

    // store incoming radar measurement
    VectorXd z = meas_package.raw_measurements_;

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, dim);

    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < n_sig_; i++) {

      //residual
      VectorXd z_diff = Zsig.col(i) - z_pred;
      //angle normalization
      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

      // state difference
      VectorXd x_diff = X_sig_pred_.col(i) - x_;
      //angle normalization
      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    //calculate Kalman gain K;
    MatrixXd Si = S.inverse();
    MatrixXd K = Tc * Si;


    //update state mean and covariance matrix
    VectorXd y = z - z_pred;
    x_ += K * (y);
    P_ -= K * S * K.transpose();

    // calculate NIS for radar
    nis_radar_ = y.transpose() * Si * y;
  }
}
