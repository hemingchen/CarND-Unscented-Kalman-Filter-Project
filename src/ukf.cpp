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

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.55;

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
  // Program controls
  is_initialized_ = false;

  // State and measurement dims
  n_x_ = 5;
  n_aug_ = 7;
  n_z_laser_ = 2;
  n_z_radar_ = 3;

  // States and augmented states
  x_ = VectorXd(n_x_);
  x_aug_ = VectorXd(n_aug_);

  // Spreading parameter
  lambda_ = 3 - n_aug_;

  // State covariance
  P_ = MatrixXd::Identity(n_x_, n_x_);

  // Sigma points prediction
  Xsig_ = MatrixXd(n_x_, 2 * n_x_ + 1);
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1); // Ignore and remove the noises in x_aug_

  // Sigma points in measurement space
  Zsig_laser_ = MatrixXd(n_z_laser_, 2 * n_aug_ + 1);
  Zsig_radar_ = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);

  // Predicted measurement
  z_pred_laser_ = VectorXd(n_z_laser_);
  z_pred_radar_ = VectorXd(n_z_radar_);

  // Measurement matrix for laser
  H_laser_ = MatrixXd(n_z_laser_, n_x_);
  H_laser_ << 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0;

  // Measurement covariance matrix
  R_laser_ = MatrixXd(n_z_laser_, n_z_laser_);
  R_laser_ << std_laspx_ * std_laspx_, 0,
      0, std_laspy_ * std_laspy_;

  // Measurement covariance matrix S
  S_laser_ = MatrixXd(n_z_laser_, n_z_laser_);
  S_radar_ = MatrixXd(n_z_radar_, n_z_radar_);

  // Sigma point weights
  weights_ = VectorXd(2 * n_aug_ + 1);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    if (i == 0) {
      weights_(i) = lambda_ / (lambda_ + n_aug_);
    } else {
      weights_(i) = 0.5 / (lambda_ + n_aug_);
    }
  }
  cout << "weights_:" << endl << weights_ << endl;

  // Minimum px and py
  p_min_ = 0.0001;

  // Step counter
  step_ = 1;
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

  cout << "\n ################################## processing measurements with UKF at step " << step_
       << " ##################################" << endl;
  step_++;


  /*****************************************************************************
   *  Switch sensor update on/off
   ****************************************************************************/
  if ((!use_laser_) && (meas_package.sensor_type_ == MeasurementPackage::LASER)) {
    cout << "Lidar not used, skipping laser measurements" << endl;
    return;
  }

  if ((!use_radar_) && (meas_package.sensor_type_ == MeasurementPackage::RADAR)) {
    cout << "Radar not used, skipping radar measurements" << endl;
    return;
  }


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "UKF: " << endl;
    x_ = VectorXd(n_x_);
    x_ << 1, 1, 1, 1, 1;

    cout << "initiating UKF..." << endl;
    float px;
    float py;
    float vx;
    float vy;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];
      px = rho * cos(phi);
      py = rho * sin(phi);
      vx = rho_dot * cos(phi);
      vy = rho_dot * sin(phi);
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      px = meas_package.raw_measurements_[0];
      py = meas_package.raw_measurements_[1];
      vx = 0;
      vy = 0;
    }

    // Avoid division by 0
    if (fabs(px) < p_min_) {
      px = p_min_;
      cout << "init px too small" << endl;
    }

    if (fabs(py) < p_min_) {
      py = p_min_;
      cout << "init py too small" << endl;
    }

    // Init states
    x_(0) = px;
    x_(1) = py;
    x_(2) = sqrt(vx * vx + vy * vy);
    cout << "x_:" << endl << x_ << endl;

    // Update timestamp counter
    previous_timestamp_ = meas_package.timestamp_;
    cout << "previous_timestamp_:" << previous_timestamp_ << endl;

    // Done initializing, no need to predict or update
    is_initialized_ = true;

    cout << "UKF init done" << endl;

    return;
  }


  /*****************************************************************************
   *  Predict
   ****************************************************************************/
  cout << "Starting prediction..." << endl;

  // Compute the time elapsed between the current and previous measurements
  cout << "previous timestamp: " << previous_timestamp_ << endl;
  cout << "current  timestamp: " << meas_package.timestamp_ << endl;

  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;  //dt - expressed in seconds
  cout << "dt = " << dt << endl;

  previous_timestamp_ = meas_package.timestamp_;

  Prediction(dt);

  cout << "Prediction done" << endl;


  /*****************************************************************************
   *  Update
   ****************************************************************************/
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
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

  /*****************************************************************************
   *  Compute augmented sigma points
   ****************************************************************************/
  cout << "Computing augmented sigma points..." << endl;

  // Augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // Set augmented mean state
  x_aug_.head(n_x_) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;
  cout << "x_aug:" << endl << x_aug_ << endl;

  // Set augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;
  cout << "P_aug:" << endl << P_aug << endl;

  // Get square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // Get augmented sigma points
  Xsig_aug_.col(0) = x_aug_;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug_.col(i + 1) = x_aug_ + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug_.col(i + 1 + n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  cout << "Xsig_aug_:" << endl << Xsig_aug_ << endl;

  cout << "Augmented sigma points done" << endl;


  /*****************************************************************************
   *  Sigma points prediction
   ****************************************************************************/
  cout << "Predicting sigma points..." << endl;

  // Predict sigma points
  Xsig_pred_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd Xsig_aug_i = Xsig_aug_.col(i);
    VectorXd Xsig_i = Xsig_aug_i.head(n_x_);

    const float px = Xsig_aug_i(0);
    const float py = Xsig_aug_i(1);
    const float v = Xsig_aug_i(2);
    const float psi = Xsig_aug_i(3);
    const float psi_dot = Xsig_aug_i(4);
    const float nu_a = Xsig_aug_i(5);
    const float nu_psi_ddot = Xsig_aug_i(6);

    VectorXd x_incr = VectorXd(n_x_);
    VectorXd x_noise = VectorXd(n_x_);

    // Avoid division by zero
    if (fabs(psi_dot) > 0.001) {
      x_incr << v / psi_dot * (sin(psi + psi_dot * delta_t) - sin(psi)),
          v / psi_dot * (-cos(psi + psi_dot * delta_t) + cos(psi)),
          0,
          psi_dot * delta_t,
          0;
      x_noise << 0.5 * delta_t * delta_t * cos(psi) * nu_a,
          0.5 * delta_t * delta_t * sin(psi) * nu_a,
          delta_t * nu_a,
          0.5 * delta_t * delta_t * nu_psi_ddot,
          delta_t * nu_psi_ddot;
    } else {
      x_incr << v * cos(psi) * delta_t,
          v * sin(psi) * delta_t,
          0,
          psi_dot * delta_t,
          0;
      x_noise << 0.5 * delta_t * delta_t * cos(psi) * nu_a,
          0.5 * delta_t * delta_t * sin(psi) * nu_a,
          delta_t * nu_a,
          0.5 * delta_t * delta_t * nu_psi_ddot,
          delta_t * nu_psi_ddot;
    }

    // Get predicted signma points
    Xsig_pred_.col(i) = Xsig_i + x_incr + x_noise;
  }
  cout << "Xsig_pred_:" << endl << Xsig_pred_ << endl;

  cout << "Sigma point prediction done" << endl;


  /*****************************************************************************
   *  State mean and covariance prediction
   ****************************************************************************/
  cout << "Predicting state mean and covariance..." << endl;

  // Predict state mean
  x_ = Xsig_pred_ * weights_;
  cout << "Predicted state mean x_:" << endl << x_ << endl;

  // Predict state covariance matrix
  P_.fill(0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // Angle normalilzation
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
  cout << "Predicted covariance matrix P_:" << endl << P_ << endl;

  cout << "State mean and covariance prediction done" << endl;
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

  cout << "--------------- Updating with lidar measurements ---------------" << endl;


  /*****************************************************************************
   *  Predict lidar measurements
   ****************************************************************************/
  z_pred_laser_ = H_laser_ * x_;

  VectorXd z_laser;
  z_laser = meas_package.raw_measurements_;
  cout << "Actual laser measurements z_laser:" << endl << z_laser << endl;

  VectorXd z_diff;
  z_diff = z_laser - z_pred_laser_;
  cout << "Lidar z_diff (before normalization):" << endl << z_diff << endl;

  // Angle normalization
  NormalizeAngle(z_diff(1));
  cout << "Lidar z_diff (after normalization):" << endl << z_diff << endl;


  /*****************************************************************************
   *  KF update with lidar measurements
   ****************************************************************************/
  MatrixXd H_laser_t = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * H_laser_t + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * H_laser_t;
  MatrixXd K_laser = PHt * Si;

  // new estimate
  x_ = x_ + K_laser * z_diff;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K_laser * H_laser_) * P_;
  cout << "Updated state mean x_:" << endl << x_ << endl;
  cout << "Updated covariance matrix P_:" << endl << P_ << endl;

  cout << "laser update done" << endl;
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

  cout << "--------------- Updating with radar measurements ---------------" << endl;


  /*****************************************************************************
   *  Predict radar measurements
   ****************************************************************************/
  // Transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd Xsig_pred_i = Xsig_pred_.col(i);

    float px = Xsig_pred_i(0);
    float py = Xsig_pred_i(1);
    float v = Xsig_pred_i(2);
    float psi = Xsig_pred_i(3);
    float psi_dot = Xsig_pred_i(4);

    VectorXd z_i = VectorXd(n_z_radar_);

    // Avoid division by 0
    if (fabs(px) < p_min_ || fabs(py) < p_min_) {
      px = p_min_;
      py = p_min_;
    }

    z_i << sqrt(px * px + py * py),
        atan2(py, px),
        (px * cos(psi) * v + py * sin(psi) * v) / sqrt(px * px + py * py);

    Zsig_radar_.col(i) = z_i;
  }
  cout << "Zsig_radar_:" << endl << Zsig_radar_ << endl;

  // Get predicted measurement
  z_pred_radar_ = Zsig_radar_ * weights_;
  cout << "predicted measurement z_pred_radar_:" << endl << z_pred_radar_ << endl;


  /*****************************************************************************
   *  Calculate innovation covariance matrix S
   ****************************************************************************/
  // Measurement noise matrix
  MatrixXd R_radar = MatrixXd(n_z_radar_, n_z_radar_);
  R_radar.fill(0.0);
  R_radar(0, 0) = std_radr_ * std_radr_;
  R_radar(1, 1) = std_radphi_ * std_radphi_;
  R_radar(2, 2) = std_radrd_ * std_radrd_;

  S_radar_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // Residual
    VectorXd z_diff = Zsig_radar_.col(i) - z_pred_radar_;

    // Angle normalization
    NormalizeAngle(z_diff(1));

    S_radar_ += weights_(i) * z_diff * z_diff.transpose();
  }
  S_radar_ += R_radar;
  cout << "predicted covariance matrix S_radar_:" << endl << S_radar_ << endl;


  /*****************************************************************************
   *  UKF update with radar measurements
   ****************************************************************************/
  // Calculate cross correlation matrix
  MatrixXd Tc_radar = MatrixXd(n_x_, n_z_radar_);
  Tc_radar.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // Residual
    VectorXd z_diff = Zsig_radar_.col(i) - z_pred_radar_;

    // Angle normalization
    NormalizeAngle(z_diff(1));

    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // Angle normalization
    NormalizeAngle(x_diff(3));

    Tc_radar += weights_(i) * x_diff * z_diff.transpose();
  }
  cout << "Cross correlation matrix Tc_radar:" << endl << Tc_radar << endl;

  // Kalman gain
  MatrixXd K_radar = Tc_radar * S_radar_.inverse();
  cout << "Kalman gain K_radar:" << endl << K_radar << endl;

  // Residual
  VectorXd z_radar;
  z_radar = meas_package.raw_measurements_;
  cout << "Actual radar measurements z_radar:" << endl << z_radar << endl;

  VectorXd z_diff;
  z_diff = z_radar - z_pred_radar_;
  cout << "Radar z_diff (before normalization):" << endl << z_diff << endl;

  // Angle normalization
  NormalizeAngle(z_diff(1));

  cout << "Radar z_diff (after normalization):" << endl << z_diff << endl;

  // Update state mean and covariance matrix
  x_ = x_ + K_radar * z_diff;
  P_ = P_ - K_radar * S_radar_ * K_radar.transpose();
  cout << "Updated state mean x_:" << endl << x_ << endl;
  cout << "Updated covariance matrix P_:" << endl << P_ << endl;

  cout << "Radar update done" << endl;
}

void UKF::NormalizeAngle(double &phi) {
  phi = atan2(sin(phi), cos(phi));
}
