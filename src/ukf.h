#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;
  VectorXd x_aug_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* sigma points matrix
  MatrixXd Xsig_;
  MatrixXd Xsig_aug_;
  MatrixXd Xsig_pred_;

  ///* measurement vectors
  VectorXd z_pred_laser_;
  VectorXd z_pred_radar_;

  ///* sigma points in measurement space
  MatrixXd Zsig_laser_;
  MatrixXd Zsig_radar_;

  // Measurement covariance matrix S
  MatrixXd S_laser_;
  MatrixXd S_radar_;

  // Measurement matrix
  MatrixXd H_laser_;

  // Measurement covariance matrix
  MatrixXd R_laser_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Measurement dimensions
  int n_z_laser_;
  int n_z_radar_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* Previous timestamp
  long long previous_timestamp_;

  ///# Minimum px and py
  double p_min_;

  ///# Step counter
  int step_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  void NormalizeAngle(double &phi);
};

#endif /* UKF_H */
