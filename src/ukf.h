#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
 public:
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

  /**
   * Calculate the augmented sigma points for prediction
   */
  void AugmentedSigmaPoints();

  /**
   * Use the augmented sigma points to produce predicted sigma points 
   *
   * @param delta_t time delta between two measurements in us
   */
  void SigmaPointPrediction(double delta_t);

  /**
   * Use the predicted sigma points to determine the new predicted state mean
   * and covariance
   */
  void PredictMeanAndCovariance();

  /**
   * @brief Normalize angle
   *
   * @param angle [out] normalized angle
   */
  void NormalizeAngle(double& angle, bool enabled = true);

  /**
   * Calculate the mean and covariance using sigma points
   *
   * @param sigma_points [in] sigma points
   * @param mean  [out] mean vector
   * @param cov  [out] covariance matrix
   */
  void SigmaPointsToMeanCov(MatrixXd& sigma_points, VectorXd& mean,
    MatrixXd& cov, bool normalize_angle = false);

  /**
   * Pass in the measurement package and use it to update the state mean and
   * covariance. The function hanldes the update differently based on the sensor
   * used to generate the measurement package.
   *
   * @param meas_package sensor generated data package
   */
  void Update(MeasurementPackage meas_package);

  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // state vector predicted
  Eigen::VectorXd x_pred_;

  // state covariance matrix predicted
  Eigen::MatrixXd P_pred_;

  // predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred_;

  // generated sigma points matrix
  Eigen::MatrixXd Xsig_gen_;

  // generated augmented sigma points matrix
  Eigen::MatrixXd Xsig_aug_;

  // time when the state is true, in us
  long long time_us_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  // Radar measurement noise standard deviation radius in m
  double std_radr_;

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  // Weights of sigma points
  Eigen::VectorXd weights_;

  // State dimension
  int n_x_;

  // Augmented state dimension
  int n_aug_;

  // Sigma point spreading parameter
  double lambda_;
};

#endif  // UKF_H