#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  //px, py, v, phi, phi_dot
  n_x_ = 5;

  // added the process noise from std_a_ and std_yawdd_
  n_aug_ = 7;

  // spreading paramter
  lambda_ = 3 - n_aug_;

  //init weights
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_/(lambda_ + n_aug_);
  for (int i = 1; i < 2*n_aug_ + 1; ++i)
  {
      weights_(i) = 0.5/(lambda_ + n_aug_);
  }

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
      if (!is_initialized_)
    {
      if (meas_package.sensor_type_ == MeasurementPackage::LASER)
      {
        x_ << meas_package.raw_measurements_[0], //px
              meas_package.raw_measurements_[1], //py
              0, //v
              0, //yaw
              0; //yaw_dot
      }
      else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
      {
        double rho = meas_package.raw_measurements_[0];
        double phi = meas_package.raw_measurements_[1];
        x_ << rho*cosf(phi), //px
              rho*sinf(phi), //py
              0, //v
              0, //yaw
              0; //yaw_dot
      }

      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;
      return;
    }
}

void UKF::AugmentedSigmaPoints()
{
  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // fill x_aug
  x_aug.fill(0.0);
  x_aug.head(x_.rows()) = x_;

  // fill P_aug
  P_aug.fill(0.0);
  P_aug.topLeftCorner(P_.rows(), P_.cols()) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //calculate the square root of P
  MatrixXd P_aug_sqrt = P_aug.llt().matrixL();

  //the signma points are a combination of and x_aug and calculated points from P_aug
  MatrixXd sigma_points = pow(lambda_ + n_aug_, 0.5) * P_aug_sqrt.array();

  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; ++i)
  {
      Xsig_aug.col(1 + i) = x_aug + sigma_points.col(i);
      Xsig_aug.col(i + 1+ n_aug_) = x_aug - sigma_points.col(i);
  }
  Xsig_aug_ = Xsig_aug;
}

void UKF::SigmaPointPrediction(double delta_t)
{
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // predict sigma points
  for (int i = 0; i< 2*n_aug_+1; ++i) {
    // extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  Xsig_pred_ = Xsig_pred;
}

void UKF::PredictMeanAndCovariance()
{
  // predict state mean
  for (int i = 0 ; i < 2*n_aug_+1; ++i)
  {
      x_pred_ += Xsig_pred_.col(i) * weights_(i);
  }
  // predict state covariance matrix
  for (int i = 0; i < 2*n_aug_+1; ++i)
  {
      P_pred_ += weights_(i) * (Xsig_pred_.col(i) - x_pred_) * (Xsig_pred_.col(i) - x_pred_).transpose();
  }
}
void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  AugmentedSigmaPoints();
  SigmaPointPrediction(delta_t);
  PredictMeanAndCovariance();
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}