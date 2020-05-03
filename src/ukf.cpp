#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

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
  std_a_ = 1.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;//M_PI/20.0;

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
  NIS_laser_ = 0.0;
  NIS_radar_ = 0.0;
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

      //P_?
      P_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rhodot = meas_package.raw_measurements_(2);
      double x = rho * cos(phi);
      double y = rho * sin(phi);
      double vx = rhodot * cos(phi);
      double vy = rhodot * sin(phi);
      double v = sqrt(vx * vx + vy * vy);
      x_ << x, y, v, rho, rhodot;
      
      //state covariance matrix
      //***** values can be tuned *****
      P_ << 1000, 0, 0, 0, 0,
            0, 1000, 0, 0, 0,
            0, 0, 1000, 0, 0,
            0, 0, 0, 1000, 0,
            0, 0, 0, 0, 1000;
    }
  }
  else
  {
    Prediction((meas_package.timestamp_ - time_us_)/1000000.0);
    Update(meas_package);
  }
  time_us_ = meas_package.timestamp_;
  is_initialized_ = true;
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

void UKF::SigmaPointsToMeanCov(MatrixXd& sigma_points, VectorXd& mean,
  MatrixXd& cov, bool normalize_angle)
{
  // calculate mean
  for (int i = 0 ; i < 2*n_aug_+1; ++i)
  {
      mean +=  weights_(i) * sigma_points.col(i);
  }
  // calculate covariance matrix
  for (int i = 0; i < 2*n_aug_+1; ++i)
  {
    VectorXd mean_diff = sigma_points.col(i) - mean;
    if (normalize_angle)
    {
      NormalizeAngle(mean_diff(3));
    }
    cov += weights_(i) * mean_diff * mean_diff.transpose();
  }

}

void UKF::PredictMeanAndCovariance()
{
  // create vector for predicted state
  VectorXd x_pred = VectorXd(n_x_);
  // create covariance matrix for prediction
  MatrixXd P_pred = MatrixXd(n_x_, n_x_);

  //default values
  x_pred.fill(0.0);
  P_pred.fill(0.0);

  //predict state mean and covariance
  SigmaPointsToMeanCov(Xsig_pred_, x_pred, P_pred, true);

  x_ = x_pred;
  P_ = P_pred;
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

  // measurement dimension
  int n_z = 2;

  // measurement
  VectorXd z = VectorXd(n_z);

  MatrixXd Zsig_pred = MatrixXd(n_z, 2 * n_aug_ + 1);

  // predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  // predicted measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  // measurement noise matrix
  MatrixXd R = MatrixXd(n_z, n_z);

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // get the predicted measurement sigma points
  for (int i = 0; i < 2*n_aug_+1; ++i)
  {
      Zsig_pred(0, i) = Xsig_pred_(0, i);
      Zsig_pred(1, i) = Xsig_pred_(1, i);
  }

  // calculate mean predicted measurement and covariance
  z_pred.fill(0.0);
  R.fill(0.0);

  R(0,0) = std_laspx_*std_laspx_;
  R(1,1) = std_laspy_*std_laspy_;

  S = R;
  SigmaPointsToMeanCov(Zsig_pred, z_pred, S);

  //get measurement
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; ++i)
  {
    VectorXd z_diff = Zsig_pred.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngle(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }


  // residual
  VectorXd z_diff = z - z_pred;
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

  // calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // update state mean and covariance matrix
  x_ += K * z_diff;
  P_ = P_ - (K * S * K.transpose());
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  //measurement dimension
  int n_z = 3;

  // measurement
  VectorXd z = VectorXd(n_z);

  MatrixXd Zsig_pred = MatrixXd(n_z, 2 * n_aug_ + 1);

  // predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  // predicted measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  // measurement noise matrix
  MatrixXd R = MatrixXd(n_z, n_z);

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // get the predicted measurement sigma points
  for (int i = 0; i < 2*n_aug_+1; ++i)
  {
      double px = Xsig_pred_(0, i);
      double py = Xsig_pred_(1, i);
      double v = Xsig_pred_(2, i);
      double yaw = Xsig_pred_(3, i);

      Zsig_pred(0, i) = sqrt(px*px + py*py); //rho
      Zsig_pred(1, i) = atan2(py,px); //phi
      Zsig_pred(2, i) = v*(px*cos(yaw) + py*sin(yaw))/Zsig_pred(0, i); //rho_dot
  }

  // calculate mean predicted measurement and covariance
  z_pred.fill(0.0);
  R.fill(0.0);

  R(0,0) = std_radr_*std_radr_;
  R(1,1) = std_radphi_*std_radphi_;
  R(2,2) = std_radrd_*std_radrd_;
  S = R;
  SigmaPointsToMeanCov(Zsig_pred, z_pred, S);

  //get measurement
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1],
    meas_package.raw_measurements_[2];

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; ++i)
  {
    VectorXd z_diff = Zsig_pred.col(i) - z_pred;
    NormalizeAngle(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngle(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }


  // residual
  VectorXd z_diff = z - z_pred;
  NormalizeAngle(z_diff(1));
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

  // calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // update state mean and covariance matrix
  x_ += K * z_diff;
  P_ = P_ - (K * S * K.transpose());
}

void UKF::Update(MeasurementPackage meas_package)
{
  if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    // update for Lidar
    UpdateLidar(meas_package);
  }
  else if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    // update for Radar
    UpdateRadar(meas_package);
  }
}

void UKF::NormalizeAngle(double& angle, bool enabled)
{
  if (enabled)
  {
    // angle normalization
    while (angle> M_PI) angle-=2.*M_PI;
    while (angle<-M_PI) angle+=2.*M_PI;
  }
}