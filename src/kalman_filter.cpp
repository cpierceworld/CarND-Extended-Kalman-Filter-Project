#include "kalman_filter.h"
#include <iostream>
#include <cmath>

#define PI 3.14159265358979323846

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
	this->x_ = this->F_ * this->x_;
	this->P_ = this->F_ * this->P_ * this->F_.transpose() + this->Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
	MatrixXd I = MatrixXd::Identity(4, 4);
	MatrixXd Ht = this->H_.transpose();
	MatrixXd S = this->H_ * this->P_ * Ht + this->R_;
	MatrixXd Si = S.inverse();

	//Calculating K
	MatrixXd K =  this->P_ * Ht * Si;

	//Calculating innovation
	VectorXd y = z - this->H_ * this->x_;

	//Update State
	this->x_ = this->x_ + (K * y);

	//Update Error Covariance
	this->P_ = (I - K * this->H_) * this->P_ + this->Q_;;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	MatrixXd I = MatrixXd::Identity(4, 4);
	MatrixXd Ht = this->H_.transpose();
	MatrixXd S = this->H_ * this->P_ * Ht + this->R_;
	MatrixXd Si = S.inverse();

	//Calculating K
	MatrixXd K =  this->P_ * Ht * Si;

	//recover state parameters
	double px = this->x_(0);
	double py = this->x_(1);
	double vx = this->x_(2);
	double vy = this->x_(3);

	//Convert current state to "radar"
	double rho = sqrt(px*px + py*py);

	if(fabs(rho) < 0.00001) {
		// avoid division by zero
		rho = 0.00001;
	}

	double phi = atan2(py,px);
	double rho_dot = (px * vx + py * vy) / rho;
	
	VectorXd x_radar = VectorXd(3);
	x_radar(0) = rho;
	x_radar(1) = phi;
	x_radar(2) = rho_dot;

	//Calculating innovation
	VectorXd y = z - x_radar;

	// normalize angle: -PI <= phi <= PI
	if(y(1) < (-1 * PI)) {
		y(1) += 2 * PI;
	}
	else if (y(1) > PI) {
		y(1) -= 2 * PI;
	}

	//Update State
	this->x_ = this->x_ + (K * y);

	//Update Error Covariance
	this->P_ = (I - K * this->H_) * this->P_ + this->Q_;
}
