/* Copyright 2018, Jesse Chen
* This is a EPnP implementation with Eigen only.
*
*/
#ifndef EPNP_EIGEN_H_
#define EPNP_EIGEN_H_

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Dense>
#include <Eigen/Geometry>

struct DistPattern {
  int a;
  int b;
};

class EPnPEigen {
 public:
  EPnPEigen(Eigen::MatrixXd& points3d, Eigen::MatrixXd& points2d, Eigen::Matrix3d& K);
  ~EPnPEigen(){}

  void computePose(Eigen::Matrix3d &R, Eigen::Vector3d &t);

 private:
  void chooseControlPoints(void);
  void computeBaryCentricCoordinates(void);
  void calculateM(Eigen::MatrixXd& M);
  void computeL6x10(const Eigen::MatrixXd& U, Eigen::MatrixXd& L6x10); 
  void computeRho(Eigen::VectorXd& rho);
  void findBetasApprox1(const Eigen::MatrixXd& L6x10, const Eigen::VectorXd& rho, double* betas);
  void findBetasApprox2(const Eigen::MatrixXd& L6x10, const Eigen::VectorXd& rho, double* betas);
  void findBetasApprox3(const Eigen::MatrixXd& L6x10, const Eigen::VectorXd& rho, double* betas);
  void computeGaussNewtonJacobian(const Eigen::MatrixXd& L6x10, double betas[4], Eigen::MatrixXd& jacobian);
  void computeResiduals(const Eigen::MatrixXd& U, double betas[4], Eigen::VectorXd& residuals);
  void doGaussNewtonOptimization(const Eigen::MatrixXd& U, const Eigen::MatrixXd& L6x10, double betas[4]);
  void computeControlPointsUnderCameraCoord(const Eigen::MatrixXd& U, double betas[4]);
  void computeReferencePointsUnderCameraCoord(void);
  void solveForSign(void);
  void estimateRt(Eigen::Matrix3d& R, Eigen::Vector3d& t);
  double computeRt(const Eigen::MatrixXd&U, double betas[4], Eigen::Matrix3d& R, Eigen::Vector3d& t);
  double reprojectionError(Eigen::Matrix3d& R, Eigen::Vector3d& t);

  Eigen::MatrixXd reference_3d_points_;
  Eigen::MatrixXd reference_2d_points_;
  Eigen::MatrixXd reference_3d_points_camera_coord_;
  Eigen::MatrixXd control_3d_points_;
  Eigen::MatrixXd control_3d_points_camera_coord_;
  Eigen::MatrixXd bary_centric_coord_;
  int reference_points_count_;
  double fu_, fv_, uc_, vc_;
};


#endif