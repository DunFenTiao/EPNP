/* Copyright 2018, Jesse Chen
* This is a EPnP implementation with Eigen only.
*
*/
/*
COMPLETE BY 
JIAYAOMA 2019.4.17
*/
#include "eigenEPNP.h"
#include <gtest/gtest.h>

#include <iostream>
using namespace std;
using namespace Eigen;


EPnPEigen::EPnPEigen(Eigen::MatrixXd& points3d, Eigen::MatrixXd& points2d, Eigen::Matrix3d& K) {
  //初始化3D参考点pwi,2D参考点，点个数N
  reference_3d_points_ = points3d;
  reference_2d_points_ = points2d;  
  reference_points_count_ = reference_3d_points_.rows();
  //初始化4个控制点pci
  control_3d_points_ = Eigen::MatrixXd::Zero(4, 3);
  control_3d_points_camera_coord_ = Eigen::MatrixXd::Zero(4, 3);
  bary_centric_coord_ = Eigen::MatrixXd::Zero(reference_points_count_, 4);
  reference_3d_points_camera_coord_ = Eigen::MatrixXd::Zero(reference_points_count_, 3);

  fu_ = K(0, 0);
  fv_ = K(1, 1);
  uc_ = K(0, 2);
  vc_ = K(1, 2);
  

}


void EPnPEigen::chooseControlPoints(void){
  double lambda;
  Eigen::VectorXd eigvec;
  //第一个控制点为重心，cw1= sum{pwi} / N
  Eigen::MatrixXd pointsSum = reference_3d_points_.colwise().sum();
  pointsSum = pointsSum/reference_points_count_;
  control_3d_points_.row(0) = pointsSum;
  //去中心，A = pwi-cw1
  Eigen::MatrixXd centroidMat = pointsSum.replicate(reference_points_count_, 1);
  Eigen::MatrixXd PW0 = reference_3d_points_ - centroidMat;
  //计算ATA特征值
  Eigen::MatrixXd PW0t = PW0;
  PW0t.transposeInPlace();
  Eigen::MatrixXd PW0tPW0 = PW0t * PW0;


  //计算剩下三个控制点
  // 1. 用特征向量，特征根来求
  /*
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(PW0tPW0);
  Eigen::VectorXd eigenval = es.eigenvalues();
  Eigen::VectorXd k = (eigenval/reference_points_count_).cwiseSqrt();
  int sign_value[3] = {1, -1, -1};

  for (int i = 2; i >= 0; i--){
    lambda = k(i);
    eigvec = es.eigenvectors().col(i);
    cout << "i" << i << endl;
    cout << "lamda" << lambda << endl;
    cout << "eigvec" << eigvec << endl;
    cout << "eigvec.transpose()" << eigvec.transpose() << endl;
    control_3d_points_.row(3-i) = control_3d_points_.row(0) + sign_value[i] * lambda * eigvec.transpose();
  }
  */
  
  // 2. 用奇异值分解来求
  double k;
  JacobiSVD<Eigen::MatrixXd> svd(PW0tPW0, ComputeThinU | ComputeThinV );
  Matrix3d V = svd.matrixV(), U = svd.matrixU();
  Matrix3d S = U.inverse() * PW0tPW0 * V.transpose().inverse();
  MatrixXd UT= U.transpose();
  //cout << U << endl;
  //cout << V << endl;
  
  //计算剩下三个控制点  cwj = cw1 + sqrt(lamdba{c,j-1}) * v{c,j-1}
  for (int i = 1; i < 4; i++){
    for(int j =0 ; j < 3 ; j++){
      k = sqrt(S(i-1,i-1)/reference_points_count_);
      control_3d_points_(i,j) = control_3d_points_(0,j)+ k * UT(i-1,j);
    }
  }
  
  //cout << "control_3d_points_" << control_3d_points_ << endl;

}


void EPnPEigen::computeBaryCentricCoordinates(){
  Eigen::MatrixXd CC(3,3);
  // 其他三个控制点都减去重心控制点
  for (int i = 0; i < 3; i++){
    CC.row(i) = control_3d_points_.row(i+1) - control_3d_points_.row(0);
  }	
  CC.transposeInPlace();

  Eigen::MatrixXd CC_inv = CC.inverse();
  // barycentric 是 4*N维矩阵,记录权重参数 1 = ai0+ai1+ai2+ai3
  //[0,1,2,3,4,5,6] n个参考点,参考点到所有控制点 = 其他控制点到重心 * 参考点到重心
  //[a]到重心的ai0
  //[1]到第一个控制点ai1
  //[2]到第一个控制点ai2
  //[3]到第一个控制点ai3
  Eigen::MatrixXd pt_3d_diff_mat(1,3);
  Eigen::Vector3d pt_3d_diff_vec;
  double a;
  for (int i = 0; i < reference_points_count_; i++){
    pt_3d_diff_mat = reference_3d_points_.row(i) - control_3d_points_.row(0);
    pt_3d_diff_vec = pt_3d_diff_mat.transpose();
    pt_3d_diff_vec = CC_inv * pt_3d_diff_vec;
    a = 1.0 - pt_3d_diff_vec.sum();
    bary_centric_coord_(i,0) = a;
    bary_centric_coord_.block(i, 1, 1, 3) = pt_3d_diff_vec.transpose();
  }
}


void EPnPEigen::calculateM(Eigen::MatrixXd& M){
  /*求M矩阵，mx =0 其中m为2n*12维，X为12*1维度 
  \sum^4_{j=1}a_{ij}f_ux^c_j+a_{ij}(u_c-u_i)z^c_j=0,\\
  \sum^4_{j=1}a_{ij}f_vy^c_j+a_{ij}(v_c-v_i)z^c_j=0,
  [ aij*fu,   0  ,aij*uci
      0   ,aij*fv,aij*vci ] * [x,y,z]T =0
  */
  double uci, vci, barycentric;
  for (int i = 0; i < reference_points_count_; i++){
  	uci = uc_ - reference_2d_points_(i, 0);
    vci = vc_ - reference_2d_points_(i, 1);
    for (int j = 0; j < 4; j++){
      barycentric = bary_centric_coord_(i,j);
      M(2*i, 3*j    ) = barycentric*fu_;
      M(2*i, 3*j + 1) = 0;
      M(2*i, 3*j + 2) = barycentric*uci;

      M(2*i + 1, 3*j    ) = 0;
      M(2*i + 1, 3*j + 1) = barycentric*fv_;
      M(2*i + 1, 3*j + 2) = barycentric*vci;
    }  	
  }

}


void EPnPEigen::computeL6x10(const Eigen::MatrixXd& U, Eigen::MatrixXd& L6x10){
  //MTM为12*2n *2n *12，特征根维度取N=4时
  // ||cci-ccj||^2 = ||cwi-cwj||^2
  // L6*10 * beta 10*1 = p 6*1
  Eigen::MatrixXd V = U.block(0, 0, 12, 4);
  Eigen::MatrixXd DiffMat = Eigen::MatrixXd::Zero(18, 4);
  DistPattern diff_pattern[6] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
  
  for (int i = 0; i < 6; i++){
    DiffMat.block(3*i, 0, 3, 4) = V.block(3 * diff_pattern[i].a, 0, 3, 4) - V.block(3 * diff_pattern[i].b, 0, 3, 4);  	
  } 	

  Eigen::Vector3d v1, v2, v3, v4;
  for (int i = 0; i < 6; i++){
  	v1 = DiffMat.block(3*i, 0, 3, 1);
  	v2 = DiffMat.block(3*i, 1, 3, 1);
  	v3 = DiffMat.block(3*i, 2, 3, 1);
  	v4 = DiffMat.block(3*i, 3, 3, 1);

  	L6x10.block(i, 0, 1, 10) << v1.dot(v1), 2*v1.dot(v2), v2.dot(v2), 2*v1.dot(v3), 2*v2.dot(v3), 
  	                            v3.dot(v3), 2*v1.dot(v4), 2*v2.dot(v4), 2*v3.dot(v4), v4.dot(v4);
  }
}


void EPnPEigen::computeRho(Eigen::VectorXd& rho){
  Eigen::Vector3d control_point_a, control_point_b, control_point_diff;
  DistPattern diff_pattern[6] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};

  for (int i = 0; i < 6; i++){
    control_point_a << control_3d_points_(diff_pattern[i].a, 0), control_3d_points_(diff_pattern[i].a, 1), control_3d_points_(diff_pattern[i].a, 2);
    control_point_b << control_3d_points_(diff_pattern[i].b, 0), control_3d_points_(diff_pattern[i].b, 1), control_3d_points_(diff_pattern[i].b, 2);
    control_point_diff = control_point_a - control_point_b;

    rho(i) = control_point_diff.dot(control_point_diff);
  }
}

// betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
// betas_approx_1 = [B11 B12     B13         B14]
// for N = 4
void EPnPEigen::findBetasApprox1(const Eigen::MatrixXd& L6x10, const Eigen::VectorXd& rho, double* betas){
  Eigen::MatrixXd L6x4(6, 4);

  L6x4.block(0, 0, 6, 1) = L6x10.block(0, 0, 6, 1);
  L6x4.block(0, 1, 6, 1) = L6x10.block(0, 1, 6, 1);
  L6x4.block(0, 2, 6, 1) = L6x10.block(0, 3, 6, 1);
  L6x4.block(0, 3, 6, 1) = L6x10.block(0, 6, 6, 1);

  //Eigen::VectorXd B = L6x4.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rho);
  Eigen::VectorXd B = L6x4.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rho);
  
  if (B(0) < 0) {
    betas[0] = sqrt(-B(0));
    betas[1] = -B(1) / betas[0];
    betas[2] = -B(2) / betas[0];
    betas[3] = -B(3) / betas[0];
  } else {
    betas[0] = sqrt(B(0));
    betas[1] = B(1) / betas[0];
    betas[2] = B(2) / betas[0];
    betas[3] = B(3) / betas[0];
  }

}

// betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
// betas_approx_2 = [B11 B12 B22                            ]              N=2
void EPnPEigen::findBetasApprox2(const Eigen::MatrixXd& L6x10, const Eigen::VectorXd& rho, double* betas){
  Eigen::MatrixXd L6x3(6, 3);

  L6x3.block(0, 0, 6, 1) = L6x10.block(0, 0, 6, 1);
  L6x3.block(0, 1, 6, 1) = L6x10.block(0, 1, 6, 1);
  L6x3.block(0, 2, 6, 1) = L6x10.block(0, 2, 6, 1);
  

  //Eigen::VectorXd B = L6x4.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rho);
  Eigen::VectorXd B = L6x3.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rho);
  
  if (B(0) < 0) {
    betas[0] = sqrt(-B(0));
    betas[1] = (B(2)<0)?sqrt(-(B(2))):0.0;
   
  } else {
    betas[0] = sqrt(B(0));
    betas[1] = (B(2)>0)?sqrt((B(2))):0.0;
  }
  if (B[1] < 0) betas[0] = -betas[0];
  betas[2] = 0.0;
  betas[3] = 0.0;
}

// betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
// betas_approx_3 = [B11 B12 B22 B13 B23                    ]               N=3
void EPnPEigen::findBetasApprox3(const Eigen::MatrixXd& L6x10, const Eigen::VectorXd& rho, double* betas){
  Eigen::MatrixXd L6x5(6, 5);

  L6x5.block(0, 0, 6, 1) = L6x10.block(0, 0, 6, 1);
  L6x5.block(0, 1, 6, 1) = L6x10.block(0, 1, 6, 1);
  L6x5.block(0, 2, 6, 1) = L6x10.block(0, 2, 6, 1);
  L6x5.block(0, 3, 6, 1) = L6x10.block(0, 3, 6, 1);
  L6x5.block(0, 4, 6, 1) = L6x10.block(0, 4, 6, 1);
  

  //Eigen::VectorXd B = L6x4.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rho);
  Eigen::VectorXd B = L6x5.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rho);
  
  if (B(0) < 0) {
    betas[0] = sqrt(-B(0));
    betas[1] = (B(2)<0)?sqrt(-(B(2))):0.0;
   
  } else {
    betas[0] = sqrt(B(0));
    betas[1] = (B(2)>0)?sqrt((B(2))):0.0;
  }
  if (B[1] < 0) betas[0] = -betas[0];
  betas[2] = B(3) / betas[0];
  betas[3] = 0.0;
}

void EPnPEigen::computeResiduals(const Eigen::MatrixXd& U, double betas[4], Eigen::VectorXd& residuals){
  Eigen::MatrixXd V = U.block(0, 0, 12, 4);
  DistPattern diff_pattern[6] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
  Eigen::VectorXd CC(12, 1);
  Eigen::Vector3d Ca, Cb;
  Eigen::MatrixXd Wa, Wb;
  Eigen::Vector3d Vwa, Vwb;

  CC = betas[0] * V.block(0, 0, 12, 1) + betas[1] * V.block(0, 1, 12, 1) + betas[2] * V.block(0, 2, 12, 1) + betas[3] * V.block(0, 3, 12, 1);

  for (int i = 0; i < 6; i++){
    Ca = CC.block(3*diff_pattern[i].a, 0, 3, 1);
    Cb = CC.block(3*diff_pattern[i].b, 0, 3, 1);
    Wa = control_3d_points_.block(diff_pattern[i].a, 0, 1, 3);
    Wb = control_3d_points_.block(diff_pattern[i].b, 0, 1, 3);

    Ca = Ca - Cb;
    Cb = Ca;
    double d1 = Ca.dot(Cb);

    Wa = Wa - Wb;
    Wa.transposeInPlace();
    Vwa = Wa;
    Vwb = Vwa;
    double d2 = Vwa.dot(Vwb);

    residuals(i) = d1 - d2;  	
  }
}


void EPnPEigen::computeGaussNewtonJacobian(const Eigen::MatrixXd& L6x10, double betas[4], Eigen::MatrixXd& jacobian){
  Eigen::MatrixXd L2J(10, 4);

  L2J << 2*betas[0],          0,          0,          0,  // 0
           betas[1],   betas[0],          0,          0,  // 1
                  0, 2*betas[1],          0,          0,  // 2
           betas[2],          0,   betas[0],          0,  // 3
                  0,   betas[2],   betas[1],          0,  // 4
                  0,          0, 2*betas[2],          0,  // 5
           betas[3],          0,          0,   betas[0],  // 6
                  0,   betas[3],          0,   betas[1],  // 7
                  0,          0,   betas[3],   betas[2],  // 8
                  0,          0,          0, 2*betas[3];  // 9

  jacobian = L6x10 * L2J;
}


void EPnPEigen::doGaussNewtonOptimization(const Eigen::MatrixXd& U, const Eigen::MatrixXd& L6x10, double betas[4]){
  const int iterations_number = 5;
  Eigen::MatrixXd jacobian(6, 4);
  Eigen::VectorXd residuals(6, 1);
  Eigen::MatrixXd JtJ(4, 4), JtJ_inv(4, 4);
  //Eigen::Map<Eigen::Vector4d> Vb(betas);
  Eigen::Vector4d Vb;
  Eigen::Vector4d jacobian_res;

  Vb << betas[0], betas[1], betas[2], betas[3];
  for (int i = 0; i < iterations_number; i++){
    computeGaussNewtonJacobian(L6x10, betas, jacobian);
    computeResiduals(U, betas, residuals);

    JtJ = jacobian.transpose() * jacobian;
    JtJ_inv = JtJ.inverse(); 
    jacobian_res = jacobian.transpose() * residuals;

    Vb = Vb - JtJ_inv * jacobian_res;

    betas[0] = Vb(0);
    betas[1] = Vb(1);
    betas[2] = Vb(2);
    betas[3] = Vb(3);
  }
  
}


void EPnPEigen::computeControlPointsUnderCameraCoord(const Eigen::MatrixXd& U, double betas[4]){
  Eigen::MatrixXd V = U.block(0, 0, 12, 4);
  Eigen::MatrixXd control_3d_points_camera_coord_vector(12, 1);

  control_3d_points_camera_coord_vector = betas[0] * V.block(0, 0, 12, 1) + betas[1] * V.block(0, 1, 12, 1) +
                                          betas[2] * V.block(0, 2, 12, 1) + betas[3] * V.block(0, 3, 12, 1);
  
  for (int i = 0; i < 4; i++){
    control_3d_points_camera_coord_.block(i, 0, 1, 3) << control_3d_points_camera_coord_vector(3*i), 
                                                         control_3d_points_camera_coord_vector(3*i + 1),
                                                         control_3d_points_camera_coord_vector(3*i + 2);
  }
}


void EPnPEigen::computeReferencePointsUnderCameraCoord(void){
  reference_3d_points_camera_coord_ = bary_centric_coord_ * control_3d_points_camera_coord_;
}


void EPnPEigen::solveForSign(void){
  if (reference_3d_points_camera_coord_(0,2) < 0){
    control_3d_points_camera_coord_ = -1 * control_3d_points_camera_coord_;
    reference_3d_points_camera_coord_ = -1 * reference_3d_points_camera_coord_;
  }
}


void EPnPEigen::estimateRt(Eigen::Matrix3d& R, Eigen::Vector3d& t){
  Eigen::MatrixXd pointsSum = reference_3d_points_.colwise().sum();
  pointsSum = pointsSum/reference_points_count_;
  Eigen::Vector3d P0w = pointsSum.transpose();

  Eigen::MatrixXd centroidMat = pointsSum.replicate(reference_points_count_, 1);
  Eigen::MatrixXd Piw = reference_3d_points_ - centroidMat;
  
  pointsSum = reference_3d_points_camera_coord_.colwise().sum();
  pointsSum = pointsSum/reference_points_count_;
  Eigen::Vector3d P0c = pointsSum.transpose();
  centroidMat = pointsSum.replicate(reference_points_count_, 1);
  Eigen::MatrixXd Pic = reference_3d_points_camera_coord_ - centroidMat;

  Eigen::Matrix3d M = Pic.transpose() * Piw;

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::MatrixXd U = svd.matrixU();
  Eigen::MatrixXd V = svd.matrixV();

  R = U*V.transpose();
  double detR = R.determinant();

  if (detR < 0){
    R(2,0) = -R(2,0);
    R(2,1) = -R(2,1);
    R(2,2) = -R(2,2);    
  }

  t = P0c - R*P0w;  

}


double EPnPEigen::computeRt(const Eigen::MatrixXd&U, double betas[4], Eigen::Matrix3d& R, Eigen::Vector3d& t){
  computeControlPointsUnderCameraCoord(U, betas);
  computeReferencePointsUnderCameraCoord();
  solveForSign();

  estimateRt(R, t);
  return reprojectionError(R, t);
}

double EPnPEigen::reprojectionError(Eigen::Matrix3d& R, Eigen::Vector3d& t)
{
  double sum2 = 0.0;

  for(int i = 0; i < reference_points_count_; i++) {

    
    double Xc = R.row(0).dot(reference_3d_points_.row(i)) + t(0);                                // pws经外参（R和t）变为pcs
    double Yc = R.row(1).dot(reference_3d_points_.row(i)) + t(1);
    double inv_Zc = 1.0 / R.row(2).dot(reference_3d_points_.row(i)) + t(2);
    double ue = uc_ + fu_ * Xc * inv_Zc;                               // pcs经内参变为uv
    double ve = vc_ + fv_ * Yc * inv_Zc;
    double u = reference_2d_points_(i,0), v = reference_2d_points_(i,1);

    sum2 += sqrt( (u - ue) * (u - ue) + (v - ve) * (v - ve) );       // 计算与真实2d点之间的误差
  }

  return sum2 / reference_points_count_;
}


void EPnPEigen::computePose(Eigen::Matrix3d &R, Eigen::Vector3d &t){
  // 选择4个控制点（质点+3个主轴方向的单位向量）
  chooseControlPoints();
  // 根据4个控制点计算所有3d空间点的阿尔法系数
  // pi = ai_0*cws0 + ai_1*cws1 + ai_2*cws2 + ai_3*cws3
  computeBaryCentricCoordinates();	
  // Mx=0 M为2n*12的矩阵，x为4个控制点在相机坐标系下的xyz值，为12*1矩阵
  Eigen::MatrixXd M(2*reference_points_count_, 12);
  M = Eigen::MatrixXd::Zero(2*reference_points_count_, 12);
  calculateM(M);
  
  Eigen::MatrixXd MtM = M.transpose() * M;
  // 也可以对M^tM进行SVD分解，得到的特征向量既是M的右奇异向量
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(MtM);
  Eigen::VectorXd eigenval = es.eigenvalues();
  Eigen::MatrixXd eigvector = es.eigenvectors();

  // L * betas = rho
  // rho 一直为 6*1 的矩阵，记录着4个控制点之间各自的距离; N=4;L为6*10的矩阵
  Eigen::MatrixXd L6x10 = Eigen::MatrixXd::Zero(6, 10);
  Eigen::VectorXd rho(6, 1);
  computeL6x10(eigvector, L6x10);
  computeRho(rho);

  /* 算出4个betas的初值 */
  double betas[4][4],rep_errors[4];;
  
  Eigen::Matrix3d RR = Eigen::Matrix3d::Zero(3,3);
  Eigen::Vector3d tt = Eigen::Vector3d::Zero(3,1);
  vector<Eigen::Matrix3d> Rs(4,RR);
  vector<Eigen::Vector3d> ts(4,tt);
  
  // （N=4: B1 B2 B3 B4）
  findBetasApprox1(L6x10, rho, betas[1]);
  printf("betas[1] = %f, %f, %f, %f\n", betas[1][0], betas[1][1], betas[1][2], betas[1][3]);
  doGaussNewtonOptimization(eigvector, L6x10, betas[1]);
  printf("betas[1] = %f, %f, %f, %f\n", betas[1][0], betas[1][1], betas[1][2], betas[1][3]);
  rep_errors[1]=computeRt(eigvector, betas[1], Rs[1], ts[1]);

  // （N=2: B1 B2 (B3B4为0)）
  findBetasApprox2(L6x10, rho, betas[2]);
  printf("betas[2] = %f, %f, %f, %f\n", betas[2][0], betas[2][1], betas[2][2], betas[2][3]);
  doGaussNewtonOptimization(eigvector, L6x10, betas[2]);
  rep_errors[2]=computeRt(eigvector, betas[2], Rs[2], ts[2]);
  printf("betas[2] = %f, %f, %f, %f\n", betas[2][0], betas[2][1], betas[2][2], betas[2][3]);
  
  // （N=3: B1 B2 B3 (B4为0)）
  findBetasApprox3(L6x10, rho, betas[3]);
  printf("betas[3] = %f, %f, %f, %f\n", betas[3][0], betas[3][1], betas[3][2], betas[3][3]);
  doGaussNewtonOptimization(eigvector, L6x10, betas[3]);
  rep_errors[3]=computeRt(eigvector, betas[3], Rs[3], ts[3]);
  printf("betas[3] = %f, %f, %f, %f\n", betas[3][0], betas[3][1], betas[3][2], betas[3][3]);
  

  // 选出重投影误差最小的R和t
  int N = 1;
  if (rep_errors[2] < rep_errors[1]) { N = 2;}                              
  if (rep_errors[3] < rep_errors[N]) { N = 3;}
  cout << "rep_errors[1]: " << rep_errors[1] << endl;
  cout << "rep_errors[2]: " << rep_errors[2] << endl;
  cout << "rep_errors[3]: " << rep_errors[3] << endl;
  cout << "choose N= " << N << endl;

  R = Rs[N];
  t = ts[N];
}