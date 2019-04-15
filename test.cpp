#include <iostream>
#include <opencv2/opencv.hpp>

#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/calib3d/calib3d.hpp>


#include "eigenEPNP.h"

using namespace cv;
using namespace std;
using namespace Eigen;

void find_feature_matches ( const Mat& img_1, const Mat& img_2,
                            std::vector<KeyPoint>& keypoints_1,
                            std::vector<KeyPoint>& keypoints_2,
                            std::vector< DMatch >& matches );

Point2d pixel2cam ( const Point2d& p, const Mat& K );

int main(int argc, char** argv){
	cout << "Hello EPNP" << endl;
    Mat image;
    //image = imread("/Users/apple/project/EPNP/1.png", 1);
    //image = imread("../data/1.png", 1);
    //namedWindow("Display Image", WINDOW_AUTOSIZE);
    //imshow("Display Image", image);
    //waitKey(0);
    
    //-- 读图像2张
    Mat img_1 = imread ( "../data/1.png", CV_LOAD_IMAGE_COLOR );
    Mat img_2 = imread ( "../data/2.png", CV_LOAD_IMAGE_COLOR );

 
    // 初始化
    vector<KeyPoint> keypoints_1, keypoints_2;
    vector<DMatch> matches;
    find_feature_matches ( img_1, img_2, keypoints_1, keypoints_2, matches );
    cout<<"一共找到了"<<matches.size() <<"组匹配点"<<endl;

    // 建立3D点
    Mat d1 = imread ( "../data/1_depth.png", CV_LOAD_IMAGE_UNCHANGED );       // 深度图为16位无符号数，单通道图像
    Mat K = ( Mat_<double> ( 3,3 ) << 520.9, 0, 325.1, 0, 521.0, 249.7, 0, 0, 1 );
    vector<Point3f> pts_3d;
    vector<Point2f> pts_2d;
    for ( DMatch m:matches )
    {
        ushort d = d1.ptr<unsigned short> (int ( keypoints_1[m.queryIdx].pt.y )) [ int ( keypoints_1[m.queryIdx].pt.x ) ];
        if ( d == 0 )   // bad depth
            continue;
        float dd = d/5000.0;
        Point2d p1 = pixel2cam ( keypoints_1[m.queryIdx].pt, K );
        pts_3d.push_back ( Point3f ( p1.x*dd, p1.y*dd, dd ) );
        pts_2d.push_back ( keypoints_2[m.trainIdx].pt );
    }

    cout<<"3d-2d pairs: "<<pts_3d.size() <<endl;

    // 用opencv 自带方法来解PNP
    cout << "---using opencv solvePNP---" << endl;
    Mat r, t;
    solvePnP ( pts_3d, pts_2d, K, Mat(), r, t, false,cv::SOLVEPNP_EPNP); // 调用OpenCV 的 PnP 求解，可选择EPNP，DLS等方法
    Mat R;
    cv::Rodrigues ( r, R ); // r为旋转向量形式，用Rodrigues公式转换为矩阵
    
    cout << "R1: " << R << endl;;
    cout << "t1:" << t << endl;
    
    // 用eigen 求EPNP
    cout << "---using egien to solve EPNP---" << endl;
    
    // 输入转化
	  Eigen::MatrixXd points3d(pts_3d.size(),3);
    Eigen::MatrixXd points2d(pts_2d.size(),2);
   
    for(int i = 0; i< pts_3d.size(); i++){
      points3d(i,0) = pts_3d[i].x;
      points3d(i,1) = pts_3d[i].y;
      points3d(i,2) = pts_3d[i].z;
    }
   

    for(int i = 0; i< pts_2d.size(); i++){
      points2d(i,0) = pts_2d[i].x;
      points2d(i,1) = pts_2d[i].y; 
    }
    
    Eigen::Matrix3d KK;
    KK << 520.9, 0, 325.1, 0, 521.0, 249.7, 0, 0, 1;
    
    //EPnPEigen(Eigen::MatrixXd& points3d, Eigen::MatrixXd& points2d, Eigen::Matrix3d& K);
    EPnPEigen test(points3d, points2d, KK);
    test.computePose();
   
     
    

    
     //    p1w-c1w
    //PWO = ...
    //    pn2-c1w
    // PWO*PWO^T的特征值和特征向量
     
     //Sum aijfuxcj + aij*(uc-ui)zcj = 0;
    //Sum aijfvycj + aij*(vc-vi)zcj = 0;
     //Mx = 0


    //根据旋转矩阵求坐标旋转角
    // theta_x = atan2(R.at<double>(2, 1), R.at<double>(2, 2));
    // theta_y = atan2(-R.at<double>(2, 0),
    // sqrt(R.at<double>(2, 1)*R.at<double>(2, 1) + R.at<double>(2, 2)*R.at<double>(2, 2)));
    // theta_z = atan2(R.at<double>(1, 0), R.at<double>(0, 0));

    //将弧度转化为角度
    // theta_x = theta_x * (180 / PI);
    //theta_y = theta_y * (180 / PI);
    // theta_z = theta_z * (180 / PI);



    return 0;
}


 
 

void find_feature_matches ( const Mat& img_1, const Mat& img_2,
                            std::vector<KeyPoint>& keypoints_1,
                            std::vector<KeyPoint>& keypoints_2,
                            std::vector< DMatch >& matches )
{
    //-- 初始化
    Mat descriptors_1, descriptors_2;
    // used in OpenCV3
    Ptr<FeatureDetector> detector = ORB::create();
    Ptr<DescriptorExtractor> descriptor = ORB::create();
    // use this if you are in OpenCV2
    // Ptr<FeatureDetector> detector = FeatureDetector::create ( "ORB" );
    // Ptr<DescriptorExtractor> descriptor = DescriptorExtractor::create ( "ORB" );
    Ptr<DescriptorMatcher> matcher  = DescriptorMatcher::create ( "BruteForce-Hamming" );
    //-- 第一步:检测 Oriented FAST 角点位置
    detector->detect ( img_1,keypoints_1 );
    detector->detect ( img_2,keypoints_2 );

    //-- 第二步:根据角点位置计算 BRIEF 描述子
    descriptor->compute ( img_1, keypoints_1, descriptors_1 );
    descriptor->compute ( img_2, keypoints_2, descriptors_2 );

    //-- 第三步:对两幅图像中的BRIEF描述子进行匹配，使用 Hamming 距离
    vector<DMatch> match;
    // BFMatcher matcher ( NORM_HAMMING );
    matcher->match ( descriptors_1, descriptors_2, match );

    //-- 第四步:匹配点对筛选
    double min_dist=10000, max_dist=0;

    //找出所有匹配之间的最小距离和最大距离, 即是最相似的和最不相似的两组点之间的距离
    for ( int i = 0; i < descriptors_1.rows; i++ )
    {
        double dist = match[i].distance;
        if ( dist < min_dist ) min_dist = dist;
        if ( dist > max_dist ) max_dist = dist;
    }

    printf ( "-- Max dist : %f \n", max_dist );
    printf ( "-- Min dist : %f \n", min_dist );

    //当描述子之间的距离大于两倍的最小距离时,即认为匹配有误.但有时候最小距离会非常小,设置一个经验值30作为下限.
    for ( int i = 0; i < descriptors_1.rows; i++ )
    {
        if ( match[i].distance <= max ( 2*min_dist, 30.0 ) )
        {
            matches.push_back ( match[i] );
        }
    }

    //show
    //Mat img_goodmatch;
    //drawMatches(img_1, keypoints_1, img_2, keypoints_2, matches, img_goodmatch);
    //imshow("good matches", img_goodmatch);
    //waitKey(0);
}

Point2d pixel2cam ( const Point2d& p, const Mat& K )
{
    return Point2d
           (
               ( p.x - K.at<double> ( 0,2 ) ) / K.at<double> ( 0,0 ),
               ( p.y - K.at<double> ( 1,2 ) ) / K.at<double> ( 1,1 )
           );
}
