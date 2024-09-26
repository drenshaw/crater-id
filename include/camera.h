#ifndef CAMERA_H
#define CAMERA_H

#include <eigen3/Eigen/Dense>
#include <opencv2/viz/types.hpp>

#define CAMERA_INTRINSIC_PARAM 5

class Camera {
  public:
    Camera( const double dx,
            const double dy,
            const double up,
            const double vp,
            const double skew,
            const cv::Size2i& image_size,
            const Eigen::Quaterniond& attitude=Eigen::Quaterniond::Identity(),
            const Eigen::Vector3d& position=Eigen::Vector3d::Zero());
    Camera( const Eigen::Matrix3d& intrinsic, 
                const cv::Size2i& image_size,
            const Eigen::Quaterniond& attitude=Eigen::Quaterniond::Identity(),
            const Eigen::Vector3d& position=Eigen::Vector3d::Zero());
    Camera( const std::vector<double>& intrinsics,
                const cv::Size2i& image_size,
            const Eigen::Quaterniond& attitude=Eigen::Quaterniond::Identity(),
            const Eigen::Vector3d& position=Eigen::Vector3d::Zero());
    Camera( const Eigen::Matrix3d& intrinsic, 
                const cv::Size2i& image_size,
            const Eigen::Matrix4d& extrinsic);
    Camera( const Eigen::Matrix3d& intrinsic, 
                const cv::Size2i& image_size,
            const Eigen::Matrix3d& attitude,
            const Eigen::Vector3d& translation);
    Camera( const Eigen::Matrix3d& intrinsic,
                const cv::Size2i& image_size,
            const Eigen::Vector3d& position=Eigen::Vector3d::Zero(),
            const Eigen::Quaterniond& attitude=Eigen::Quaterniond::Identity());
    // Camera( const Eigen::Matrix3d& intrinsic,
                // const cv::Size2i& image_size,
    //         const Eigen::Vector3d& position=Eigen::Vector3d::Zero(),
    //         const Eigen::Matrix3d& attitude=Eigen::Matrix3d::Identity());
    Eigen::Vector3d getPosition() const;
    Eigen::Vector3d getLocation() const;
    Eigen::Quaterniond getAttitude() const;
    void getAttitude(Eigen::Quaterniond& attitude) const;
    void getAttitude(Eigen::Matrix3d& attitude) const;
    void getAttitude(Eigen::Vector3d& attitude) const;
    Eigen::Quaterniond getAttitudeQuaternion() const ;
    Eigen::Matrix3d getAttitudeMatrix() const ;
    Eigen::Vector3d getAttitudeEuler(Eigen::Index a0, Eigen::Index a1, Eigen::Index a2) const ;
    Eigen::Matrix3d getIntrinsicMatrix() const ;
    Eigen::Matrix3d getInverseIntrinsicMatrix() const;
    Eigen::Matrix4d getHomogeneousExtrinsicMatrix() const ;
    Eigen::MatrixXd getExtrinsicMatrix() const;
    Eigen::Matrix4d getHomogeneousProjectionMatrix() const ;
    Eigen::MatrixXd getProjectionMatrix() const;
    std::array<double, CAMERA_INTRINSIC_PARAM> getIntrinsicParams() const ;
    void getIntrinsicParams(std::array<double, CAMERA_INTRINSIC_PARAM>& params) const;
    Eigen::Vector3d getPointWrtCameraFrame(const Eigen::Vector3d& pt) const;
    bool isInFrontOfCamera(const Eigen::Vector3d& pt) const;
    bool pointUVinImage(const Eigen::Vector2d& pt_uv) const;
    bool isInCameraFrame(const Eigen::Vector3d& pt) const ;
    bool isInCameraFrame(const Eigen::Vector3d& pt, Eigen::Vector2d& pt_pxl) const;
    void projectXYZtoImage(const Eigen::Vector3d& pt, Eigen::Vector2d& pt_pxl) const ;
    Eigen::Vector2d projectXYZtoImage(const Eigen::Vector3d& pt) const ;
    bool setAttitude(const Eigen::Quaterniond& orientation);
    bool setAttitude(const Eigen::Matrix3d& orientation);
    bool setAttitude(const Eigen::Vector3d& orientation);
    bool setIntrinsicMatrix(const Eigen::Matrix3d& intrinsic);
    bool setExtrinsicMatrix(const Eigen::Matrix3d& extrinsic);
    // bool setProjectionMatrix(const Eigen::Matrix3d& projection_matrix);
    // Camera(const Eigen::Matrix4d& projection_matrix);

  private:
    double dx_;
    double dy_;
    double skew_;
    double up_;
    double vp_;
    cv::Size2i image_size_;
    Eigen::Matrix3d intrinsic_matrix_;
    Eigen::Vector3d position_;
    Eigen::Quaterniond attitude_;
    //  Eigen::Matrix4d extrinsic_matrix_;
    //  Eigen::Matrix4d extrinsic_matrix4_;
    // Eigen::Matrix4d projection_matrix_;
};

#endif