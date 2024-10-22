#pragma once

#include <eigen3/Eigen/Dense>
#include <opencv2/viz/types.hpp>

#define CAMERA_INTRINSIC_PARAM 5

class Camera {
  public:
    Camera();
    Camera( const double dx,
            const double dy,
            const double up,
            const double vp,
            const double skew,
            const cv::Size2i& image_size,
            const Eigen::Isometry3d& transform);
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
    Camera( const Eigen::Matrix3d& intrinsic, 
            const cv::Size2i& image_size,
            const Eigen::Isometry3d& attitude=Eigen::Isometry3d::Identity());
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
    Eigen::Isometry3d getTransformation() const;
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
    Eigen::Isometry3d getHomogeneousExtrinsicMatrix() const ;
    Eigen::MatrixXd getExtrinsicMatrix() const;
    Eigen::Affine3d getHomogeneousProjectionMatrix() const ;
    Eigen::MatrixXd getProjectionMatrix() const;
    cv::Mat getBlankCameraImage() const;
    void resetImage(cv::Mat& image) const;
    std::array<double, CAMERA_INTRINSIC_PARAM> getIntrinsicParams() const ;
    void getIntrinsicParams(std::array<double, CAMERA_INTRINSIC_PARAM>& params) const;
    double getImageWidth() const;
    double getImageHeight() const;
    Eigen::Vector2d getImageMidpoint() const;
    void getImageMidpoint(Eigen::Vector2d& img_midpoint) const;
    cv::Size2i getImageSize() const;
    double getFovX() const;
    double getFovXDeg() const;
    double getFovY() const;
    double getFovYDeg() const;
    bool isInFrontOfCamera(const Eigen::Vector3d& pt) const;
    bool isInCameraFrame(const Eigen::Vector3d& pt) const ;
    bool isInCameraFrame(const Eigen::Vector3d& pt, Eigen::Vector2d& pt_pxl) const;
    bool isInCameraFrame(const Eigen::Vector2d& pt_pxl) const;
    void world2Camera(const Eigen::Vector3d& pt, Eigen::Vector3d& pt_cam) const ;
    void world2Pixel(const Eigen::Vector3d& pt, Eigen::Vector2d& pt_pxl) const ;

    void world2Pixel(const Eigen::Vector3d& pt, cv::Point2d& pt_pxl) const;
    Eigen::Vector3d world2Camera(const Eigen::Vector3d& pt) const;
    Eigen::Vector2d world2Pixel(const Eigen::Vector3d& pt) const ;
    void setAttitude(const Eigen::Quaterniond& orientation);
    void setAttitude(const Eigen::Matrix3d& orientation);
    void setAttitude(const Eigen::AngleAxisd& orientation);
    void setLocation(const Eigen::Vector3d& location);
    void setPosition(const Eigen::Vector3d& location);
    void setIntrinsicMatrix(const Eigen::Matrix3d& intrinsic);
    void setExtrinsicMatrix(const Eigen::Matrix3d& extrinsic);
    // void setProjectionMatrix(const Eigen::Matrix3d& projection_matrix);
    // Camera(const Eigen::Matrix4d& projection_matrix);
    void moveCamera(const Eigen::Isometry3d& transform);
    void rotate(const Eigen::Quaterniond& rotation);
    void rotate(const Eigen::Matrix3d& rotation);
    void rotate(const Eigen::AngleAxisd& rotation);
    void move(const Eigen::Vector3d& translation);
    void move(const Eigen::Translation3d& translation);
    void moveX(const double x_offset);
    void moveY(const double y_offset);
    void moveZ(const double z_offset);
    void moveRelative(const Eigen::Vector3d& translation);
    void moveRelative(const Eigen::Translation3d& translation);
    void rotateRelative(const Eigen::Quaterniond& rotation);
    void pointTo(const Eigen::Vector3d& point, const Eigen::Vector3d& up_axis);
    void pointTo(const double lat, const double lon, const Eigen::Vector3d& up_axis);
    void moveTo(const Eigen::Vector3d& point);
    void moveTo(const Eigen::Translation3d& translation);
    void resetCameraState();
    // Eigen::Transform<double, 3, Eigen::Isometry>
    Eigen::Matrix3d projectQuadric(const Eigen::Matrix4d& quadric) const;
    Eigen::Matrix3d projectQuadricToLocus(const Eigen::Matrix4d& quadric_locus) const;
    Eigen::Matrix3d getImagePlaneLocus(const Eigen::Matrix3d& image_locus) const;

  private:
    double dx_;
    double dy_;
    double skew_;
    double up_;
    double vp_;
    cv::Size2i image_size_;
    Eigen::Matrix3d intrinsic_matrix_;
    // Eigen::Vector3d position_;
    // Eigen::Quaterniond attitude_;
    bool isInPixelArray(const Eigen::Vector2d& pt_uv) const;
    Eigen::Isometry3d state_;

    friend std::ostream& operator<<(std::ostream& os, const Camera& cam);
    friend std::ostream& operator<<(std::ostream& os, const Camera* cam);
};


/********************************************/
void getCameraIntrinsicParams(const Eigen::Matrix3d& cam_intrinsic_mtx,
                              std::array<double, CAMERA_INTRINSIC_PARAM>& params);
Eigen::Matrix3d getCameraInverseIntrinsicMatrix(const Eigen::Matrix3d& camera_intrinsic_mtx);
double getCameraFovX(const Eigen::Matrix3d& camera_intrinsic_mtx, const cv::Size2i& image_size);
double getCameraFovXDeg(const Eigen::Matrix3d &camera_intrinsic_mtx, const cv::Size2i &image_size);
double getCameraFovY(const Eigen::Matrix3d& camera_intrinsic_mtx, const cv::Size2i& image_size);
double getCameraFovYDeg(const Eigen::Matrix3d &camera_intrinsic_mtx, const cv::Size2i &image_size);
Eigen::Matrix3d getAttitudeTransformBetweenPoints(const Eigen::Vector3d& camera_position, const Eigen::Vector3d& desired_location);
Eigen::Matrix3d lookAt( const Eigen::Vector3d& camera_position, 
                        const Eigen::Vector3d& desired_location, 
                        const Eigen::Vector3d& up_vector);
Eigen::Matrix3d lookAt( const Eigen::Vector3d& camera_position, 
                        const Eigen::Vector3d& desired_location);
Eigen::Matrix3d lookAt( const Eigen::Isometry3d& transform, 
                        const Eigen::Vector3d& desired_location, 
                        const Eigen::Vector3d& up_vector);                        

bool isInImage(const Eigen::Vector2d& pt_uv, const cv::Size2i image_size);
bool isInImage(const cv::Point& pt_uv, const cv::Size2i image_size);
bool isInImage(const cv::Point& pt_uv, const cv::MatSize image_size);
