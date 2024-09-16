#ifndef CAMERA_H
#define CAMERA_H

#include <eigen3/Eigen/Dense>

class Camera {
    public:
        Camera(const double dx,
               const double dy,
               const double up,
               const double vp,
               const double skew=0);
        // Camera(const double dx,
        //        const double dy,
        //        const double up,
        //        const double vp);
        Camera(const std::vector<double>& intrinsics);
        Camera(const Eigen::Matrix3d& intrinsic);
        Camera(const Eigen::Matrix3d& intrinsic, 
               const Eigen::MatrixXd& extrinsic);
        Camera(const Eigen::Matrix3d& intrinsic, 
               const Eigen::Matrix3d& rotation,
               const Eigen::Vector3d& translation);
        Camera(const Eigen::Matrix3d& intrinsic, 
               const Eigen::Quaterniond& rotation,
               const Eigen::Vector3d& translation);
        Eigen::Vector3d getPosition();
        Eigen::Vector3d getLocation();
        Eigen::Matrix3d getAttitude();
        void get_attitude(Eigen::Quaterniond& attitude);
        Eigen::Matrix3d getIntrinsicMatrix();
        Eigen::MatrixXd getExtrinsicMatrix();
        Eigen::MatrixXd getProjectionMatrix();
        bool getIntrinsicParams(std::vector<double>&);
        Eigen::Vector2d projectXYZtoImage();
        bool setAttitude(const Eigen::Matrix3d& orientation);
        bool setIntrinsicMatrix(const Eigen::Matrix3d& intrinsic);
        bool setExtrinsicMatrix(const Eigen::Matrix3d& extrinsic);
        Eigen::Matrix3d getAttitudeQuaternion();
        // bool setProjectionMatrix(const Eigen::Matrix3d& projection_matrix);
        // Camera(const Eigen::MatrixXd& projection_matrix);

    private:
        double dx_;
        double dy_;
        double skew_;
        double up_;
        double vp_;
        Eigen::Matrix3d intrinsic_matrix_;
        Eigen::Vector3d position_;
        Eigen::Matrix3d attitude_;
       //  Eigen::MatrixXd extrinsic_matrix_;
       //  Eigen::Matrix4d extrinsic_matrix4_;
        Eigen::MatrixXd projection_matrix_;
        void makeCamera(const double dx_,
                        const double dy_,
                        const double up_,
                        const double vp_,
                        const double skew_=0);
        void makeCamera(const std::vector<double>& intrinsics);
        void makeCamera(const Eigen::Matrix3d& intrinsic);
        void makeCamera(const Eigen::Matrix3d& intrinsic, 
                        const Eigen::MatrixXd& extrinsic);
        void makeCamera(const Eigen::Matrix3d& intrinsic, 
                        const Eigen::Matrix3d& rotation,
                        const Eigen::Vector3d& translation);
        void makeCamera(const Eigen::Matrix3d& intrinsic, 
                        const Eigen::Quaterniond& rotation,
                        const Eigen::Vector3d& translation);
       void makeCamera(const Eigen::Matrix3d& intrinsic, 
                   const Eigen::Quaterniond& rotation,
                   const Eigen::Vector3d& translation, bool yes);
        // void makeCamera(const Eigen::MatrixXd& projection_matrix);
};

#endif