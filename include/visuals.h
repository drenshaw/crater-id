#pragma once
#include <opencv2/viz/types.hpp>

#include "conics.h"
#include "camera.h"
#include "quadrics.h"

template <typename Iterator>
    void advance_all (Iterator & iterator) {
        ++iterator;
    }
template <typename Iterator, typename ... Iterators>
    void advance_all (Iterator & iterator, Iterators& ... iterators) {
        ++iterator;
        advance_all(iterators...);
    } 
template <typename Function, typename Iterator, typename ... Iterators>
    Function zip (Function func, Iterator begin, 
            Iterator end, 
            Iterators ... iterators)
    {
        for(;begin != end; ++begin, advance_all(iterators...))
            func(*begin, *(iterators)... );
        //could also make this a tuple
        return func;
    }


namespace viz {

const std::vector<cv::Scalar> CV_colors = {
  cv::viz::Color::amethyst(),
  cv::viz::Color::apricot(),
  cv::viz::Color::azure(),
  // cv::viz::Color::black(),
  cv::viz::Color::bluberry(),
  cv::viz::Color::blue(),
  cv::viz::Color::brown(),
  cv::viz::Color::celestial_blue(),
  cv::viz::Color::chartreuse(),
  cv::viz::Color::cherry(),
  cv::viz::Color::cyan(),
  cv::viz::Color::gold(),
  // cv::viz::Color::gray(), // used for the Moon color
  cv::viz::Color::green(),
  cv::viz::Color::indigo(),
  cv::viz::Color::lime(),
  cv::viz::Color::magenta(),
  cv::viz::Color::maroon(),
  cv::viz::Color::mlab(),
  cv::viz::Color::navy(),
  // cv::viz::Color::not_set(),
  cv::viz::Color::olive(),
  cv::viz::Color::orange(),
  cv::viz::Color::orange_red(),
  cv::viz::Color::pink(),
  cv::viz::Color::purple(),
  cv::viz::Color::raspberry(),
  cv::viz::Color::red(),
  cv::viz::Color::rose(),
  cv::viz::Color::silver(),
  cv::viz::Color::teal(),
  cv::viz::Color::turquoise(),
  cv::viz::Color::violet(),
  cv::viz::Color::white(),
  cv::viz::Color::yellow(),
};

void drawEllipse( cv::Mat& image, const Conic& conic, const cv::Scalar& color=cv::Scalar(0, 255, 255));
void drawEllipse( const Conic& conic, const cv::Scalar& color=cv::Scalar(0, 255, 255));
void drawEllipses(cv::Mat& image, const std::vector<Conic>& conics, const std::vector<cv::Scalar>& colors={});
void drawEllipses(cv::Mat& image, const Camera& camera, const std::vector<Quadric>& quadrics, const std::vector<cv::Scalar>& colors={});
void drawEllipses(const std::vector<Conic>& conics, const std::vector<cv::Scalar>& colors={});
void drawLine(cv::Mat& image, const Eigen::Vector3d& my_line, const std::string& text, const cv::Scalar& my_color);
void drawLine(const Eigen::Vector3d& my_line, const std::string& text, const cv::Scalar& my_color);
void drawLines(cv::Mat& image, const std::vector<Eigen::Vector3d>& lines, const std::vector<std::string>& text, const std::vector<cv::Scalar>& colors);
void drawLines(const std::vector<Eigen::Vector3d>& lines, const std::string& text, const std::vector<cv::Scalar>& colors);
void drawPoint(cv::Mat& image, const Eigen::Vector2d& point, const cv::Scalar& color);
void drawPoints(cv::Mat& image, const std::vector<Eigen::Vector2d>& points, const std::vector<cv::Scalar>& colors);
void drawPoints(cv::Mat& image, const std::vector<std::vector<Eigen::Vector2d> >& points, const std::vector<cv::Scalar>& colors);
void getSlopeInterceptFromStandard(const Eigen::Vector3d& my_line, double& slope, double& intercept);
bool getEndpointsFromLine(const cv::Mat& image, const Eigen::Vector3d& my_line, cv::Point2l& start_pt, cv::Point2l& end_pt);
void get3dAxes( const Camera& cam, 
                Eigen::Vector2d& origin, Eigen::Vector2d& x_axis, 
                Eigen::Vector2d& y_axis, Eigen::Vector2d& z_axis);
void draw3dAxes(cv::Mat& image, const Camera& cam);
void interactiveZoom(cv::Mat& image);
void drawNoisyPoints( const Camera& cam, const std::vector<Conic>& conics,
                      const std::vector<std::vector<Eigen::Vector2d> >& noisy_pts);


// /* VTK */
// #include <vtk-9.1/vtkActor.h>
// void Other();
// void Sphere();
// void Cone();
// void Ellipsoid();
// void Cylinder();
// void HyperboloidOneSheet();
// void HyperboloidTwoSheets();
// void HyperbolicParaboloid();
// void EllipticParaboloid();

// void PlotFunction(vtkQuadric* quadric, double value);
} // namespace
