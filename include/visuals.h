#ifndef VISUALS_H
#define VISUALS_H

#include <opencv2/core/core.hpp> 
#include <opencv2/imgproc.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include <opencv2/viz/types.hpp>
#include <eigen3/Eigen/Dense>

#include "conics.h"

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
  cv::viz::Color::gray(),
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

void plotline(cv::Mat& image, const Eigen::Vector3d& my_line, const cv::Scalar my_color);
void plotEllipse( cv::Mat& image, const Conic& conic, const cv::Scalar& color=cv::Scalar(0, 255, 255));
void plotEllipse( const Conic& conic, const cv::Scalar& color=cv::Scalar(0, 255, 255));
void plotEllipses(cv::Mat& image, const std::vector<Conic> conics, const std::vector<cv::Scalar> colors={});
void plotEllipses(const std::vector<Conic> conics, const std::vector<cv::Scalar> colors={});

/* VTK */
// #include <cstdlib>
// #include <vtkActor.h>
// #include <vtkCamera.h>
// #include <vtkContourFilter.h>
// #include <vtkImageData.h>
// #include <vtkNamedColors.h>
// #include <vtkNew.h>
// #include <vtkOutlineFilter.h>
// #include <vtkPolyDataMapper.h>
// #include <vtkProperty.h>
// #include <vtkQuadric.h>
// #include <vtkRenderWindow.h>
// #include <vtkRenderWindowInteractor.h>
// #include <vtkRenderer.h>
// #include <vtkSampleFunction.h>
// #include <cstdlib>
// #include <vtkActor.h>
// #include <vtkCamera.h>
// #include <vtkContourFilter.h>
// #include <vtkImageData.h>
// #include <vtkNamedColors.h>
// #include <vtkNew.h>
// #include <vtkOutlineFilter.h>
// #include <vtkPolyDataMapper.h>
// #include <vtkProperty.h>
// #include <vtkQuadric.h>
// #include <vtkRenderWindow.h>
// #include <vtkRenderWindowInteractor.h>
// #include <vtkRenderer.h>
// #include <vtkSampleFunction.h>


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


#endif