#include "visuals.h"
#include "camera.h"
#include "opencv2/core/types.hpp"
#include "opencv2/viz/types.hpp"
#include <iostream>
#include <stdexcept>
#include <cmath>

#include <opencv2/core/core.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include <opencv2/imgproc.hpp> 
#include <eigen3/Eigen/Dense>

namespace viz {

void drawEllipse(cv::Mat& image, const Conic& conic, const cv::Scalar& color) {
  if (!image.data) { 
    std::cout << "Could not open or "
              << "find the image\n"; 
    return; 
  }
  cv::Size2i img_size(image.rows, image.cols);
  cv::Point ellipse_center;
  conic.getCenter(ellipse_center);
  if(!isInImage(ellipse_center, image.size)) {
    std::cout << "Ellipse is not in the image: " << ellipse_center << std::endl;
    return;
  }
  // Drawing the ellipse 
  cv::ellipse(image, ellipse_center, 
              conic.getSize(), conic.getAngleDeg(), 
              0, 360, 
              color, -1, cv::LINE_AA); 

  // int font = cv::FONT_HERSHEY_SIMPLEX;
  // float font_scale = 1; // scale factor from base size
  // int thickness = 3; //in pixels
  // cv::Scalar rect_color = cv::viz::Color::orange();
  // cv::Scalar text_color = cv::viz::Color::azure();
  // cv::Point ll_offset(-10, -35);
  // cv::Point ur_offset( 50,  20);
  // cv::rectangle (image, center+ll_offset, center+ur_offset, rect_color, cv::FILLED, cv::LINE_8, 0);
  // cv::putText(image, std::to_string(conic.getID()), center, font, font_scale,
  //             text_color, thickness, cv::LINE_AA, false);
}

void drawEllipse(const Conic& conic, const cv::Scalar& color) {
  cv::Mat image(1024, 1296, CV_8UC3, 
                cv::Scalar(50, 50, 50)); 
  drawEllipse(image, conic, color);
  // Showing image inside a window 
  cv::imshow("Ellipse", image); 
  cv::waitKey(0); 
}

void drawEllipses(cv::Mat& image, const std::vector<Conic>& conics, const std::vector<cv::Scalar>& colors) {
  assert(conics.size() <= colors.size());
  for (auto conic_it = conics.begin(); conic_it != conics.end(); ++conic_it) {
    int index = std::distance(conics.begin(), conic_it);
    drawEllipse(image, *conic_it, CV_colors.at(index));
  }
}

void drawEllipses(const std::vector<Conic>& conics, const std::vector<cv::Scalar>& colors) {
  cv::Mat image(1024, 1296, CV_8UC3, 
                cv::Scalar(50, 50, 50)); 
  drawEllipses(image, conics, colors);
  // Showing image inside a window 
  cv::imshow("Ellipses", image); 
  cv::waitKey(0); 
}

void drawLine(cv::Mat& image, const Eigen::Vector3d& line, const std::string& text, const cv::Scalar& color) {
  cv::Point2l start_pt, end_pt;
  cv::Size image_size(image.cols,image.rows);
  if(!getEndpointsFromLine(image, line, start_pt, end_pt)) {
    return;
  }
  int thickness = 2; //in pixels
  int lineType = cv::LINE_AA;
  int shift = 0;
  cv::line(image, start_pt, end_pt, color, thickness, lineType, shift);
  int font = cv::FONT_HERSHEY_SIMPLEX;
  float font_scale = 1; // scale factor from base size
  // int thickness = 3; //in pixels
  cv::Point2l center = (start_pt + end_pt) / 2;
  cv::Scalar rect_color = cv::viz::Color::orange();
  cv::Scalar text_color = cv::viz::Color::azure();
  cv::Point2l ll_offset(-10, -35);
  cv::Point2l ur_offset( 50,  20);
  cv::rectangle (image, center+ll_offset, center+ur_offset, rect_color, cv::FILLED, cv::LINE_8, 0);
  cv::putText(image, text, center, font, font_scale,
              text_color, thickness, lineType, false);
}

void drawLine(const Eigen::Vector3d& line, const std::string text, const cv::Scalar& color) {
  cv::Mat image(500, 500, CV_8UC3, 
                cv::Scalar(255, 255, 255)); 
  drawLine(image, line, text, color);
  // Showing image inside a window 
  cv::imshow("Line", image); 
  cv::waitKey(0); 
}

void drawLines(cv::Mat& image, const std::vector<Eigen::Vector3d>& lines, const std::vector<std::string>& text, const std::vector<cv::Scalar>& colors) {
  assert(lines.size() <= colors.size());
  for (auto conic_it = lines.begin(); conic_it != lines.end(); ++conic_it) {
    int index = std::distance(lines.begin(), conic_it);
    drawLine(image, *conic_it, text.at(index), colors.at(index));
  }
}

void drawLines(const std::vector<Eigen::Vector3d>& lines, const std::vector<std::string>& text, const std::vector<cv::Scalar>& colors) {
  cv::Mat image(500, 500, CV_8UC3, 
                cv::Scalar(255, 255, 255)); 
  drawLines(image, lines, text, colors);
  // Showing image inside a window 
  cv::imshow("Lines", image); 
  cv::waitKey(0); 
}

void drawPoint(cv::Mat& image, const Eigen::Vector2d& point, const cv::Scalar& color) {

    cv::Point2d pt_cv = {point(0), point(1)};
    cv::drawMarker(image, pt_cv, color);
}

void drawPoints(cv::Mat& image, const std::vector<Eigen::Vector2d>& points, const std::vector<cv::Scalar>& colors) {
  assert(points.size() <= colors.size());
  for (auto pt_it = points.begin(); pt_it != points.end(); ++pt_it) {
    int index = std::distance(points.begin(), pt_it);
    drawPoint(image, *pt_it, viz::CV_colors.at(index));
  }
}

void getSlopeInterceptFromStandard(const Eigen::Vector3d& my_line, double& slope, double& intercept) {
  double A = my_line(0);
  double B = my_line(1);
  double C = my_line(2);
  if(A == 0 && B == 0) {
    slope     = std::numeric_limits<double>::signaling_NaN();
    intercept = std::numeric_limits<double>::signaling_NaN();
    throw std::runtime_error("Division by zero");
  }
  assert(A != 0 || B != 0);
  if(std::abs(B) > 1e-5) {
    slope =     -A/B;
    intercept = -C/B;
  }
  else {
    std::cerr << "Line is vertical/nearly vertical: "
              << A << "x + " << B << "y + " << C << " = 0\n";
    slope =     std::numeric_limits<double>::infinity();
    intercept = -C/A;
  }
}

bool getEndpointsFromLine(const cv::Mat& image, const Eigen::Vector3d& my_line, cv::Point2l& start_pt, cv::Point2l& end_pt) {
  double slope, intercept, x0, y0, x1, y1;
  try {
    getSlopeInterceptFromStandard(my_line, slope, intercept);
  }
  catch (const std::runtime_error& e) {
    std::cerr << "Exception: " << e.what() << std::endl;
    return false;
  }
  cv::Size image_size(image.cols,image.rows);
  if(!std::isinf(slope)) {
    x0 = 0;
    x1 = image.cols;
    y0 = slope*x0 + intercept;
    y1 = slope*x1 + intercept;
  }
  else {
    x0 = intercept;
    x1 = intercept;
    y0 = 0;
    y1 = image.rows;
  }
  start_pt.x = x0;
  start_pt.y = y0;
  end_pt.x = x1;
  end_pt.y = y1;
  bool is_visible = cv::clipLine(image_size, start_pt, end_pt);
  if (!is_visible) {
    std::cerr << "Line is not visible: l = (" << slope << ")x + " << intercept << std::endl << my_line << std::endl;
    return false;
  }
  return true;
}

void get3dAxes( const Camera& cam, 
                Eigen::Vector2d& origin, Eigen::Vector2d& x_axis, 
                Eigen::Vector2d& y_axis, Eigen::Vector2d& z_axis) {
  Eigen::Vector3d x, y, z, o;
  o = Eigen::Vector3d::Zero();
  x = R_MOON/2*Eigen::Vector3d::UnitX();
  y = R_MOON/2*Eigen::Vector3d::UnitY();
  z = R_MOON/2*Eigen::Vector3d::UnitZ();
  cam.world2Pixel(o, origin);
  cam.world2Pixel(x, x_axis);
  cam.world2Pixel(y, y_axis);
  cam.world2Pixel(z, z_axis);
}

void draw3dAxes(cv::Mat& image, const Camera& cam) {
  Eigen::Vector2d origin, x_axis, y_axis, z_axis;
  viz::get3dAxes(cam, origin, x_axis, y_axis, z_axis);
  cv::Point2d o(origin(0), origin(1));
  cv::Point2d x(x_axis(0), x_axis(1));
  cv::Point2d y(y_axis(0), y_axis(1));
  cv::Point2d z(z_axis(0), z_axis(1));
  
  cv::circle(image, o, 8, cv::viz::Color::black());
  cv::circle(image, o, 5, cv::viz::Color::white());
  cv::arrowedLine(image, o, x, cv::viz::Color::red());
  cv::arrowedLine(image, o, y, cv::viz::Color::green());
  cv::arrowedLine(image, o, z, cv::viz::Color::blue());
}

// void onMouse(int event, int x, int y, int, void*) { if (event == cv::EVENT_LBUTTONDOWN) { dragging = true; prev_pt = cv::Point(x, y); } else if (event == cv::EVENT_LBUTTONUP) { dragging = false; } else if (event == cv::EVENT_MOUSEMOVE && dragging) { cv::Point delta = cv::Point(x, y) - prev_pt; roi.x -= delta.x / scale; roi.y -= delta.y / scale; prev_pt = cv::Point(x, y); cv::resize(img, temp_img, cv::Size(), scale, scale); cv::imshow("Interactive Zoom", temp_img(roi)); } else if (event == cv::EVENT_MOUSEWHEEL) { if (cv::getMouseWheelDelta(event) > 0) { scale *= 1.1; } else { scale *= 0.9; } roi = cv::Rect2f(roi.x, roi.y, img.cols / scale, img.rows / scale); cv::resize(img, temp_img, cv::Size(), scale, scale); cv::imshow("Interactive Zoom", temp_img(roi)); } }

void interactiveZoom(cv::Mat& image) {

    namedWindow("Image: Press 'x' to close", cv::WINDOW_NORMAL);
    imshow("Image: Press 'x' to close", image);

    int zoomFactor = 1;
    int x = 0;
    int y = 0;

    while (true) {
        int key = cv::waitKey(1);

        if (key == 'x') {
            break;
        } else if (key == 'q') {
            zoomFactor += 1;
        } else if (key == 'e' && zoomFactor > 1) {
            zoomFactor -= 1;
        } else if (key == 'a') {
            x -= 10; 
        } else if (key == 'd') {
            x += 10;
        } else if (key == 'w') {
            y -= 10;
        } else if (key == 's') {
            y += 10;
        }

        // Make sure the zoom center is within the image bounds
        x = std::max(0, std::min(x, image.cols - 1));
        y = std::max(0, std::min(y, image.rows - 1));

        // Calculate the zoom region
        int zoomWidth = image.cols / zoomFactor;
        int zoomHeight = image.rows / zoomFactor;
        int x1 = std::max(0, x - zoomWidth / 2);
        int y1 = std::max(0, y - zoomHeight / 2);
        int x2 = std::min(image.cols, x1 + zoomWidth);
        int y2 = std::min(image.rows, y1 + zoomHeight);

        // Extract and resize the zoom region
        cv::Mat zoomedImage = image(cv::Rect(x1, y1, x2 - x1, y2 - y1));
        resize(zoomedImage, zoomedImage, cv::Size(image.cols, image.rows), 0, 0, cv::INTER_LINEAR);

        imshow("Image: Press 'x' to close", zoomedImage);
    }

    cv::destroyAllWindows();
}

// void Sphere()
// {
//   // create the quadric function definition
//   vtkNew<vtkQuadric> quadric;
//   quadric->SetCoefficients(1, 1, 1, 0, 0, 0, 0, 0, 0, 0);
// void Sphere()
// {
//   // create the quadric function definition
//   vtkNew<vtkQuadric> quadric;
//   quadric->SetCoefficients(1, 1, 1, 0, 0, 0, 0, 0, 0, 0);

//   // F(x,y,z) = a0*x^2 + a1*y^2 + a2*z^2 + a3*x*y + a4*y*z + a5*x*z + a6*x +
//   // a7*y + a8*z + a9 F(x,y,z) = 1*x^2 + 1*y^2 + 1*z^2

//   PlotFunction(quadric, 1.0);
// }

// void EllipticParaboloid()
// {
//   // create the quadric function definition
//   vtkNew<vtkQuadric> quadric;
//   quadric->SetCoefficients(1, 1, 0, 0, 0, 0, 0, 0, -1, 0);

//   // F(x,y,z) = a0*x^2 + a1*y^2 + a2*z^2 + a3*x*y + a4*y*z + a5*x*z + a6*x +
//   // a7*y + a8*z + a9 F(x,y,z) = 1*x^2 + 1*y^2

//   PlotFunction(quadric, 10.0);
// }

// void HyperbolicParaboloid()
// {
//   // create the quadric function definition
//   vtkNew<vtkQuadric> quadric;
//   quadric->SetCoefficients(1, -1, 0, 0, 0, 0, 0, 0, 0, 0);

//   // F(x,y,z) = a0*x^2 + a1*y^2 + a2*z^2 + a3*x*y + a4*y*z + a5*x*z + a6*x +
//   // a7*y + a8*z + a9 F(x,y,z) = 1*x^2 - 1*y^2

//   PlotFunction(quadric, 10.0);
// }

// void Cylinder()
// {
//   // create the quadric function definition
//   vtkNew<vtkQuadric> quadric;
//   quadric->SetCoefficients(1, 1, 0, 0, 0, 0, 0, 0, 0, 0);

//   // F(x,y,z) = a0*x^2 + a1*y^2 + a2*z^2 + a3*x*y + a4*y*z + a5*x*z + a6*x +
//   // a7*y + a8*z + a9 F(x,y,z) = 1*x^2 + 1*y^2

//   PlotFunction(quadric, 1.0);
// }

// void HyperboloidOneSheet()
// {
//   // create the quadric function definition
//   vtkNew<vtkQuadric> quadric;
//   quadric->SetCoefficients(1, 1, -1, 0, 0, 0, 0, 0, 0, 0);

//   // F(x,y,z) = a0*x^2 + a1*y^2 + a2*z^2 + a3*x*y + a4*y*z + a5*x*z + a6*x +
//   // a7*y + a8*z + a9 F(x,y,z) = 1*x^2 + 1*y^2

//   PlotFunction(quadric, 1.0);
// }

// void HyperboloidTwoSheets()
// {
//   // create the quadric function definition
//   vtkNew<vtkQuadric> quadric;
//   quadric->SetCoefficients(1, 1, -1, 0, 0, 0, 0, 0, 0, 0);

//   // F(x,y,z) = a0*x^2 + a1*y^2 + a2*z^2 + a3*x*y + a4*y*z + a5*x*z + a6*x +
//   // a7*y + a8*z + a9 F(x,y,z) = 1*x^2 + 1*y^2

//   PlotFunction(quadric, -1.0);
// }

// void Ellipsoid()
// {
//   // create the quadric function definition
//   vtkNew<vtkQuadric> quadric;
//   quadric->SetCoefficients(1, 1, 2, 0, 0, 0, 0, 0, 0, 0);

//   // F(x,y,z) = a0*x^2 + a1*y^2 + a2*z^2 + a3*x*y + a4*y*z + a5*x*z + a6*x +
//   // a7*y + a8*z + a9 F(x,y,z) = 1*x^2 + 1*y^2 + 1*z^2

//   PlotFunction(quadric, -1.0);
// }

// void Cone()
// {
//   // create the quadric function definition
//   vtkNew<vtkQuadric> quadric;
//   quadric->SetCoefficients(1, 1, -1, 0, 0, 0, 0, 0, 0, 0);

//   // F(x,y,z) = a0*x^2 + a1*y^2 + a2*z^2 + a3*x*y + a4*y*z + a5*x*z + a6*x +
//   // a7*y + a8*z + a9 F(x,y,z) = 1*x^2 + 1*y^2 - 1*z^2
//   PlotFunction(quadric, 0.0);
// }

// void Other()
// {
//   // create the quadric function definition
//   vtkNew<vtkQuadric> quadric;
//   quadric->SetCoefficients(.5, 1, .2, 0, 0.1, 0, 0, .2, 0, 0);

//   // F(x,y,z) = a0*x^2 + a1*y^2 + a2*z^2 + a3*x*y + a4*y*z + a5*x*z + a6*x +
//   // a7*y + a8*z + a9 F(x,y,z) = 0.5*x^2 + 1*y^2 + 0.2*z^2 + 0*x*y + 0.1*y*z +
//   // 0*x*z + 0*x + 0.2*y + 0*z + 0
//   PlotFunction(quadric, 1.0);
// }

// void PlotFunction(vtkQuadric* quadric, double value)
// {

//   vtkNew<vtkNamedColors> colors;

//   // sample the quadric function
//   vtkNew<vtkSampleFunction> sample;
//   sample->SetSampleDimensions(50, 50, 50);
//   sample->SetImplicitFunction(quadric);
//   // double xmin = 0, xmax=1, ymin=0, ymax=1, zmin=0, zmax=1;
//   double bounds[6]{-100, 110, -100, 100, -100, 100};
//   sample->SetModelBounds(bounds);

//   // Create five surfaces F(x,y,z) = constant between range specified
//   /*
//   vtkContourFilter *contours = vtkContourFilter::New();
//   contours->SetInput(sample->GetOutput());
//   contours->GenerateValues(5, 0.0, 1.2);
//   */

//   // create the 0 isosurface
//   vtkNew<vtkContourFilter> contours;
//   contours->SetInputConnection(sample->GetOutputPort());
//   contours->GenerateValues(1, value, value);

//   // map the contours to graphical primitives
//   vtkNew<vtkPolyDataMapper> contourMapper;
//   contourMapper->SetInputConnection(contours->GetOutputPort());
//   contourMapper->SetScalarRange(0.0, 1.2);

//   // create an actor for the contours
//   vtkNew<vtkActor> contourActor;
//   contourActor->SetMapper(contourMapper);

//   // -- create a box around the function to indicate the sampling volume --

//   // create outline
//   vtkNew<vtkOutlineFilter> outline;
//   outline->SetInputConnection(sample->GetOutputPort());

//   // map it to graphics primitives
//   vtkNew<vtkPolyDataMapper> outlineMapper;
//   outlineMapper->SetInputConnection(outline->GetOutputPort());

//   // create an actor for it
//   vtkNew<vtkActor> outlineActor;
//   outlineActor->SetMapper(outlineMapper);
//   outlineActor->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());

//   // setup the window
//   vtkNew<vtkRenderer> ren1;
//   vtkNew<vtkRenderWindow> renWin;
//   renWin->AddRenderer(ren1);
//   renWin->SetWindowName("DisplayQuadricSurfaces");

//   vtkNew<vtkRenderWindowInteractor> iren;
//   iren->SetRenderWindow(renWin);

//   // add the actors to the scene
//   ren1->AddActor(contourActor);
//   ren1->AddActor(outlineActor);
//   ren1->SetBackground(colors->GetColor3d("AliceBlue").GetData());

//   // render and interact
//   renWin->Render();
//   ren1->GetActiveCamera()->Azimuth(-55);
//   ren1->GetActiveCamera()->Elevation(15);
//   iren->Start();
//   std::cout << "What?!\n";
// }

} // namespace
