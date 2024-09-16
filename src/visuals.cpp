#include "visuals.h"

namespace viz {

void plotEllipse(cv::Mat& image, const Conic& conic, const cv::Scalar& color) {
  if (!image.data) { 
    std::cout << "Could not open or "
              << "find the image\n"; 
    return; 
  }
  cv::Point center;
  conic.GetCenter(center);
  // Drawing the ellipse 
  cv::ellipse(image, center, 
              conic.GetSize(), conic.GetAngle(), 
              0, 360, 
              color, -1, cv::LINE_AA); 

  int font = cv::FONT_HERSHEY_SIMPLEX;
  float font_scale = 1; // scale factor from base size
  int thickness = 1; //in pixels
  cv::Scalar text_neg(255, 255, 255);
  cv::Scalar text_color = text_neg - color;
  cv::Point offset(10, 5);
  cv::putText(image, std::to_string(conic.GetID()), center+offset, font, font_scale,
              text_color, thickness, true);
}

void plotEllipse(const Conic& conic, const cv::Scalar& color) {
  cv::Mat image(500, 500, CV_8UC3, 
                cv::Scalar(255, 255, 255)); 
  plotEllipse(image, conic, color);
  // Showing image inside a window 
  cv::imshow("Output", image); 
  cv::waitKey(0); 
}

void plotEllipses(cv::Mat& image, const std::vector<Conic> conics, const std::vector<cv::Scalar> colors) {
  assert(conics.size() <= colors.size());
  for (auto conic_it = conics.begin(); conic_it != conics.end(); ++conic_it) {
    int index = std::distance(conics.begin(), conic_it);
    plotEllipse(image, *conic_it, CV_colors.at(index));
  }
}

void plotEllipses(const std::vector<Conic> conics, const std::vector<cv::Scalar> colors) {
  cv::Mat image(500, 500, CV_8UC3, 
                cv::Scalar(255, 255, 255)); 
  plotEllipses(image, conics, colors);
  // Showing image inside a window 
  cv::imshow("Output", image); 
  cv::waitKey(0); 
}

void plotline(cv::Mat& image, const Eigen::Vector3d& my_line, const cv::Scalar my_color) {
  cv::Point start_pt, end_pt;
  cv::Size image_size(image.rows,image.cols);
  double slope, intercept, x0, y0, x1, y1;
  if(std::abs(my_line(1)) > 1e-5) {
    x0 = 0;
    x1 = image.cols;
    slope =     -my_line(0)/my_line(1);
    intercept = -my_line(2)/my_line(1);
    y0 = slope*x0 + intercept;
    y1 = slope*x1 + intercept;
  }
  else {
    std::cerr << "Line is (nearly) vertical:\n" << my_line << std::endl;
    slope =     -my_line(1)/my_line(0);
    intercept = -my_line(2)/my_line(0);
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
    return;
  std::cerr << "Line eq: l = (" << slope << ")x + " << intercept << " | " << start_pt << ", " << end_pt << std::endl;
  }
  int thickness = 1; //in pixels
  int lineType = cv::LINE_AA;
  int shift = 0;
  cv::line(image, start_pt, end_pt, my_color, thickness, lineType, shift);
  
}

void plotline(const Eigen::Vector3d& line, const cv::Scalar& color) {
  cv::Mat image(500, 500, CV_8UC3, 
                cv::Scalar(255, 255, 255)); 
  plotline(image, line, color);
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
