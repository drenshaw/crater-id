#ifndef VISUALS_H
#define VISUALS_H

/* VTK */
#include <cstdlib>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkContourFilter.h>
#include <vtkImageData.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkOutlineFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkQuadric.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSampleFunction.h>


namespace VIS {
void Other();
void Sphere();
void Cone();
void Ellipsoid();
void Cylinder();
void HyperboloidOneSheet();
void HyperboloidTwoSheets();
void HyperbolicParaboloid();
void EllipticParaboloid();

void PlotFunction(vtkQuadric* quadric, double value);
} // namespace


#endif