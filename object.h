#ifndef OBJECT_H
#define OBJECT_H

#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkXMLDataReader.h>
#include <vtkPolyData.h>
#include <vtkLight.h>
#include <vtkProperty.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>

class Object
{
public:
    Object();
    ~Object(){}
    void SetInput(vtkSmartPointer<vtkPolyData> t_model);
    vtkSmartPointer<vtkPolyData> GetOutput(){return model;}
    vtkSmartPointer<vtkActor> GetActor(){return actor;}
protected:
    vtkSmartPointer<vtkPolyData> model;
    vtkSmartPointer<vtkPolyDataMapper> mapper;
    vtkSmartPointer<vtkActor> actor;
};

#endif // OBJECT_H
