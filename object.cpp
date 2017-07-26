#include "object.h"
#include <iostream>
Object::Object()
{
    model = vtkSmartPointer<vtkPolyData>::New();
    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    actor = vtkSmartPointer<vtkActor>::New();
    actor->GetProperty()->SetInterpolationToGouraud();
    actor->SetMapper(mapper);
}
void Object::SetInput(vtkSmartPointer<vtkPolyData> t_model)
{
    model->DeepCopy(t_model);
    //model = t_model;
    std::cout<<"Points:"<<model->GetNumberOfPoints()<<endl;
    mapper->SetInputData(model);
    mapper->Update();
}

