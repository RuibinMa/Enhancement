#include "data.h"
#include <iostream>
#include <vtkProperty.h>

Data::Data()
{
    m_model = vtkSmartPointer<vtkPolyData>::New();

    m_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();

    m_actor = vtkSmartPointer<vtkActor>::New();
    m_actor->SetMapper(m_mapper);
    m_actor->GetProperty()->SetInterpolationToGouraud();
    m_actor->GetProperty()->SetOpacity(0.5);

    m_light = vtkSmartPointer<vtkLight>::New();
    m_light->SetLightTypeToHeadlight();
    m_light->SetAttenuationValues(0, 0.01, 0);
    m_light->SetPositional(true);
    m_light->SetConeAngle(180);
    m_light->SetDiffuseColor(10,10,10);
    m_light->SetAmbientColor(10,10,10);
    m_light->SetSpecularColor(10,10,10);

    m_renderer = vtkSmartPointer<vtkRenderer>::New();
    m_renderer->AddActor(m_actor);
}
Data::~Data()
{
}
void Data::SetInput(vtkSmartPointer<vtkPolyData> model)
{
    m_model = model;
}

void Data::rendermodel()
{
    m_mapper->SetInputData(m_model);
    m_renderer->ResetCamera();
}
void Data::addlight()
{
    m_renderer->AddLight(m_light);
}
