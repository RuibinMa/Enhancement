#ifndef DATA_H
#define DATA_H

#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkXMLDataReader.h>
#include <vtkPolyData.h>
#include <vtkLight.h>


class Data
{
public:
    Data();
    ~Data();
    vtkSmartPointer<vtkRenderer> getrenderer() {return m_renderer;}
    void rendermodel();
    void addlight();
    void SetInput(vtkSmartPointer<vtkPolyData> model);
private:
    vtkSmartPointer<vtkPolyData> m_model;
    vtkSmartPointer<vtkPolyDataMapper> m_mapper;
    vtkSmartPointer<vtkActor> m_actor;

    vtkSmartPointer<vtkRenderer> m_renderer;
    vtkSmartPointer<vtkLight> m_light;
};

#endif // DATA_H
