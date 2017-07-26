#ifndef GENESYNDATA_H
#define GENESYNDATA_H

#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkLight.h>
#include <vtkLightActor.h>
#include <vtkWindowToImageFilter.h>

class Genesyndata
{
public:
    Genesyndata();
    ~Genesyndata();
    vtkSmartPointer<vtkRenderer> getrenderer() {return m_renderer;}
    void addlight();

    void rendermodel(vtkSmartPointer<vtkPolyData> t_model);

    void loadcamerapath(vtkSmartPointer<vtkPolyData> t_model);

    void get_z_values(vtkRenderWindow*t_renderwin);

    int get_num_cams(){return m_numcams;}


    void updatecamera();
    void setconstantlight(int value);
    void setlinearlight(int value);
    void setqudraticlight(int value);
    void setambient(int value);
    void setdiffuse(int value);
    void setspecular(int value);
    void randomcampos(float x,float y,float z, float elevation,float azimuth);


private:
     vtkSmartPointer<vtkSphereSource> m_sphereSource;

    vtkSmartPointer<vtkLight> m_light;
    vtkSmartPointer<vtkPolyDataMapper> m_Mapper;
    vtkSmartPointer<vtkActor> m_Actor;
    vtkSmartPointer<vtkRenderer> m_renderer;
    vtkSmartPointer<vtkPolyData> m_camerapath;
    vtkSmartPointer<vtkLightActor> m_lightActor;

    int counter;
    int totalcount;
    int m_numcams;


};

#endif // GENESYNDATA_H
