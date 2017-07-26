#ifndef RENDERMANAGER_H
#define RENDERMANAGER_H
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>

#include <vtkSmartPointer.h>
#include <vtkLight.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>

class RenderManager
{
private:
    vtkSmartPointer<vtkRenderer> render;
    vtkSmartPointer<vtkLight> light;
public:
    RenderManager();
    ~RenderManager(){}
    void renderModel(vtkSmartPointer<vtkActor> t_actor);

    void addlight();
    void setconstantlight(int value);
    void setlinearlight(int value);
    void setqudraticlight(int value);
    void setambient(int value);
    void setdiffuse(int value);
    void setspecular(int value);
    vtkSmartPointer<vtkLight> GetLight(){return light;}

    vtkSmartPointer<vtkRenderer> GetRender(){return render;}
};

#endif // RENDERMANAGER_H
