#include "rendermanager.h"

RenderManager::RenderManager()
{
    render = vtkSmartPointer<vtkRenderer>::New();

    light = vtkSmartPointer<vtkLight>::New();
    light->SetLightTypeToHeadlight();
    light->SetIntensity(1);
    light->SetAttenuationValues(0, 0.5, 0);
    light->SetPositional(true);
    light->SetConeAngle(60);
    light->SetDiffuseColor(50, 50, 50);
    light->SetAmbientColor(30, 30, 30);
    light->SetSpecularColor(30, 30, 30);
    //light->SetColor(0.1,0.1,0.1);
    //light->SetLightTypeToSceneLight();
}

void RenderManager::renderModel(vtkSmartPointer<vtkActor> t_actor)
{
    render->AddActor(t_actor);
    render->ResetCamera();
}

void RenderManager::addlight()
{
    render->AddLight(light);
}
void RenderManager::setconstantlight(int value)
{
    light->SetAttenuationValues(float(value)/100,0,0);
}

void RenderManager::setlinearlight(int value)
{
    light->SetAttenuationValues(0,float(value)/100,0);
}

void RenderManager::setqudraticlight(int value)
{
    light->SetAttenuationValues(0,0,float(value)/100);
}

void RenderManager::setambient(int value)
{
    light->SetAmbientColor(value,value,value);
}

void RenderManager::setdiffuse(int value)
{
    light->SetDiffuseColor(value,value,value);
}

void RenderManager::setspecular(int value)
{
    light->SetSpecularColor(value,value,value);
}
