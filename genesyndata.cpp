#include "genesyndata.h"
#include <iostream>
#include <vtkCenterOfMass.h>
#include <vtkProperty.h>
#include <vtkPolyDataNormals.h>
#include <vtkLightCollection.h>

#include <vtkBMPWriter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkImageShiftScale.h>
#include <vtkPNGWriter.h>
#include <fstream>
#include <vtkMatrix4x4.h>
#include <vtkFloatArray.h>
#include <string>
#include <sstream>





Genesyndata::Genesyndata()
{

//    m_sphereSource =
//      vtkSmartPointer<vtkSphereSource>::New();
    //Mapper
    m_Mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
//      m_Mapper->SetInputConnection(m_sphereSource->GetOutputPort());

    //Actor
    m_Actor =
        vtkSmartPointer<vtkActor>::New();
    m_Actor->SetMapper(m_Mapper);
    m_Actor->GetProperty()->SetInterpolationToGouraud();

    //light
    m_light = vtkSmartPointer<vtkLight>::New();
    m_light->SetLightTypeToHeadlight();
    m_light->SetAttenuationValues(0,0.2,0);
    m_light->SetPositional(true); // required for vtkLightActor below
    m_light->SetConeAngle(180);
    m_light->SetDiffuseColor(100,100,100);
    m_light->SetAmbientColor(100,100,100);
    m_light->SetSpecularColor(100,100,100);
    m_light->SetLightTypeToHeadlight();

    //Renderer
    m_renderer =
        vtkSmartPointer<vtkRenderer>::New();

    m_renderer->AddActor(m_Actor);


    counter=0;
    totalcount=0;

}

//==================
//Render 3D model
//==================
void Genesyndata::rendermodel(vtkSmartPointer<vtkPolyData> t_model)
{


    // calculate normals
     vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
     normalGenerator->SetInputData(t_model);
     normalGenerator->ComputePointNormalsOn();
     normalGenerator->ComputeCellNormalsOn();
     normalGenerator->Update();

     t_model = normalGenerator->GetOutput();


    m_Mapper->SetInputData(t_model);

    m_renderer->ResetCamera();


}


//==================
//Store camera path
//==================
void Genesyndata::loadcamerapath(vtkSmartPointer<vtkPolyData> t_model)
{

    //m_renderer->ResetCamera();
    m_camerapath=t_model;
    m_numcams=t_model->GetNumberOfPoints();
}


//==================
//Set current camera position and look at
//==================
void Genesyndata::updatecamera()
{
    //Get inital camera position

    double p[3];
    m_camerapath->GetPoint(counter,p);


    vtkSmartPointer<vtkCamera> t_camera=m_renderer->GetActiveCamera();
    t_camera->SetPosition(p[0],p[1],p[2]);

    m_camerapath->GetPoint(counter+1,p);  

    t_camera->SetFocalPoint(p[0],p[1],p[2]);
    counter++;
}


void Genesyndata::randomcampos(float x,float y,float z, float elevation,float azimuth)
{
    double p[3];
    m_camerapath->GetPoint(counter-1,p);
    vtkSmartPointer<vtkCamera> t_camera=m_renderer->GetActiveCamera();


    t_camera->SetPosition(p[0]+x,p[1]+y,p[2]+z);

    m_camerapath->GetPoint(counter,p);

    t_camera->SetFocalPoint(p[0]+x,p[1]+y,p[2]+z);

    t_camera->Elevation(elevation);
    t_camera->Azimuth(azimuth);

    totalcount+=1;
}


//===============
// Obtain the depth map observed at current camera view
//===============

void Genesyndata::get_z_values(vtkRenderWindow*t_renderwin)
{


    vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
      vtkSmartPointer<vtkWindowToImageFilter>::New();
    windowToImageFilter->SetInput(t_renderwin);
    windowToImageFilter->SetMagnification(1); //set the resolution of the output image (3 times the current resolution of vtk render windowtotalcount)
    windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
    windowToImageFilter->ReadFrontBufferOff(); // read from the back buffer
    windowToImageFilter->Update();

    vtkSmartPointer<vtkPNGWriter> writer =
      vtkSmartPointer<vtkPNGWriter>::New();

    std::ostringstream oss;

    int t_count=totalcount;
    oss<<"frame"<<setfill('0')<<setw(5)<<t_count<<".png";

    std::string filename=oss.str();

    writer->SetFileName( ("./images/"+filename).c_str());
    writer->SetInputConnection(windowToImageFilter->GetOutputPort());
    writer->Write();

    vtkSmartPointer<vtkCamera> t_camera=m_renderer->GetActiveCamera();


    double* cliprange=t_camera->GetClippingRange();

    double nearclip=cliprange[0];
    double farclip=cliprange[1];


    int * t_sz=t_renderwin->GetSize();
    int t_width=t_sz[0];
    int t_height=t_sz[1];


    float *zvalues = new float[t_width * t_height];

    t_renderwin->GetZbufferData(0,0,t_width-1,t_height-1,zvalues);

    //rescale z from [0,1] to [-1,1]
    for (int i=0;i<t_width*t_height;i++){
        zvalues[i]=zvalues[i]*2-1;
        zvalues[i]=2*nearclip*farclip/(farclip+nearclip-zvalues[i]*(farclip-nearclip));
    }

    std::string depthfilename="./depth/"+filename+".bin";
    ofstream myFile(depthfilename.c_str(),ios::out | ios::binary);
    myFile.write((char*)zvalues,sizeof(float)*t_width*t_height);

    delete zvalues;

}


void Genesyndata::addlight()
{
    m_renderer->AddLight(m_light);
}

void Genesyndata::setconstantlight(int value)
{
    m_light->SetAttenuationValues(float(value)/100,0,0);

}

void Genesyndata::setlinearlight(int value)
{
    m_light->SetAttenuationValues(0,float(value)/100,0);

}

void Genesyndata::setqudraticlight(int value)
{
    m_light->SetAttenuationValues(0,0,float(value)/100);

}

void Genesyndata::setambient(int value)
{
    m_light->SetAmbientColor(value,value,value);

}

void Genesyndata::setdiffuse(int value)
{
    m_light->SetDiffuseColor(value,value,value);

}

void Genesyndata::setspecular(int value)
{
    m_light->SetSpecularColor(value,value,value);

}



Genesyndata::~Genesyndata(){}

