#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <unistd.h>
#include <cstring>
#include <string>

namespace
{
void CallbackFunction (vtkObject* caller, long unsigned int eventId,
                        void* clientData, void* callData );
void CallbackFunction_inverse (vtkObject* caller, long unsigned int eventId,
                        void* clientData, void* callData );
}
namespace
{
void CallbackFunction (vtkObject* caller,
                       long unsigned int vtkNotUsed(eventId),
                       void* clientData,
                       void* vtkNotUsed(callData) )
{
  vtkImageTracerWidget* tracerWidget = static_cast<vtkImageTracerWidget*>(caller);
  MainWindow* t_window = static_cast<MainWindow*>(clientData);

  vtkSmartPointer<vtkPolyData> path =
    vtkSmartPointer<vtkPolyData>::New();

  tracerWidget->GetPath(path);

  if(path->GetNumberOfPoints() > 0)
  {
      double p[3];
      path->GetPoint(0, p);
      vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
      pointLocator->SetDataSet(t_window->GetData(0)->GetOutput());
      pointLocator->BuildLocator();
      int id = pointLocator->FindClosestPoint(p);
      std::cout<<id<<" :"<<p[0]<<" "<<p[1]<<" "<<p[2]<<endl;

      double p1[3], p2[3];
      t_window->GetData(0)->GetOutput()->GetPoint(id, p1);
      t_window->GetData(1)->GetOutput()->GetPoint(id, p2);
      vtkSmartPointer<vtkPoints> points1 = vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkPoints> points2 = vtkSmartPointer<vtkPoints>::New();
      points1->DeepCopy(t_window->GetTracerMark(0)->GetOutput()->GetPoints());
      points2->DeepCopy(t_window->GetTracerMark(1)->GetOutput()->GetPoints());
      points1->InsertNextPoint(p1);
      points2->InsertNextPoint(p2);

      t_window->GetTracerMark(0)->SetColor(0,0,1);
      t_window->GetTracerMark(1)->SetColor(0,0,1);
      t_window->GetTracerMark(0)->InputPoints(points1);
      t_window->GetTracerMark(1)->InputPoints(points2);
      t_window->GetRenderManager(0)->GetRender()->AddActor(t_window->GetTracerMark(0)->GetActor());
      //t_window->GetRenderManager(0)->renderModel(t_window->GetTracerMark(0)->GetActor());
      t_window->GetRenderManager(1)->GetRender()->AddActor(t_window->GetTracerMark(1)->GetActor());
      //t_window->GetRenderManager(1)->renderModel(t_window->GetTracerMark(1)->GetActor());

      QVTKWidget* left = t_window->findChild<QVTKWidget*>("left");
      left->GetRenderWindow()->Render();
      QVTKWidget* right = t_window->findChild<QVTKWidget*>("right");
      right->GetRenderWindow()->Render();
  }
}
void CallbackFunction_inverse(vtkObject* caller,
                       long unsigned int vtkNotUsed(eventId),
                       void* clientData,
                       void* vtkNotUsed(callData) )
{
  vtkImageTracerWidget* tracerWidget = static_cast<vtkImageTracerWidget*>(caller);
  MainWindow* t_window = static_cast<MainWindow*>(clientData);

  vtkSmartPointer<vtkPolyData> path =
    vtkSmartPointer<vtkPolyData>::New();

  tracerWidget->GetPath(path);

  if(path->GetNumberOfPoints() > 0)
  {
      double p[3];
      path->GetPoint(0, p);
      vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
      pointLocator->SetDataSet(t_window->GetData(1)->GetOutput());
      pointLocator->BuildLocator();
      int id = pointLocator->FindClosestPoint(p);

      double p1[3], p2[3];
      t_window->GetData(0)->GetOutput()->GetPoint(id, p1);
      t_window->GetData(1)->GetOutput()->GetPoint(id, p2);
      vtkSmartPointer<vtkPoints> points1 = vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkPoints> points2 = vtkSmartPointer<vtkPoints>::New();
      points1->DeepCopy(t_window->GetTracerMark(2)->GetOutput()->GetPoints());
      points2->DeepCopy(t_window->GetTracerMark(3)->GetOutput()->GetPoints());
      points1->InsertNextPoint(p1);
      points2->InsertNextPoint(p2);

      t_window->GetTracerMark(2)->SetColor(0,1,0);
      t_window->GetTracerMark(3)->SetColor(0,1,0);
      t_window->GetTracerMark(2)->InputPoints(points1);
      t_window->GetTracerMark(3)->InputPoints(points2);
      t_window->GetRenderManager(0)->GetRender()->AddActor(t_window->GetTracerMark(2)->GetActor());
      //t_window->GetRenderManager(0)->renderModel(t_window->GetTracerMark(2)->GetActor());
      t_window->GetRenderManager(1)->GetRender()->AddActor(t_window->GetTracerMark(3)->GetActor());
      //t_window->GetRenderManager(1)->renderModel(t_window->GetTracerMark(3)->GetActor());

      QVTKWidget* left = t_window->findChild<QVTKWidget*>("left");
      left->GetRenderWindow()->Render();
      QVTKWidget* right = t_window->findChild<QVTKWidget*>("right");
      right->GetRenderWindow()->Render();
  }
}
} // end anonymous namespace

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    m_filemanager = new FileManager;
    m_centerline = new Centerline;
    m_colon = new Colon;
    m_colon_new = new Colon;
    m_rendermanager = new RenderManager;
    m_rendermanager_right = new RenderManager;

    tracer = vtkSmartPointer<vtkImageTracerWidget>::New();
    tracer_inverse = vtkSmartPointer<vtkImageTracerWidget>::New();
    m_tracermark_1 = new TracerMarker;
    m_tracermark_2 = new TracerMarker;
    m_tracermark_inverse_1 = new TracerMarker;
    m_tracermark_inverse_2 = new TracerMarker;

    m_oldcenterline = vtkSmartPointer<vtkPoints>::New();
    m_newcenterline = vtkSmartPointer<vtkPoints>::New();
    m_oldplanenormals = vtkSmartPointer<vtkDoubleArray>::New();
    m_newplanenormals = vtkSmartPointer<vtkDoubleArray>::New();
    m_sweepingplane_l = new SweepingPlane;
    m_sweepingplane_r = new SweepingPlane;

    QVTKWidget* left = this->findChild<QVTKWidget*>("left");
    left->GetRenderWindow()->AddRenderer(m_rendermanager->GetRender());
    QVTKWidget* right = this->findChild<QVTKWidget*>("right");
    right->GetRenderWindow()->AddRenderer(m_rendermanager_right->GetRender());

    QCheckBox* checkbox_l = this->findChild<QCheckBox*>("tracer");
    checkbox_l->setDisabled(true);
    QCheckBox* checkbox_r = this->findChild<QCheckBox*>("tracer_r");
    checkbox_r->setDisabled(true);
    twoWindowReady = 0;

    QScrollBar* scrollbar = this->findChild<QScrollBar*>("horizontalScrollBar");
    scrollbar->setMaximum(10000);
    scrollbar->setDisabled(true);
}

MainWindow::~MainWindow()
{
    delete ui;
    delete m_filemanager;
    delete m_centerline;
    delete m_colon;
    delete m_colon_new;
    delete m_rendermanager;
    delete m_rendermanager_right;
    delete m_tracermark_1;
    delete m_tracermark_2;
    delete m_tracermark_inverse_1;
    delete m_tracermark_inverse_2;
    delete m_sweepingplane_l;
    delete m_sweepingplane_r;
}
void MainWindow::addlight()
{
    m_rendermanager->GetRender()->RemoveAllLights();
    m_rendermanager->addlight();
}

// load colon surface
void MainWindow::on_actionNew_file_triggered()
{
    //Open file dialog
    QString filePath = QFileDialog::getOpenFileName(
                this, tr("Open File"), "",
                tr("3Dmodels (*.stl *.vtp *.ply *.obj *.off)"));
    if(filePath.isEmpty()) return;

    m_filemanager->LoadNewFile(filePath);
    m_colon->Object::SetInput(m_filemanager->getfile());

    m_rendermanager->renderModel(m_colon->GetActor());
    QVTKWidget* left = this->findChild<QVTKWidget*>("left");

    tracer->GetLineProperty()->SetLineWidth(5);
    tracer->SetInteractor(left->GetRenderWindow()->GetInteractor());
    tracer->SetViewProp(m_colon->GetActor());
    tracer->AutoCloseOn();

    left->GetRenderWindow()->Render();
    vtkSmartPointer<vtkCallbackCommand> callback =
            vtkSmartPointer<vtkCallbackCommand>::New();
    callback->SetCallback(CallbackFunction);
    callback->SetClientData(this);
    tracer->AddObserver(vtkCommand::EndInteractionEvent, callback);

    QVTKWidget* right = this->findChild<QVTKWidget*>("right");
    right->GetRenderWindow()->Render();

    if(twoWindowReady == 0)
        twoWindowReady++;
    else
    {
        twoWindowReady = 2;
        QCheckBox* checkbox_l = this->findChild<QCheckBox*>("tracer");
        checkbox_l->setEnabled(true);
        QCheckBox* checkbox_r = this->findChild<QCheckBox*>("tracer_r");
        checkbox_r->setEnabled(true);
    }
}

// centerline and colon deformation are done in this function
void MainWindow::on_actionLoad_Centerline_triggered()
{
    //m_centerline->ConnectTwoContoursTest(m_rendermanager, m_filemanager);
    QString filePath = QFileDialog::getOpenFileName(
                this, tr("Open File"),"",
                tr("Centerline File (*.vtp)"));
    if(filePath.isEmpty()) return;
    m_filemanager->LoadNewFile(filePath);
    m_centerline->Object::SetInput(m_filemanager->getfile());
    // Uniform Sampling
    // m_centerline->UniformSample(400);
    // Gaussian Smoothing
    //m_centerline->SmoothCenterline(3);
    //m_filemanager->SaveFile(m_centerline->GetOutput(), "SmoothedCenterline.vtp");
    // Centerline-Driven Colon Deformation

    // get the origin point

    double p[3];
    m_centerline->GetOutput()->GetPoint(0, p);
    std::cout<<p[0]<<" "<<p[1]<<" "<<p[2]<<endl;


    vtkSmartPointer<vtkPolyData> newColonPoly = vtkSmartPointer<vtkPolyData>::New();
    //newColonPoly = m_centerline->EliminateTorsion(m_rendermanager,m_rendermanager_right, m_colon->GetOutput(), m_filemanager);
    m_colon_new->Object::SetInput(newColonPoly);

    m_filemanager->SaveFile(m_centerline->GetOutput(), "ModifiedCenterline.vtp");

    m_rendermanager->renderModel(m_centerline->GetActor());
    m_rendermanager_right->renderModel(m_colon_new->GetActor());
    QVTKWidget* right = this->findChild<QVTKWidget*>("right");

    tracer_inverse->GetLineProperty()->SetLineWidth(5);
    tracer_inverse->SetInteractor(right->GetRenderWindow()->GetInteractor());
    tracer_inverse->SetViewProp(m_colon_new->GetActor());

    right->GetRenderWindow()->Render();
    vtkSmartPointer<vtkCallbackCommand> callback =
            vtkSmartPointer<vtkCallbackCommand>::New();
    callback->SetCallback(CallbackFunction_inverse);
    callback->SetClientData(this);
    tracer_inverse->AddObserver(vtkCommand::EndInteractionEvent, callback);


    QVTKWidget* left = this->findChild<QVTKWidget*>("left");
    left->GetRenderWindow()->Render();
}

// Give two points on colon surface, find the geodesic path between them.
// This is used to find the slitting line on the surface
vtkSmartPointer<vtkPolyData> MainWindow::GeodesicPath(vtkSmartPointer<vtkPolyData> points)
{
    vtkSmartPointer<vtkPolyData> path = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> path_temp = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
    vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
    vtkSmartPointer<vtkDijkstraGraphGeodesicPath> dijkstra = vtkSmartPointer<vtkDijkstraGraphGeodesicPath>::New();

    dijkstra->SetInputData(m_colon->GetOutput());

    vtkSmartPointer<vtkPolyData> colon_temp = vtkSmartPointer<vtkPolyData>::New();
    colon_temp = m_colon->GetOutput();
    vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
    pointLocator->SetDataSet(m_colon->GetOutput());
    pointLocator->BuildLocator();
    for(vtkIdType i = 0; i<points->GetNumberOfPoints()-1; i++)
    {
        vtkIdType StartVertex, EndVertex;
        double closestPoint1[3], closestPoint2[3];
        double point1[3];
        points->GetPoint(i, point1);
        double point2[3];
        points->GetPoint(i+1, point2);

        StartVertex = pointLocator->FindClosestPoint(point1);
        EndVertex = pointLocator->FindClosestPoint(point2);
        colon_temp->GetPoint(StartVertex, closestPoint1);
        colon_temp->GetPoint(EndVertex, closestPoint2);

        //std::cout<<point1[0]<<" "<<point1[1]<<" "<<point1[2]<<std::endl;
        //std::cout<<point2[0]<<" "<<point2[1]<<" "<<point2[2]<<std::endl;
        //std::cout<<closestPoint1[0]<<" "<<closestPoint1[1]<<" "<<closestPoint1[2]<<std::endl;
        //std::cout<<closestPoint2[0]<<" "<<closestPoint2[1]<<" "<<closestPoint2[2]<<std::endl;
        //std::cout<<"from  "<<StartVertex<<" to  "<<EndVertex<<std::endl;

        dijkstra->SetStartVertex(StartVertex);
        dijkstra->SetEndVertex(EndVertex);
        dijkstra->Update();
        path_temp->ShallowCopy(dijkstra->GetOutput());

        appendFilter->RemoveAllInputs();
        appendFilter->AddInputData(path);
        appendFilter->AddInputData(path_temp);
        appendFilter->Update();

        cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
        cleanFilter->Update();

        path->DeepCopy(cleanFilter->GetOutput());
    }
    return path;
}

// Upsample the path by factor 10. Currently unused.
vtkSmartPointer<vtkPolyData> MainWindow::Upsampling(vtkSmartPointer<vtkPolyData> path)
{
    vtkSmartPointer<vtkCardinalSpline> spline = vtkSmartPointer<vtkCardinalSpline>::New();
    vtkSmartPointer<vtkSplineFilter> splineFilter = vtkSmartPointer<vtkSplineFilter>::New();
    splineFilter->SetInputData(path);
    splineFilter->SetNumberOfSubdivisions(path->GetNumberOfPoints() * 10);
    splineFilter->SetSpline(spline);
    splineFilter->Update();
    return splineFilter->GetOutput();
}

void MainWindow::on_action_Deform_Colon_triggered(bool test)
{
    test = false;
    double factor = 2.0, adjust = 0.75;
    //double r0 = 18.5793;
    double r0 = 1.01247;

    double k, b;
    k=0.5;
    k=1;

    b = r0*(1-k);
    double aver = 0;
    //double origin[3] = {667.6, 491.462, -213.051};
    double origin[3] = {0.00419734, 0.443927, -0.148548};
    double normal[3] = {0,0,1}, d1[3] = {1,0,0}, d2[3] = {0,1,0};
    double overlapangle = 45.0 / 180.0 * 3.1415926;
    double normal_r[3]; normal_r[0] = 0; normal_r[1] = -sin(overlapangle); normal_r[2] = cos(overlapangle);
    double normal_l[3]; normal_l[0] = 0; normal_l[1] = -sin(overlapangle); normal_l[2] = -cos(overlapangle);

    if(test)
    {
        // configuration

        r0 = 2.79012; adjust = 0.9; b = r0*(1-k);
        double f1[3] = {7.904, -4.674, 11.167}, f2[3] = {-5.976, 3.461, -3.272}, f3[3] = {-0.826, -1.501, 10.889};
        origin[0] = f2[0];
        origin[1] = f2[1];
        origin[2] = f2[2];

        /* // visualize the fiducial points
        vtkSmartPointer<vtkPoints> fiducials = vtkSmartPointer<vtkPoints>::New();
        fiducials->InsertNextPoint(f1);
        fiducials->InsertNextPoint(f2);
        fiducials->InsertNextPoint(f3);
        vtkSmartPointer<vtkPolyData> fiducialspoly = vtkSmartPointer<vtkPolyData>::New();
        fiducialspoly->SetPoints(fiducials);
        vtkSmartPointer<vtkVertexGlyphFilter> vertexFiducial = vtkSmartPointer<vtkVertexGlyphFilter>::New();
        vertexFiducial->SetInputData(fiducialspoly);
        vertexFiducial->Update();
        vtkSmartPointer<vtkPolyDataMapper> fiducialsMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        fiducialsMapper->SetInputConnection(vertexFiducial->GetOutputPort());
        fiducialsMapper->Update();
        vtkSmartPointer<vtkActor> fiducialActor = vtkSmartPointer<vtkActor>::New();
        fiducialActor->SetMapper(fiducialsMapper);
        fiducialActor->GetProperty()->SetPointSize(5);
        m_rendermanager->renderModel(fiducialActor);
        */

        double pd1[3], pd2[3], pd3[3], vpar[3], v2[3];
        double projection;
        vtkMath::Subtract(f1, f2, pd1);
        vtkMath::Normalize(pd1);
        vtkMath::Subtract(f3, f2, v2);
        projection = vtkMath::Dot(pd1, v2);
        vpar[0] = projection * pd1[0];
        vpar[1] = projection * pd1[1];
        vpar[2] = projection * pd1[2];
        vtkMath::Subtract(v2, vpar, pd2);

        //vtkMath::MultiplyScalar(pd2, -1);

        vtkMath::Normalize(pd2);
        vtkMath::Cross(pd1, pd2, pd3);
        vtkMath::Normalize(pd3);

        vtkSmartPointer<vtkPoints> newpoints = vtkSmartPointer<vtkPoints>::New();
        double pave[3] = {0,0,0};
        for(int i = 0; i < m_colon->GetOutput()->GetNumberOfPoints(); i++)
        {
            double p[3], v[3];
            m_colon->GetOutput()->GetPoint(i, p);
            vtkMath::Subtract(p, f2, v);
            double x, y, z;
            x = vtkMath::Dot(v, pd1);
            y = vtkMath::Dot(v, pd2);
            z = vtkMath::Dot(v, pd3);
            double vx[3], vy[3], vz[3], temp1[3], temp2[3], newp[3];
            vx[0] = d1[0]; vx[1] = d1[1]; vx[2] = d1[2];
            vy[0] = d2[0]; vy[1] = d2[1]; vy[2] = d2[2];
            vz[0] = normal[0]; vz[1] = normal[1]; vz[2] = normal[2];
            vtkMath::MultiplyScalar(vx, x);
            vtkMath::MultiplyScalar(vy, y);
            vtkMath::MultiplyScalar(vz, z);
            vtkMath::Add(origin, vx, temp1);
            vtkMath::Add(temp1, vy, temp2);
            vtkMath::Add(temp2, vz, newp);
            newpoints->InsertNextPoint(newp);
            vtkMath::Add(pave, newp, pave);
        }
        m_colon->GetOutput()->SetPoints(newpoints);
        m_rendermanager->renderModel(m_colon->GetActor());

        double position[3];
        vtkMath::MultiplyScalar(pave, 1.0/newpoints->GetNumberOfPoints());
        pave[0] += 2;
        position[0] = 40*0+pave[0];
        position[1] = 40*-1+pave[1];
        position[2] = 40*-0.5+pave[2];
        vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
        camera->SetPosition(position);
        camera->SetFocalPoint(pave);
        m_rendermanager_right->GetRender()->SetActiveCamera(camera);
        m_rendermanager_right->GetRender()->SetBackground(0, 0.5, 3);
    }

    vtkSmartPointer<vtkPlane> clipPlane = vtkSmartPointer<vtkPlane>::New();
    clipPlane->SetOrigin(origin);
    clipPlane->SetNormal(normal);
    vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
    clipper->SetClipFunction(clipPlane);
    clipper->SetInputData(m_colon->GetOutput());
    clipper->Update();
    vtkSmartPointer<vtkPolyData> leftpart = vtkSmartPointer<vtkPolyData>::New();
    leftpart->DeepCopy(clipper->GetOutput());
    vtkSmartPointer<vtkPoints> leftnewpoints = vtkSmartPointer<vtkPoints>::New();
    for(int i=0; i < leftpart->GetNumberOfPoints(); i++)
    {
        double p[3], v[3];
        leftpart->GetPoint(i, p);
        vtkMath::Subtract(p, origin, v);
        double x = vtkMath::Dot(d1, v);
        double z = vtkMath::Dot(normal, v);
        if(z < 0) z = -z;
        double y = vtkMath::Dot(d2, v);
        double r = sqrt(z*z + y*y);
        double angle = acos(y/r);
        //double newr = r * factor * exp(adjust *((r0 - r)/r0));
        double newr = factor * (r * k + b);
        double newangle = angle / factor;

        double newy = newr * cos(newangle);
        double newz = newr * sin(newangle);

        double vx[3], vy[3], vz[3], temp1[3], temp2[3], newp[3];
        vx[0] = d1[0]; vx[1] = d1[1]; vx[2] = d1[2];
        vy[0] = d2[0]; vy[1] = d2[1]; vy[2] = d2[2];
        vz[0] = normal[0]; vz[1] = normal[1]; vz[2] = normal[2];
        vtkMath::MultiplyScalar(vx, x);
        vtkMath::MultiplyScalar(vy, newy);
        vtkMath::MultiplyScalar(vz, newz);
        vtkMath::Add(origin, vx, temp1);
        vtkMath::Add(temp1, vy, temp2);
        vtkMath::Add(temp2, vz, newp);
        leftnewpoints->InsertNextPoint(newp);
        aver += r;
            //std::cout<<newp[0]<<" "<<newp[1]<<" "<<newp[2]<<endl;
    }
    // cut from left part for right overlap
    vtkSmartPointer<vtkPolyData> rightoverlap = vtkSmartPointer<vtkPolyData>::New();
    clipPlane->SetNormal(normal_l);
    clipper->RemoveAllInputs();
    clipper->SetInputData(leftpart);
    clipper->Update();
    rightoverlap->DeepCopy(clipper->GetOutput());
    vtkSmartPointer<vtkPoints> rightoverlapnewpoints = vtkSmartPointer<vtkPoints>::New();
    vtkMath::MultiplyScalar(normal, -1);
    for(int i=0; i < rightoverlap->GetNumberOfPoints(); i++)
    {
        double p[3], v[3];
        rightoverlap->GetPoint(i, p);
        vtkMath::Subtract(p, origin, v);
        double x = vtkMath::Dot(d1, v);
        double z = vtkMath::Dot(normal, v);
        if(z < 0) z = -z;
        double y = vtkMath::Dot(d2, v);
        double r = sqrt(z*z + y*y);
        double angle = 2*3.1415926-acos(y/r);
        std::cout<<angle/3.1415926*180<<endl;
        //double newr = r * factor * exp(adjust *((r0 - r)/r0));
        double newr = factor * (r * k + b);
        double newangle = angle / factor;

        double newy = newr * cos(newangle);
        double newz = newr * sin(newangle);

        double vx[3], vy[3], vz[3], temp1[3], temp2[3], newp[3];
        vx[0] = d1[0]; vx[1] = d1[1]; vx[2] = d1[2];
        vy[0] = d2[0]; vy[1] = d2[1]; vy[2] = d2[2];
        vz[0] = normal[0]; vz[1] = normal[1]; vz[2] = normal[2];
        vtkMath::MultiplyScalar(vx, x);
        vtkMath::MultiplyScalar(vy, newy);
        vtkMath::MultiplyScalar(vz, newz);
        vtkMath::Add(origin, vx, temp1);
        vtkMath::Add(temp1, vy, temp2);
        vtkMath::Add(temp2, vz, newp);
        rightoverlapnewpoints->InsertNextPoint(newp);
        //aver += r;
            //std::cout<<newp[0]<<" "<<newp[1]<<" "<<newp[2]<<endl;
    }
    rightoverlap->SetPoints(rightoverlapnewpoints);
    vtkSmartPointer<vtkPolyDataMapper> rightoverlapmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    rightoverlapmapper->SetInputData(rightoverlap);
    rightoverlapmapper->Update();
    vtkSmartPointer<vtkActor> rightoverlapactor = vtkSmartPointer<vtkActor>::New();
    rightoverlapactor->SetMapper(rightoverlapmapper);


    leftpart->SetPoints(leftnewpoints);


    vtkSmartPointer<vtkPolyData> rightpart = vtkSmartPointer<vtkPolyData>::New();
    clipPlane->SetNormal(normal);
    clipper->RemoveAllInputs();
    clipper->SetInputData(m_colon->GetOutput());
    clipper->SetClipFunction(clipPlane);
    clipper->Update();
    rightpart->DeepCopy(clipper->GetOutput());
    vtkSmartPointer<vtkPoints> rightnewpoints = vtkSmartPointer<vtkPoints>::New();
    for(int i=0; i < rightpart->GetNumberOfPoints(); i++)
    {
        double p[3], v[3];
        rightpart->GetPoint(i, p);
        vtkMath::Subtract(p, origin, v);
        double x = vtkMath::Dot(d1, v);
        double z = vtkMath::Dot(normal, v);
        if(z < 0) z = -z;
        double y = vtkMath::Dot(d2, v);
        double r = sqrt(z*z + y*y);
        double angle = acos(y/r);
        //double newr = r * factor * exp(adjust *((r0 - r)/r0));
        double newr = factor * (r * k + b);
        double newangle = angle / factor;

        double newy = newr * cos(newangle);
        double newz = newr * sin(newangle);

        double vx[3], vy[3], vz[3], temp1[3], temp2[3], newp[3];
        vx[0] = d1[0]; vx[1] = d1[1]; vx[2] = d1[2];
        vy[0] = d2[0]; vy[1] = d2[1]; vy[2] = d2[2];
        vz[0] = normal[0]; vz[1] = normal[1]; vz[2] = normal[2];
        vtkMath::MultiplyScalar(vx, x);
        vtkMath::MultiplyScalar(vy, newy);
        vtkMath::MultiplyScalar(vz, newz);
        vtkMath::Add(origin, vx, temp1);
        vtkMath::Add(temp1, vy, temp2);
        vtkMath::Add(temp2, vz, newp);
        rightnewpoints->InsertNextPoint(newp);
        aver += r;
            //std::cout<<newp[0]<<" "<<newp[1]<<" "<<newp[2]<<endl;
    }
    rightpart->SetPoints(rightnewpoints);
    aver = aver / (leftpart->GetNumberOfPoints() + rightpart->GetNumberOfPoints());
    std::cout<<"Average radius is "<<aver<<endl;

    // one to one mapping used for interaction
    vtkSmartPointer<vtkPoints> onetoonepoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPolyData> onetoone = vtkSmartPointer<vtkPolyData>::New();
    onetoone->DeepCopy(m_colon->GetOutput());
    for(vtkIdType i=0; i<m_colon->GetOutput()->GetNumberOfPoints(); i++)
    {
        double p[3], v[3];
        m_colon->GetOutput()->GetPoint(i, p);
        vtkMath::Subtract(p, origin, v);
        if(vtkMath::Dot(v, normal) > 0) // left points
        {
            double x = vtkMath::Dot(d1, v);
            double z = vtkMath::Dot(normal, v);
            if(z < 0) z = -z;
            double y = vtkMath::Dot(d2, v);
            double r = sqrt(z*z + y*y);
            double angle = acos(y/r);
            //double newr = r * factor * exp(adjust *((r0 - r)/r0));
            double newr = factor * (r * k + b);
            double newangle = angle / factor;

            double newy = newr * cos(newangle);
            double newz = newr * sin(newangle);

            double vx[3], vy[3], vz[3], temp1[3], temp2[3], newp[3];
            vx[0] = d1[0]; vx[1] = d1[1]; vx[2] = d1[2];
            vy[0] = d2[0]; vy[1] = d2[1]; vy[2] = d2[2];
            vz[0] = normal[0]; vz[1] = normal[1]; vz[2] = normal[2];
            vtkMath::MultiplyScalar(vx, x);
            vtkMath::MultiplyScalar(vy, newy);
            vtkMath::MultiplyScalar(vz, newz);
            vtkMath::Add(origin, vx, temp1);
            vtkMath::Add(temp1, vy, temp2);
            vtkMath::Add(temp2, vz, newp);
            onetoonepoints->InsertNextPoint(newp);
        }
        else // right points
        {
            vtkMath::MultiplyScalar(normal, -1);
            double x = vtkMath::Dot(d1, v);
            double z = vtkMath::Dot(normal, v);
            if(z < 0) z = -z;
            double y = vtkMath::Dot(d2, v);
            double r = sqrt(z*z + y*y);
            double angle = acos(y/r);
            //double newr = r * factor * exp(adjust *((r0 - r)/r0));
            double newr = factor * (r * k + b);
            double newangle = angle / factor;

            double newy = newr * cos(newangle);
            double newz = newr * sin(newangle);

            double vx[3], vy[3], vz[3], temp1[3], temp2[3], newp[3];
            vx[0] = d1[0]; vx[1] = d1[1]; vx[2] = d1[2];
            vy[0] = d2[0]; vy[1] = d2[1]; vy[2] = d2[2];
            vz[0] = normal[0]; vz[1] = normal[1]; vz[2] = normal[2];
            vtkMath::MultiplyScalar(vx, x);
            vtkMath::MultiplyScalar(vy, newy);
            vtkMath::MultiplyScalar(vz, newz);
            vtkMath::Add(origin, vx, temp1);
            vtkMath::Add(temp1, vy, temp2);
            vtkMath::Add(temp2, vz, newp);
            onetoonepoints->InsertNextPoint(newp);
            vtkMath::MultiplyScalar(normal, -1);
        }
    }
    onetoone->SetPoints(onetoonepoints);

    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
    appendFilter->AddInputData(leftpart);
    appendFilter->AddInputData(rightpart);
    appendFilter->Update();
    vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
    cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
    cleanFilter->Update();
    vtkSmartPointer<vtkPolyData> bended = vtkSmartPointer<vtkPolyData>::New();
    bended->DeepCopy(cleanFilter->GetOutput());

    vtkSmartPointer<vtkPolyDataMapper> selectedMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    selectedMapper->SetInputData(bended);
    //selectedMapper->SetInputData(onetoone);
    vtkSmartPointer<vtkActor> selectedActor = vtkSmartPointer<vtkActor>::New();
    selectedActor->SetMapper(selectedMapper);

    vtkSmartPointer<vtkFeatureEdges> featureEdges =
            vtkSmartPointer<vtkFeatureEdges>::New();
    featureEdges->SetInputData(bended);
    featureEdges->BoundaryEdgesOn();
    featureEdges->FeatureEdgesOff();
    featureEdges->ManifoldEdgesOff();
    featureEdges->NonManifoldEdgesOff();
    featureEdges->Update();

    // Visualize
    vtkSmartPointer<vtkPolyDataMapper> edgeMapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
    edgeMapper->SetInputConnection(featureEdges->GetOutputPort());
    vtkSmartPointer<vtkActor> edgeActor =
            vtkSmartPointer<vtkActor>::New();
    edgeActor->SetMapper(edgeMapper);
    edgeActor->GetProperty()->SetColor(255,0,0);

    m_colon_new->Object::SetInput(onetoone);

    QVTKWidget* right = this->findChild<QVTKWidget*>("right");
    if(test)
    {
        m_rendermanager_right->GetRender()->AddActor(selectedActor);
        //m_rendermanager_right->GetRender()->AddActor(edgeActor);
        m_rendermanager_right->GetRender()->AddActor(rightoverlapactor);
        vtkSmartPointer<vtkWindowToImageFilter> screenshot = vtkSmartPointer<vtkWindowToImageFilter>::New();
        screenshot->SetInput(right->GetRenderWindow());
        screenshot->SetMagnification(2);
        screenshot->SetInputBufferTypeToRGBA();
        screenshot->ReadFrontBufferOff();
        screenshot->Update();
        vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
        char filename[50];
        sprintf(filename, "%.2lf.png", factor);
        std::cout<<filename<<endl;
        writer->SetFileName(filename);
        writer->SetInputConnection(screenshot->GetOutputPort());
        writer->Write();
    }
    else
        m_rendermanager_right->renderModel(selectedActor);


    tracer_inverse->GetLineProperty()->SetLineWidth(5);
    tracer_inverse->SetInteractor(right->GetRenderWindow()->GetInteractor());
    tracer_inverse->SetViewProp(selectedActor);

    right->GetRenderWindow()->Render();
    vtkSmartPointer<vtkCallbackCommand> callback =
            vtkSmartPointer<vtkCallbackCommand>::New();
    callback->SetCallback(CallbackFunction_inverse);
    callback->SetClientData(this);
    tracer_inverse->AddObserver(vtkCommand::EndInteractionEvent, callback);


    QVTKWidget* left = this->findChild<QVTKWidget*>("left");
    left->GetRenderWindow()->Render();

    //m_showselectedwindow.show();
    //m_showselectedwindow.RenderSelected(selectedActor);

    m_filemanager->SaveFile(bended, "openedcolon.off", test);

    //m_showselectedwindow.RenderSelected(edgeActor);
    //m_showselectedwindow.GetRenderManager().GetRender()->SetBackground(0.1, 0.6, 1);
}

void MainWindow::on_tracer_toggled(bool checked)
{
    if(checked)
    {
        tracer->On();
    }
    else
    {
        tracer->Off();

        m_rendermanager->GetRender()->RemoveActor(m_tracermark_1->GetActor());
        m_rendermanager_right->GetRender()->RemoveActor(m_tracermark_2->GetActor());
        vtkSmartPointer<vtkPoints> empty = vtkSmartPointer<vtkPoints>::New();
        m_tracermark_1->InputPoints(empty);
        m_tracermark_2->InputPoints(empty);

        QVTKWidget* left = this->findChild<QVTKWidget*>("left");
        left->GetRenderWindow()->Render();
        QVTKWidget* right = this->findChild<QVTKWidget*>("right");
        right->GetRenderWindow()->Render();
    }
}

void MainWindow::on_action_Load_Deformed_Surface_triggered()
{
    //Open file dialog
    QString filePath = QFileDialog::getOpenFileName(
                this, tr("Open File"), "",
                tr("3Dmodels (*.stl *.vtp *.ply *.obj *.off)"));
    if(filePath.isEmpty()) return;

    m_filemanager->LoadNewFile(filePath);
    m_colon_new->Object::SetInput(m_filemanager->getfile());

    m_rendermanager_right->renderModel(m_colon_new->GetActor());
    QVTKWidget* right = this->findChild<QVTKWidget*>("right");

    tracer_inverse->GetLineProperty()->SetLineWidth(5);
    tracer_inverse->SetInteractor(right->GetRenderWindow()->GetInteractor());
    tracer_inverse->SetViewProp(m_colon_new->GetActor());

    right->GetRenderWindow()->Render();
    vtkSmartPointer<vtkCallbackCommand> callback =
            vtkSmartPointer<vtkCallbackCommand>::New();
    callback->SetCallback(CallbackFunction_inverse);
    callback->SetClientData(this);
    tracer_inverse->AddObserver(vtkCommand::EndInteractionEvent, callback);


    QVTKWidget* left = this->findChild<QVTKWidget*>("left");
    left->GetRenderWindow()->Render();

    if(twoWindowReady == 0)
        twoWindowReady++;
    else
    {
        twoWindowReady = 2;
        QCheckBox* checkbox_l = this->findChild<QCheckBox*>("tracer");
        checkbox_l->setEnabled(true);
        QCheckBox* checkbox_r = this->findChild<QCheckBox*>("tracer_r");
        checkbox_r->setEnabled(true);
    }
}

void MainWindow::on_tracer_r_toggled(bool checked)
{
    if(checked)
    {
        tracer_inverse->On();
    }
    else
    {
        tracer_inverse->Off();

        m_rendermanager->GetRender()->RemoveActor(m_tracermark_inverse_1->GetActor());
        m_rendermanager_right->GetRender()->RemoveActor(m_tracermark_inverse_2->GetActor());
        vtkSmartPointer<vtkPoints> empty = vtkSmartPointer<vtkPoints>::New();
        m_tracermark_inverse_1->InputPoints(empty);
        m_tracermark_inverse_2->InputPoints(empty);

        QVTKWidget* left = this->findChild<QVTKWidget*>("left");
        left->GetRenderWindow()->Render();
        QVTKWidget* right = this->findChild<QVTKWidget*>("right");
        right->GetRenderWindow()->Render();
    }
}

void MainWindow::on_horizontalScrollBar_sliderMoved(int position)
{
    std::cout<<position<<endl;
    vtkSmartPointer<vtkRegularPolygonSource> polygonSource_l = vtkSmartPointer<vtkRegularPolygonSource>::New();

    double p_l[3], n_l[3];
    m_oldcenterline->GetPoint(position, p_l);
    m_oldplanenormals->GetTuple(position, n_l);
    polygonSource_l->SetNumberOfSides(50);
    polygonSource_l->SetRadius(2); // 40 for colon
    polygonSource_l->SetCenter(p_l);
    polygonSource_l->SetNormal(n_l);
    polygonSource_l->Update();

    m_sweepingplane_l->InputData(polygonSource_l->GetOutput());
    QVTKWidget* left = this->findChild<QVTKWidget*>("left");
    left->GetRenderWindow()->Render();

    vtkSmartPointer<vtkRegularPolygonSource> polygonSource_r = vtkSmartPointer<vtkRegularPolygonSource>::New();

    double p_r[3], n_r[3];
    m_newcenterline->GetPoint(position, p_r);
    m_newplanenormals->GetTuple(position, n_r);
    polygonSource_r->SetNumberOfSides(50);
    polygonSource_r->SetRadius(2); // 40 for colon
    polygonSource_r->SetCenter(p_r);
    polygonSource_r->SetNormal(n_r);
    polygonSource_r->Update();

    m_sweepingplane_r->InputData(polygonSource_r->GetOutput());
    QVTKWidget* right = this->findChild<QVTKWidget*>("right");
    right->GetRenderWindow()->Render();
}

void MainWindow::on_action_Load_Centerlines_triggered()
{
    QString filePath = QFileDialog::getOpenFileName(
                this, tr("Open File"),"",
                tr("Centerline Config (*.txt)"));
    if(filePath.isEmpty()) return;
    ifstream file;
    file.open(filePath.toStdString().c_str());
    int N = 0;
    file>>N;
    std::cout<<N<<" points"<<endl;
    vtkSmartPointer<vtkPoints> emptypoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkDoubleArray> emptydouble = vtkSmartPointer<vtkDoubleArray>::New();
    emptydouble->SetNumberOfComponents(3);
    m_oldcenterline->DeepCopy(emptypoints);
    m_newcenterline->DeepCopy(emptypoints);
    m_oldplanenormals->DeepCopy(emptydouble);
    m_newplanenormals->DeepCopy(emptydouble);

    for(int i=0; i<N; i++)
    {
        double pold[3], pnew[3], nold[3], nnew[3];
        file >> pold[0] >> pold[1] >> pold[2];
        file >> nold[0] >> nold[1] >> nold[2];
        file >> pnew[0] >> pnew[1] >> pnew[2];
        file >> nnew[0] >> nnew[1] >> nnew[2];
        m_oldcenterline->InsertNextPoint(pold);
        m_newcenterline->InsertNextPoint(pnew);
        m_oldplanenormals->InsertNextTuple(nold);
        m_newplanenormals->InsertNextTuple(nnew);
    }
    VisualizePoints(m_oldcenterline, 1, 0, 0, 1, m_rendermanager);
    VisualizePoints(m_newcenterline, 1, 0, 0, 1, m_rendermanager_right);


    QScrollBar* scrollbar = this->findChild<QScrollBar*>("horizontalScrollBar");
    scrollbar->setMaximum(N-1);
    scrollbar->setEnabled(true);

    m_rendermanager->renderModel(m_sweepingplane_l->GetActor());
    QVTKWidget* left = this->findChild<QVTKWidget*>("left");
    left->GetRenderWindow()->Render();

    m_rendermanager_right->renderModel(m_sweepingplane_r->GetActor());
    QVTKWidget* right = this->findChild<QVTKWidget*>("right");
    right->GetRenderWindow()->Render();
}

void MainWindow::on_action_Compute_Height_triggered()
{
    vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoothFilter->SetInputData(m_colon->GetOutput());
    smoothFilter->SetNumberOfIterations(30);
    smoothFilter->SetRelaxationFactor(0.2);
    //smoothFilter->FeatureEdgeSmoothingOff();
    smoothFilter->BoundarySmoothingOn();
    smoothFilter->Update();

    vtkSmartPointer<vtkPolyData> smoothed = vtkSmartPointer<vtkPolyData>::New();
    smoothed->DeepCopy(smoothFilter->GetOutput());
    std::cout<<smoothed->GetNumberOfPoints()<<endl;

    vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
    normalGenerator->SetInputData(smoothed);
    normalGenerator->ComputePointNormalsOn();
    normalGenerator->ComputeCellNormalsOff();
    normalGenerator->ConsistencyOn();
    normalGenerator->SplittingOff();
    normalGenerator->Update();
    smoothed = normalGenerator->GetOutput();

    vtkFloatArray* normalDataFloat =
      vtkFloatArray::SafeDownCast(smoothed->GetPointData()->GetArray("Normals"));

    if(normalDataFloat)
    {
      int nc = normalDataFloat->GetNumberOfTuples();
      std::cout << "There are " << nc
              << " components in normalDataFloat" << std::endl;
    }



    unsigned char red[3] = {0,255,0};
    unsigned char blue[3] = {0,0,255};
    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");

    double minVal = INFINITY;
    double maxVal = -INFINITY;

    vtkSmartPointer<vtkPolyData> enhancement = vtkSmartPointer<vtkPolyData>::New();
    enhancement->DeepCopy(m_colon->GetOutput());
    vtkSmartPointer<vtkPoints> enhancedpoints = vtkSmartPointer<vtkPoints>::New();
    enhancedpoints->DeepCopy(enhancement->GetPoints());

    vtkSmartPointer<vtkDoubleArray> Heights = vtkSmartPointer<vtkDoubleArray>::New();
    for(int i = 0; i < smoothed->GetNumberOfPoints(); i++){
        double newp[3];
        double oldp[3];
        m_colon->GetOutput()->GetPoint(i, oldp);
        smoothed->GetPoint(i, newp);

        double v[3];
        vtkMath::Subtract(newp, oldp, v);

        float normal[3];
        normalDataFloat->GetTypedTuple(i, normal);

        if((v[0]*(double)normal[0] + v[1]*(double)normal[1] + v[2]*(double)normal[2]) < 0){
            double p[3];
            enhancement->GetPoint(i, p);
            double vv[3];
            vv[0] = v[0]; vv[1] = v[1]; vv[2] = v[2];
            vtkMath::MultiplyScalar(vv, -1);
            vtkMath::Add(p, vv, p);
            enhancedpoints->SetPoint(i, p);
        }

        double val = vtkMath::Norm(v) * ((v[0]*(double)normal[0] + v[1]*(double)normal[1] + v[2]*(double)normal[2]) >=0 ? 1 : -1);

        minVal = val < minVal ? val : minVal;
        maxVal = val > maxVal ? val : maxVal;

        Heights->InsertNextValue(val);
    }

    enhancement->SetPoints(enhancedpoints);

    for(int i = 0; i < Heights->GetNumberOfValues(); i++){
        double r = (Heights->GetValue(i) - minVal) / (maxVal - minVal);
        unsigned char color[3];
        color[0] = (unsigned char)(r * (double)blue[0] + (1-r) * (double)red[0]);
        color[1] = (unsigned char)(r * (double)blue[1] + (1-r) * (double)red[1]);
        color[2] = (unsigned char)(r * (double)blue[2] + (1-r) * (double)red[2]);

        //std::cout<<color[0]<<" "<<color[1]<<" "<<color[2]<<endl;
        colors->InsertNextTypedTuple(color);
    }


    //m_colon->GetOutput()->GetPointData()->SetScalars(colors);

    vtkSmartPointer<vtkPolyData> heightmap = vtkSmartPointer<vtkPolyData>::New();
    heightmap->DeepCopy(m_colon->GetOutput());
    heightmap->GetPointData()->SetScalars(colors);

    //m_filemanager->SaveFile(heightmap, "heightmap.ply");


    vtkSmartPointer<vtkPolyDataMapper> Mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    Mapper->SetInputData(heightmap);
    //Mapper->SetInputData(enhancement);
    Mapper->Update();

    vtkSmartPointer<vtkActor> Actor = vtkSmartPointer<vtkActor>::New();
    Actor->SetMapper(Mapper);

    m_rendermanager_right->renderModel(Actor);
    QVTKWidget* right = this->findChild<QVTKWidget*>("right");
    right->GetRenderWindow()->Render();
}

/*
vtkSmartPointer<vtkDoubleArray> MainWindow::ComputeNormals(vtkSmartPointer<vtkPolyData> mesh){

    int count = 0;
    int seed = 0;
    vtkSmartPointer<vtkIdList> currentlevel = vtkSmartPointer<vtkIdList>::New();
    // init
    currentlevel->InsertNextId(seed);

    while(currentlevel->GetNumberOfIds() > 0)
    {
        count += currentlevel->GetNumberOfIds();
        //std::cout<<"processed "<<count<<" points"<<endl;
        vtkSmartPointer<vtkIdList> nextlevel = vtkSmartPointer<vtkIdList>::New();
        for(vtkIdType n = 0; n < currentlevel->GetNumberOfIds(); n++)
        {
            vtkIdType pointid = currentlevel->GetId(n);
            vtkSmartPointer<vtkIdList> connectedVertices = vtkSmartPointer<vtkIdList>::New();
            connectedVertices = GetConnectedVertices(t_colon, pointid);



        }
        currentlevel->Reset();
        currentlevel->DeepCopy(nextlevel);
    } // while(1)
}
*/

void MainWindow::on_action_Angular_Missing_Ratio_triggered()
{
    m_centerline->GetOutput()->Reset();
    QString filePath = QFileDialog::getOpenFileName(
                this, tr("Open File"),"",
                tr("Centerline File (*.vtp)"));
    std::cout<<filePath.toStdString().c_str()<<endl;
    const char * path = filePath.toStdString().c_str();
    int len = strlen(path);
    std::cout<<len<<endl;
    int start = -1, end = -1;
    for(int i = len-1; i >= 0; i--){
        if(path[i] >= '0' && path[i] <= '9'){
            if(end < 0){
                end = i;
            }
            else{
                // do nothing
            }
        }
        else{
            if(end < 0){
                // do nothing
            }
            else{
                start = i + 1;
                break;
            }
        }
    }
    std::string casename;
    for(int i = start; i <= end; i++){
        casename += path[i];
    }
    std::cout<<casename<<endl;



    if(filePath.isEmpty()) return;
    m_filemanager->LoadNewFile(filePath);
    m_centerline->Object::SetInput(m_filemanager->getfile());
    // Uniform Sampling
    m_centerline->UniformSample(100);


    m_centerline->ComputeAngularMissing(m_rendermanager, m_rendermanager_right, m_colon->GetOutput(), m_filemanager, casename);

    m_rendermanager->renderModel(m_centerline->GetActor());
    QVTKWidget* left = this->findChild<QVTKWidget*>("left");
    left->GetRenderWindow()->Render();
}

void MainWindow::on_action_Maunally_Create_Centerline_triggered()
{
    int resolution = 200;
    vtkSmartPointer<vtkParametricSpline> spline = vtkSmartPointer<vtkParametricSpline>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(-8.290, 15.469, 44.268);
    points->InsertNextPoint(-5.426, 5.429, 29.437);
    points->InsertNextPoint(-4.140, -5.792, 12.248);
    spline->SetPoints(points);

    vtkSmartPointer<vtkParametricFunctionSource> functionSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
    functionSource->SetParametricFunction(spline);
    functionSource->Update();
    functionSource->SetUResolution(resolution);
    functionSource->Update();

    m_filemanager->SaveFile(functionSource->GetOutput(), "ManualCenterline.vtp");
    std::cout<<"ManualCenterline.vtp created ("<<functionSource->GetOutput()->GetNumberOfPoints()<<" points)"<<std::endl;
}
