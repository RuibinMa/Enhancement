#include "showselectedwindow.h"
#include "ui_showselectedwindow.h"

ShowSelectedWindow::ShowSelectedWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::ShowSelectedWindow)
{
    ui->setupUi(this);

    QVTKWidget* widget_showselected = this->findChild<QVTKWidget*>("showselectedwindow");
    widget_showselected->GetRenderWindow()->AddRenderer(m_showselectedrender.GetRender());
}

ShowSelectedWindow::~ShowSelectedWindow()
{
    delete ui;
}

void ShowSelectedWindow::RenderSelected(vtkSmartPointer<vtkActor> t_actor)
{
    m_showselectedrender.renderModel(t_actor);
    QVTKWidget* widget_showselected = this->findChild<QVTKWidget*>("showselectedwindow");
    widget_showselected->GetRenderWindow()->Render();
}
