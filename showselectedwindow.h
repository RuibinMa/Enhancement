#ifndef SHOWSELECTEDWINDOW_H
#define SHOWSELECTEDWINDOW_H

#include <QMainWindow>
#include "rendermanager.h"
#include <QVTKWidget.h>
#include <vtkActor.h>
#include <vtkSmartPointer.h>

namespace Ui {
class ShowSelectedWindow;
}

class ShowSelectedWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit ShowSelectedWindow(QWidget *parent = 0);
    ~ShowSelectedWindow();
    RenderManager GetRenderManager(){return m_showselectedrender;}
    void RenderSelected(vtkSmartPointer<vtkActor> t_actor);

private:
    Ui::ShowSelectedWindow *ui;
    RenderManager m_showselectedrender;
};

#endif // SHOWSELECTEDWINDOW_H
