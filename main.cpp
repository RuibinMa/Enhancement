#include <vtkAutoInit.h>
#include "mainwindow.h"
#include <QApplication>
#include <QVTKWidget.h>

VTK_MODULE_INIT(vtkRenderingOpenGL2)

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}
