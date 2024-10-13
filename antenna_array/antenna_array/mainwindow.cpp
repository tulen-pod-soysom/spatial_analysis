#include "mainwindow.h"

#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->widget->set_array_range(11,11);

    auto& p = ui->widget_2;

    const int N = 256;

    m = new QCPColorMap(p->xAxis,p->yAxis);
    m->data()->setSize(N, N);
    // m->data()->setRange(QCPRange(-180,180), QCPRange(90,-90));
    m->data()->setRange(QCPRange(-1000,1000),QCPRange(-1000,1000));


    s = new QCPColorScale(p);
    p->plotLayout()->addElement(0, 1, s);
    s->setType(QCPAxis::atRight);
    m->setColorScale(s);
    m->setGradient(QCPColorGradient::gpHot);
    m->rescaleDataRange();
    p->rescaleAxes();
    p->replot();
}

MainWindow::~MainWindow()
{
    delete ui;
}


QVector<double> linspace(double left, double right, size_t size)
{
    QVector<double> x(size);
    for (int i = 0; i < size; ++i) {
        x[i] = left + i*(right-left)/ double(size);
    }
    return x;
}

void MainWindow::on_pushButton_2_clicked()
{
    array.create_grid_configuration_from(ui->widget->array);

    auto& p = ui->widget_2;
    const int N = 256;


    auto z = get_directional_diagram_projection(array,N,1000);

    for (auto i = 0; i < N; ++i) {
        for (auto j = 0; j <N; ++j) {
            m->data()->setCell(i, j, z[i*N+j]);
        }
    }

    m->rescaleDataRange(true);
    auto r = m->dataRange();
    r.lower = 0;
    m->setDataRange(r);
    p->rescaleAxes();
    ui->widget_2->replot();
}


void MainWindow::on_pushButton_clicked()
{
    ui->pushButton_2->setEnabled(true);
}


void MainWindow::on_pushButton_4_clicked()
{
    ui->widget->array_set_empty();
    ui->widget->repaint();
}


void MainWindow::on_pushButton_3_clicked()
{
    ui->widget->array_set_all();
    ui->widget->repaint();
}

