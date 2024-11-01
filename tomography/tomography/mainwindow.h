#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFileDialog>
#include <QMessageBox>
#include <armadillo>
#include "qcustomplot.h"

#include "image_io.hpp"
#include "tomography.hpp"

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_pushButton_clicked();


    void on_pushButton_2_clicked();

    void on_horizontalSlider_sliderMoved(int position);

private:
    Ui::MainWindow *ui;

    QString filename;
    void start();

    void get_projections();
    void get_sinogram();;

    template<typename Image>
    QImage create_QImage(Image img, unsigned w, unsigned h)
    {
        QImage image(w,h,QImage::Format_Grayscale8);
        for (auto i = 0; i < w; ++i)
            for (auto j =0 ;j < h; ++j)
            {
                auto c = qRgb(img(i,j),img(i,j),img(i,j));
                image.setPixel(i,j,c);
            }

        return image;
    }

    arma::mat object;
};
#endif // MAINWINDOW_H
