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

    void on_pushButton_3_clicked();

    void on_pushButton_4_clicked();

private:
    Ui::MainWindow *ui;

    QString filename;
    void start();

    void get_projections();
    void get_sinogram();

    template<typename Image>
    QImage create_QImage(Image& img, unsigned w, unsigned h)
    {
        QImage image(w,h,QImage::Format_Grayscale8);
        for (auto i = 0; i < w; ++i)
            for (auto j = 0 ; j < h; ++j)
            {
                int v = img(i,j);
                v = std::clamp(v,0,255);

                QRgb c = qRgb(v,v,v);

                image.setPixel(i,j,c);
            }

        return image;
    }

    template<typename Image>
    QImage create_QImage_from_complex(Image& img, unsigned w, unsigned h)
    {
        double max = std::abs(*std::max_element(std::begin(img),std::end(img),[](auto& a, auto& b){return std::abs(a) < std::abs(b);}));

        QImage image(w,h,QImage::Format_Grayscale8);
        for (auto i = 0; i < w; ++i)
            for (auto j = 0 ; j < h; ++j)
            {
                double v = std::sqrt(img(i,j).real()*img(i,j).real() + img(i,j).imag()*img(i,j).imag());
                int v_ = v / max * 255;

                QRgb c = qRgb(v_,v_,v_);

                image.setPixel(i,j,c);
            }

        return image;
    }

    arma::mat object;
    arma::mat sinogram;
    arma::cx_mat spectre;
    arma::mat restored_object;

};
#endif // MAINWINDOW_H
