#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    ui->widget->addGraph();
    ui->widget->graph(0)->valueAxis()->setLabel("Коэффициент поглощения");
    ui->widget->graph(0)->keyAxis()->setLabel("сдвиг t, пикс");

    ui->object_image->setStyleSheet("QLabel { background-color : rgb(192,192,192) ;}");
    ui->label_2->setStyleSheet("QLabel { background-color : rgb(192,192,192) ;}");
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    QFileDialog dialog;

    dialog.setFileMode(QFileDialog::ExistingFile);
    dialog.setAcceptMode(QFileDialog::AcceptOpen);

    filename = dialog.getOpenFileName();

    if (!filename.isEmpty())
    start();
}

void MainWindow::start(){
    object = gray_image_from_file<arma::mat>(filename.toStdString());

    auto image = create_QImage(object,object.n_rows,object.n_cols);

    ui->object_image->setPixmap(QPixmap::fromImage(image));
}

void MainWindow::get_sinogram(){

    double x_shift_max = sqrt(object.n_rows*object.n_rows/4. + object.n_cols*object.n_cols/4.);
    arma::mat m(180,2*x_shift_max);

    for (auto i = 0 ; i < 180; ++i)
    {
        for (auto j = 0; j < m.n_cols; ++j)
        {
            double x_shift =  double (j) - double(m.n_cols)/2.0;
            m(i,j) = get_projection(object,object.n_rows,object.n_cols,
                                     i*3.14159265358979323/180.0,
                                    x_shift);
        }
    }

    m /= *std::max_element(m.begin(),m.end()) / 255.0;

    sinogram = m;
    auto image = create_QImage(m,m.n_rows,m.n_cols);
    ui->label_2->setPixmap(QPixmap::fromImage(image));
}


void MainWindow::on_pushButton_2_clicked()
{
    get_sinogram();
}


void MainWindow::on_horizontalSlider_sliderMoved(int position)
{
    double angle = position * 3.14159265358979323 / 180.0;
    QVector<double> y(object.n_rows);
    QVector<double> x(object.n_rows);
    for (auto i =0; i < object.n_rows; ++i)
    {
        // y[i] = -log10(1 + get_projection(object,object.n_rows,object.n_cols,angle,(double)i - double(object.n_rows)/2.0));
        y[i] = get_projection(object,object.n_rows,object.n_cols,angle,(double)i - double(object.n_rows)/2.0);
        x[i] = i;
    }

    auto g = ui->widget->graph(0);

    g->setData(x,y,true);
    g->rescaleAxes(true);
    // g->valueAxis()->setScaleType(QCPAxis::stLogarithmic);
    ui->widget->replot();

    ui->label->setText(QString::number(position) + " градусов: ");
}


void MainWindow::on_pushButton_3_clicked()
{
    auto w = sinogram.n_rows;
    auto h = sinogram.n_cols;

    sinogram = bilinear_interpolation<arma::mat>(sinogram,w,h,w,closest_pow_2(h));
    h = closest_pow_2(h);

    spectre = arma::cx_mat(h,h);
    sinogram_to_spectre(sinogram,w,h,spectre);


    auto image = create_QImage_from_complex(spectre,spectre.n_rows,spectre.n_cols);
    ui->label_3->setPixmap(QPixmap::fromImage(image));

}


void MainWindow::on_pushButton_4_clicked()
{
    auto s = spectre;
    // arma::mat s(spectre.n_rows,spectre.n_cols);

    for (auto i = 0; i < s.n_rows/2;++i)
    for (auto j = 0; j < s.n_cols/2;++j)
    {
        std::swap(s(i,j),s(i+s.n_rows/2,j+s.n_cols/2));
        std::swap(s(i+s.n_rows/2,j),s(i,j+s.n_cols/2));
    }

    s = bilinear_interpolation<arma::cx_mat>(s,s.n_rows,s.n_cols,closest_pow_2(s.n_rows),closest_pow_2(s.n_cols));


    auto s__ = fft_2d<arma::cx_mat,arma::cx_mat>(s,s.n_rows,s.n_cols,true);
    // arma::cx_mat s__ = arma::ifft2(s);


    for (auto i = 0; i < s__.n_rows/2;++i)
    for (auto j = 0; j < s__.n_cols/2;++j)
    {
        std::swap(s__(i,j),s__(i+s__.n_rows/2,j+s__.n_cols/2));
        std::swap(s__(i+s__.n_rows/2,j),s__(i,j+s__.n_cols/2));
    }

    restored_object = arma::mat(s__.n_rows,s__.n_cols);


    double max = std::abs(*std::max_element(s__.begin(),s__.end(),[](auto a, auto b){return std::abs(a) < std::abs(b);}));

    std::transform(s__.begin(),s__.end(),restored_object.begin(),[&](auto& v){return double(std::abs(v)/max * 255);});



    auto image_1 = create_QImage_from_complex(s,s.n_rows,s.n_cols);
    ui->label_3->setPixmap(QPixmap::fromImage(image_1));

    auto image = create_QImage(restored_object,restored_object.n_rows,restored_object.n_cols);
    ui->label_4->setPixmap(QPixmap::fromImage(image));
}

