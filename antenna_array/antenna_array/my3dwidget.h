#ifndef MY3DWIDGET_H
#define MY3DWIDGET_H

#include <QObject>
#include <QWidget>
#include <glm/glm.hpp>
#include <boost/numeric/ublas/matrix.hpp>

constexpr double PI = 3.1415265358979323;

template <typename Real = double>
auto from_sphere_coords(Real r, Real phi, Real theta = PI/2.0)
{
    Real x = r * cos(phi) * sin(theta);
    Real y = r * sin(phi) * sin(theta);
    Real z = r * cos(theta);
    return std::make_tuple(x,y,z);
}




class My3DWidget : public QWidget
{
    Q_OBJECT
public:
    explicit My3DWidget(QWidget *parent = nullptr);


signals:

    // QWidget interface
protected:
    void paintEvent(QPaintEvent *event);


public:
    boost::numeric::ublas::matrix<glm::vec3> mesh;

    template <typename InputIt>
    void set_data(InputIt begin, InputIt end, unsigned matrix_width)
    {
        mesh.resize((end - begin)/matrix_width,matrix_width);

        for (auto i =0; i < mesh.size1(); ++i)
        {
            for (auto j = 0; j < mesh.size2();++j)
            {
                // auto& p = mesh.at_element(i,j);
                auto& p = *(begin + i*mesh.size1() + j);
                auto [x,y,z] = from_sphere_coords(p.x,p.y,p.z);

                mesh.at_element(i,j) = glm::vec3{x,y,z};
            }
        }
    }
};

#endif // MY3DWIDGET_H
