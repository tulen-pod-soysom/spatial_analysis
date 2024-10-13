#include "my3dwidget.h"
#include "qpainter.h"
#include <glm/ext/matrix_clip_space.hpp>
#include <glm/ext/matrix_transform.hpp>

My3DWidget::My3DWidget(QWidget *parent)
    : QWidget{parent}
{

}

void My3DWidget::paintEvent(QPaintEvent *event){
    QPainter painter(this);

    painter.fillRect(this->rect(),QBrush(Qt::GlobalColor::black));

    if ((mesh.size1() == 0) || (mesh.size2() == 0)) return;

    auto m = mesh.size1();
    auto n = mesh.size2();

    float w = rect().width();
    float h = rect().height();

    glm::mat4 mat(1.0);
    mat = glm::rotate(mat,glm::radians(45.0f),glm::vec3(0.0,0.0,1.0));
    mat *= glm::ortho(-3.f,3.f,-3.f,3.f,0.1f,100.f);
    mat *= glm::mat4(
        w,0,0,0,
        0,h,0,0,
        0,0,1,0,
        0,0,0,1);

    for (auto i = 0; i < m - 1; ++i)
        for (auto j =0 ; j < n - 1; ++j)
        {
            auto& p1 = mesh(i+0,j+0);
            auto& p2 = mesh(i+1,j+0);
            auto& p3 = mesh(i+1,j+1);
            auto& p4 = mesh(i+0,j+1);

            auto p1_ = mat * glm::vec4(p1,1.0);
            auto p2_ = mat * glm::vec4(p2,1.0);
            auto p3_ = mat * glm::vec4(p3,1.0);
            auto p4_ = mat * glm::vec4(p4,1.0);

            QPointF poly[4];
            poly[0] = {p1_.x,p1_.y};
            poly[1] = {p2_.x,p2_.y};
            poly[2] = {p3_.x,p3_.y};
            poly[3] = {p4_.x,p4_.y};

            painter.drawPolygon(poly,4);
        }

}
