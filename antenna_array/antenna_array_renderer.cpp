#include "antenna_array_renderer.h"
#include "qobject.h"
#include "qpoint.h"
#include "qsize.h"
#include "qtransform.h"
#include <QMouseEvent>
#include <QPainter>
#include <cmath>

template <typename T, typename T_, typename T__>
inline bool in_interval(T elem, T_ left, T__ right)
{
    return (elem >= left) && (elem < right);
}




antenna_array_renderer::antenna_array_renderer(QWidget *parent)
    : QWidget{parent}
{}

unsigned int antenna_array_renderer::array_width(){if (array.empty()) return 0; else return array[0].size();}

unsigned int antenna_array_renderer::array_height(){return array.size();}

bool antenna_array_renderer::array_empty(){return array.empty() || array[0].empty();}

void antenna_array_renderer::set_array_range(unsigned int height, unsigned int width)
{
    array = std::vector<std::vector<bool>>(height,std::vector<bool>(width,1));
}

void antenna_array_renderer::toggle_antenna(int i, int j)
{
    if (in_interval(i,0,array_height()) && in_interval(j,0,array_width()))
    {
        array[i][j] = !array[i][j];
    }
}

void antenna_array_renderer::array_set_empty()
{
    for (int i = 0; i < array.size(); ++i)
        for (int j =0 ; j < array[0].size(); ++j)
            array[i][j] = false;
}

void antenna_array_renderer::array_set_all()
{
    for (int i = 0; i < array.size(); ++i)
        for (int j =0 ; j < array[0].size(); ++j)
            array[i][j] = true;
}

void draw_antenna(QPainter& painter, QPointF center, double width, double height)
{
    QPointF l1 = center + QPointF(-width/2.,height/2.);
    QPointF r1 = center + QPointF(width/2.,height/2.);
    QPointF up = center + QPointF(0,height/2.);
    QPointF down = center + QPointF(0,-height/2.);

    painter.drawLine(QLineF(l1,center));
    painter.drawLine(QLineF(r1,center));
    painter.drawLine(QLineF(up,down));
}
void draw_antenna(QPainter& painter, QPointF center)
{
    draw_antenna(painter,center,0.7,0.7);
}

void draw_table(QPainter& painter, unsigned w, unsigned h)
{
    painter.fillRect(QRectF(QPointF(0,0),QSizeF(w,h)),QBrush(QColor(255,255,255)));
    for (auto j = 1; j < w; ++j) {
        painter.drawLine(QLineF(QPointF(j,h),QPointF(j,0)));
    }
    for (auto i = 1; i < h; ++i) {
        painter.drawLine(QLineF(QPointF(w,i),QPointF(0,i)));
    }
}


void antenna_array_renderer::paintEvent(QPaintEvent *event)
{
    if (array_empty()) return;

    QPainter painter(this);
    QRect r = this->rect();

    double w = array_width();
    double h = array_height();

    transform.reset();
    transform.scale(r.width()/w,-r.height()/h);
    transform.translate(0,-h);

    painter.setTransform(transform);
    painter.setPen(QPen(QColor(0,0,0),0));

    draw_table(painter, w, h);
    for (auto i = 0; i < w; i++) {
        for (auto j = 0; j < h; j++) {
            if (array[j][i])
            draw_antenna(painter,QPointF{i + 0.5,j + 0.5});
        }
    }
}

void antenna_array_renderer::mousePressEvent(QMouseEvent *event)
{
    auto it = transform.inverted();
    auto p = it.map(QPointF(event->pos()));

    int i = std::floor(p.y());
    int j = std::floor(p.x());

    toggle_antenna(i,j);

    repaint();
    last_pressed_point = {i,j};
}

void antenna_array_renderer::mouseMoveEvent(QMouseEvent *event)
{
    if (event->buttons() & Qt::LeftButton)
    {
        auto it = transform.inverted();
        auto p = it.map(QPointF(event->pos()));

        int i = std::floor(p.y());
        int j = std::floor(p.x());

        static int prev_i = i;
        static int prev_j = j;

        if ((prev_i == i) && (prev_j == j))
        {return;}
        if (QPoint(i,j) == last_pressed_point)
        {return;}

        prev_i = i; prev_j = j;

        toggle_antenna(i,j);

        repaint();
    }
}
