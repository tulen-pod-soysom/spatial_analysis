#ifndef ANTENNA_ARRAY_RENDERER_H
#define ANTENNA_ARRAY_RENDERER_H

#include <QWidget>

class antenna_array_renderer : public QWidget
{
    Q_OBJECT
public:
    explicit antenna_array_renderer(QWidget *parent = nullptr);

    std::vector<std::vector<bool>> array;

    unsigned array_width();
    unsigned array_height();
    bool array_empty();
    void set_array_range(unsigned height, unsigned width);
    void toggle_antenna(int i, int j);

    void array_set_empty();
    void array_set_all();

signals:

    // QWidget interface
protected:
    void paintEvent(QPaintEvent *event);

    // QWidget interface
protected:
    void mousePressEvent(QMouseEvent *event);


    QTransform transform;

    // QWidget interface
protected:
    void mouseMoveEvent(QMouseEvent *event);
private:

    QPoint last_pressed_point;
};

#endif // ANTENNA_ARRAY_RENDERER_H
