/***************************************************************************/
/* Copyright(C) 2021


The authors of

Reliable Feature-Line Driven Quad-Remeshing
Siggraph 2021


 All rights reserved.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
****************************************************************************/

/* LZ change
- header QGLWidget -> QOpenGLWidget
- class QGLWidget -> QOpenGLWidget
*/ 

#ifndef GLWIDGET_H
#define GLWIDGET_H


#include <QOpenGLWidget>
#include <wrap/gui/trackball.h>
#include <AntTweakBar.h>
#include <QDir>
#include "bvh.h"

class GLWidget : public QOpenGLWidget
{
    Q_OBJECT

public:
    bool hasToPick;

    GLWidget(QWidget *parent = 0);

    ~GLWidget(){};


    QSize sizeHint() const {
      return QSize(1400, 1200);
  };

protected:
    //void DrawPolyMesh();

    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseDoubleClickEvent ( QMouseEvent * event );
    void keyPressEvent(QKeyEvent *event);
    void keyReleaseEvent(QKeyEvent *event);
    void wheelEvent(QWheelEvent *event);

};

#endif
