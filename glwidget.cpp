#include "glwidget.h"
#include <QtGui>
#include "QDebug"
//#include <GL/freeglut.h>

GLWidget::GLWidget(QWidget *parent) :
    QGLWidget(parent)
{
	m_image = NULL;
	m_image2 = NULL;
	m_u = 500;
	m_v = 500;
	msDown.setX( -1 );
}

void GLWidget::initializeGL()
{
    glClearColor(0.5,0.5,0.5,1);

}

//opengl (0,0) is the left botton corner by the glViewPort definition
void GLWidget::paintGL()
{
	if( m_image != NULL ){
		glViewport(0, 0, m_u, m_v);
		glClear( GL_COLOR_BUFFER_BIT );

		glColor3f(1,0,0);
		glRasterPos2f( -1,-1);
		glDrawPixels( m_u, m_v, GL_RGBA, GL_UNSIGNED_BYTE, m_image );
	}

	if( msDown.x() != -1 ){
		float xDown = (msDown.x()/(float)m_u)*2-1.0;
		float yDown = (msDown.y()/(float)m_v)*2-1.0;
		float xNow = (msNow.x()/(float)m_u)*2-1.0;
		float yNow = (msNow.y()/(float)m_v)*2-1.0;
		glBegin(GL_LINE_LOOP);              
			glColor3f(1.0f, 0.0f, 0.0f); 
			glVertex2f(xDown, yDown);
			glVertex2f(xNow, yDown);
			glVertex2f(xNow, yNow);
			glVertex2f(xDown, yNow);
		glEnd();
	}

	if( m_image2 != NULL ){
		glViewport(0, 200, m_u, m_v);

		glColor3f(1,0,0);
		glRasterPos2f( -1,-1);
		glDrawPixels( m_u, m_v, GL_RGBA, GL_UNSIGNED_BYTE, m_image2 );
	}
}

void GLWidget::resizeGL( int w, int h )
{
	
}

void GLWidget::setImagePointer( unsigned char* image )
{
    m_image = image;
}

void GLWidget::setImagePointer2( unsigned char* image2 )
{
    m_image2 = image2;
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
    msDown = event->pos();
	msDown.setY( GLWIDGET_HEIGHT - msDown.y() );
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    msNow = event->pos();
	msNow.setY( GLWIDGET_HEIGHT - msNow.y() );
	repaint();
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
    msNow = event->pos();
	msNow.setY( GLWIDGET_HEIGHT - msNow.y() );
	if( msNow.x() == msDown.x() && msNow.y() == msDown.y() )
		msDown.setX( -1 );
	repaint();
}

void GLWidget::setImageSize( int w, int h )
{
	m_u = w;
	m_v = h;
}

void GLWidget::getMsSltInfo( int &msDownX, int &msDownY, int &msRlsX, int &msRlsY )
{
	msDownX = msDown.x();
	msDownY = msDown.y();
	msRlsX = msNow.x();
	msRlsY = msNow.y();
}