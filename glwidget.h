#ifndef GLWIDGET_H
#define GLWIDGET_H

#define GLWIDGET_WIDTH 530
#define GLWIDGET_HEIGHT 520

#include <QGLWidget>

class GLWidget : public QGLWidget
{
    Q_OBJECT
public:
    explicit GLWidget(QWidget *parent = 0);

    void initializeGL();
    void paintGL();
    void resizeGL( int w, int h );

	//the setting in QmouseEvent, the left-upper coner is (0,0)
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);

	void setImagePointer( unsigned char* image );
	void setImagePointer2( unsigned char* image2 );
	void setImageSize( int w, int h );

	void getMsSltInfo( int &msDownX, int &msDownY, int &msRlsX, int &msRlsY );
	

private:
	unsigned char* m_image;
	unsigned char* m_image2;
	int m_u;
	int m_v;
	QPoint msDown;
	QPoint msNow;
};

#endif // GLWIDGET_H
