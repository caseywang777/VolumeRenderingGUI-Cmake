#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QtGui>
#include <QtCore>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

protected:
    void paintEvent(QPaintEvent *e);
    void mouseDoubleClickEvent(QMouseEvent *me);
    void mouseMoveEvent(QMouseEvent *me);
    void mousePressEvent(QMouseEvent *me);
    void mouseReleaseEvent( QMouseEvent *me );

private slots:
    void on_loadFileButton_clicked();


    void on_horizontalSlider0_sliderMoved(int position);

    void on_horizontalSlider1_sliderMoved(int position);

    void on_horizontalSlider2_sliderMoved(int position);

    void on_horizontalSlider3_sliderMoved(int position);

    void on_horizontalSlider4_sliderMoved(int position);

    void on_horizontalSlider5_sliderMoved(int position);

    void on_horizontalSlider6_sliderMoved(int position);

    void on_horizontalSlider7_sliderMoved(int position);

    void on_horizontalSlider8_sliderMoved(int position);

    void on_horizontalSlider9_sliderMoved(int position);


    void on_renderingButton_clicked();

    void on_startButton_clicked();

    void on_exitButton_clicked();

    void on_pushButtonDisQuery_clicked();

    void on_checkBoxDisQuery_clicked();

    void on_testButton_clicked();

private:
    Ui::MainWindow *ui;
	void outputMessage( QString opt );
	void computeVisibilityHistogram( float* monoVisHist, float* opaTF, int Kgroups, int Depth, int Bins );
};

#endif // MAINWINDOW_H
