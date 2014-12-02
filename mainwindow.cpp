#define CommFILE "C:\\GravityLabDataSet\\VolumeRenderingGUIComm.txt" //communication file

#define BINS 128

#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "stdio.h"
#include "stdlib.h"

#include "QFileDialog"
#include "QMessageBox"
#include <QKeyEvent>
#include <Qt>
#include <QDebug>
#include <QMessageBox>
#include <QFont>
#include <vector>
#include <QInputDialog>

#include "rendererThread.h"
#include "glwidget.h"
#include "QColorDialog"

#include "Image.h"


//extern "C" double runCudaKernel( float* visHist, int K, int D, int B, float* otf, int* rawHistogramRays );

GLWidget* glWg;

QString g_filename;
#define NumCtrlPt 9
#define NumSaCtrlPt 9
float g_ctrlPtR[NumCtrlPt];
float g_ctrlPtG[NumCtrlPt];
float g_ctrlPtB[NumCtrlPt];
float g_ctrlPtA[NumCtrlPt];
float g_ctrlPtV[NumCtrlPt];
int g_tfpX=50;//transfer function point X
int g_tfpY=700;
int g_tfpW=500;
int g_tfpH=-200;

int g_sltPoint = -1;
float g_distribution[BINS];
float g_distribution2[BINS];

rendererThread* rt;

Image* img = NULL;
QString opt = "";
int msgCounter = 1;
vector<Ray*> m_rays;

//for visibility histogram estimation
float ctfr[128];
float ctfg[128];
float ctfb[128];
float ctfa[128];

int K = 128;
int D = 126;
int B = 128;

int m_dim1 = 512;
int m_dim2 = 126;
int m_dim3 = 126;

float *rawData;
unsigned char* m_color;

float trueVisHist[128];

//char rayClusterFile[500] = "C:\\Users\\VisKC\\Dropbox\\00_ImportanceNow\\Research\\WindowsCode\\VisibilityHistogram\\isabel_k128_d100_b128.bin";
//char rayRawFile[500] = "C:\\Users\\VisKC\\Dropbox\\00_ImportanceNow\\Research\\WindowsCode\\VisibilityHistogram\\isabel_500_500_100.bin";
char rayClusterFile[500] = "C:\\Users\\VisKC\\Dropbox\\00_ImportanceNow\\Research\\WindowsCode\\VisibilityHistogram\\plume_k128_d126_b128.bin";
char rayRawFile[500] = "C:\\Users\\VisKC\\Dropbox\\00_ImportanceNow\\Research\\WindowsCode\\VisibilityHistogram\\plume_512_126_126.bin";
int* rawHistogramRays;

//simulated Annealing
float temperature = 1000;
float innerIter = 10;
bool saFirstFlag = true;
int l = 0;
int lmax = 1000;
float tfV[NumSaCtrlPt];
float tfA[NumSaCtrlPt];
float targetVH[128];
float opcTF[128];
float saDensity = 1.0; //volume density in sa, to adjust transfer function
float randomBound = 1.00;	//perturbe TF bound
float cost1 = 10000;



void loadClusterVisibilityHistogramData()
{
	//test from matlab
	int size = K * D * B;
	int* arr = (int*) malloc( sizeof( int) * size );
	FILE* fp  = fopen( rayClusterFile, "rb" );
	fread(arr, sizeof(int), size, fp);
	fclose( fp );
	rawHistogramRays = arr;	//rawHistogramRay keeps the data to load into gpu



	//load cluster ray-histogram to evaluate visibility histogram
	int segments = 0;
	float* iof = (float*) malloc( sizeof(float) * D );
	float* monoVisHist = (float*) malloc( sizeof(float) * B );
	memset( monoVisHist, 0, sizeof(float)*B );
	Histogram hist(B);
	float* input = (float*)malloc( sizeof(float) * D * B );
	for( int k = 0; k< K; k++ ){
		int base = k * D * B;
		int slcSample = 0;
		for( int i = base; i < base + B; i++ )
			slcSample += arr[i];
		for( int i = base; i < base + D * B; i++ ){
			input[i - base] = arr[i];
		}

		Ray* r = new Ray( input, D, B, slcSample );
		m_rays.push_back( r );
		float rayError, rayCnt;
		int seg = m_rays[k]->produceMonoRepresentation( 0.001, 1, rayError, rayCnt );
	}

	//load raw data
	rawData = (float*) malloc( sizeof(float) * m_dim1 *m_dim2 *m_dim3 );
	fp = fopen( rayRawFile, "rb" );
	fread( rawData, sizeof(float), m_dim1 *m_dim2 *m_dim3, fp );
	fclose( fp );

	 m_color = (unsigned char*) malloc( m_dim1 *m_dim2 * 4 * sizeof(unsigned char) );

	 glWg->setImageSize( m_dim2, m_dim1 );
}


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    glWg = new GLWidget();
    ui->gridLayout->addWidget(glWg,0,0,1,1);

    //init transfer function control point
    float v=0;
    for( int i=0; i<NumCtrlPt; i++ ){
        g_ctrlPtV[i] = v;
        v+= 1/(float)(NumCtrlPt-1);
    }
    g_ctrlPtR[0] = 0.00; g_ctrlPtG[0] = 0.00; g_ctrlPtB[0] = 0.50; g_ctrlPtA[0] = 1.00;
    g_ctrlPtR[1] = 0.00; g_ctrlPtG[1] = 0.00; g_ctrlPtB[1] = 1.00; g_ctrlPtA[1] = 1.00;
    g_ctrlPtR[2] = 0.00; g_ctrlPtG[2] = 0.50; g_ctrlPtB[2] = 1.00; g_ctrlPtA[2] = 1.00;
    g_ctrlPtR[3] = 0.00; g_ctrlPtG[3] = 1.00; g_ctrlPtB[3] = 1.00; g_ctrlPtA[3] = 1.00;
    g_ctrlPtR[4] = 0.50; g_ctrlPtG[4] = 1.00; g_ctrlPtB[4] = 0.50; g_ctrlPtA[4] = 1.00;
    g_ctrlPtR[5] = 1.00; g_ctrlPtG[5] = 1.00; g_ctrlPtB[5] = 0.00; g_ctrlPtA[5] = 1.00;
    g_ctrlPtR[6] = 1.00; g_ctrlPtG[6] = 0.50; g_ctrlPtB[6] = 0.00; g_ctrlPtA[6] = 1.00;
    g_ctrlPtR[7] = 1.00; g_ctrlPtG[7] = 0.00; g_ctrlPtB[7] = 0.00; g_ctrlPtA[7] = 1.00;
    g_ctrlPtR[8] = 0.50; g_ctrlPtG[8] = 0.00; g_ctrlPtB[8] = 0.00; g_ctrlPtA[8] = 1.00;
    //end - init transfer function

	loadClusterVisibilityHistogramData();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_startButton_clicked()
{
    on_renderingButton_clicked();
	
	//load inc/dec file
	QByteArray byteArray = g_filename.toLocal8Bit();
    char *c = byteArray.data();
	img = new Image( c , ui->resxLineEdit->text().toInt()-1, 
							ui->resyLineEdit->text().toInt()-1,  
							ui->reszLineEdit->text().toInt(),
							BINS,
							ui->samplePerUnitLineEdit->text().toInt() );
	outputMessage( "Loaded" );
	glWg->setImageSize( ui->resxLineEdit->text().toInt()-1, ui->resyLineEdit->text().toInt()-1 );
	
}

void MainWindow::on_exitButton_clicked()
{
    glWg->updateGL();
    //rt->terminate();
    //exit(0);
}

void MainWindow::on_loadFileButton_clicked()
{
    g_filename = QFileDialog::getOpenFileName(
                    this,
                    tr("Open File"),
                    "C:\\GravityLabDataSet\\myData\\",
                    "All files (*.*)"
                );

    ui->fileNameLabel->setText( g_filename );
}


void MainWindow::on_horizontalSlider0_sliderMoved(int position)
{
    ui->lcdNumber0->display( position );
}

void MainWindow::on_horizontalSlider1_sliderMoved(int position)
{
    ui->lcdNumber1->display( position );
}

void MainWindow::on_horizontalSlider2_sliderMoved(int position)
{
    ui->lcdNumber2->display( position );
}

void MainWindow::on_horizontalSlider3_sliderMoved(int position)
{
    ui->lcdNumber3->display( position );
}

void MainWindow::on_horizontalSlider4_sliderMoved(int position)
{
    ui->lcdNumber4->display( position );
}

void MainWindow::on_horizontalSlider5_sliderMoved(int position)
{
    ui->lcdNumber5->display( position );
}

void MainWindow::on_horizontalSlider6_sliderMoved(int position)
{
    ui->lcdNumber6->display( position );
}

void MainWindow::on_horizontalSlider7_sliderMoved(int position)
{
    ui->lcdNumber7->display( position );
}

void MainWindow::on_horizontalSlider8_sliderMoved(int position)
{
    ui->lcdNumber8->display( position );
}

void MainWindow::on_horizontalSlider9_sliderMoved(int position)
{
    ui->lcdNumber9->display( position/1000.0 );

	//only use emphsize, run this block
	if( ui->renderingRadioButton2->isChecked() ){
		on_renderingButton_clicked();
	}
}

void local2GlocalCoordinate( float &x, float &y, float ix, float iy )
{
    x = ix*g_tfpW + g_tfpX;
    y = iy*g_tfpH + g_tfpY;
}

void global2LocalCoordinate( float &x, float &y, float ix, float iy )
{
    x = (ix - g_tfpX)/(float)g_tfpW;
    y = (iy - g_tfpY)/(float)g_tfpH;
}

float interpolation1D( float p, float ps, float pt, float vs, float vt)
{
    float ratio = (p-ps)/(pt-ps);
    return ((vt-vs)*ratio + vs);
}

void MainWindow::paintEvent(QPaintEvent *e)
{
    QPainter painter(this);

    QPen pointPen;
    pointPen.setColor(Qt::red);
    pointPen.setWidth(20);

    QPen linePen;
    linePen.setColor(Qt::blue);
    linePen.setWidth(10);

    QPen bluePen;
    bluePen.setColor(Qt::blue);
    bluePen.setWidth(2);

    QPen blackPen5;
    blackPen5.setColor(Qt::black);
    blackPen5.setWidth(5);

    QPen blackPen1;
    blackPen1.setColor(Qt::black);
    blackPen1.setWidth(1);

    QPen textPen;
    textPen.setColor(Qt::black);
    textPen.setWidth(4);

    QFont font =painter.font();
    font.setPointSize(10);

    //draw axis line
    painter.setPen(blackPen5);
    painter.drawLine( QPointF((float)g_tfpX, (float)g_tfpY), QPointF((float)(g_tfpX+g_tfpW), (float)(g_tfpY)) );
    painter.drawLine( QPointF((float)g_tfpX, (float)g_tfpY), QPointF((float)(g_tfpX), (float)(g_tfpY+g_tfpH)) );
    for( int i=0; i<NumCtrlPt; i++ ){
        painter.setPen(blackPen1);
        painter.drawLine( QPointF((float)g_tfpX+g_tfpW*(i/(float)(NumCtrlPt-1)), (float)g_tfpY), QPointF((float)(g_tfpX+g_tfpW*(i/(float)(NumCtrlPt-1))), (float)(g_tfpY+g_tfpH)) );
        painter.setFont(font);
        QString str;
        str.sprintf("%4.3f", i/(float)(NumCtrlPt-1));
        painter.drawText( QPoint( (float)g_tfpX+g_tfpW*(i/(float)(NumCtrlPt-1)) -15 , (float)g_tfpY+15 ), str);
    }
    for( int i=0; i<5; i++ ){
        painter.setPen(blackPen1);
        painter.drawLine( QPointF((float)g_tfpX, (float)g_tfpY+g_tfpH*(i/(float)(5-1))), QPointF((float)(g_tfpX+g_tfpW), (float)(g_tfpY+g_tfpH*(i/(float)(5-1)))) );
        painter.setFont(font);
        QString str;
        str.sprintf("%3.2f", i/(float)(5-1));
        painter.drawText( QPoint( (float)(g_tfpX)-30 , (float)g_tfpY+g_tfpH*(i/(float)(5-1)) +5 ), str);
    }


    //draw control point and line
    painter.setPen(bluePen);
    QPointF points[NumCtrlPt];
    for( int i=0; i<NumCtrlPt; i++ ){
        float x,y;
        local2GlocalCoordinate( x, y, g_ctrlPtV[i], g_ctrlPtA[i] );
        points[i].setX( x );
        points[i].setY( y );
    }


    for( int i=0; i<NumCtrlPt; i++ ){
        //draw Line
        float step = (1.0/50.0);

        if( i!=NumCtrlPt-1){
            for( float j=step; j<=1; j+=step){
                linePen.setColor( QColor((g_ctrlPtR[i]+(g_ctrlPtR[i+1]-g_ctrlPtR[i])*j)*255,
                                         (g_ctrlPtG[i]+(g_ctrlPtG[i+1]-g_ctrlPtG[i])*j)*255,
                                         (g_ctrlPtB[i]+(g_ctrlPtB[i+1]-g_ctrlPtB[i])*j)*255  ) );
                linePen.setWidth(10);
                painter.setPen(linePen);
                painter.drawLine( QPointF(points[i].x()+(points[i+1].x()-points[i].x())*(j-step), points[i].y()+(points[i+1].y()-points[i].y())*(j-step))
                        ,QPointF(points[i].x()+(points[i+1].x()-points[i].x())*(j), points[i].y()+(points[i+1].y()-points[i].y())*(j)) );
            }
        }
    }

    for( int i=0; i<NumCtrlPt; i++ ){
        //draw Point
        pointPen.setColor( QColor(g_ctrlPtR[i]*255, g_ctrlPtG[i]*255, g_ctrlPtB[i]*255 ) );
        painter.setPen(pointPen);
        painter.drawPoint( QPointF( points[i].x(), points[i].y() ));
    }

    if( g_sltPoint != -1 ){
        painter.setPen(blackPen1);
        painter.setFont(font);
        QString str;
        str.sprintf("(R:%3.2f\n G:%3.2f\n B:%3.2f\n A:%3.2f\n V:%4.3f)",
                    g_ctrlPtR[g_sltPoint],
                    g_ctrlPtG[g_sltPoint],
                    g_ctrlPtB[g_sltPoint],
                    g_ctrlPtA[g_sltPoint],
                    g_ctrlPtV[g_sltPoint]);
        painter.drawText( QPoint( points[g_sltPoint].x()-50 , points[g_sltPoint].y()-10 ), str);
    }

	//update control point
	float tf[1000];
	int tfCnt=0;
	for( int i=0; i<NumCtrlPt; i++ ){
		float r, g, b, a, v;
		r = g_ctrlPtR[i];
		g = g_ctrlPtG[i];
		b = g_ctrlPtB[i];
		a = g_ctrlPtA[i];
		v = g_ctrlPtV[i];
		tf[tfCnt++] = r;
		tf[tfCnt++] = g;
		tf[tfCnt++] = b;
		tf[tfCnt++] = a;
		tf[tfCnt++] = v;
	}
	if( img!=NULL ){
		img->makeTransferFunction(NumCtrlPt, tf, 0);
		img->makeTransferFunction(NumCtrlPt, tf, 1);
	}


    //draw distribution
    int bins = 128;
    float prob=0;

    //find the max probablity to make the PDF normalize by the max probability
    float maxProb = 0;
    for( int i=0; i<bins-1; i++ ){
        if( maxProb < g_distribution[i] ){
            maxProb = g_distribution[i];
        }
    }
    if( maxProb == 0 )maxProb=1;

    for( int i=0; i<bins-1; i++ ){
        linePen.setColor(QColor(0,255,0));
        linePen.setWidth(2);
        painter.setPen(linePen);
        float x1, y1;
        float x2, y2;
        if(  ui->checkBoxDisQuery->isChecked() ){
            prob+=g_distribution[i];
            local2GlocalCoordinate( x1, y1, i/(float)(bins-1), prob );
            local2GlocalCoordinate( x2, y2, (i+1)/(float)(bins-1), prob+g_distribution[i+1] );
        }else{
            local2GlocalCoordinate( x1, y1, i/(float)(bins-1), g_distribution[i]/maxProb );
            local2GlocalCoordinate( x2, y2, (i+1)/(float)(bins-1), g_distribution[i+1]/maxProb );
        }
        painter.drawLine( QPointF(x1,y1),QPointF(x2,y2) );
    }

    //draw distribution2
    bins = 128;
    prob=0;

    //find the max probablity to make the PDF normalize by the max probability
    maxProb = 0;
    for( int i=0; i<bins-1; i++ ){
        if( maxProb < g_distribution2[i] ){
            maxProb = g_distribution2[i];
        }
    }
    if( maxProb == 0 )maxProb=1;

	//calculate transfer function
	int cpNow = 0;
	for( int i=0; i<128; i++ ){
		float vxl = i/127.0;
		//current vxl > vxl in current control point, move to next control point
		if( vxl > g_ctrlPtV[cpNow+1] )cpNow++;
		ctfr[i] = interpolation1D(vxl, g_ctrlPtV[cpNow], g_ctrlPtV[cpNow+1],
										g_ctrlPtR[cpNow], g_ctrlPtR[cpNow+1]);//r
		ctfg[i] = interpolation1D(vxl, g_ctrlPtV[cpNow], g_ctrlPtV[cpNow+1],
										g_ctrlPtG[cpNow], g_ctrlPtG[cpNow+1]);//g
		ctfb[i] = interpolation1D(vxl, g_ctrlPtV[cpNow], g_ctrlPtV[cpNow+1],
										g_ctrlPtB[cpNow], g_ctrlPtB[cpNow+1]);//b
		ctfa[i] = interpolation1D(vxl, g_ctrlPtV[cpNow], g_ctrlPtV[cpNow+1],
										g_ctrlPtA[cpNow], g_ctrlPtA[cpNow+1]) * ui->lineEdit2->text().toFloat();//a
	}

	////calculate visibility histogram
	//float monoVisHist[128];
	//float dataFreqHist[128];
	//memset( monoVisHist, 0 , sizeof(float)*128 );
	//memset( dataFreqHist, 0 , sizeof(float)*128 );

	//clock_t start_time, end_time;
	//float total_time = 0;
	//start_time = clock(); 

	//#pragma omp parallel for
	//for( int k = 0; k< K; k++ ){
	//	float alpha = 0;
	//	float iof[128];
	//	Histogram hist( bins );
	//	Histogram dataHist( bins );
	//	for( int d = 0; d < D; d++ ){
	//		iof[d] = alpha;
	//		m_rays[k]->getHistByIdm( &hist, d, d );
	//		float sumAlpha = 0;
	//		float sumCount = 0;
	//		for( int b = 0; b < B; b ++ ){
	//			sumCount += hist.getBin( b );
	//			sumAlpha += hist.getBin( b ) * ctfa[b];
	//		}
	//		float a = sumAlpha / sumCount;
	//		alpha += ( 1.0 - alpha ) * a;
	//	}

	//	#pragma omp critical
	//	{
	//		for( int d = 0; d < D; d++ ){
	//			m_rays[k]->getHistByIdm( &hist, d, d );
	//			for( int b = 0; b < B; b ++ ){
	//				dataFreqHist[b] += hist.getBin( b );
	//				monoVisHist[b] += (1.0 - iof[d]) * hist.getBin( b ) * ctfa[b];
	//			}
	//		}
	//	}
	//}

	//end_time = clock();
	//total_time = (float)(end_time - start_time)/CLOCKS_PER_SEC;
	//FILE* fp = fopen( "timetest.txt", "w" );
	//fprintf(fp, "%f \n", total_time);
	//fclose( fp );


	//float maxMono = 0;
	//float maxTrue = 0;

	//float trueVisHistShow[128];

	//for( int i=0; i<128; i++ )
	//	trueVisHistShow[i] = trueVisHist[i];

	//for( int i=0; i<128; i++ ){
	//	if( dataFreqHist[i] != 0 ){
	//		monoVisHist[i] /= dataFreqHist[i];
	//		trueVisHistShow[i] /= dataFreqHist[i];
	//	}
	//}
	//for( int i=0; i<128; i++ ){
	//	if( monoVisHist[i] > maxMono )
	//		maxMono = monoVisHist[i];
	//	if( trueVisHistShow[i] > maxTrue )
	//		maxTrue = trueVisHistShow[i];
	//}
	//if( maxMono == 0 ) maxMono = 1;
	//if( maxTrue == 0 ) maxTrue = 1;
	//for( int i=0; i<128; i++ ){
	//	monoVisHist[i] /= maxMono;
	//	trueVisHistShow[i] /= maxTrue;
	//}

	////use function
	//computeVisibilityHistogram( monoVisHist, ctfa, K, D, B );


	float monoVisHist[128];//it is transfer function now
	intepolateTransferFunctionFromControlPoint( monoVisHist, 128, tfV, tfA );
	float trueVisHistShow[128];
	computeVisibilityHistogram( trueVisHistShow, monoVisHist, K, D, B );

	for( int i=0; i<bins-1; i++ ){
        linePen.setColor(QColor(0,0,255));
        linePen.setWidth(2);
        painter.setPen(linePen);
        float x1, y1;
        float x2, y2;
		local2GlocalCoordinate( x1, y1, i/(float)(bins-1), monoVisHist[i] );
        local2GlocalCoordinate( x2, y2, (i+1)/(float)(bins-1), monoVisHist[i+1] );
		painter.drawLine( QPointF(x1,y1),QPointF(x2,y2) );
        if(  ui->checkBoxDisQuery->isChecked() ){
            /*prob+=g_distribution2[i];
            local2GlocalCoordinate( x1, y1, i/(float)(bins-1), prob );
            local2GlocalCoordinate( x2, y2, (i+1)/(float)(bins-1), prob+g_distribution2[i+1] );*/
			local2GlocalCoordinate( x1, y1, i/(float)(bins-1), trueVisHistShow[i] );
            local2GlocalCoordinate( x2, y2, (i+1)/(float)(bins-1), trueVisHistShow[i+1] );
        }
    }

	if(  ui->checkBoxDisQuery->isChecked() ){
		for( int i=0; i<bins-1; i++ ){
			linePen.setColor(QColor(255,0,0));
			linePen.setWidth(2);
			painter.setPen(linePen);
			float x1, y1;
			float x2, y2;
	
			local2GlocalCoordinate( x1, y1, i/(float)(bins-1), trueVisHistShow[i] );
			local2GlocalCoordinate( x2, y2, (i+1)/(float)(bins-1), trueVisHistShow[i+1] );
			painter.drawLine( QPointF(x1,y1),QPointF(x2,y2) );

        }
    }


    //keeping updating status
    if( ui->checkBoxkeepingUpdata->isChecked() ){
        on_renderingButton_clicked();
    }

    update();
}

void MainWindow::outputMessage( QString msg )
{
	QString counterStr= "(" + QString::number( msgCounter ++ ) + ") ";
	opt = counterStr + msg + ("\n") + opt;
	ui->outputDisplay->setText( opt );
}

int matchPoint( float x, float y )
{
    float threshold = 8;

    int i = 0;
    for( i=0; i<NumCtrlPt; i++ ){
        float gx,gy;
        local2GlocalCoordinate( gx, gy, g_ctrlPtV[i], g_ctrlPtA[i] );
        float xDis = fabs( x - gx );
        float yDis = fabs( y - gy );
        float dis = sqrt( xDis*xDis + yDis*yDis );
        if( dis < threshold )return i;
    }

    return -1;
}

void MainWindow::mouseDoubleClickEvent(QMouseEvent *me)
{
    int x = me->x();
    int y = me->y();

    if( me->button() == Qt::RightButton ){
        int sltPoint = matchPoint(x,y);
        if( sltPoint != -1 ){
            qDebug()<<QString::number(x) + "   " +QString::number(y);
            QColor color = QColorDialog::getColor( QColor(g_ctrlPtR[sltPoint]*255, g_ctrlPtG[sltPoint]*255, g_ctrlPtB[sltPoint]*255), this );
            int r, g, b, a;
            color.getRgb( &r, &g, &b, &a );
            g_ctrlPtR[sltPoint] = r/255.0;
            g_ctrlPtG[sltPoint] = g/255.0;
            g_ctrlPtB[sltPoint] = b/255.0;
        }
    }

    repaint();
}


void MainWindow::mousePressEvent(QMouseEvent *me)
{
    int x = me->x();
    int y = me->y();

    g_sltPoint = matchPoint(x,y);
    repaint();
}

void MainWindow::mouseReleaseEvent(QMouseEvent *me)
{
    int x = me->x();
    int y = me->y();


    g_sltPoint = -1;
    repaint();
}

void MainWindow::mouseMoveEvent(QMouseEvent *me)
{
    int x = me->x();
    int y = me->y();

    if( g_sltPoint!=-1 ){
        float lx, ly;
        global2LocalCoordinate( lx, ly, x, y );
        if( ly < 0 )ly=0;
        if( ly > 1 )ly=1;

        if( g_sltPoint == 0 )lx=0;
        if( g_sltPoint == NumCtrlPt-1 )lx=1;

        if( g_sltPoint!=NumCtrlPt-1 && lx >= g_ctrlPtV[g_sltPoint+1] )lx = g_ctrlPtV[g_sltPoint+1];
        if( g_sltPoint!=0 && lx <= g_ctrlPtV[g_sltPoint-1] )lx = g_ctrlPtV[g_sltPoint-1];

        g_ctrlPtV[g_sltPoint] = lx;
        g_ctrlPtA[g_sltPoint] = ly;
    }

    repaint();
}

void MainWindow::computeVisibilityHistogram( float* monoVisHist, float* opaTF, int Kgroups, int Depth, int Bins )
{
	//calculate visibility histogram
	float dataFreqHist[128];
	memset( monoVisHist, 0 , sizeof(float)*128 );
	memset( dataFreqHist, 0 , sizeof(float)*128 );

	#pragma omp parallel for
	for( int k = 0; k< Kgroups; k++ ){
		float alpha = 0;
		float iof[128];
		Histogram hist( 128 );
		Histogram dataHist( 128 );
		for( int d = 0; d < Depth; d++ ){
			iof[d] = alpha;
			m_rays[k]->getHistByIdm( &hist, d, d );
			float sumAlpha = 0;
			float sumCount = 0;
			for( int b = 0; b < B; b ++ ){
				sumCount += hist.getBin( b );
				sumAlpha += hist.getBin( b ) * opaTF[b];
			}
			float a = sumAlpha / sumCount;
			alpha += ( 1.0 - alpha ) * a;
		}

		#pragma omp critical
		{
			for( int d = 0; d < Depth; d++ ){
				m_rays[k]->getHistByIdm( &hist, d, d );
				for( int b = 0; b < Bins; b ++ ){
					dataFreqHist[b] += hist.getBin( b );
					monoVisHist[b] += (1.0 - iof[d]) * hist.getBin( b ) * opaTF[b];
				}
			}
		}
	}


	float maxMono = 0;

	for( int i=0; i<128; i++ ){
		if( dataFreqHist[i] != 0 ){
			monoVisHist[i] /= dataFreqHist[i];
		}
	}
	for( int i=0; i<128; i++ ){
		if( monoVisHist[i] > maxMono )
			maxMono = monoVisHist[i];
	}
	if( maxMono == 0 ) maxMono = 1;
	for( int i=0; i<128; i++ ){
		monoVisHist[i] /= maxMono;
	}
}

void MainWindow::on_renderingButton_clicked()
{
	//if( img == NULL ) return;	//if img has not been constructed, do not run this function

	//choose rendering technique
    int renderingTech = 0;
    if( ui->renderingRadioButton0->isChecked() )renderingTech = 0;
    else if( ui->renderingRadioButton1->isChecked() )renderingTech = 1;
    else if( ui->renderingRadioButton2->isChecked() )renderingTech = 2;
    else if( ui->renderingRadioButton3->isChecked() )renderingTech = 3;
    else if( ui->renderingRadioButton4->isChecked() )renderingTech = 4;
    else if( ui->renderingRadioButton5->isChecked() )renderingTech = 5;
    else if( ui->renderingRadioButton6->isChecked() )renderingTech = 6;
    else if( ui->renderingRadioButton7->isChecked() )renderingTech = 7;
    else if( ui->renderingRadioButton8->isChecked() )renderingTech = 8;
    else if( ui->renderingRadioButton9->isChecked() )renderingTech = 9;

	//slider bar values
	int sldBar0 = ui->horizontalSlider0->value();
	int sldBar1 = ui->horizontalSlider1->value();
	int sldBar2 = ui->horizontalSlider2->value();
	int sldBar3 = ui->horizontalSlider3->value();
	int sldBar4 = ui->horizontalSlider4->value();
	int sldBar5 = ui->horizontalSlider5->value();
	int sldBar6 = ui->horizontalSlider6->value();
	int sldBar7 = ui->horizontalSlider7->value();
	int sldBar8 = ui->horizontalSlider8->value();
	int sldBar9 = ui->horizontalSlider9->value();

	//line Edit value
	float lineEdit0 = ui->lineEdit0->text().toFloat();
	float lineEdit1 = ui->lineEdit1->text().toFloat();
	float lineEdit2 = ui->lineEdit2->text().toFloat();
	float lineEdit3 = ui->lineEdit3->text().toFloat();
	float lineEdit4 = ui->lineEdit4->text().toFloat();
	float lineEdit5 = ui->lineEdit5->text().toFloat();
	float lineEdit6 = ui->lineEdit6->text().toFloat();
	float lineEdit7 = ui->lineEdit7->text().toFloat();
	float lineEdit8 = ui->lineEdit8->text().toFloat();
	float lineEdit9 = ui->lineEdit9->text().toFloat();

	/*img->setUIParameters( sldBar0, sldBar1, sldBar2, sldBar3, sldBar4, 
						  sldBar5, sldBar6, sldBar7, sldBar8, sldBar9, 
						  lineEdit0, lineEdit1, lineEdit2, lineEdit3, lineEdit4, 
						  lineEdit5, lineEdit6, lineEdit7, lineEdit8, lineEdit9 );*/

	if( renderingTech == 0 ){

	}else if( renderingTech == 1 ){
		memset( trueVisHist, 0 , sizeof(float)*128 );
		int m_u = m_dim1;
		int m_v = m_dim2;
		int m_d = m_dim3;
		float bgColor = 0;

		clock_t start_time, end_time;
		float total_time = 0;
		start_time = clock(); 

		#pragma omp parallel for
        for( int v = 0; v< m_v; v+=1 ){
            for( int u = 0; u< m_u; u+=1 ){
				float cr = 0, cg = 0, cb = 0, ca = 0;

				int base = ( v * m_u + u ) * m_d;
				float visHist[128];
				memset( visHist, 0 , sizeof(float)*128 );
				for( int d = 0; d< m_d; d++ ){
                    float vxl = rawData[ base + d ];
					int vxlIdx = (int)(vxl * 127.0);

					float r, g, b, a;
					r = ctfr[vxlIdx];
					g = ctfg[vxlIdx];
					b = ctfb[vxlIdx];
					//a = ctfa[vxlIdx];
					a = opcTF[vxlIdx];
                    
                    cr += (r*a) * (1.0 - ca);
                    cg += (g*a) * (1.0 - ca);
                    cb += (b*a) * (1.0 - ca);
                    ca += a * (1.0 - ca);
					visHist[vxlIdx] += a * (1.0 - ca);
                    //if( ca > 0.95 )break;
                }
				cr += (bgColor) * (1.0 - ca);
                cg += (bgColor) * (1.0 - ca);
                cb += (bgColor) * (1.0 - ca);
                unsigned char optR = (unsigned char)(int)floor(cr * 255.0);
                unsigned char optG = (unsigned char)(int)floor(cg * 255.0);
                unsigned char optB = (unsigned char)(int)floor(cb * 255.0);
				
				m_color[(v*m_u+u)*4+0] = optR;
                m_color[(v*m_u+u)*4+1] = optG;
                m_color[(v*m_u+u)*4+2] = optB;
                //printf("%d %d %d\n", optR, optG, optB);

				#pragma omp critical
				{
					for( int i=0; i<128; i++ )
						trueVisHist[i] += visHist[i];
				}
            }//end u
        }
		end_time = clock();
		total_time = (float)(end_time - start_time)/CLOCKS_PER_SEC;
		FILE* fp = fopen( "timetest2.txt", "w" );
		fprintf(fp, "%f \n", total_time);
		fclose( fp );
		glWg->setImagePointer( m_color );
		glWg->updateGL();
	}else if( renderingTech == 2 ){//simulated annealing
		float cost;
		float tmpTfA[NumSaCtrlPt];
		float tmpOpcTF[128];

		//if this is first run of SA, some init for SA
		if( saFirstFlag ){
			//init tfV, tfA
			//for( int i= 0; i < NumCtrlPt; i ++ ){
			//	tfV[i] = g_ctrlPtV[i];
			//	tfA[i] = g_ctrlPtA[i] * saDensity; //0.x is the density of volume
			//}
			for( int i= 0; i < NumSaCtrlPt; i ++ ){
				tfV[i] = ( 1 / (float)(NumSaCtrlPt - 1) ) * i;
				tfA[i] = 0.5;//g_ctrlPtA[i] * saDensity; //0.x is the density of volume
			}
			cost1 = 100000;
			//compute the target VH
			intepolateTransferFunctionFromControlPoint( targetVH, 128, g_ctrlPtV, g_ctrlPtA );
			saFirstFlag = false;
		}

		//the following is one run of SA, some monte carlo run
		for( int i = 0; i < innerIter; i ++ ){
			//perturbe the trasfer function (pertube one control point)
			int randomCP = rand() % NumSaCtrlPt;
			float randomAlpha = ( rand() / (float)RAND_MAX ) * randomBound;
			for( int i = 0; i < NumSaCtrlPt; i ++ )tmpTfA[i] = tfA[i];
			tmpTfA[ randomCP ] = randomAlpha;
			if( targetVH[(int) ((( 1 / (float)(NumSaCtrlPt - 1) ) * randomCP ) * 127.0)] < 0.001 )tmpTfA[ randomCP ] = 0;

			//the cost of perturbe transfer function
			intepolateTransferFunctionFromControlPoint( tmpOpcTF, 128, tfV, tmpTfA );
			float tmpVisHist[128];
			computeVisibilityHistogram( tmpVisHist, tmpOpcTF, K, D, B );
			float cost2 = computeEMD( 128, tmpVisHist, targetVH );

			//cost between, perturbe TF and current TF
			float e = cost2 - cost1;
			if( e < 0 || 
				( exp( -e/temperature ) > rand()/(float)RAND_MAX ) ){
					for( int i = 0; i < NumSaCtrlPt; i ++ )
						tfA[i] = tmpTfA[i];

				//the cost of new transfer function
				intepolateTransferFunctionFromControlPoint( opcTF, 128, tfV, tfA );
				float visHist[128];
				computeVisibilityHistogram( visHist, opcTF, K, D, B );
				cost1 = computeEMD( 128, visHist, targetVH );
			}//if
		}//for
		l ++;

		//compute current transfer function cost
		intepolateTransferFunctionFromControlPoint( opcTF, 128, tfV, tfA );
		float visHist[128];
		computeVisibilityHistogram( visHist, opcTF, K, D, B );
		cost = computeEMD( 128, visHist, targetVH );

		//lower the temperature
		temperature *= 0.95;

		//output cost and temperature
		QString str;
		str =  QString::number(cost, 'g', 5) + "  " +  QString::number(temperature, 'g', 5);
		outputMessage( str );

	}else if( renderingTech == 3 ){
		FILE* fp = fopen("time.txt", "a+" );
		clock_t start_time, end_time;
		  float total_time = 0;
		  start_time = clock(); /* mircosecond */
		  img->multiResolutionEntropyEstimationRendering( sldBar8/99.0 );
		  end_time = clock();

		   /* CLOCKS_PER_SEC is defined at time.h */
			total_time = (float)(end_time - start_time)/CLOCKS_PER_SEC;

			fprintf(fp, "Time : %f sec \n", total_time);
			fclose(fp);
		
		glWg->setImagePointer( img->getImgColorPointer() );
		//glWg->setImagePointer2( img->getImgColorPointer2() );
		glWg->updateGL();
		img->outputPPMImage( "img.ppm" );
	}else if( renderingTech == 4 ){
		//img->statisticRenderingEntropy(8);
		//img->entropyBaseOpacityRendering( lineEdit3, 0, lineEdit0, sldBar8/99.0, lineEdit4, lineEdit5 );
		img->entropyBaseSampleRendering(  );
		img->outputPPMImage( "img.ppm" );
		//img->visibilityBaseSampleRendering();
		//img->selfInformationBaseSampleRendering();
		glWg->setImagePointer( img->getImgColorPointer() );
		glWg->updateGL();
	}else if( renderingTech == 5 ){
		img->colorRenderingBySegmentMean(lineEdit3, lineEdit0, sldBar8/99.0, sldBar3, lineEdit9, lineEdit8 );		//use for feature extraction
		//img->colorRenderingByHistogramInfomation();
		//img->visibilityBaseMultiResolutionRendering(lineEdit3, lineEdit0, sldBar8/99.0, sldBar3, lineEdit9, lineEdit8 );
		//img->multiResolutionColorRendering(122);
		glWg->setImagePointer( img->getImgColorPointer() );
		//glWg->setImagePointer2( img->getImgColorPointer2() );
		glWg->updateGL();
		img->outputPPMImage( "img.ppm" );
		img->outputPPMImage2( "img2.ppm" );

		Histogram optVisHist( 128 );
		img->calculateVisibilityHistogram( &optVisHist, lineEdit0, 1 );
		
	}else if( renderingTech == 6 ){
		for( int i=0; i<10; i++ ){
			img->sampleRendering( lineEdit3, sldBar3, lineEdit0, sldBar8/99.0 );
			glWg->setImagePointer( img->getImgColorPointer() );
			glWg->updateGL();
		}
	}else if( renderingTech == 7 ){
		//img->isosurfaceRenderingAO( lineEdit1 );
		img->fitStorageMultiResolution();
		glWg->setImagePointer( img->getImgColorPointer() );
		glWg->setImagePointer2( img->getImgColorPointer2() );
		glWg->updateGL();
	}else if( renderingTech == 8 ){
		//img->colorRenderingAO( lineEdit3, lineEdit1);
		img->summerizeSampleRendering( lineEdit3, sldBar3, lineEdit0, sldBar8/99.0 );
		glWg->setImagePointer( img->getImgColorPointer() );
		glWg->updateGL();
	}else if( renderingTech == 9 ){
		//test cuda
		float visHist[128];
		float otf[128];
		intepolateTransferFunctionFromControlPoint( otf, 128, g_ctrlPtV, g_ctrlPtA );
   //     runCudaKernel(visHist, K, D, B, otf, rawHistogramRays );
	}

}


void MainWindow::on_pushButtonDisQuery_clicked()
{
	//int msDownX, msDownY, msRlsX, msRlsY;
	//glWg->getMsSltInfo( msDownX, msDownY, msRlsX, msRlsY );
	//int tmp = msDownY;
	//msDownY = msRlsY;
	//msRlsY =  tmp;

	//Histogram hist(BINS);
	//img->calculateDistributionQuery(  &hist, msDownX, msDownY, msRlsX, msRlsY,
	//	ui->horizontalSlider6->value()/99.0, ui->horizontalSlider7->value()/99.0 ); 

	////Histogram hist( 128 );
	////img->getRayTrend( &hist, 128, msDownX, msDownY );

	//for( int i=0; i<BINS; i++ ){
	//	g_distribution[i] = hist.getBin(i);
	//}
	
	Histogram hist( 128 );
	img->calculateVisibilityHistogram( &hist, ui->lineEdit0->text().toFloat(), ui->lineEdit3->text().toFloat() );
	for( int i=0; i<BINS; i++ ){
		g_distribution[i] = hist.getBin(i);
	}

	/*FILE* fp = fopen( "tmp.txt", "w" );
	for( int i=0; i< 128 ; i++ ){
		fprintf( fp, "%d %f\n", i, hist.getBin(i) );
	}
	fclose( fp );*/
}

void MainWindow::on_checkBoxDisQuery_clicked()
{
    repaint();
}

void MainWindow::on_testButton_clicked()
{

}
