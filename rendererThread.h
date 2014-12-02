#ifndef RENDERERTHREAD_H
#define RENDERERTHREAD_H

#define RendererPath "C:\\Users\\VisKC\\Dropbox\\00_ImportanceNow\\Research\\WindowsCode\\volumeRender\\bin\\volumeRender.exe"

#include "stdio.h"
#include "stdlib.h"
#include <QThread>

class rendererThread : public QThread
{
    private:

    public:
        rendererThread()
        {


        }

        void run()
        {
            system(RendererPath);
        }

};

#endif // RENDERERTHREAD_H
