project(VolumeRenderingGUI)

cmake_minimum_required(VERSION 2.6)

Set(CMAKE_PREFIX_PATH "C:/Program Files (x86)/Windows Kits/8.0/Lib/win8/um/x64")

if(POLICY CMP0020)
  cmake_policy(SET CMP0020 NEW)
endif()

# Set your files and resources here
set( Srcs mainwindow.cpp main.cpp glwidget.cpp)

set( Hdrs rendererThread.h RegularScalarData.h Ray.h myUtility.h mainwindow.h IntegralHistogram.h Image.h Histogram.h glwidget.h)

set( UIs mainwindow.ui )

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)


set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

find_package(Qt5Widgets)
qt5_wrap_ui(UI_Srcs ${UIs})
qt5_add_resources(QtApp_RCC_SRCS ${Rcs})
add_executable(VolumeRenderingGUI ${Srcs} ${Hdrs} ${SHARED_SOURCES} ${QtApp_RCC_SRCS} ${UI_Srcs})#)# ${MOC_Hdrs})
qt5_use_modules(VolumeRenderingGUI OpenGL GUI)
target_link_libraries(VolumeRenderingGUI  Qt5::Widgets opengl32.lib ${SHARED_LIBRARIES}) 
