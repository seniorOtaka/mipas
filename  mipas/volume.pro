TEMPLATE      = vclib
CONFIG       += thread
DEFINES			 += MITK_DLL

INCLUDEPATH  += ../core/ ../main
INCLUDEPATH  += 

HEADERS = mitkdefines.h	\
volumewidget.h	\
volumewindow.h

SOURCES = 
volumewidget.cpp	\
volumewindow.cpp

win32:DEFINES -= NOMINMAX
LIBS				 += 

DESTDIR			  = ../../bin/
TARGET        = $$qtLibraryTarget(volume)
INSTALLS		 += ../../bin/
