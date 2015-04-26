QMAKE_CXXFLAGS += -static 
QMAKE_MAKEFILE = C5Makefile
DESTDIR = .
CONFIG += qt thread 
TEMPLATE = app
TARGET = C5
DEPENDPATH += .
INCLUDEPATH += . \
	QCodeEdit \
	QCodeEdit/document \
	/home/gio/CoCoALib-0.99/include  
LIBS += -L/home/gio/CoCoALib-0.99/lib -L/home/gio/cocoa-hg/QCodeEdit -lboost_thread-mt -lQtXml -lpthread -lzip -lcocoa -lgmp.dll -Wl,-whole-archive -lqcodeedit -Wl,-no-whole-archive
DEFINES += C5IDE
MOC_DIR = 
UI_DIR =
OBJECTS_DIR = Debug
QT += xml

# Input
HEADERS += AST.H \
           C5.H \
           CompilationDate.H \
           Interpreter.H \
           Lexer.H \
           LineProviders.H \
           Parser.H \
	   OnlineHelp.H
FORMS += Console.ui MainWindow.ui SourceEditor.ui Debugger.ui
SOURCES += AST.C \
           BuiltInFunctions.C \
           C5.C \
           CompilationDate.C \
           Interpreter.C \
           Lexer.C \
           LineProviders.C \
           Main.C \
           Parser.C \
	   OnlineHelp.C
RESOURCES += C5.qrc

