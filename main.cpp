#include "mainwindow.h"
#include <QApplication>
#include <QLocale>
#include <QTranslator>
#include "primenumber.h"
#include <iostream>
int main(int argc, char *argv[])
{

    QApplication a(argc, argv);
    PrimeNumer prime;
    //prime.miller_rabin();
    std::vector<uint32_t>t1={123456789};
    std::vector<uint32_t>t2={123456789};
    std::vector<uint32_t>t3={1000000007};
    std::vector<uint32_t>res=prime.quickExp(t1,t2,t3);
    std::cout<<res.size()<<std::endl;
    for(uint32_t &x:res) std::cout<<x<<' ';
    // QTranslator translator;
    // const QStringList uiLanguages = QLocale::system().uiLanguages();
    // for (const QString &locale : uiLanguages) {
    //     const QString baseName = "untitled_" + QLocale(locale).name();
    //     if (translator.load(":/i18n/" + baseName)) {
    //         a.installTranslator(&translator);
    //         break;
    //     }
    // }
    // MainWindow w;
    // w.show();
    // return a.exec();
    return 0;
}
