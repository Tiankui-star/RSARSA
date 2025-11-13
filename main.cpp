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
    std::vector<uint32_t>t1={2};
    std::vector<uint32_t>t2={1000};
    std::vector<uint32_t>t3={1000001,122,22,556,21,15,166,244,46};
    std::vector<uint32_t>mu=prime.compute_mu(t3);
    std::vector<uint32_t> res = prime.quickExp(t1, t2, t3,mu);
    for (uint32_t &x : res)
        std::cout << x << ' ';
    // if(prime.miller_rabin()){
    //     for (uint32_t &x : prime.prime1)
    //         std::cout << x << ' ';
    // }
    // else {
    //     std::cout<<"not found";
    // }


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
