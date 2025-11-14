#include "mainwindow.h"
#include <QApplication>
#include <QLocale>
#include <QTranslator>
#include "primenumber.h"
#include <iostream>
int main(int argc, char *argv[])
{
    //QApplication a(argc, argv);
     PrimeNumer prime;
    // std::vector<uint32_t>t1={2};
    // std::vector<uint32_t>t2={1000,26161,1616,15151,11,6161,616,16,161,61,61,61,1};
    // std::vector<uint32_t>t3={1000001,122,22,556,21,15,166,244,4615151,5111,65515111,1552,5116,5226,25662,55225,2533,226,294,941,191,91,19,441,551,5512225,111555,1144,114,112};
    // std::vector<uint32_t>mu=prime.compute_mu(t3);
    // //std::cout<<1;
    // std::vector<uint32_t> res = prime.quickExp(t1, t2, t3,mu);
    // for (uint32_t &x : res)
    //     std::cout << x << ' ';
    if(prime.miller_rabin()){
        for (uint32_t &x : prime.prime1)
            std::cout << x << ',';
    }
    else {
        std::cout<<"not found";
    }


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
