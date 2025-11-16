#include "mainwindow.h"
#include <QApplication>
#include <QLocale>
#include <QTranslator>
#include "primenumber.h"
#include <iostream>
#include"olfunct.h"
int main(int argc, char *argv[])
{
    //QApplication a(argc, argv);
      PrimeNumer prime;
    // if(prime.miller_rabin()){
    //     for (uint32_t &x : prime.prime1)
    //         std::cout << x << ',';
    //     std::cout<<std::endl;
    //     for (uint32_t &x : prime.p1_n_1)
    //         std::cout << x << ',';
    //     std::cout<<std::endl;
    //     for (uint32_t &x : prime.prime2)
    //         std::cout << x << ',';
    //     for (uint32_t &x : prime.p2_n_1)
    //         std::cout << x << ',';
    // }
      olfunct test(prime);
      // for(auto t:test.mu) std::cout<<t<<',';
      // std::cout<<std::endl;
      for(auto t:test.olfunction) std::cout<<t<<',';
      std::cout<<std::endl;
      test.solve();
      for(auto t:test.d) std::cout<<t<<',';


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
