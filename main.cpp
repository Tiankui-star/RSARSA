#include "mainwindow.h"
#include <QApplication>
#include <QLocale>
#include <QTranslator>
#include <QProcess>
#include <QDebug>
#include<iostream>
#include"primenumber.h"
#include <QCoreApplication>
int main(int argc, char *argv[])
{
    //prime=new PrimeNumer();
    //prime->miller_rabin();
    QApplication a(argc, argv);


     /*std::vector<uint32_t> textnum={263356,5959595,126515165,19651561,191119};
     std::cout<<"yuanwen"<<std::endl;





      olfunct test(prime);
      test.solve();

      std::cout<<"siyao"<<std::endl;
      for(auto &t:prime.prime1) std::cout<<t<<',';
      std::cout<<std::endl;

      textnum=getmod(textnum,test.pube,test.mu);
       // textnum=(prime.quickExp(textnum,test.pube,test.mu,test.mumu));
      std::cout<<"miwen"<<std::endl;
      for(auto &t:textnum) std::cout<<t<<' ';
      std::cout<<std::endl;
      std::cout<<"jieji"<<std::endl;

      std::vector<uint32_t>rsu=getmod(textnum,test.d,test.mu);;
      for(auto &t:rsu) std::cout<<t<<*/' ';





    QTranslator translator;
    const QStringList uiLanguages = QLocale::system().uiLanguages();
    for (const QString &locale : uiLanguages) {
        const QString baseName = "untitled_" + QLocale(locale).name();
        if (translator.load(":/i18n/" + baseName)) {
            a.installTranslator(&translator);
            break;
        }
    }
    MainWindow w;
    w.show();
   return a.exec();
    // return 0;

}



