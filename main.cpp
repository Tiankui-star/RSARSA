#include "mainwindow.h"
#include <QApplication>
#include <QLocale>
#include <QTranslator>
#include "primenumber.h"
#include <iostream>
#include <QProcess>
#include <QDebug>
#include"olfunct.h"
#include <QCoreApplication>
std::vector<uint32_t> getmod(std::vector<uint32_t>&base,std::vector<uint32_t>&exp,std::vector<uint32_t>&mod);
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
      PrimeNumer prime;
     std::vector<uint32_t> textnum={263356,5959595,126515165,19651561,191119};
     std::cout<<"yuanwen"<<std::endl;



    // prime.miller_rabin();

      olfunct test(prime);
      test.solve();
      for(auto &t:textnum) std::cout<<t<<',';
      std::cout<<std::endl;
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
      for(auto &t:rsu) std::cout<<t<<' ';
      std::cout<<std::endl;

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

}
std::vector<uint32_t> getmod(std::vector<uint32_t>&base,std::vector<uint32_t>&exp,std::vector<uint32_t>&mod){
    QProcess process;
    QString program = "C:/Users/Gdell/AppData/Local/Programs/Python/Python312/python.exe";
    QStringList arguments;
    QString textnumStr1;
    for(size_t i=0; i<base.size(); ++i) {
        textnumStr1 += QString::number(base[i]);
        if(i != base.size()-1) textnumStr1 += ",";
    }
    QString textnumStr2;
    for(size_t i=0; i<exp.size(); ++i) {
        textnumStr2 += QString::number(exp[i]);
        if(i != exp.size()-1) textnumStr2 += ",";
    }
    QString textnumStr3;
    for(size_t i=0; i<mod.size(); ++i) {
        textnumStr3 += QString::number(mod[i]);
        if(i != mod.size()-1) textnumStr3 += ",";
    }


    arguments << "D:\\QTProject\\RSA\\learning.py" << textnumStr1 << textnumStr2<<textnumStr3<<"jiemi";
    process.start(program, arguments);
    if (!process.waitForStarted()) {
        qDebug() << "Failed to start Python!";

    }
    process.waitForFinished();
    QByteArray output = process.readAllStandardOutput();
    std::vector<uint32_t>rsu;
    int i=1;
    for(;i<=output.size()-4;i++){

        uint32_t val=0;
        while(i<=output.size()-4&&output[i]!=','){
            val=val*10+output[i]-'0';
            i++;
        }

        rsu.push_back(val);
        i++;
        //val=0;
    }
    return rsu;

}


