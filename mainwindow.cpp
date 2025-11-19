#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>
#include <QMessageBox>
#include "primenumber.h"
#include "olfunct.h"
#include<vector>
#include<QChar>
#include<QProcess>
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
    ,prime()
{
    ui->setupUi(this);

    // 设置输出框为只读（根据实际组件名修改）

     ui->plainTextEdit->setReadOnly(true);  // 如果是QPlainTextEdit且未重命名
    //prime.miller_rabin();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_encrypt_clicked()
{
    // 根据实际组件名获取输入
    QString inputText = ui->textEdit_input->toPlainText();  // 如果是QTextEdit


    if(inputText.isEmpty()) {
        QMessageBox::warning(this, "警告", "请输入要加密的文本!");
        return;
    }
    //qDebug()<<inputText.size();
    // 调用加密函数
    std::vector<uint32_t> encryptedText = encryptText(inputText);

    // 根据实际组件名设置输出
    QStringList decList;
    for (uint32_t value : encryptedText) {
        decList.append(QString::number(value));  // 直接转换为十进制字符串
        decList.append(QString(','));
    }

    ui->textEdit_output->setPlainText(decList.join(" "));
}

void MainWindow::on_pushButton_decrypt_clicked()
{
    QString inputText = ui->inputText->toPlainText();

    if(inputText.isEmpty()) {
        QMessageBox::warning(this, "警告", "请输入要解密的文本!");
        return;
    }

    std::vector<uint32_t> decryptedText = decryptText(inputText);

    QStringList decList;
    for (uint32_t value : decryptedText) {
        decList.append(QString::number(value));  // 直接转换为十进制字符串
        decList.append(QString(','));
    }

    ui->plainTextEdit->setPlainText(decList.join(" "));
}
std::vector<uint32_t> MainWindow::getmod(std::vector<uint32_t>&base,std::vector<uint32_t>&exp,std::vector<uint32_t>&mod){
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
std::vector<uint32_t> MainWindow::encryptText(const QString &plainText)
{
    // 这里实现加密逻辑
    std::vector<uint32_t>plain;
    int i=0;
    for(;i<plainText.size();i++){
        uint32_t val=0;
        while(i<plainText.size()&&plainText[i]!=','){
            val=val*10+(plainText[i].digitValue());
            i++;
        }
        plain.push_back(val);
    }
    olfunct test(prime);

    std::vector<uint32_t> secret=getmod(plain,test.pube,test.mu);
    return secret;
}




std::vector<uint32_t> MainWindow::decryptText(const QString &cipherText)
{
    // 这里实现解密逻辑
    qDebug()<<cipherText[1];
    olfunct test(prime);
    std::vector<uint32_t>cplain;
    int i=0;
    for(;i<cipherText.size();i++){
        uint32_t val=0;
        while(i<cipherText.size()&&cipherText[i]!=','){
            val=val*10+(cipherText[i].digitValue());
            i++;
        }
        cplain.push_back(val);
        qDebug()<<val<<' ';
    }

    std::vector<uint32_t>rsu=getmod(cplain,test.d,test.mu);
    return rsu;
}
