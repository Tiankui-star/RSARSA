#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include<vector>
#include<QObject>
#include<QThread>
#include"primenumber.h"
#include"olfunct.h"
QT_BEGIN_NAMESPACE

namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_pushButton_encrypt_clicked();
    void on_pushButton_decrypt_clicked();
    void on_pushButton_getpk_clicked();
    void on_pushButton_getmod_clicked();

private:
    Ui::MainWindow *ui;

    PrimeNumer *prime;
    olfunct test;
    std::vector<uint32_t> getmod(std::vector<uint32_t>&base,std::vector<uint32_t>&exp,std::vector<uint32_t>&mod);
    std::vector<uint32_t>  encryptText(const QString &plainText);
    std::vector<uint32_t> decryptText(const QString &cipherText);
    std::vector<uint32_t>getpk();
    std::vector<uint32_t>getmod();
};
#endif
