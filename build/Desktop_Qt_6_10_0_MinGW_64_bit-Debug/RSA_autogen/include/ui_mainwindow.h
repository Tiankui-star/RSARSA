/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 6.10.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPlainTextEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralwidget;
    QHBoxLayout *horizontalLayout;
    QGridLayout *gridLayout;
    QPlainTextEdit *plainTextEdit;
    QPlainTextEdit *textEdit_output;
    QLabel *label_2;
    QLabel *label_4;
    QLabel *label_3;
    QPlainTextEdit *inputText;
    QPushButton *pushButton_encrypt;
    QPushButton *pushButton_decrypt;
    QLabel *label;
    QPlainTextEdit *textEdit_input;
    QPushButton *pushButton_getpk;
    QMenuBar *menubar;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName("MainWindow");
        MainWindow->resize(940, 457);
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName("centralwidget");
        horizontalLayout = new QHBoxLayout(centralwidget);
        horizontalLayout->setObjectName("horizontalLayout");
        gridLayout = new QGridLayout();
        gridLayout->setObjectName("gridLayout");
        gridLayout->setSizeConstraint(QLayout::SizeConstraint::SetNoConstraint);
        gridLayout->setHorizontalSpacing(0);
        gridLayout->setVerticalSpacing(7);
        plainTextEdit = new QPlainTextEdit(centralwidget);
        plainTextEdit->setObjectName("plainTextEdit");

        gridLayout->addWidget(plainTextEdit, 2, 3, 1, 2);

        textEdit_output = new QPlainTextEdit(centralwidget);
        textEdit_output->setObjectName("textEdit_output");

        gridLayout->addWidget(textEdit_output, 2, 1, 1, 1);

        label_2 = new QLabel(centralwidget);
        label_2->setObjectName("label_2");

        gridLayout->addWidget(label_2, 1, 0, 1, 1);

        label_4 = new QLabel(centralwidget);
        label_4->setObjectName("label_4");

        gridLayout->addWidget(label_4, 1, 2, 1, 1);

        label_3 = new QLabel(centralwidget);
        label_3->setObjectName("label_3");

        gridLayout->addWidget(label_3, 1, 1, 1, 1);

        inputText = new QPlainTextEdit(centralwidget);
        inputText->setObjectName("inputText");

        gridLayout->addWidget(inputText, 2, 2, 1, 1);

        pushButton_encrypt = new QPushButton(centralwidget);
        pushButton_encrypt->setObjectName("pushButton_encrypt");

        gridLayout->addWidget(pushButton_encrypt, 0, 0, 1, 2);

        pushButton_decrypt = new QPushButton(centralwidget);
        pushButton_decrypt->setObjectName("pushButton_decrypt");

        gridLayout->addWidget(pushButton_decrypt, 0, 2, 1, 3);

        label = new QLabel(centralwidget);
        label->setObjectName("label");

        gridLayout->addWidget(label, 1, 3, 1, 2);

        textEdit_input = new QPlainTextEdit(centralwidget);
        textEdit_input->setObjectName("textEdit_input");

        gridLayout->addWidget(textEdit_input, 2, 0, 1, 1);

        pushButton_getpk = new QPushButton(centralwidget);
        pushButton_getpk->setObjectName("pushButton_getpk");

        gridLayout->addWidget(pushButton_getpk, 3, 2, 1, 1);


        horizontalLayout->addLayout(gridLayout);

        MainWindow->setCentralWidget(centralwidget);
        menubar = new QMenuBar(MainWindow);
        menubar->setObjectName("menubar");
        menubar->setGeometry(QRect(0, 0, 940, 25));
        MainWindow->setMenuBar(menubar);
        statusbar = new QStatusBar(MainWindow);
        statusbar->setObjectName("statusbar");
        MainWindow->setStatusBar(statusbar);
        QWidget::setTabOrder(textEdit_input, plainTextEdit);
        QWidget::setTabOrder(plainTextEdit, textEdit_output);
        QWidget::setTabOrder(textEdit_output, pushButton_decrypt);
        QWidget::setTabOrder(pushButton_decrypt, pushButton_encrypt);
        QWidget::setTabOrder(pushButton_encrypt, inputText);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QCoreApplication::translate("MainWindow", "MainWindow", nullptr));
        label_2->setText(QCoreApplication::translate("MainWindow", "<html><head/><body><p align=\"center\">    \346\230\216\346\226\207</p></body></html>", nullptr));
        label_4->setText(QCoreApplication::translate("MainWindow", "<html><head/><body><p align=\"center\">\345\257\206\346\226\207\350\276\223\345\205\245</p></body></html>", nullptr));
        label_3->setText(QCoreApplication::translate("MainWindow", "<html><head/><body><p align=\"center\">\345\257\206\346\226\207\350\276\223\345\207\272</p></body></html>", nullptr));
        pushButton_encrypt->setText(QCoreApplication::translate("MainWindow", "\345\212\240\345\257\206", nullptr));
        pushButton_decrypt->setText(QCoreApplication::translate("MainWindow", "\350\247\243\345\257\206", nullptr));
        label->setText(QCoreApplication::translate("MainWindow", "<html><head/><body><p align=\"center\">\350\247\243\345\257\206</p></body></html>", nullptr));
        pushButton_getpk->setText(QCoreApplication::translate("MainWindow", "\346\237\245\347\234\213\347\247\201\351\222\245", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
