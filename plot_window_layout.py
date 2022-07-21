# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'plot_window.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog2(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(845, 699)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("chart_scatter_plot_icon_135788.ico"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        Dialog.setWindowIcon(icon)
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.scrollArea = QtWidgets.QScrollArea(Dialog)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setObjectName("scrollArea")
        self.scrollAreaWidgetContents = QtWidgets.QWidget()
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, -4, 808, 713))
        self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")
        self.formLayout = QtWidgets.QFormLayout(self.scrollAreaWidgetContents)
        self.formLayout.setObjectName("formLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label_2 = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        font = QtGui.QFont()
        font.setPointSize(9)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout.addWidget(self.label_2)
        self.lineEditDataForPlot = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)
        self.lineEditDataForPlot.setMaximumSize(QtCore.QSize(200, 16777215))
        font = QtGui.QFont()
        font.setPointSize(9)
        self.lineEditDataForPlot.setFont(font)
        self.lineEditDataForPlot.setObjectName("lineEditDataForPlot")
        self.horizontalLayout.addWidget(self.lineEditDataForPlot)
        self.pushButtonSearchFroPlot = QtWidgets.QPushButton(self.scrollAreaWidgetContents)
        self.pushButtonSearchFroPlot.setMaximumSize(QtCore.QSize(100, 27))
        self.pushButtonSearchFroPlot.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.pushButtonSearchFroPlot.setObjectName("pushButtonSearchFroPlot")
        self.horizontalLayout.addWidget(self.pushButtonSearchFroPlot)
        self.formLayout.setLayout(0, QtWidgets.QFormLayout.LabelRole, self.horizontalLayout)
        self.tableWidgetSelectedDataForPlot = QtWidgets.QTableWidget(self.scrollAreaWidgetContents)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tableWidgetSelectedDataForPlot.sizePolicy().hasHeightForWidth())
        self.tableWidgetSelectedDataForPlot.setSizePolicy(sizePolicy)
        self.tableWidgetSelectedDataForPlot.setMinimumSize(QtCore.QSize(700, 250))
        self.tableWidgetSelectedDataForPlot.setMaximumSize(QtCore.QSize(800, 290))
        self.tableWidgetSelectedDataForPlot.setStyleSheet("background-color: rgb(212, 214, 255);")
        self.tableWidgetSelectedDataForPlot.setObjectName("tableWidgetSelectedDataForPlot")
        self.tableWidgetSelectedDataForPlot.setColumnCount(0)
        self.tableWidgetSelectedDataForPlot.setRowCount(0)
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.tableWidgetSelectedDataForPlot)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.pushButtonAddData = QtWidgets.QPushButton(self.scrollAreaWidgetContents)
        self.pushButtonAddData.setMaximumSize(QtCore.QSize(100, 27))
        self.pushButtonAddData.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.pushButtonAddData.setObjectName("pushButtonAddData")
        self.horizontalLayout_2.addWidget(self.pushButtonAddData)
        self.pushButtonClear = QtWidgets.QPushButton(self.scrollAreaWidgetContents)
        self.pushButtonClear.setMaximumSize(QtCore.QSize(100, 27))
        self.pushButtonClear.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.pushButtonClear.setObjectName("pushButtonClear")
        self.horizontalLayout_2.addWidget(self.pushButtonClear)
        self.formLayout.setLayout(2, QtWidgets.QFormLayout.LabelRole, self.horizontalLayout_2)
        self.line = QtWidgets.QFrame(self.scrollAreaWidgetContents)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.line)
        self.listWidget = QtWidgets.QListWidget(self.scrollAreaWidgetContents)
        self.listWidget.setMinimumSize(QtCore.QSize(400, 150))
        self.listWidget.setMaximumSize(QtCore.QSize(400, 500))
        self.listWidget.setObjectName("listWidget")
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.listWidget)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.pushButtonGeneratePlot = QtWidgets.QPushButton(self.scrollAreaWidgetContents)
        self.pushButtonGeneratePlot.setMaximumSize(QtCore.QSize(500, 27))
        self.pushButtonGeneratePlot.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.pushButtonGeneratePlot.setStyleSheet("font-weight: bold;\n"
"color: rgb(85, 85, 255);")
        self.pushButtonGeneratePlot.setObjectName("pushButtonGeneratePlot")
        self.horizontalLayout_3.addWidget(self.pushButtonGeneratePlot)
        self.pushButton = QtWidgets.QPushButton(self.scrollAreaWidgetContents)
        self.pushButton.setMaximumSize(QtCore.QSize(100, 27))
        self.pushButton.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.pushButton.setStyleSheet("font-weight: bold;\n"
"color: rgb(255, 0, 0);")
        self.pushButton.setObjectName("pushButton")
        self.horizontalLayout_3.addWidget(self.pushButton)
        self.formLayout.setLayout(5, QtWidgets.QFormLayout.LabelRole, self.horizontalLayout_3)
        self.label = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        font = QtGui.QFont()
        font.setPointSize(9)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.formLayout.setWidget(6, QtWidgets.QFormLayout.LabelRole, self.label)
        self.labelPlot = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        self.labelPlot.setText("")
        self.labelPlot.setObjectName("labelPlot")
        self.formLayout.setWidget(7, QtWidgets.QFormLayout.LabelRole, self.labelPlot)
        self.widget = QtWidgets.QWidget(self.scrollAreaWidgetContents)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.widget.sizePolicy().hasHeightForWidth())
        self.widget.setSizePolicy(sizePolicy)
        self.widget.setObjectName("widget")
        self.formLayout.setWidget(8, QtWidgets.QFormLayout.LabelRole, self.widget)
        self.label_3 = QtWidgets.QLabel(self.scrollAreaWidgetContents)
        font = QtGui.QFont()
        font.setPointSize(9)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.formLayout.setWidget(9, QtWidgets.QFormLayout.LabelRole, self.label_3)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.pushButtonPasteElements = QtWidgets.QPushButton(self.scrollAreaWidgetContents)
        self.pushButtonPasteElements.setMaximumSize(QtCore.QSize(150, 27))
        self.pushButtonPasteElements.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.pushButtonPasteElements.setObjectName("pushButtonPasteElements")
        self.horizontalLayout_4.addWidget(self.pushButtonPasteElements)
        self.pushButtonAddToWorkspace = QtWidgets.QPushButton(self.scrollAreaWidgetContents)
        self.pushButtonAddToWorkspace.setMinimumSize(QtCore.QSize(150, 27))
        self.pushButtonAddToWorkspace.setMaximumSize(QtCore.QSize(150, 27))
        self.pushButtonAddToWorkspace.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.pushButtonAddToWorkspace.setObjectName("pushButtonAddToWorkspace")
        self.horizontalLayout_4.addWidget(self.pushButtonAddToWorkspace)
        self.formLayout.setLayout(11, QtWidgets.QFormLayout.LabelRole, self.horizontalLayout_4)
        self.textEdit = QtWidgets.QTextEdit(self.scrollAreaWidgetContents)
        self.textEdit.setMaximumSize(QtCore.QSize(600, 120))
        self.textEdit.setObjectName("textEdit")
        self.formLayout.setWidget(10, QtWidgets.QFormLayout.SpanningRole, self.textEdit)
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.verticalLayout_4.addWidget(self.scrollArea)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Plot creator"))
        self.label_2.setText(_translate("Dialog", "Type in organism name or BioProject id:"))
        self.pushButtonSearchFroPlot.setText(_translate("Dialog", "Search"))
        self.pushButtonAddData.setText(_translate("Dialog", "Add data"))
        self.pushButtonClear.setText(_translate("Dialog", "Clear"))
        self.pushButtonGeneratePlot.setText(_translate("Dialog", "Generate plot"))
        self.pushButton.setText(_translate("Dialog", "Delete data"))
        self.label.setText(_translate("Dialog", "Plot:"))
        self.label_3.setText(_translate("Dialog", "Notes:"))
        self.pushButtonPasteElements.setText(_translate("Dialog", "Paste elements from list"))
        self.pushButtonAddToWorkspace.setText(_translate("Dialog", "Add to workspace"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = QtWidgets.QDialog()
    ui = Ui_Dialog2()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())