# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'downloading_dialog.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(483, 282)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("database_theapplication_3365.ico"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        Dialog.setWindowIcon(icon)
        self.layoutWidget = QtWidgets.QWidget(Dialog)
        self.layoutWidget.setGeometry(QtCore.QRect(20, 10, 441, 251))
        self.layoutWidget.setObjectName("layoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.layoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        self.pushButtonStartDownload = QtWidgets.QPushButton(self.layoutWidget)
        font = QtGui.QFont()
        font.setPointSize(9)
        self.pushButtonStartDownload.setFont(font)
        self.pushButtonStartDownload.setObjectName("pushButtonStartDownload")
        self.gridLayout.addWidget(self.pushButtonStartDownload, 4, 1, 1, 1)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QtWidgets.QLabel(self.layoutWidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.gridLayout.addLayout(self.horizontalLayout, 1, 0, 1, 2)
        self.labelDownloadedInfo = QtWidgets.QLabel(self.layoutWidget)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.labelDownloadedInfo.setFont(font)
        self.labelDownloadedInfo.setText("")
        self.labelDownloadedInfo.setObjectName("labelDownloadedInfo")
        self.gridLayout.addWidget(self.labelDownloadedInfo, 5, 0, 1, 2)
        self.labelDownloadInfo = QtWidgets.QLabel(self.layoutWidget)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.labelDownloadInfo.setFont(font)
        self.labelDownloadInfo.setText("")
        self.labelDownloadInfo.setObjectName("labelDownloadInfo")
        self.gridLayout.addWidget(self.labelDownloadInfo, 8, 0, 1, 1)
        self.comboBoxDatabases = QtWidgets.QComboBox(self.layoutWidget)
        font = QtGui.QFont()
        font.setPointSize(9)
        self.comboBoxDatabases.setFont(font)
        self.comboBoxDatabases.setObjectName("comboBoxDatabases")
        self.comboBoxDatabases.addItem("")
        self.comboBoxDatabases.addItem("")
        self.comboBoxDatabases.addItem("")
        self.comboBoxDatabases.addItem("")
        self.comboBoxDatabases.addItem("")
        self.comboBoxDatabases.addItem("")
        self.gridLayout.addWidget(self.comboBoxDatabases, 4, 0, 1, 1)
        spacerItem1 = QtWidgets.QSpacerItem(20, 20, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Maximum)
        self.gridLayout.addItem(spacerItem1, 2, 0, 1, 1)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Download database"))
        self.pushButtonStartDownload.setText(_translate("Dialog", "Download"))
        self.label.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:10pt; font-weight:600;\">Which database would you like to download?</span></p></body></html>"))
        self.comboBoxDatabases.setItemText(0, _translate("Dialog", "Overview"))
        self.comboBoxDatabases.setItemText(1, _translate("Dialog", "Organelles"))
        self.comboBoxDatabases.setItemText(2, _translate("Dialog", "Eucaryota"))
        self.comboBoxDatabases.setItemText(3, _translate("Dialog", "Procaryota"))
        self.comboBoxDatabases.setItemText(4, _translate("Dialog", "Viruses"))
        self.comboBoxDatabases.setItemText(5, _translate("Dialog", "Plasmids"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = QtWidgets.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())
