import sys

from PyQt5.QtWidgets import QDialog, QApplication, QTableWidgetItem, QTableWidget
from PyQt5.QtCore import QThread
import requests
import csv

from layout.downloading_dialog_layout import *
from main_window import *
from addresses import *

##################################################################
### MyDialog class - instantiation in application's main class.
### Responsible for downloading chosen Genome database from NCBI:
### oganelles, procaryota, eucaryota, viruses, plasmids, overview
##################################################################

class MyDialog(QDialog):
    def __init__(self):
        super().__init__()
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        self.ui.pushButtonStartDownload.clicked.connect(self.download_database_window)
        self.thread = {}
        self.adresses_radio_check = {'Overview':overview,'Organelles':organelles,
                                     'Eucaryota':eucaryotes,'Procaryota':procaryotes,
                                     'Viruses':viruses,'Plasmids':plasmids}

    @staticmethod
    def download_csv(url,file) -> None:
        """
        :param url: address of database
        :param file: specififed file name by user
        :return: saves downloaded database as csv file
        """
        req = requests.get(url)
        url_content = req.content
        csv_file = open(file, 'wb')
        csv_file.write(url_content)
        csv_file.close()
        # self.thread[1].stop()

    def download_database_window(self) -> None:
        """
        :return: dowloads chosen database
        """
        # self.thread[1] = ThreadClass(parent = None, index = 1)
        # self.thread[1].start()
        self.ui.labelDownloadedInfo.clear()
        try:
            chosen_option = self.ui.comboBoxDatabases.itemText(self.ui.comboBoxDatabases.currentIndex())
            file, _ = QFileDialog.getSaveFileName(self, "Save file",
                                                   f"{chosen_option}",
                                                   "All files (*);;CSV files (*.csv)")
            if file:
                url = self.adresses_radio_check[chosen_option]
                self.download_csv(url,file)
                self.ui.labelDownloadedInfo.setText(f'{self.ui.comboBoxDatabases.itemText(self.ui.comboBoxDatabases.currentIndex())} database'
                                                    f' successfully downloaded')
            else:
                self.thread[1].stop()
        except Exception as e:
            '''In case of lack of internet connection return error'''
            QMessageBox.critical(self,'Error',f'Something went wrong. Check your internet connection: {e}')


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    w = MyDialog()
    w.show()
    sys.exit(app.exec_())
