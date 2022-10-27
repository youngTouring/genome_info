from PyQt5.QtWidgets import QDialog, QMessageBox,QMainWindow, QInputDialog, QListWidgetItem
from PyQt5.QtGui import QPainter, QPixmap
from plot_window_layout import *
import matplotlib.pyplot as plt
from main_window import *
import seaborn as sns
import pandas as pd
import numpy as np
import sys

#################################################################
### MyPlotDialog class - instantiation in application's main class.
### Responsible for plotting data from database, previosuly loaded
### to main table.
#################################################################

class MyPlotDialog(QDialog,QMainWindow):
    def __init__(self):
        super().__init__()
        self.ui = Ui_Dialog2()
        self.ui.setupUi(self)

        self.ui.pushButtonSearchFroPlot.clicked.connect(self.search_for_plot)
        self.ui.pushButtonClear.clicked.connect(self.clear_data)
        self.ui.pushButtonGeneratePlot.clicked.connect(self.plot)
        self.ui.pushButtonAddData.clicked.connect(self.add_data)
        self.ui.pushButton.clicked.connect(self.delete_data)
        self.ui.pushButtonPasteElements.clicked.connect(self.paste_list)

        self.data_from_database = pd.DataFrame()
        self.data_for_plot= pd.DataFrame()
        self.dataNames = ['Size(Mb)', 'GC%', 'Scaffolds', 'CDS', 'Chromosomes', 'Organelles', 'Assemblies', 'Plasmids']
        self.data = pd.DataFrame()
        self.tab = []
        self.tab_concat = []
        self.header = []
        self.item = ''
        self.items = []

    def clear_data(self):
        """
        :return: clears data from widget with selected featrures from table
        """
        self.ui.tableWidgetSelectedDataForPlot.clear()
        
    def search_for_plot(self):
        """
        :return: filters data for plott by organism nama or BioProject id
        """
        try:
            find = self.ui.lineEditDataForPlot.text()
            try:
                if "BioProject" in self.data_from_database:
                    self.data = pd.DataFrame(self.data_from_database.query(
                        f'Organism_name == ["{find}"] or BioProject == ["{find}"]'))
                else:
                    self.data = pd.DataFrame(self.data_from_database.query(
                        f'Organism_name == ["{find}"]'))
                self.tab = [self.data.columns.values.tolist()] + self.data.values.tolist()
                self.item = self.data.iloc[0]['Organism_name']
                numrows, numcols = self.data.shape
                self.ui.tableWidgetSelectedDataForPlot.setColumnCount(numcols)
                self.ui.tableWidgetSelectedDataForPlot.setRowCount(numrows)
                self.ui.tableWidgetSelectedDataForPlot.setHorizontalHeaderLabels(self.data.columns)

                for row in range(numrows):
                    for column in range(numcols):
                        self.ui.tableWidgetSelectedDataForPlot.setItem(row, column, QTableWidgetItem(
                            str(self.data.iloc[row, column])))
            except Exception:
                QMessageBox.information(self,'Info','Try to search for another data')
        
        except Exception as e:
            QMessageBox.critical(self, 'Error', e)

    def add_data(self):
        """
        :return: adds data to the list
        """
        try:
            if len(self.item) > 0:
                self.ui.listWidget.clear()
                self.tab_concat.append(self.tab[1])
                self.items.append(self.item)
                for i in range(len(self.items)):
                    self.ui.listWidget.addItem(self.items[i])
            else:
                # if data was deleted or user did not searched for more return error
                QMessageBox.information(self,'Info','Please search for data again')
        except Exception as e:
            QMessageBox.critical(self,'Error',f'Something went wrong: {e}')

    def delete_data(self):
        """
        :return: deletes data from the list
        """
        try:
            if len(self.items) !=  0:
                self.tab_concat.clear()
                self.item = ''
                self.items.clear()
                self.dane_melted = None
                self.ui.listWidget.clear()
            else:
                # if there is no data to delete return info
                QMessageBox.information(self,'Info','No data to remove')
        except Exception as e:
            QMessageBox.critical(self,'Error',f'Something went wrong: {e}')

    def paste_list(self):
        """
        :return: pastes organism names from the list to describe it
        """
        if len(self.items) > 0:
            data_text = ', '.join(self.items)
            self.ui.textEdit.append(data_text)
        else:
            QMessageBox.information(self,'Info','No data to paste')

    def plot(self):
        """
        :return: creates data plot
        """
        try:
            self.header = self.tab[0]
            self.data_for_plot = pd.DataFrame(self.tab_concat,columns=self.header)
            organism_names_for_plot = self.data_for_plot['Organism_name'].tolist()
            headers_names = self.data_for_plot.columns.values.tolist()
            self.valuesToConvert = []
            for i in headers_names:
                for j in self.dataNames:
                    if i == j:
                        self.valuesToConvert.append(i)
            for i in self.valuesToConvert:
                self.data_for_plot[i] = self.data_for_plot[i].astype(float)
            data_plot = self.data_for_plot[self.valuesToConvert]
            data_plot.insert(0, 'Organism_name', organism_names_for_plot, True)
            data_plot['Size(Kb)'] = data_plot['Size(Mb)'] * 1000
            data_plot = data_plot.drop('Size(Mb)', axis=1)
            self.dane_melted = data_plot.melt("Organism_name", var_name="Organism", value_name="Value")
            plt.figure()
            self.plot = sns.barplot(data=self.dane_melted, x='Organism', y='Value', hue="Organism_name")
            self.PlotFile = 'plot.png'
            self.plot.get_figure().savefig(self.PlotFile)
            self.plot_png = QPixmap(self.PlotFile)
            self.ui.labelPlot.setPixmap(self.plot_png)
        except Exception as e:
            QMessageBox.critical(self,'Error',f'Something went wrong: {e}')

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    w = MyPlotDialog()
    w.show()
    sys.exit(app.exec_())
