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

        self.dataFromDatabase = pd.DataFrame()
        self.dataForPlot= pd.DataFrame()
        self.dataNames = ['Size(Mb)', 'GC%', 'Scaffolds', 'CDS', 'Chromosomes', 'Organelles', 'Assemblies', 'Plasmids']
        self.data = pd.DataFrame()
        self.tab = []
        self.tabConcat = []
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
            self.find = self.ui.lineEditDataForPlot.text()
            try:
                if "BioProject" in self.dataFromDatabase:
                    self.data = pd.DataFrame(self.dataFromDatabase.query(
                        f'Organism_name == ["{self.find}"] or BioProject == ["{self.find}"]'))
                else:
                    self.data = pd.DataFrame(self.dataFromDatabase.query(
                        f'Organism_name == ["{self.find}"]'))
                self.tab = [self.data.columns.values.tolist()] + self.data.values.tolist()
                self.item = self.data.iloc[0]['Organism_name']
                numrows, numcols = self.data.shape
                self.ui.tableWidgetSelectedDataForPlot.setColumnCount(numcols)
                self.ui.tableWidgetSelectedDataForPlot.setRowCount(numrows)
                self.ui.tableWidgetSelectedDataForPlot.setHorizontalHeaderLabels(self.data.columns)

                # sends search result to the table

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
                self.tabConcat.append(self.tab[1])
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
                self.tabConcat.clear()
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
            dataText = ', '.join(self.items)
            self.ui.textEdit.append(dataText)
        else:
            QMessageBox.information(self,'Info','No data to paste')

    def plot(self):
        """
        :return: creates data plot
        """
        try:
            self.header = self.tab[0]
            self.dataForPlot = pd.DataFrame(self.tabConcat,columns=self.header)
            organism_names_for_plot = self.dataForPlot['Organism_name'].tolist()
            headersNames = self.dataForPlot.columns.values.tolist()
            self.valuesToConvert = []

            # converts data to float

            for i in headersNames:
                for j in self.dataNames:
                    if i == j:
                        self.valuesToConvert.append(i)
            for i in self.valuesToConvert:
                self.dataForPlot[i] = self.dataForPlot[i].astype(float)
            dataPlot = self.dataForPlot[self.valuesToConvert]

            # dynamic dataframe creation for plot

            dataPlot.insert(0, 'Organism_name', organism_names_for_plot, True)
            dataPlot['Size(Kb)'] = dataPlot['Size(Mb)'] * 1000
            dataPlot = dataPlot.drop('Size(Mb)', axis=1)
            self.dane_melted = dataPlot.melt("Organism_name", var_name="Organism", value_name="Value")
            plt.figure()

            # shows plot in QLabel

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
    