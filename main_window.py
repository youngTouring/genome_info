from PyQt5.QtWidgets import QApplication, QTableWidgetItem, QMainWindow, QScrollArea, QFileDialog,QMessageBox, QLabel
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from Bio import Entrez, SeqIO, Seq
from main_window_layout import *
from downloading_dialog import *
import matplotlib.pyplot as plt
import reportlab.platypus
from plot_window import *
from table_model import *
from PIL import Image
import time, datetime
import pandas as pd
import subprocess
import Bio
import uuid
import sys

##############################################################################
### MyForm class - application's main class. Responsible for exposing ########
### database view, sending queries to Entrez database, saving data to ########
### file, loading data from file, and exposing plot view from MyPlotDialog ###
### class ####################################################################
##############################################################################

class MyForm(QMainWindow, QScrollArea):

    def __init__(self):
        super().__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.ui.actionDownloadDatabase.triggered.connect(self.open_downloading_dialog)
        self.ui.actionOpenDatabase.triggered.connect(self.open_file)
        self.ui.actionGenerete_plot.triggered.connect(self.open_plot_window)
        self.ui.actionGenerate_report.triggered.connect(self.generate_report)
        self.ui.actionHelp.triggered.connect(self.help)

        self.ui.pushButtonSearchInEntrez.clicked.connect(self.search_in_entrez)
        self.ui.pushButtonExportToTxt.clicked.connect(self.save_query_result_to_txt)
        self.ui.pushButtonSearchInDatabse.clicked.connect(self.search_in_database)
        self.ui.pushButtonLoadFromTxt.clicked.connect(self.load_from_txt)

        self.ui.tableWidgetSelectedData.cellClicked.connect(self.get_clicked_cell)
        self.queryResult = ''
        self.queryTranslated = ''
        self.data_loaded = False
        self.dialog = MyDialog()
        self.plotCreator = MyPlotDialog()
        self.plotCreator.ui.pushButtonAddToWorkspace.clicked.connect(self.add_chart_to_workspace)
        self.file_name_nocsv = ''
        self.thread = {}
        self.show()

    def open_downloading_dialog(self):
        '''Opens downloading dialog'''
        self.dialog.show()
        
    def open_plot_window(self):
        '''Opens plotting dialog. In case of lack of data return error'''
        try:
            self.plotCreator.dataFromDatabase = self.df_new
            self.plotCreator.show()
        except Exception as e:
            QMessageBox.critical(self,'Error','No data to plot. Load database first')

    def add_chart_to_workspace(self):
        '''Adds plot from plotting dialog to main window'''
        try:
            self.ui.labelPlottoWorkspace.setPixmap(self.plotCreator.plot_png)
            text = self.plotCreator.ui.textEdit.toPlainText()
            self.ui.textEditNoteForPlot.setText(text)
        except Exception as e:
            QMessageBox.critical(self,'Error',f'Something went wrong: {e}')

    def open_file(self):
        '''Opens .csv database file'''
        if self.plotCreator.isVisible():
            self.plotCreator.close()
        try:
            file, _ = QFileDialog.getOpenFileName(self,"Open file","","All files (*);;CSV files (*.csv)")
            if file:
                dane = pd.read_csv(file, keep_default_na=False)
                self.df_new = pd.DataFrame(dane)
                self.df_new.rename(columns = {'#Organism Name': 'Organism_name'}, inplace=True)
                self.df_new.sort_values(by = ['Organism_name'], ascending = True, inplace=True)
                self.df_new.insert(0,'Index',[i for i in range(1,self.df_new.shape[0] + 1)],True)
                self.data_loaded = True
                numrows, numcols = self.df_new.shape
                dane_tab = TabelaModel(self.df_new)
                self.ui.tableView.setModel(dane_tab)
                file_name = str(file).split('/')
                self.file_name_nocsv = file_name[-1].split('.')
                self.ui.labelDataName.setText('Data view ' + self.file_name_nocsv[0])
                self.ui.labelDatabaseSizeInfo.setText(f'Database size: {numcols - 1} columns, {numrows} rows')
                QMessageBox.information(self,'Info','Database successfully loaded!')
        except Exception as e:
            '''In case of incorrect file extension return error'''
            QMessageBox.critical(self,'Error','Chose correct file extension: csv')

    def search_in_database(self):
        '''Filters data for table by organism name or BioProject id'''
        try:
            if self.data_loaded:
                searchText = self.ui.lineEditSearchInDatabase.text()
                data_from_database = pd.DataFrame()
                try:
                    if "BioProject" in self.df_new:
                        data_from_database = pd.DataFrame(self.df_new.query(
                                f'Organism_name == ["{searchText}"] or BioProject == ["{searchText}"]'))
                    else:
                        data_from_database = pd.DataFrame(self.df_new.query(
                            f'Organism_name == ["{searchText}"]'))
                    numrows, numcols = data_from_database.shape
                    self.ui.tableWidgetSelectedData.setColumnCount(numcols)
                    self.ui.tableWidgetSelectedData.setRowCount(numrows)
                    self.ui.tableWidgetSelectedData.setHorizontalHeaderLabels(data_from_database.columns)
                    '''Send search result to the table'''
                    for row in range(numrows):
                        for column in range(numcols):
                            self.ui.tableWidgetSelectedData.setItem(row, column, QTableWidgetItem(str(data_from_database.iloc[row, column])))
                    self.ui.labelNumbeOfRows.setText('Number of results: ' + str(self.ui.tableWidgetSelectedData.rowCount()))
                except Exception:
                    QMessageBox.information(self,'Info','Try to search for another data')
        except Exception as e:
            QMessageBox.warning(self, 'Warning', e)

    def get_clicked_cell(self, row, column):
        '''When cell clicked, copy content to LineEdit below'''
        item_from_table = self.ui.tableWidgetSelectedData.item(row,column).text()
        self.ui.lineEditSelectedForQuery.setText(item_from_table)

    def search_in_entrez(self):
        '''
           Sends query to chosen Entrez databse. Possible databases:
           nucleotide, assemble, BioProject
        '''
        chosenEntrezDatabase = self.ui.comboBoxWhichDatabase.itemText(self.ui.comboBoxWhichDatabase.currentIndex()).lower()
        try:
            Entrez.email = f"{self.ui.lineEditEntrezEmail.text()}"
            query = f"{self.ui.lineEditSelectedForQuery.text()}"
            handle = Entrez.esearch(db=f'{chosenEntrezDatabase}',term=query,idtype='acc',retmax=5)
            result = Entrez.read(handle)
            handle = Entrez.efetch(db="nucleotide",id=result["IdList"],rettype="fasta")
            handle = Entrez.efetch(db="nucleotide",id=result["IdList"][0],rettype="gb")
            self.queryTranslated = SeqIO.read(handle,'genbank')
            handle = Entrez.efetch(db="nucleotide",id=result["IdList"],rettype="gb")
            self.queryResult = handle.read()
            '''Sends query result to TextEdit below'''
            self.ui.textEditQueryResult.setText(self.queryResult)
            handle.close()
        except Exception as c:
            '''In case of incorrect query return error'''
            QMessageBox.critical(self,'Error',f'Wrong query: {c}.\nTry to use another database for this query')

    def save_query_result_to_txt(self):
        '''Saves data to .txt file'''
        try:
            if len(self.queryResult) > 0:
                file,_ = QFileDialog.getSaveFileName(self, "Save file","","All files (*);;TXT files (*.txt)")
                if file:
                    saving = open(file,'w')
                    saving.write(self.queryResult)
        except Exception as e:
            '''In case of lack of data retur  error'''
            QMessageBox.information(self, 'Info',f'No data to export. Please send query\n''to receive data')

    def load_from_txt(self):
        '''Loads data from .txt file'''
        try:
            file, _ = QFileDialog.getOpenFileName(self,"Open file","","All files (*);;TXT files (*.txt)")
            if file:
                try:
                    self.ui.textEditQueryResult.clearFocus()
                    with open(file,'r') as f:
                        data = f.read()
                    self.ui.textEditQueryResult.setText(data)
                except Exception as e:
                    '''In case of incorrect file extension return error'''
                    QMessageBox.critical(self,'Error',f'Wrong file extension: {e}')
        except Exception as e:
            QMessageBox.show(self,'Info','Something went wrong. Try to load file again')

    ########################
    ### Report generation
    ########################
    def add_amino_percentage(self):
        try:
            data = []
            nucleotides = ['A','T','C','G']
            values = {}
            for nucleotide in nucleotides:
                count = self.queryTranslated.seq.count(nucleotide)
                values[nucleotide] = count
            all = sum(list(values.values()))
            for i in values:
                value = values[i]
                percent = value/all*100
                data.append(percent)
            explode = (0, 0.1, 0, 0)
            fig1, ax1 = plt.subplots()
            plot = ax1.pie(data, explode=explode, labels=nucleotides, autopct='%1.1f%%',
                    shadow=True, startangle=90)
            ax1.axis('equal')
            plt.title(f'Nucleotide percentage of {self.queryTranslated.annotations["organism"]}')
            plt.savefig('plot2.png')
            PlotFile = 'plot2.png'
            return PlotFile
        except Exception as e:
            QMessageBox.critical(self,'Error',f'No data! {e}')

    def generate_report(self):
        '''
        Generates report as a .pdf extension file. Data: report id, date&time of report generation
        entrez query results, plots: nucleotide percentage of organism from query
        and plot which comes from MyPlotDialog class
        '''
        try:
            org_name = '_'.join(self.queryTranslated.annotations["organism"].split())
            file, _ = QFileDialog.getSaveFileName(self, "Save file",f'{self.file_name_nocsv[0]}_{org_name}',
                                                  "All files (*);;PDF files (*.pdf)")
            if file:
                doc = SimpleDocTemplate(file, pagesize=letter,
                                        rightMargin=25, leftMargin=25,
                                        topMargin=25, bottomMargin=18)
                logo = 'C:/Users/miko5/Desktop/TDS/Project/icons/dna_icon_132453.ico'

                Story = []
                img = reportlab.platypus.Image(filename = logo, height = 0.5 * inch, width = 0.5 * inch, hAlign='RIGHT')

                name = 'GenomeInfo'
                report_id = f'Report ID: {uuid.uuid4()}'
                now = datetime.datetime.now()
                generation_date_and_time = now.strftime("%d/%m/%Y %H:%M:%S")
                date = f'Date and time of report generation: {generation_date_and_time}'
                result_txt = 'Result of Entrez query:'
                result_info_1 = f'ID number: {self.queryTranslated.id}'
                result_info_2 = f'Organism: {self.queryTranslated.annotations["organism"]}'
                result_info_3 = f'Sequence: {self.queryTranslated.seq[:100]}...'
                result_info_4 = f'Sequence length: {len(self.queryTranslated.seq)}'
                result_info_5 = f'Description: {self.queryTranslated.description}'
                information = [name,report_id,date,result_txt]
                results = [result_info_1,result_info_2,result_info_3,result_info_4,result_info_5]
                notes = self.ui.textEditNoteForPlot.toPlainText()

                Story.append(img)
                styles = getSampleStyleSheet()
                styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))

                for info in information:
                    Story.append(Paragraph(info, styles["Normal"]))
                    Story.append(Spacer(1, 12))

                for result in results:
                    Story.append(Spacer(1, 12))
                    Story.append(Paragraph(result, styles["Justify"]))

                Story.append(Spacer(1, 10))
                percentage_plot = self.add_amino_percentage()
                img2 = reportlab.platypus.Image(percentage_plot,4*inch,3*inch,hAlign='LEFT')
                Story.append(img2)
                Story.append(Spacer(1, 1))

                img3 = reportlab.platypus.Image('plot.png', 4 * inch, 3 * inch, hAlign='LEFT')
                Story.append(img3)
                Story.append(Spacer(1, 2))

                Story.append(Paragraph(notes, styles["Normal"]))
                doc.build(Story)
        except Exception as e:
            QMessageBox.critical(self,'Error',f'Something went wrong: {e}')

    def help(self):
        '''Open .pdf document with user guide'''
        help_document = 'Help.pdf'
        subprocess.Popen([help_document],shell=True)

#################################################################
### ThreadCLass - responsible for threading downloading
#################################################################

class ThreadClass(QtCore.QThread):

    any_signal = QtCore.pyqtSignal()
    def __init__(self,parent = None, index = 0):
        super(ThreadClass,self).__init__(parent)
        self.index = index
        self.isRunning = True

    def run(self):
        print('Starting thread', self.index)

    def stop(self):
        self.isRunning = False
        print('Stoping thread', self.index)
        self.terminate()

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    w= MyForm()
    w.show()
    sys.exit(app.exec_())