import datetime
import subprocess
import sys
import traceback
import uuid

import pandas as pd
import reportlab.platypus
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from Bio import Entrez, SeqIO, Seq
from PyQt5.QtWidgets import QApplication, QTableWidgetItem, QMainWindow, \
    QScrollArea, QFileDialog, QMessageBox, QLabel
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image

from downloading_dialog import *
from plot_window import *
from table_model import *
from layout.main_window_layout import *

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
        self.ui.actionOpenDatabase.triggered.connect(self.open_database)
        self.ui.actionGenerete_plot.triggered.connect(self.open_plot_window)
        self.ui.actionGenerate_report.triggered.connect(self.generate_report)
        self.ui.actionHelp.triggered.connect(self.help)

        self.ui.pushButtonSearchInEntrez.clicked.connect(self.search_in_entrez)
        self.ui.pushButtonExportToTxt.clicked.connect(self.save_query_result_to_txt)
        self.ui.pushButtonSearchInDatabse.clicked.connect(self.search_in_database)
        self.ui.pushButtonLoadFromTxt.clicked.connect(self.load_from_txt)

        self.ui.tableWidgetSelectedData.cellClicked.connect(self.get_clicked_cell)
        self.query_result = ''
        self.query_translated = ''
        self.data_loaded = False
        self.dialog = MyDialog()
        self.plot_creator = MyPlotDialog()
        self.plot_creator.ui.pushButtonAddToWorkspace.clicked.connect(self.add_chart_to_workspace)
        self.file_name_no_csv = ''
        self.thread = {}
        self.show()

    @staticmethod
    def help():
        """
        :return: opens .pdf document "Help.pdf" with user guide
        """
        help_document = 'Help.pdf'
        subprocess.Popen([help_document], shell=True)

    def open_downloading_dialog(self):
        """
        :return: opens downloading dialog
        """
        self.dialog.show()

    def open_plot_window(self):
        """
        :return:opens plotting dialog. In case of lack of data return error
        """
        try:
            self.plot_creator.data_from_database = self.df_new
            self.plot_creator.show()
        except Exception as e:
            QMessageBox.critical(self, 'Error', 'No data to plot. Load database first')

    def add_chart_to_workspace(self):
        """
        :return: adds plot from plotting dialog to main window
        """
        try:
            self.ui.labelPlottoWorkspace.setPixmap(self.plot_creator.plot_png)
            text = self.plot_creator.ui.textEdit.toPlainText()
            self.ui.textEditNoteForPlot.setText(text)
        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Something went wrong: {e}')

    def open_database(self):
        """
        :return: opens .csv database file
        """
        if self.plot_creator.isVisible():
            self.plot_creator.close()
        try:
            file, _ = QFileDialog.getOpenFileName(self, "Open file", "", "All files (*);;CSV files (*.csv)")
            if file:
                dane = pd.read_csv(file, keep_default_na=False)
                self.df_new = pd.DataFrame(dane)
                self.df_new.rename(columns={'#Organism Name': 'Organism_name'}, inplace=True)
                self.df_new.sort_values(by=['Organism_name'], ascending=True, inplace=True)
                self.df_new.insert(0, 'Index', [i for i in range(1, self.df_new.shape[0] + 1)], True)
                self.data_loaded = True
                num_rows, num_cols = self.df_new.shape
                data_tab = TabelaModel(self.df_new)
                self.ui.tableView.setModel(data_tab)
                file_name = str(file).split('/')
                self.file_name_no_csv = file_name[-1].split('.')
                self.ui.labelDataName.setText('Data view ' + self.file_name_no_csv[0])
                self.ui.labelDatabaseSizeInfo.setText(f'Database size: {num_cols - 1} columns, {num_rows} rows')
                QMessageBox.information(self, 'Info', 'Database successfully loaded!')
        except Exception as e:
            QMessageBox.critical(self, 'Error', 'Chose correct file extension: csv')

    def search_in_database(self):
        """
        :return: filters data for table by organism name or BioProject id
        """
        try:
            if self.data_loaded:
                search_text = self.ui.lineEditSearchInDatabase.text()
                data_from_database = pd.DataFrame()
                try:
                    if "BioProject" in self.df_new:
                        data_from_database = pd.DataFrame(self.df_new.query(
                            f'Organism_name == ["{search_text}"] or BioProject == ["{search_text}"]'))
                    else:
                        data_from_database = pd.DataFrame(self.df_new.query(
                            f'Organism_name == ["{search_text}"]'))
                    num_rows, num_cols = data_from_database.shape
                    self.ui.tableWidgetSelectedData.setColumnCount(num_cols)
                    self.ui.tableWidgetSelectedData.setRowCount(num_rows)
                    self.ui.tableWidgetSelectedData.setHorizontalHeaderLabels(data_from_database.columns)
                    for row in range(num_rows):
                        for column in range(num_cols):
                            self.ui.tableWidgetSelectedData.setItem(row, column, QTableWidgetItem(
                                str(data_from_database.iloc[row, column])))
                    self.ui.labelNumbeOfRows.setText(
                        'Number of results: ' + str(self.ui.tableWidgetSelectedData.rowCount()))
                except Exception:
                    QMessageBox.information(self, 'Info', 'Try to search for another data')
        except Exception as e:
            QMessageBox.warning(self, 'Warning', e)

    def get_clicked_cell(self, row, column):
        """
        :param row: vertical position of clicked cell
        :param column: horizontal position of clicked cell
        :return: when cell is clicked, copy content to LineEdit below
        """
        item_from_table = self.ui.tableWidgetSelectedData.item(row, column).text()
        self.ui.lineEditSelectedForQuery.setText(item_from_table)

    def search_in_entrez(self):
        """
        :return: sends query to chosen Entrez databse. Possible databases: nucleotide, assemble, BioProject
        """
        chosen_entrez_database = self.ui.comboBoxWhichDatabase.itemText(
            self.ui.comboBoxWhichDatabase.currentIndex()).lower()
        try:
            Entrez.email = f"{self.ui.lineEditEntrezEmail.text()}"
            query = f"{self.ui.lineEditSelectedForQuery.text()}"
            handle = Entrez.esearch(db=f"{chosen_entrez_database}", term=query, idtype='acc', retmax=5)
            result = Entrez.read(handle)
            handle = Entrez.efetch(db=f"{chosen_entrez_database}", id=result["IdList"], rettype="fasta")
            handle = Entrez.efetch(db=f"{chosen_entrez_database}", id=result["IdList"][0], rettype="gb")
            self.query_translated = SeqIO.read(handle, 'genbank')
            handle = Entrez.efetch(db=f"{chosen_entrez_database}", id=result["IdList"], rettype="gb")
            self.query_result = handle.read()
            self.ui.textEditQueryResult.setText(self.query_result)
            handle.close()
        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Wrong query:\n{e}\nTry to use another database for this query.')

    def save_query_result_to_txt(self):
        """
        :return: saves data to .txt file
        """
        try:
            if len(self.query_result) > 0:
                file, _ = QFileDialog.getSaveFileName(self, "Save file", "", "All files (*);;TXT files (*.txt)")
                if file:
                    saving = open(file, 'w')
                    saving.write(self.query_result)
                    saving.close()
        except Exception as e:
            QMessageBox.critical(self, 'Info', f'No data to export.\n{e}')

    def load_from_txt(self):
        """
        :return: loads data from .txt file
        """
        try:
            file, _ = QFileDialog.getOpenFileName(self, "Open file", "", "All files (*);;TXT files (*.txt)")
            if file:
                try:
                    self.ui.textEditQueryResult.clearFocus()
                    with open(file, 'r') as f:
                        data = f.read()
                    self.ui.textEditQueryResult.setText(data)
                    f.close()
                except Exception as e:
                    QMessageBox.critical(self, 'Error', f'Wrong file extension: {e}')
        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Something went wrong. Try to load file again\n{e}')

    def add_amino_percentage(self):
        """
        :return: generates plot with amino percentage from entrez query
        """
        try:
            data = []
            nucleotides = ['A', 'T', 'C', 'G']
            values = {}
            for nucleotide in nucleotides:
                count = self.query_translated.seq.count(nucleotide)
                values[nucleotide] = count
            all = sum(list(values.values()))
            for i in values:
                value = values[i]
                percent = value / all * 100
                data.append(percent)
            explode = (0, 0.1, 0, 0)
            fig1, ax1 = plt.subplots()
            plot = ax1.pie(data, explode=explode, labels=nucleotides, autopct='%1.1f%%',
                           shadow=True, startangle=90)
            ax1.axis('equal')
            plt.title(f'Nucleotide percentage of {self.query_translated.annotations["organism"]}')
            plt.savefig('plot2.png')
            plot_file = 'plot2.png'
            return plot_file
        except Exception as e:
            QMessageBox.critical(self, 'Error', f'No data! {e}')

    def generate_report(self):
        """
        :return:
        Generates report as a .pdf extension file. Data: report id, date&time of report generation,
        entrez query results, plots: nucleotide percentage of organism from query
        and plot from MyPlotDialog class.
        """
        try:
            org_name = '_'.join(self.query_translated.annotations["organism"].split())
            file, _ = QFileDialog.getSaveFileName(self, "Save file", f'{self.file_name_no_csv[0]}_{org_name}',
                                                  "All files (*);;PDF files (*.pdf)")
            if file:
                doc = SimpleDocTemplate(file, pagesize=letter,
                                        rightMargin=25, leftMargin=25,
                                        topMargin=25, bottomMargin=18)
                logo = 'C:/Users/miko5/Desktop/TDS/genome-info/icons/dna_icon_132453.ico'
                story = []
                logo_image = reportlab.platypus.Image(filename=logo, height=0.5 * inch, width=0.5 * inch,
                                                      hAlign='RIGHT')
                name = 'GenomeInfo'
                report_id = f'Report ID: {uuid.uuid4()}'
                now = datetime.datetime.now()
                generation_date_and_time = now.strftime("%d/%m/%Y %H:%M:%S")
                date = f'Date and time of report generation: {generation_date_and_time}'
                result_txt = 'Result of Entrez query:'
                result_info_id = f'ID number: {self.query_translated.id}'
                result_info_organism = f'Organism: {self.query_translated.annotations["organism"]}'
                try:
                    result_info_seq = f'Sequence: {self.query_translated.seq[:100]}...'
                except Exception as e:
                    QMessageBox.warning(self, 'Error',
                                        f'Lack of access to organism\' sequence: {traceback.format_exc()}')
                finally:
                    result_info_seq = 'None'
                result_info_seq_length = f'Sequence length: {len(self.query_translated.seq)}'
                result_info_description = f'Description: {self.query_translated.description}'
                information = [name, report_id, date, result_txt]
                results = [result_info_id, result_info_organism, result_info_seq, result_info_seq_length,
                           result_info_description]
                notes = self.ui.textEditNoteForPlot.toPlainText()

                story.append(logo_image)
                styles = getSampleStyleSheet()
                styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))

                for info in information:
                    story.append(Paragraph(info, styles["Normal"]))
                    story.append(Spacer(1, 12))

                for result in results:
                    story.append(Spacer(1, 12))
                    story.append(Paragraph(result, styles["Justify"]))

                story.append(Spacer(1, 10))
                percentage_plot = self.add_amino_percentage()
                percentage_plot_image2 = reportlab.platypus.Image(percentage_plot, 4 * inch, 3 * inch, hAlign='LEFT')
                story.append(percentage_plot_image2)
                story.append(Spacer(1, 1))

                percentage_plot_image3 = reportlab.platypus.Image('plot.png', 4 * inch, 3 * inch, hAlign='LEFT')
                story.append(percentage_plot_image3)
                story.append(Spacer(1, 2))

                story.append(Paragraph(notes, styles["Normal"]))
                doc.build(story)
        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Something went wrong: {e}')


class ThreadClass(QtCore.QThread):
    any_signal = QtCore.pyqtSignal()

    def __init__(self, parent=None, index=0):
        super(ThreadClass, self).__init__(parent)
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
    w = MyForm()
    w.show()
    sys.exit(app.exec_())
