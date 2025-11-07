#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import sys
from PyQt5.QtSql import QSqlDatabase, QSqlTableModel, QSqlQuery
from PyQt5 import QtWidgets

class MainUi(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Nuclear DataBase")
        self.resize(800,600)

        self.table_widget = QtWidgets.QTableView()
        self.table_widget.setSortingEnabled(True)
        
        self.setCentralWidget(self.table_widget)
        self.show_db()

    def show_db(self):
        con = QSqlDatabase.addDatabase("QSQLITE")
        con.setDatabaseName('./ionic_data.db')
        con.open()
        self.model = QSqlTableModel()
        self.table_widget.setModel(self.model)
        self.model.setTable('IONICDATA')
        self.model.select()

app = QtWidgets.QApplication(sys.argv)
win = MainUi()
win.show()
sys.exit(app.exec())
