import sys
import json
import os

from PyQt5.QtWidgets import (QApplication, QMainWindow, QPushButton, QFileDialog, QTabWidget,
                             QVBoxLayout, QHBoxLayout, QWidget, QPlainTextEdit, QLineEdit, QLabel, QMessageBox, QToolBar, QFormLayout, QSizePolicy)
from PyQt5.QtGui import QFont

import matplotlib
matplotlib.use('Qt5Agg')

import importlib

from airplane import *

tabs = ['class_I_weight', 'matching', 'planform', 'airfoil', 'hld', 'sar']
for m in tabs:
    globals()[m] = __import__(m)

from airfoil import *

class TextBoxWriter:
    def __init__(self, text_box):
        self.text_box = text_box

    def write(self, text):
        self.text_box.insertPlainText(text)
        self.text_box.verticalScrollBar().setValue(self.text_box.verticalScrollBar().maximum())

    def flush(self):
        pass

class Uberflugzeug(QMainWindow):
    def __init__(self):
        super().__init__()

        self.tabs = [globals()[m].Tab(self) for m in tabs]

        self.init_ui()
        self.airplane = None
        self.unsaved_changes = False

        self.setMinimumSize(1500, 700) 

        try:
            with open(os.path.join(os.path.dirname(__file__), 'airplane.json')) as f:
                data = json.load(f)
            self.airplane = Airplane(data)
            self.init_ui_for_airplane()
        except:
            pass

    def init_ui(self):
        self.setWindowTitle('Uberflugzeug')

        main_layout = QHBoxLayout()

        self.tab_widget = QTabWidget()
        main_layout.addWidget(self.tab_widget)

        self.console_text_box = QPlainTextEdit(self)
        self.console_text_box.setFont(QFont("Courier", 10))  
        self.console_text_box.setReadOnly(True)
        self.console_text_box.setLineWrapMode(QPlainTextEdit.WidgetWidth) 
        self.console_text_box.setFixedWidth(400)

        sys.stdout = TextBoxWriter(self.console_text_box)

        main_layout.addWidget(self.console_text_box)

        container = QWidget()
        container.setLayout(main_layout)
        self.setCentralWidget(container)

        # Add toolbar
        toolbar = QToolBar("File")
        self.addToolBar(toolbar)

        load_action = QPushButton("Load Airplane")
        load_action.clicked.connect(self.load_airplane)
        toolbar.addWidget(load_action)

        save_action = QPushButton("Save Airplane")
        save_action.clicked.connect(self.save_airplane)
        toolbar.addWidget(save_action)

        self.setGeometry(300, 300, 800, 600)

        for t in self.tabs:
            t.init_ui()

    def init_ui_for_airplane(self):
        for t in self.tabs:
            t.init_ui_for_airplane()

    def recalculate(self):
        if not self.airplane: return

        self.console_text_box.clear()

        try:
            print("##############################################")
            print("#          AIRPLANE INFORMATION              #")
            print("##############################################")
            for i, t in enumerate(self.tabs):
                print()
                m = globals()[tabs[i]]
                importlib.reload(m)
                m.recalculate(t)
                print()
                print("--------------")

        except AirplaneValueError as v:
            self.console_text_box.clear()
            print(v.msg)

    def load_airplane(self):
        options = QFileDialog.Options()
        filename, _ = QFileDialog.getOpenFileName(self, "Load Airplane", "", "JSON Files (*.json);;All Files (*)", options=options)
        if filename:
            with open(filename, 'r') as file:
                data = json.load(file)
                self.airplane = Airplane(data)
                self.init_ui_for_airplane()
                self.unsaved_changes = False

    def save_airplane(self):
        options = QFileDialog.Options()
        filename, _ = QFileDialog.getSaveFileName(self, "Save Airplane", "", "JSON Files (*.json);;All Files (*)", options=options)
        if filename:
            with open(filename, 'w') as file:
                data = self.airplane.loaded_from_data
                for key, _ in data.items():
                    if not key.startswith("__comment_"):
                        data[key] = getattr(self.airplane, key)
                json.dump(data, file, indent=4)
                self.unsaved_changes = False
            return True
        else:
            return False

    def closeEvent(self, event):
        if self.unsaved_changes:
            reply = QMessageBox.question(self, 'Unsaved Changes',
                                         "You have unsaved changes. Do you want to save before exiting?",
                                         QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel, QMessageBox.Cancel)
            if reply == QMessageBox.Yes:
                if self.save_airplane():
                    event.accept()
                else:
                    event.ignore()
            elif reply == QMessageBox.No:
                event.accept()
            else:
                event.ignore()
        else:
            event.accept()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = Uberflugzeug()
    window.show()
    sys.exit(app.exec_())
