from PyQt5.QtWidgets import (QApplication, QMainWindow, QPushButton, QFileDialog, QTabWidget, QFrame,
                             QVBoxLayout, QHBoxLayout, QWidget, QPlainTextEdit, QLineEdit, QLabel, QMessageBox, QToolBar, QFormLayout, QSizePolicy)
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

class Tab(QFormLayout):
    def __init__(self, uber, name):
        super().__init__()
        self.uber = uber
        self.name = name

    def create_input_box(self, label, attr):
        label = QLabel(f'{label}: ')

        initial_value = getattr(self.uber.airplane, attr)

        line_edit = QLineEdit()
        line_edit.setText(str(initial_value))
        line_edit.setMaximumWidth(100)

        def on_editing_finished():
            try:
                value = float(line_edit.text())
                setattr(self.uber.airplane, attr, value)
                self.uber.unsaved_changes = True
            except:
                pass
        line_edit.textChanged.connect(on_editing_finished)
        self.addRow(label, line_edit)

    def init_ui(self):
        tab = QWidget()
        self.uber.tab_widget.addTab(tab, self.name)

        self.tab_layout = QHBoxLayout() 
        tab.setLayout(self.tab_layout)

    def init_ui_for_airplane(self, plot_nrows=1, plot_ncols=1):
        info_tab_general_container = QWidget(self.uber)
        info_tab_general_container.setLayout(self)
        info_tab_general_container.setFixedWidth(300)
        self.tab_layout.addWidget(info_tab_general_container)

        recalc_button = QPushButton('Recalculate', self.uber)
        recalc_button.clicked.connect(lambda: self.uber.recalculate())
        self.addRow(recalc_button)

        canvas_and_toolbar = QVBoxLayout()
        
        if plot_nrows != 1 or plot_ncols != 1:
            self.fig, self.ax = plt.subplots(nrows=plot_nrows, ncols=plot_ncols, figsize=(10, 8))
        else:
            self.fig, self.ax = plt.subplots(figsize=(10, 8))

        self.canvas = FigureCanvas(self.fig)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.canvas.updateGeometry()

        canvas_and_toolbar.addWidget(self.canvas)

        toolbar = NavigationToolbar(self.canvas, self.uber)
        canvas_and_toolbar.addWidget(toolbar)

        self.tab_layout.addLayout(canvas_and_toolbar)

    def separator(self, title=None):
        line = QFrame()
        line.setFrameShape(QFrame.HLine)  
        line.setFrameShadow(QFrame.Sunken)
        line.setFixedHeight(2)
        self.addRow(line)

        if title is not None:
            label = QLabel(title)
            label.setStyleSheet("font-weight: bold")
            self.addRow(label)

