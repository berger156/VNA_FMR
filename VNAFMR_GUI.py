# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 16:55:11 2017

@author: ajb7

This simple app creates a main app window with a single button that can be clicked.
"""

import sys
from PyQt5 import QtWidgets, uic
import VNAFMR_FitFun
 
qtCreatorFile = "VNAFMR_UI.ui" # Enter UI filename here.
 
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)
 
class MyApp(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.run_button.clicked.connect(self.runFits) 
        self.fileSelect_button.clicked.connect(self.getFile)
        
    def getFile(self):
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(self,'Open file')
        self.filename_display.setText(str(fname))
        
    def runFits(self):
        #filename
        fname = str(self.filename_display.toPlainText())
        
        resultsDir = str(self.outputDir.toPlainText())
        
        #data columns
        field_col = int(self.field_col.value())
        reS21_col = int(self.reS21_col.value())
        imS21_col = int(self.imS21_col.value())
        
        if self.yellow_cols_bool.isChecked():
            print("Yellow Magnet")
            field_col = 6
            reS21_col = 7
            imS21_col = 8            
        if self.blue_cols_bool.isChecked():
            print("Blue Magnet")
            field_col = 13
            reS21_col = 18
            imS21_col = 19  
            
        #frequency range
        f_min = float(self.f_min.toPlainText())
        f_max = float(self.f_max.toPlainText())  
        
        f_min_lin = float(self.f_min_lin.toPlainText())
        f_max_lin = float(self.f_max_lin.toPlainText())
        
        #windowing
        if self.window_bool.isChecked():
            window = True
            windowSize = float(self.window_size.toPlainText())
        else:
            window = False
            windowSize = 10
           
        VNAFMR_FitFun.fitAll(f_min, f_max, f_min_lin, f_max_lin, field_col, reS21_col, imS21_col, window, windowSize, fname, resultsDir)                
                
if __name__ == "__main__":
    #this only executes if VNAFMR_GUI.py is run by itself, and not imported from another module
    app = 0 #to prevent Python kernel from crashing
    app = QtWidgets.QApplication(sys.argv)
    window = MyApp()
    window.show()
    #sys.exit(app.exec_())
    app.exec_()