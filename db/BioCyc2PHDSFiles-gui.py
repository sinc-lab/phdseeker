import tkinter as tk
from tkinter.ttk import Progressbar
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import scrolledtext

import os
import urllib3
import requests
import json

from BioCyc2PHDSFiles import BioCyc


####################################################
class CreateToolTip(object):
    """
    create a tooltip for a given widget
    """
    def __init__(self, widget, text='widget info'):
        self.waittime = 500     #miliseconds
        self.wraplength = 180   #pixels
        self.widget = widget
        self.text = text
        self.widget.bind("<Enter>", self.enter)
        self.widget.bind("<Leave>", self.leave)
        self.widget.bind("<ButtonPress>", self.leave)
        self.id = None
        self.tw = None

    def enter(self, event=None):
        self.schedule()

    def leave(self, event=None):
        self.unschedule()
        self.hidetip()

    def schedule(self):
        self.unschedule()
        self.id = self.widget.after(self.waittime, self.showtip)

    def unschedule(self):
        id = self.id
        self.id = None
        if id:
            self.widget.after_cancel(id)

    def showtip(self, event=None):
        x = y = 0
        x, y, cx, cy = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 20
        # creates a toplevel window
        self.tw = tk.Toplevel(self.widget)
        # Leaves only the label and removes the app window
        self.tw.wm_overrideredirect(True)
        self.tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(self.tw, text=self.text, justify='left',
                       background="#ffffff", relief='solid', borderwidth=1,
                       fg="#000000", padx=5,pady=5,
                       wraplength = self.wraplength)
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tw
        self.tw= None
        if tw:
            tw.destroy()
###################################################


#=============================
# Build and configure window
#=============================
current_path = os.getcwd()

window = tk.Tk()

window.title("BioCyc dataset builder")

window.geometry('650x280')


#----------------------

REACTIONS = tk.LabelFrame(window, padx=0, pady=5, text="", borderwidth=0, highlightthickness=0)
REACTIONS.pack(padx=10, pady=5, side="top", fill="both", expand=False)

COMPOUNDS = tk.LabelFrame(window, padx=0, pady=5, text="", borderwidth=0, highlightthickness=0)
COMPOUNDS.pack(padx=10, pady=5, side="top", fill="both", expand=False)

PROCESSING = tk.LabelFrame(window, padx=0, pady=5, text="", borderwidth=0, highlightthickness=0)
PROCESSING.pack(padx=10, pady=5, side="top", fill="both", expand=False)

#----------------------


#=======================
# REACTIONS
#=======================
def select_reactions_file():
    '''
    '''
    btn_select_reactions_file.configure(text="Loading...")
    btn_select_reactions_file.update()
    

    #===========================
    select_reactions_file.reactions_file = filedialog.askopenfilename(initialdir = current_path + "/db",
                                                                      title = "Select reactions file...",
                                                                      filetypes = (("DAT files","*.dat"),("all files","*.*")))
    select_reactions_file.path_reactions, select_reactions_file.filename_reactions = os.path.split(select_reactions_file.reactions_file)
    txt_select_reactions_file.configure(text='{}'.format(select_reactions_file.filename_reactions.strip()))  # SPLIT and SHOW FILENAME
    #===========================
    
    
    btn_select_reactions_file.configure(text="Load file")


#----------------------


group_reactions = tk.LabelFrame(REACTIONS, padx=15, pady=20, text=" Reactions ")
group_reactions.pack(padx=10, pady=5, fill='both')

btn_select_reactions_file = tk.Button(group_reactions, text="Load file", bg="#02B844", fg="black", command=select_reactions_file, height=1)
btn_select_reactions_file.pack(side="left", fill="both", expand=True)

txt_select_reactions_file = tk.Label(group_reactions, text='', height=1, width=60, bg='#bcbcbc', fg="black")
txt_select_reactions_file.pack(side="left", fill="both", expand=True)

btn_select_reactions_help = CreateToolTip(btn_select_reactions_file, 'Select the reactions file')


#=======================
# COMPOUNDS
#=======================
def select_compounds_file():
    '''
    '''
    
    btn_select_compounds_file.configure(text="Loading...")
    btn_select_compounds_file.update()
    

    #===========================
    select_compounds_file.compounds_file = filedialog.askopenfilename(initialdir = current_path + "/db",
                                                                      title = "Select compounds file...",
                                                                      filetypes = (("DAT files","*.dat"),("all files","*.*")))
    select_compounds_file.path_compounds, select_compounds_file.filename_compounds = os.path.split(select_compounds_file.compounds_file)
    txt_select_compounds_file.configure(text='{}'.format(select_compounds_file.filename_compounds.strip()))  # SPLIT and SHOW FILENAME
    #===========================
    
    
    btn_select_compounds_file.configure(text="Load file")



group_compounds = tk.LabelFrame(COMPOUNDS, padx=15, pady=20, text=" Compounds ")
group_compounds.pack(padx=10, pady=5, fill='both')

btn_select_compounds_file = tk.Button(group_compounds, text="Load file", bg="#02B844", fg="black", command=select_compounds_file, height=1)
btn_select_compounds_file.pack(side="left", fill="both", expand=True)

txt_select_compounds_file = tk.Label(group_compounds, text='', height=1, width=60, bg='#bcbcbc', fg="black")
txt_select_compounds_file.pack(side="left", fill="both", expand=True)

btn_select_compounds_help = CreateToolTip(btn_select_compounds_file, 'Select the compounds file')


#=======================
# PROCESSING
#=======================

#-----------------------------
def run_processing():
    
    btn_processing.configure(text="Building...")
    btn_processing.update()
    
    
    # CHECK FOR REVERSIBILITY OPTION
    use_reversibility_info = not chk_processing.instate(['selected'])
    
    if not os.path.exists('BioCyc'):
        os.makedirs('BioCyc')
        
    folder_selected = filedialog.askdirectory(initialdir = current_path + "/BioCyc",
                                              title = "Select destination folder...")
    
    #-------------------------------------------------
    
    DATA = BioCyc(select_compounds_file.compounds_file,select_reactions_file.reactions_file)

    DATA.save_DB_to_file(use_reversibility_info=use_reversibility_info)
    
    #-------------------------------------------------
    
    messagebox.showinfo('', 'Success!')
    
    btn_processing.configure(text="Build")
    
    
    

group_processing = tk.LabelFrame(PROCESSING, padx=15, pady=2, text="", borderwidth=0, highlightthickness=0)
group_processing.pack(padx=5, pady=1, side="top", fill="both", expand=False)

chk_processing = ttk.Checkbutton(group_processing, text="Use reversibility information")
chk_processing.pack(side='left', fill=tk.X, expand=False, padx=5)

chk_processing_help = CreateToolTip(chk_processing, 'Check if all reactions should be processed as reversibile')

#-------------------

btn_processing = tk.Button(group_processing, text=" Build ", bg="#02B844", fg="black", command=run_processing, height=2, width=8)
btn_processing.pack(side="right", fill=tk.X, expand=False)

btn_processing_help = CreateToolTip(btn_processing, 'Build dataset')


#---------------------


window.mainloop()
