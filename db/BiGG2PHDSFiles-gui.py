import tkinter as tk
from tkinter.ttk import Progressbar
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import scrolledtext

import os
import requests
import json

from BiGG2PHDSFiles import BiGG2PHDSFiles


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

print('\nConnecting to the database and downloading information.\n\nThis can take a few seconds...\n')


#=============================
# Build and configure window
#=============================
current_path = os.getcwd()

window = tk.Tk()

window.title("BiGG dataset builder")

window.geometry('650x680')


container = tk.LabelFrame(window, padx=15, pady=10, text="")

#----------------------

viewer = tk.LabelFrame(window, padx=15, pady=10, text=" Download BiGG model ")
viewer.pack(padx=10, pady=5, side="top", fill="both", expand=False)

viewer1 = tk.LabelFrame(viewer, padx=15, pady=10, text="")
viewer1.pack(padx=10, pady=5, side="top", fill="both", expand=False)

viewer2 = tk.LabelFrame(viewer, padx=0, pady=2, text="", borderwidth=0, highlightthickness=0)
viewer2.pack(padx=10, pady=5, side="top", fill="both", expand=False)

#----------------------


processing = tk.LabelFrame(window, padx=15, pady=5, text=" Process BiGG model ")
processing.pack(padx=10, pady=5, side="top", fill="both", expand=False)

processing1 = tk.LabelFrame(processing, padx=0, pady=2, text="", borderwidth=0, highlightthickness=0)
processing1.pack(padx=10, pady=5, side="top", fill="both", expand=False)

processing2 = tk.LabelFrame(processing, padx=0, pady=2, text="", borderwidth=0, highlightthickness=0)
processing2.pack(padx=10, pady=5, side="top", fill="both", expand=False)

#----------------------


#=========================
# VIEWER
#=========================
scrollbar = tk.Scrollbar(viewer1, orient='vertical')

listbox = tk.Listbox(viewer1,
                     selectmode='browse',
                     yscrollcommand=scrollbar.set,
                     font=('Tahoma', 10), height=26)

scrollbar.config(command=listbox.yview)
listbox.pack(side='left', fill='both', expand=True)
scrollbar.pack(side='left', fill=tk.Y)

#---------------------------------------

r = requests.get('http://bigg.ucsd.edu/api/v2/models')
models = r.json()

MODELS = []

for model in models['results']:
    
    MODELS.append((model['bigg_id'],
                   model['organism'],
                   model['metabolite_count'],
                   model['reaction_count'],
                   model['gene_count']))
    
    
MODELS = sorted(MODELS, key=lambda x: x[1])

for (code,organism,metabolites,reactions,genes) in MODELS:
    
    spacing = ''.join([' ']*(20-len(code)))
    listbox.insert(tk.END, u'{}  (ID: {})   [{} metabolites, {} reactions, {} genes]'.format(organism,code,metabolites,reactions,genes))

listbox.select_set(0)


#=========================
# DOWNLOAD BUTTON
#=========================
def download_button():
    '''
    '''
    
    btn_download.configure(text="Downloading...")
    btn_download.update()
    
    
    # SELECT ITEM
    items = map(int, listbox.curselection())
    item = [MODELS[int(item)] for item in items][0]
    bigg_id = item[0]
    
    #-------------------------------------------------------------------------
    # DOWNLOAD DATA
    #------------------
    url = 'http://bigg.ucsd.edu/api/v2/models/{}/download'.format(bigg_id)
    r = requests.get(url)
    model = r.json()
    
    filename = '{}.json'.format(bigg_id)
    
    if not os.path.exists('BiGG'):
        os.makedirs('BiGG')
    
    
    folder_selected = filedialog.askdirectory(initialdir = current_path + "/BiGG",
                                              title = "Select destination folder...")
    
    
    with open(os.path.join(folder_selected, filename), 'w') as fp:
        json.dump(model, fp, indent=4, sort_keys=True)
    
    #-------------------------------------------------------------------------
    
    messagebox.showinfo('', 'Success!')
    
    #------------------------
    
    btn_download.configure(text="Run", bg="#02B844")


btn_download = tk.Button(viewer2, text="Download", bg="#02B844", fg="black", height=1, width=10, command=download_button)
btn_download.pack(side="right", fill=None, expand=False)

btn_download_help = CreateToolTip(btn_download, 'Select an organism to download the model')


#--------------------------------


#=========================
# LOAD MODEL
#=========================
def load_model():
    '''
    Define action when button is press.
    '''
    load_model.model_file = filedialog.askopenfilename(initialdir = current_path + "/BiGG",
                                                       title = "Select BiGG model...",
                                                       filetypes = (("JSON files","*.json"),("all files","*.*")))
    
    load_model.path_model, load_model.filename_model = os.path.split(load_model.model_file)
    txt_model.configure(text='{}'.format(load_model.filename_model.strip()))  # SPLIT and SHOW FILENAME


btn_model = tk.Button(processing1, text="Load file", bg="#02B844", fg="black", command=load_model)
btn_model.pack(side="left", fill="both", expand=False)

txt_model = tk.Label(processing1, text='', height=1, width=60, bg='#bcbcbc', fg="black")
txt_model.pack(side="right", fill="both", expand=True)

btn_model_help = CreateToolTip(btn_model, 'Select a file with the model of the organism to build the dataset')



#=========================
# PROCESSING BUTTON
#=========================
def processing_button():
    '''
    '''
    
    btn_processing.configure(text="Building...")
    btn_processing.update()
    
    #-------------------------------------------------------------------------
    
    allreversible = btn_chk.instate(['selected'])
    
    BiGG2PHDSFiles(os.path.join(load_model.path_model, load_model.filename_model), allreversible)
    
    messagebox.showinfo('', 'Success!')
    
    #------------------------
    
    btn_processing.configure(text="Run", bg="#02B844")


#=========================
# BOX SELECTOR
#=========================
btn_chk = ttk.Checkbutton(processing2, text="All reversible")
btn_chk.pack(side='left', fill=tk.X, expand=False, padx=5)

btn_chk_help = CreateToolTip(btn_chk, 'Check if all reactions should be processed as reversibile')


#=========================
# PROCESSING BUTTON
#=========================
btn_processing = tk.Button(processing2, text="Build", bg="#02B844", fg="black", height=1, width=10, command=processing_button)
btn_processing.pack(side="right", fill=tk.X, expand=False)

btn_processing_help = CreateToolTip(btn_processing, 'Build dataset')


# ----------------


window.mainloop()
