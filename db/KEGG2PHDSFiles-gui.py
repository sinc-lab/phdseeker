import tkinter as tk
#from tkinter.ttk import Progressbar
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import scrolledtext

import os
import urllib3
import requests
import json

from KEGG2PHDSFiles import KEGG


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

window.title("KEGG dataset builder")

window.geometry('650x700')


#----------------------

ORGANISM = tk.LabelFrame(window, padx=15, pady=10, text="")
ORGANISM.pack(padx=10, pady=5, side="top", fill="both", expand=False)

PATHWAYS = tk.LabelFrame(window, padx=15, pady=10, text="")
PATHWAYS.pack(padx=10, pady=5, side="top", fill="both", expand=False)

PROCESSING = tk.LabelFrame(window, padx=15, pady=10, text="")
PROCESSING.pack(padx=10, pady=5, side="top", fill="both", expand=False)

#----------------------


#=======================
# ORGANISM
#=======================
def select_organism():
    '''
    '''
    btn_organism.configure(text="Getting pathways...")
    btn_organism.update()
    
    listbox_pathways.delete(0,tk.END)
    
    
    items = map(int, listbox_organism.curselection())
    item = [ORGANISMS[int(item)] for item in items][0]
    organism_id = item[0]
    
    select_organism.organism_id = organism_id  # STORE VALUE
    
    lines = query.request('GET','http://rest.kegg.jp/list/pathway/{}'.format(organism_id))

    lines = lines.data.decode('UTF-8')

    lines = lines.strip().split('\n')

    PATHWAYS = []

    for line in lines:
        
        code,data = line.split('\t')
        
        code = code.split(':')[-1]
        
        name = data.split(' - ')[0].strip()
        
        PATHWAYS.append((code,name))

    #===========================

    PATHWAYS = sorted(PATHWAYS, key=lambda x: x[1])
    
    select_organism.PATHWAYS = PATHWAYS  # STORE VALUES

    for (code,name) in PATHWAYS:
        listbox_pathways.insert(tk.END, '{} -- (ID: {})'.format(name, code))
    
    listbox_pathways.select_set(0)
    
    #---------------------------
    # AGREGO AL CAMPO DE TEXTO
    #---------------------------
    _code,_name = item
    txt_organism.configure(text='{} -- (ID: {})'.format(_name, _code))
    #---------------------------
    
    btn_organism.configure(text="Selected", bg='orange')

#----------------------

group_organism = tk.LabelFrame(ORGANISM, padx=15, pady=10, text=" Available Organisms ")
group_organism.pack(padx=10, pady=5, fill='both')

#----------------------

scrollbar_organism = tk.Scrollbar(group_organism, orient='vertical')

listbox_organism = tk.Listbox(group_organism,
                              selectmode='browse',
                              yscrollcommand=scrollbar_organism.set,
                              font=('Tahoma', 10), height=10)

scrollbar_organism.config(command=listbox_organism.yview)
listbox_organism.pack(side='left', fill='both', expand=True)
scrollbar_organism.pack(side='left', fill=tk.Y)

#---------------------------------------

query = urllib3.PoolManager()

lines = query.request('GET','http://rest.kegg.jp/list/genome')

lines = lines.data.decode('UTF-8')

lines = lines.strip().split('\n')

ORGANISMS = []

for line in lines:
        
    if ('; ' in line):
        data,name = line.split('; ')
        
        name = name.strip()
        
        code = data.split(',')[0]
        code = code.split('\t')[-1]
        
        ORGANISMS.append((code,name))

#===========================

ORGANISMS = sorted(ORGANISMS, key=lambda x: x[1])
ORGANISMS.insert(0,('rn', 'All available reactions'))


for (code,name) in ORGANISMS:
    listbox_organism.insert(tk.END, '{} -- (ID: {})'.format(name, code))

listbox_organism.select_set(0)

#---------------------------------

group_organism2 = tk.LabelFrame(ORGANISM, padx=10, pady=10, text="", borderwidth=0, highlightthickness=0)
group_organism2.pack(side="right", padx=0, pady=0, fill="both", expand=True)

btn_organism = tk.Button(group_organism2, text="Select", bg="#02B844", fg="black", command=select_organism)
btn_organism.pack(side="right", fill="both", expand=False)
btn_organism_help = CreateToolTip(btn_organism, 'Select the organism to build dataset')

txt_organism = tk.Label(group_organism2, text='', height=1, bg='#bcbcbc', fg="black")
txt_organism.pack(side="right", fill="both", expand=True)


#=======================
# PATHWAY
#=======================
def select_pathways():
    '''
    '''
    
    btn_pathways.configure(text="Storing information...")
    btn_pathways.update()
    
    pathways = map(int, listbox_pathways.curselection())
    pathways_id = [select_organism.PATHWAYS[int(pathway)][0] for pathway in pathways]
    
    select_pathways.all = False
    if len(select_organism.PATHWAYS) == len(pathways_id):
        select_pathways.all = True
    
    select_pathways.pathways_id = pathways_id  # STORE VALUE
    
    btn_pathways.configure(text="Selected", bg='orange')

#---------------------


group_pathways1 = tk.LabelFrame(PATHWAYS, padx=15, pady=10, text=" Available Pathways ")
group_pathways1.pack(padx=10, pady=5, fill='both')

scrollbar_pathways = tk.Scrollbar(group_pathways1, orient='vertical')

listbox_pathways = tk.Listbox(group_pathways1,
                              selectmode='extended',
                              yscrollcommand=scrollbar_pathways.set,
                              font=('Tahoma', 10), height=10)

scrollbar_pathways.config(command=listbox_pathways.yview)
listbox_pathways.pack(side='left', fill='both', expand=True)
scrollbar_pathways.pack(side='left', fill=tk.Y)

#---------------------------------------

group_pathways2 = tk.LabelFrame(PATHWAYS, padx=15, pady=10, text="", borderwidth=0, highlightthickness=0)
group_pathways2.pack(side="right", padx=0, pady=0)

btn_pathways = tk.Button(group_pathways2, text="Select", bg="#02B844", fg="black", command=select_pathways)
btn_pathways.pack(side="right", fill="both", expand=True)

btn_pathways_help = CreateToolTip(btn_pathways, 'Select pathways to build the dataset. "Ctrl" and "Shift" keys can be used to select multiple options.')



#=======================
# PROCESSING
#=======================
def load_dict():
    '''
    '''
    
    load_dict.dict_file = filedialog.askopenfilename(initialdir = current_path + "/db",
                                                     title = "Select dictionary ...",
                                                     filetypes = (("JSON files","*.json"),("all files","*.*")))
    path_dict, filename_dict = os.path.split(load_dict.dict_file)
    txt_dict.configure(text='{}'.format(filename_dict.strip()))  # SPLIT and SHOW FILENAME


load_dict.dict_file = None

#-----------------------------
def run_processing():
    
    value = False
    if (load_dict.dict_file is None):
        value = messagebox.askyesno(title='Confirm action!', message="You haven't specified a dictionary. Do you want to continue and build one? Please note that this process may take some time. However, this dictionary can be used later for other datasets.")
    
    if (value == True) or (load_dict.dict_file is not None):
        
        btn_processing.configure(text="Building...", bg='orange')
        btn_processing.update()
        
        window.wm_state('iconic')
        
        
        # CHECK FOR REVERSIBILITY OPTION
        use_reversibility = chk_processing.instate(['selected'])
        
        if not os.path.exists('KEGG'):
            os.makedirs('KEGG')
            
        folder_selected = filedialog.askdirectory(initialdir = current_path + "/KEGG",
                                                title = "Select destination folder...")
        
        #-------------------------------------------------
        
        DATA = KEGG(usedict=load_dict.dict_file)
        
        DATA.destination_folder = folder_selected
        
        selected_pathways = ','.join(select_pathways.pathways_id) if (select_pathways.all == False) else 'all'
        
        DATA.download_organism_data(organism=select_organism.organism_id,  # STORED VALUE
                                    pathways=selected_pathways,
                                    use_reversibility_info=use_reversibility
                                    )

        DATA.save_organism_to_file()
        #-------------------------------------------------
        
        messagebox.showinfo('', 'Success!')
        
        btn_processing.configure(text="Build", bg="#02B844")
        btn_organism.configure(text="Select", bg="#02B844")
        btn_pathways.configure(text="Select", bg="#02B844")
    
    
    

group_processing1 = tk.LabelFrame(PROCESSING, padx=5, pady=2, text="", borderwidth=0, highlightthickness=0)
group_processing1.pack(padx=5, pady=1, side="top", fill="both", expand=False)


btn_dict = tk.Button(group_processing1, text="Select dictionary", bg="#02B844", fg="black", command=load_dict, height=1)
btn_dict.pack(side="left", fill="both", expand=False, padx=0, pady=0)
btn_dict_help = CreateToolTip(btn_dict, 'Select the JSON file with information to translate database codes to compound names (optional file)')

txt_dict = tk.Label(group_processing1, text='', height=1, width=60, bg='#bcbcbc', fg="black")
txt_dict.pack(side="right", fill="both", expand=True)

#---------------------------------------

group_processing2 = tk.LabelFrame(PROCESSING, padx=5, pady=2, text="", borderwidth=0, highlightthickness=0)
group_processing2.pack(padx=5, pady=1, side="top", fill="both", expand=False)

chk_processing = ttk.Checkbutton(group_processing2, text="Use reversibility information")
chk_processing.pack(side='left', fill=tk.X, expand=False, padx=5)

chk_processing_help = CreateToolTip(chk_processing, 'Check if all reactions should be processed as reversibile')


btn_processing = tk.Button(group_processing2, text="Build", bg="#02B844", fg="black", command=run_processing)
btn_processing.pack(side="right", fill=tk.X, expand=False)

btn_processing_help = CreateToolTip(btn_processing, 'Build dataset')


#---------------------


window.mainloop()
