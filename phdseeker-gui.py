import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import scrolledtext

import multiprocessing


import os
import sys

import phdseeker

import webbrowser
import yaml


####################################################
class CreateToolTip(object):
    """
    create a tooltip for a given widget
    """
    def __init__(self, widget, text='widget info'):
        self.waittime = 500     #miliseconds
        self.wraplength = 280   #pixels
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

window.title("PhDSeeker app")

window.geometry('1400x675')

group = tk.LabelFrame(window, padx=5, pady=5, text=" SETTINGS ")
group.pack(padx=10, pady=10, side="left", fill="both", expand=False)


RIGHT = tk.LabelFrame(window, padx=0, pady=0, text="", borderwidth=0, highlightthickness=0)
RIGHT.pack(padx=0, pady=5, side="left", fill="both", expand=False)


results = tk.LabelFrame(RIGHT, padx=15, pady=10, text=" RESULTS ")
results.pack(padx=10, pady=5, side="top", fill="both", expand=False)


viewer = tk.LabelFrame(RIGHT, padx=15, pady=10, text=" VIEWER ")
viewer.pack(padx=10, pady=5, side="top", fill="both", expand=True)



#============================================
# OPEN REACTIONS FILE
#======================
def load_reactions():
    '''
    Define action when button is press.
    '''
    btn21.filename = filedialog.askopenfilename(initialdir = current_path + "/db",
                                                title = "Select reactions file...",
                                                filetypes = (("Text files","*.txt"),("all files","*.*")))
    
    path_reactions, filename_reactions = os.path.split(btn21.filename)
    
    txt21.configure(text='{}'.format(filename_reactions.strip()))  # SPLIT and SHOW FILENAME
#============================================


#============================================
# OPEN ENZYMES FILE
#======================
def load_enzymes():
    '''
    Define action when button is press.
    '''
    btn22.filename = filedialog.askopenfilename(initialdir = current_path + "/db",
                                                title = "Select enzymes file...",
                                                filetypes = (("Text files","*.txt"),("all files","*.*")))
    
    path_enzymes, filename_enzymes = os.path.split(btn22.filename)
    
    txt22.configure(text='{}'.format(filename_enzymes.strip()))  # SPLIT and SHOW FILENAME
#============================================


#============================================
# OPEN COMPOUNDS FILE
#======================
def load_compounds():
    '''
    Define action when button is press.
    '''
    btn23.filename = filedialog.askopenfilename(initialdir = current_path + "/db",
                                                title = "Select compounds file...",
                                                filetypes = (("Text files","*.txt"),("all files","*.*")))
    
    path_compounds, filename_compounds = os.path.split(btn23.filename)
    
    txt23.configure(text='{}'.format(filename_compounds.strip()))  # SPLIT and SHOW FILENAME
#============================================


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#=====================================================
def ok_button_abundant():
    '''
    '''
    
    #---------------
    # SELECT ITEM
    #---------------
    items = map(int, create_window_abundant.listbox_abundant.curselection())
    ok_button_abundant.compounds = [create_window_abundant.COMPOUNDS[int(item)][0].strip() for item in items]
    
    ok_button_abundant.compounds = sorted(ok_button_abundant.compounds)
    
    txt31.configure(text='{}'.format(', '.join(ok_button_abundant.compounds)))  # SPLIT and SHOW FILENAME
    
    create_window_abundant.window_abundant.destroy()

ok_button_abundant.compounds = []
#=======================================================



#=======================================================
# BUTTON
#==================
def create_window_abundant():
    
    if (btn21.filename == None):
        messagebox.showerror(title=None,message='No reaction file specified yet.')
    
    else:
    
        create_window_abundant.window_abundant = tk.Toplevel(window)
        
        create_window_abundant.window_abundant.title("Select abundant compounds")
        create_window_abundant.window_abundant.geometry('600x675')
        
        viewer_abundant = tk.LabelFrame(create_window_abundant.window_abundant, padx=15, pady=10, text=" Abundant compounds ")
        viewer_abundant.pack(padx=10, pady=5, side="top", fill="both", expand=False)

        viewer1_abundant = tk.LabelFrame(viewer_abundant, padx=15, pady=10, text="")
        viewer1_abundant.pack(padx=10, pady=5, side="top", fill="both", expand=False)
        
        scrollbar_abundant = tk.Scrollbar(viewer_abundant, orient='vertical')
        create_window_abundant.listbox_abundant = tk.Listbox(viewer_abundant,
                                                             selectmode='extended',
                                                             yscrollcommand=scrollbar_abundant.set,
                                                             font=('Tahoma', 10), height=36)
        
        scrollbar_abundant.config(command=create_window_abundant.listbox_abundant.yview)
        create_window_abundant.listbox_abundant.pack(side='left', fill='both', expand=True)
        scrollbar_abundant.pack(side='left', fill=tk.Y)
        
        
        #--------------------------------------------
        # READ REACTIONS FILE AND GETTING COMPOUNDS
        #--------------------------------------------
        create_window_abundant.COMPOUNDS = []
        
        #············································
        
        if (btn23.filename is not None):
            with open(btn23.filename,'r') as fp:
                lines = fp.read()
                lines = lines.strip().split('\n')
                
            create_window_abundant.NAMES = dict()
            for line in lines:
                code,name = line.split(': ')
                create_window_abundant.NAMES[code] = name
            
        else:
            create_window_abundant.NAMES = dict()
        #············································
        
        with open(btn21.filename,'r') as fp:
            reactions = fp.readlines()
        
        for reaction in reactions:
            
            R,equation = reaction.split(': ')
            
            if ' <=> ' in equation:
                
                left,right = equation.split(' <=> ')
                
                reverse_reaction = R + ': ' + right.strip() + ' --> ' + left.strip() + '\n'
                reactions.append(reverse_reaction)
                
                left = left.replace(' + ', ' ').replace('\n', '')
                right = right.replace(' + ', ' ').replace('\n', '')
                
            else:
                equation = equation.replace(' + ', ' ')
                equation = equation.replace('\n', '')
                left,right = equation.split(' --> ')
                
                Substrates = left.split(' ') # [Coeficiente - Compuesto - Coeficiente - Compuesto]
                Products = right.split(' ')  # [Coeficiente - Compuesto - Coeficiente - Compuesto]
                
                
                for idx in range(0,len(Substrates),2):
                    create_window_abundant.COMPOUNDS.append(Substrates[idx+1])
                
                for idx in range(0,len(Products),2):
                    create_window_abundant.COMPOUNDS.append(Products[idx+1])
        
        
        create_window_abundant.COMPOUNDS = list(set(create_window_abundant.COMPOUNDS))
        
        create_window_abundant.COMPOUNDS = [(compound,create_window_abundant.NAMES.get(compound,None)) for compound in create_window_abundant.COMPOUNDS]
        
        
        if len(create_window_abundant.NAMES) > 0:
            
            create_window_abundant.COMPOUNDS = sorted(create_window_abundant.COMPOUNDS, key=lambda x: x[1])
            
            for compound in create_window_abundant.COMPOUNDS:
                create_window_abundant.listbox_abundant.insert(tk.END, '  {} ({})'.format(compound[1],compound[0]))
            
        else:
            create_window_abundant.COMPOUNDS = sorted(create_window_abundant.COMPOUNDS, key=lambda x: x[0])
            
            for compound in create_window_abundant.COMPOUNDS:
                create_window_abundant.listbox_abundant.insert(tk.END, '  {}'.format(compound[0]))
        
        
        create_window_abundant.listbox_abundant.select_set(0)
        
        #---------------------------------------
        ok_button_abundant.compounds = None
        
        btn_abundant_ok=tk.Button(create_window_abundant.window_abundant, text="Ok", bg="#02B844", fg="black", command=ok_button_abundant, height=1, width=10)
        btn_abundant_ok.pack(side="top", fill=None, expand=False)

#-----------------------    

create_window_abundant.window_abundant = None
create_window_abundant.COMPOUNDS = []
create_window_abundant.NAMES = dict()

#btn_abundant_run=tk.Button(group, text="Select", bg="#02B844", fg="black", command=create_window_abundant, height=1, width=10)
#btn_abundant_run.pack(side="bottom", fill=None, expand=False)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#=====================================================
def ok_button_relate():
    '''
    '''
    
    #---------------
    # SELECT ITEM
    #---------------
    items = map(int, create_window_relate.listbox_relate.curselection())
    ok_button_relate.compounds = [create_window_relate.COMPOUNDS[int(item)][0].strip() for item in items]
    
    ok_button_relate.compounds = sorted(ok_button_relate.compounds)
    
    txt32.configure(text='{}'.format(', '.join(ok_button_relate.compounds)))  # SPLIT and SHOW FILENAME
    
    create_window_relate.window_relate.destroy()

ok_button_relate.compounds = []
#=======================================================



#=======================================================
# BUTTON
#==================
def create_window_relate():
    
    if (btn21.filename == None):
        messagebox.showerror(title=None, message='No reaction file specified yet.')
    
    else:
    
        create_window_relate.window_relate = tk.Toplevel(window)
        
        create_window_relate.window_relate.title("Select compounds to relate")
        create_window_relate.window_relate.geometry('600x675')
        
        viewer_relate = tk.LabelFrame(create_window_relate.window_relate, padx=15, pady=10, text=" Compounds to relate ")
        viewer_relate.pack(padx=10, pady=5, side="top", fill="both", expand=False)

        viewer1_relate = tk.LabelFrame(viewer_relate, padx=15, pady=10, text="")
        viewer1_relate.pack(padx=10, pady=5, side="top", fill="both", expand=False)
        
        scrollbar_relate = tk.Scrollbar(viewer_relate, orient='vertical')
        create_window_relate.listbox_relate = tk.Listbox(viewer_relate,
                                                             selectmode='extended',
                                                             yscrollcommand=scrollbar_relate.set,
                                                             font=('Tahoma', 10), height=36)
        
        scrollbar_relate.config(command=create_window_relate.listbox_relate.yview)
        create_window_relate.listbox_relate.pack(side='left', fill='both', expand=True)
        scrollbar_relate.pack(side='left', fill=tk.Y)
        
        
        #--------------------------------------------
        # READ REACTIONS FILE AND GETTING COMPOUNDS
        #--------------------------------------------
        if (create_window_abundant.COMPOUNDS == None):
            messagebox.showwarning(title=None,message='Abundant compounds should be defined first.')
        
        else:
            
            create_window_relate.COMPOUNDS = create_window_abundant.COMPOUNDS[:]
            
            for compound in create_window_relate.COMPOUNDS:
                
                if compound[1] is None:
                    create_window_relate.listbox_relate.insert(tk.END, '  {}'.format(compound[0]))
                    
                else:
                    create_window_relate.listbox_relate.insert(tk.END, '  {} ({})'.format(compound[1],compound[0]))
            
            
            #if create_window_relate.COMPOUNDS[0][1] is None:
                #create_window_abundant.COMPOUNDS = sorted(create_window_abundant.COMPOUNDS, key=lambda x: x[0])
            #else:
                #create_window_abundant.COMPOUNDS = sorted(create_window_abundant.COMPOUNDS, key=lambda x: x[1])
            
            
            create_window_relate.listbox_relate.select_set(0)
            
            #---------------------------------------
            ok_button_relate.compounds = None
            
            btn_relate_ok=tk.Button(create_window_relate.window_relate, text="Ok", bg="#02B844", fg="black", command=ok_button_relate, height=1, width=10)
            btn_relate_ok.pack(side="top", fill=None, expand=False)

#-----------------------    

create_window_relate.window_relate = None
create_window_relate.COMPOUNDS = []
create_window_relate.NAMES = dict()

#btn_relate_run=tk.Button(group, text="Select", bg="#02B844", fg="black", command=create_window_relate, height=1, width=10)
#btn_relate_run.pack(side="bottom", fill=None, expand=False)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#=====================================================
def ok_button_initial():
    '''
    '''
    
    #---------------
    # SELECT ITEM
    #---------------
    items = map(int, create_window_initial.listbox_initial.curselection())
    ok_button_initial.compounds = [create_window_initial.COMPOUNDS[int(item)][0].strip() for item in items]
    
    txt33.configure(text='{}'.format(', '.join(ok_button_initial.compounds)))  # SPLIT and SHOW FILENAME
    
    create_window_initial.window_initial.destroy()

ok_button_initial.compounds = []
#=======================================================



#=======================================================
# BUTTON
#==================
def create_window_initial():
    
    if (btn21.filename == None):
        messagebox.showerror(title=None,message='No reaction file specified yet.')
    
    else:
    
        create_window_initial.window_initial = tk.Toplevel(window)
        
        create_window_initial.window_initial.title("Select initial compounds")
        create_window_initial.window_initial.geometry('600x675')
        
        viewer_initial = tk.LabelFrame(create_window_initial.window_initial, padx=15, pady=10, text=" Abundant compounds ")
        viewer_initial.pack(padx=10, pady=5, side="top", fill="both", expand=False)

        viewer1_initial = tk.LabelFrame(viewer_initial, padx=15, pady=10, text="")
        viewer1_initial.pack(padx=10, pady=5, side="top", fill="both", expand=False)
        
        scrollbar_initial = tk.Scrollbar(viewer_initial, orient='vertical')
        create_window_initial.listbox_initial = tk.Listbox(viewer_initial,
                                                             selectmode='extended',
                                                             yscrollcommand=scrollbar_initial.set,
                                                             font=('Tahoma', 10), height=36)
        
        scrollbar_initial.config(command=create_window_initial.listbox_initial.yview)
        create_window_initial.listbox_initial.pack(side='left', fill='both', expand=True)
        scrollbar_initial.pack(side='left', fill=tk.Y)
        
        
        #--------------------------------------------
        # READ REACTIONS FILE AND GETTING COMPOUNDS
        #--------------------------------------------
        if (ok_button_relate.compounds == None):
            messagebox.showwarning(title=None,message='Compounds to relate should be defined first.')
        
        else:
            
            create_window_initial.COMPOUNDS = [(compound, create_window_abundant.NAMES.get(compound,None)) for compound in ok_button_relate.compounds]
            
            for compound in create_window_initial.COMPOUNDS:
                
                if compound[1] is None:
                    create_window_initial.listbox_initial.insert(tk.END, '  {}'.format(compound[0]))
                    
                else:
                    create_window_initial.listbox_initial.insert(tk.END, '  {} ({})'.format(compound[1],compound[0]))
            
            
            create_window_initial.listbox_initial.select_set(0)
            
            #---------------------------------------
            ok_button_initial.compounds = None
            
            btn_initial_ok=tk.Button(create_window_initial.window_initial, text="Ok", bg="#02B844", fg="black", command=ok_button_initial, height=1, width=10)
            btn_initial_ok.pack(side="top", fill=None, expand=False)

#-----------------------    

create_window_initial.window_initial = None
create_window_initial.COMPOUNDS = []
create_window_initial.NAMES = dict()

#btn_initial_run=tk.Button(group, text="Select", bg="#02B844", fg="black", command=create_window_initial, height=1, width=10)
#btn_initial_run.pack(side="bottom", fill=None, expand=False)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@






#############################################################################
G = tk.LabelFrame(group, padx=15, pady=10, text="", borderwidth=0, highlightthickness=0)
G.pack(padx=10, pady=5)

#----------------------------------

G1 = tk.LabelFrame(G, padx=15, pady=10, text=" ALGORITHM ", width=70)
G1.pack(padx=10, pady=5, side="left")


#----------
# NCORES
#----------
ncores = tk.Scale(G1, from_=1, to=multiprocessing.cpu_count(), digits=1, resolution=1, label="Number of cores", orient=tk.HORIZONTAL)
ncores.pack(side="top", fill="both", expand=True, pady=5)
ncores_help = CreateToolTip(ncores, 'Number of cores used for parallelizing the process.')


#-----------------------------
# NUMBER OF ANTS PER ANTHILL
#-----------------------------
nants = tk.Scale(G1, from_=5, to=100, digits=1, resolution=1, label="Number of ants", orient=tk.HORIZONTAL)
nants.pack(side="top", fill="both", expand=True, pady=5)
nants_help = CreateToolTip(nants, 'Number of ants used in the design process. Every ant searhing for an independet solution on each iteration.\nRecomended value: 10.')


#--------------------
# EVAPORATION RATE
#--------------------
rho = tk.Scale(G1, from_=0.01, to=2.0, digits=4, resolution=0.01, label="Evaporation rate", orient=tk.HORIZONTAL)
rho.pack(side="top", fill="both", expand=True, pady=5)
rho_help = CreateToolTip(rho, 'Amount of information from solutions preserved between iterations.\nRecomended value: 0.01.')


#------------------------------
# MAXIMUM NUMBER OF ITERATIONS
#------------------------------
maxIterations = tk.Scale(G1, from_=10, to=1000, digits=2, resolution=10, label="Max. iterations", orient=tk.HORIZONTAL)
maxIterations.pack(side="top", fill="both", expand=True, pady=5)
maxIterations_help = CreateToolTip(maxIterations, 'Maximum number of iterations performed by the algorithm until finish.')


#------------------------------
# ITERATIONS WITHOUT CHANGES
#------------------------------
IterationsWithoutChanges = tk.Scale(G1, from_=10, to=1000, digits=2, resolution=10, label="Iterations without changes", orient=tk.HORIZONTAL)
IterationsWithoutChanges.pack(side="top", fill="both", expand=True, pady=5)
IterationsWithoutChanges_help = CreateToolTip(IterationsWithoutChanges, 'Maximum number of iterations that the algorithm runs without changes in the best solution.\nRecomended value: 100.')


#------------------------------
# ITERATIONS WITH ALIGNED ANTS
#------------------------------
IterationsWithAlignedAnts = tk.Scale(G1, from_=10, to=100, digits=2, resolution=10, label="Iterations with aligned ants", orient=tk.HORIZONTAL)
IterationsWithAlignedAnts.pack(side="top", fill="both", expand=True, pady=5)
IterationsWithAlignedAnts_help = CreateToolTip(IterationsWithAlignedAnts, 'Minimum number of consecutive iterations for which all ants find the same solution.\nRecomended value: 10.')


#------------------------------
# STRICT INITIALIZATION
#------------------------------
chk11 = ttk.Checkbutton(G1, text="Strict initialization")
chk11.pack(side="top", fill="both", expand=True, pady=6)
chk11_help = CreateToolTip(chk11, 'Indicates if solutions should be found using only available compounds or additional compounds can be automatically added.\nRecomended value: Not checked.')


#------------------------------
# ALLOW EXTERNAL COMPOUNDS
#------------------------------
chk12 = ttk.Checkbutton(G1, text="Allow external compounds")
chk12.pack(side="top", fill="both", expand=True, pady=6)
chk12_help = CreateToolTip(chk12, 'Indicates if compounds that are not synthesized by any reaction should be seen as available ones.\nRecomended value: Checked.')


#------------------------------
# VERBOSE
#------------------------------
chk13 = ttk.Checkbutton(G1, text="Verbose")
chk13.pack(side="top", fill="both", expand=True, pady=6)
chk13_help = CreateToolTip(chk13, 'Indicates if results by iteration should be shown in the console.\nRecomended value: Checked.')

#----------------------------------

GA = tk.LabelFrame(G, padx=15, pady=10, text="", borderwidth=0, highlightthickness=0)
GA.pack(padx=0, pady=0, side="left")

G2 = tk.LabelFrame(GA, padx=15, pady=10, text=" DATASET ")
G2.pack(padx=10, pady=5, side="top")

#<<<<<<<<<<<<<<<<<

G21 = tk.LabelFrame(G2, padx=15, pady=10, text=" Reactions ")
G21.pack(padx=10, pady=5, side="top", fill="both")

btn21 = tk.Button(G21, text="Load file", bg="#02B844", fg="black", command=load_reactions, height=1, width=10)
btn21.pack(side="left", fill="both", expand=True)
btn21.filename = None

txt21 = tk.Label(G21, text='', height=1, width=40, bg='#bcbcbc', fg="black")
txt21.pack(side="right", fill="both", expand=True)
txt21_help = CreateToolTip(btn21, "Select a file with the reaction's dataset")

#<<<<<<<<<<<<<<<<<

G22 = tk.LabelFrame(G2, padx=15, pady=10, text=" Enzymes ")
G22.pack(padx=10, pady=5, side="top", fill="both")

btn22 = tk.Button(G22, text="Load file", bg="#02B844", fg="black", command=load_enzymes, height=1, width=10)
btn22.pack(side="left", fill="both", expand=True)
btn22.filename = None

txt22 = tk.Label(G22, text='', height=1, width=40, bg='#bcbcbc', fg="black")
txt22.pack(side="right", fill="both", expand=True)
txt22_help = CreateToolTip(btn22, "Select a file with the enzymes catalyzing each reaction")

#<<<<<<<<<<<<<<<<<

G23 = tk.LabelFrame(G2, padx=15, pady=10, text=" Compound Names ")
G23.pack(padx=10, pady=5, side="top", fill="both")

btn23 = tk.Button(G23, text="Load file", bg="#02B844", fg="black", command=load_compounds, height=1, width=10)
btn23.pack(side="left", fill="both", expand=True)
btn23.filename = None

txt23 = tk.Label(G23, text='', height=1, width=40, bg='#bcbcbc', fg="black")
txt23.pack(side="right", fill="both", expand=True)
txt23_help = CreateToolTip(btn23, "Select a file with the name for each compound identifier")

#----------------------------------

G3 = tk.LabelFrame(GA, padx=15, pady=10, text=" COMPOUNDS ")
G3.pack(padx=10, pady=5, side="top")


G31 = tk.LabelFrame(G3, padx=15, pady=10, text=" Abundant ")
G31.pack(padx=10, pady=5, side="top")
btn31 = tk.Button(G31, text="Select", bg="#02B844", fg="black", command=create_window_abundant, height=1, width=10)
btn31.pack(side="left", fill="both", expand=True)
btn31.filename = None
txt31 = tk.Label(G31, text='', height=1, width=40, bg='#bcbcbc', fg="black")
txt31.pack(side="right", fill="both", expand=True)
txt31_help = CreateToolTip(btn31, "Select a file with the name for each compound identifier")


G32 = tk.LabelFrame(G3, padx=15, pady=10, text=" Relate ")
G32.pack(padx=10, pady=5, side="top")
btn32 = tk.Button(G32, text="Select", bg="#02B844", fg="black", command=create_window_relate, height=1, width=10)
btn32.pack(side="left", fill="both", expand=True)
btn32.filename = None
txt32 = tk.Label(G32, text='', height=1, width=40, bg='#bcbcbc', fg="black")
txt32.pack(side="right", fill="both", expand=True)
txt32_help = CreateToolTip(btn32, "Select a file with the name for each compound identifier")


G33 = tk.LabelFrame(G3, padx=15, pady=10, text=" Initial ")
G33.pack(padx=10, pady=5, side="top")
btn33 = tk.Button(G33, text="Select", bg="#02B844", fg="black", command=create_window_initial, height=1, width=10)
btn33.pack(side="left", fill="both", expand=True)
btn33.filename = None
txt33 = tk.Label(G33, text='', height=1, width=40, bg='#bcbcbc', fg="black")
txt33.pack(side="right", fill="both", expand=True)
txt33_help = CreateToolTip(btn33, "Select a file with the name for each compound identifier")

#----------------------------------

#############################################################################




#==================
# RUN PhDSeeker
#==================
def run_button():
    '''
    Define action when button is press.
    '''
    
    txt_viewer.delete('1.0', tk.END)  # clear screen
    
    if not ((btn21.filename != None) and  ok_button_abundant.compounds and ok_button_relate.compounds and ok_button_initial.compounds):
        messagebox.showwarning(title=None,message='Some data was not definet.\nCheck configuration and run again.')
    
    else:
        
        btn_run.configure(text="Running...")
        btn_run.update()
        
        #-------------------------------------------------------------------------
        
        #filename = os.path.join(load_settings.path_settings, load_settings.filename_settings)
        
        temporal_settings = dict()
        temporal_settings['NCORES'] = ncores.get()
        temporal_settings['Nants'] = nants.get()
        temporal_settings['rho'] = rho.get()
        temporal_settings['maxIterations'] = maxIterations.get()
        temporal_settings['IterationsWithoutChanges'] = IterationsWithoutChanges.get()
        temporal_settings['IterationsWithAlignedAnts'] = IterationsWithAlignedAnts.get()
        
        temporal_settings['StrictInitialization'] = chk11.instate(['selected'])
        temporal_settings['AllowExternalCompounds'] = chk12.instate(['selected'])
        temporal_settings['Verbose'] = chk13.instate(['selected'])
        
        temporal_settings['REACTIONS'] = btn21.filename  # with full path
        temporal_settings['ENZYMES'] = btn22.filename    # with full path
        temporal_settings['COMPNAMES'] = btn23.filename  # with full path
        
        
        #*******************************
        temporal_settings['COMPOUNDS'] = dict()  # {'abundant':[], 'relate': []}
        
        #-------------------
        # ABUNDANT COMPOUNDS
        #-------------------
        temporal_settings['COMPOUNDS']['abundant'] = ok_button_abundant.compounds[:]
        
        
        temporal_settings['COMPOUNDS']['relate'] = []  # relate --> [{'compound': C00103, 'initial': True}, {'compound':C00631, 'initial': False'}]
        #-------------------
        # RELATE COMPOUNDS
        #-------------------
        relate = ok_button_relate.compounds[:]
        
        #-------------------
        # INITIALS
        #-------------------
        initials = ok_button_initial.compounds[:]
        
        for compound in relate:
            
            if compound in initials:
                
                temporal_settings['COMPOUNDS']['relate'].append({'compound': compound, 'initial': True})
                
            else:
                
                temporal_settings['COMPOUNDS']['relate'].append({'compound': compound, 'initial': False})
        
        #*******************************
        
        
        #-----------------
        
        with open(os.path.join(current_path, 'config/SETTINGS_temporal.yaml'), 'w') as fp:
            yaml.dump(temporal_settings, fp, default_flow_style=False)
            
        #-------------------------------------------------------------------------
        
        output = phdseeker.PhDSeeker(os.path.join(current_path, 'config/SETTINGS_temporal.yaml'))
        
        filename_summary,filename_history,filename_html = output
        
        
        # SELECT FOLDER
        if not os.path.exists('out'):
            os.makedirs('out')
        
        folder_selected = filedialog.askdirectory(initialdir = current_path + "/out",
                                                  title = "Select destination folder...")
        
        
        run_button.filename_summary = os.path.join(folder_selected,filename_summary)
        run_button.filename_history = os.path.join(folder_selected,filename_history)
        run_button.filename_html = os.path.join(folder_selected,filename_html)
        
        
        #------------------------
        with open(os.path.join(folder_selected,filename_history), 'r') as fp:
            data = fp.read()
            N = data.count('PATHWAY')
        #------------------------
        
        messagebox.showinfo('Execution finished!', '{} solutions found'.format(N))
        
        #------------------------
        
        # DELETE TEMPORAL FILE
        os.remove(os.path.join(current_path, 'config/SETTINGS_temporal.yaml'))
        
        btn_run.configure(text="Run", bg="#02B844")
        btn_run.update()


run_button.filename_summary = None
run_button.filename_history = None
run_button.filename_html = None


btn_run=tk.Button(group, text="Run", bg="#02B844", fg="black", command=run_button, height=1, width=10)
btn_run.pack(side="bottom", fill=None, expand=False)

txt_run_help = CreateToolTip(btn_run, 'Start searching')



#======================
# SHOW SUMMARY
#======================
txt_viewer = tk.scrolledtext.ScrolledText(viewer, width=30, height=50)
txt_viewer.pack(side="bottom", fill="both", expand=True)


#=================================================
def open_file_summary():
    
    filename = run_button.filename_summary
    
    if filename is not None:
        
        txt_viewer.delete('1.0', tk.END)  # clear screen
        
        with open(filename, 'r') as fp:
            content = fp.read()
            txt_viewer.insert(tk.INSERT, content)
#=================================================

#=================================================
def open_file_history():
    
    filename = run_button.filename_history
    
    if filename is not None:
        
        txt_viewer.delete('1.0', tk.END)  # clear screen
        
        with open(filename, 'r') as fp:
            content = fp.read()
            txt_viewer.insert(tk.INSERT, content)
#=================================================


#=================================================
def open_file_html():
    
    filename = run_button.filename_html
    
    if filename is not None:
        webbrowser.open_new(filename)
#=================================================



lf_link1 = tk.LabelFrame(results, padx=10, pady=2, text="", borderwidth=0, highlightthickness=0)
lf_link1.pack(padx=10, pady=5, side="left", fill="both", expand=False)
btn_link1=tk.Button(lf_link1, text="Summary", bg="#02B844", fg="black", command=open_file_summary, height=1, width=10)
btn_link1.pack(side=None, fill=None, expand=False)
txt_link1_help = CreateToolTip(btn_link1, 'Summary of the performed search')


lf_link2 = tk.LabelFrame(results, padx=10, pady=2, text="", borderwidth=0, highlightthickness=0)
lf_link2.pack(padx=10, pady=5, side="left", fill="both", expand=False)
btn_link2=tk.Button(lf_link2, text="History", bg="#02B844", fg="black", command=open_file_history, height=1, width=10)
btn_link2.pack(side=None, fill=None, expand=False)
txt_link2_help = CreateToolTip(btn_link2, 'Pathways found during the search')


lf_link3 = tk.LabelFrame(results, padx=10, pady=2, text="", borderwidth=0, highlightthickness=0)
lf_link3.pack(padx=10, pady=5, side="left", fill="both", expand=False)
btn_link3=tk.Button(lf_link3, text="Pathways", bg="#02B844", fg="black", command=open_file_html, height=1, width=10)
btn_link3.pack(side=None, fill=None, expand=False)
txt_link3_help = CreateToolTip(btn_link3, 'Interactive visualization of pathways found during the search')


#----------------------


window.mainloop() 
