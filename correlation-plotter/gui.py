# Use Tkinter for python 2, tkinter for python 3
import tkinter as tk
from tkinter import ttk

class Main(tk.Frame):
    def __init__(self, parent, controller, *args, **kwargs):
        "main page"
        # super().__init__(self, parent)
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text='Main Page')
        self.grid(column=0, row=0, sticky=(tk.N, tk.W, tk.E, tk.S))
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        # configf = tk.filedialog.FileDialog

        configf = tk.StringVar()
        config_entry = ttk.Entry(self, width=7, textvariable=configf)
        config_entry.grid(column=2, row=1, sticky=(tk.W, tk.E))

        ttk.Label(self, textvariable=configf).grid(column=3, row=1, sticky=tk.W)

        for child in self.winfo_children():
            child.grid_configure(padx=5, pady=5)



class PCA(tk.Frame):
    pass



class App(tk.Tk):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.title('CorrelationPlotter GUI')
        self.geometry('500x300')

        container = tk.Frame(self)

        main_frame = Main(container, self)
