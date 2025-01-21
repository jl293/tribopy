import tkinter as tk
from tkinter import ttk

class TabbedApplication(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("Tribology Chapters")
        self.geometry("600x400")
        
        # Center the window on the screen
        self.eval('tk::PlaceWindow . center')

        # Create a tab control
        self.tab_control = ttk.Notebook(self)

        # Define chapter names
        chapters = [
            "Fundamental Concepts",
            "Solid Materials",
            "Surface Roughness",
            "Non-Conformal Contact",
            "Liquid Properties",
            "Lubricant Composition",
            "Lubricant Characterization",
            "Lubricant of Conformal Contacts",
            "Lubrication of Non-Conformal Contacts",
            "Dry and Mixed Friction",
            "Wear",
            "Measuring Friction and Wear",
            "Nano- and Biotribology"
        ]

        # Create tabs for each chapter
        for chapter in chapters:
            tab = ttk.Frame(self.tab_control)
            self.tab_control.add(tab, text=chapter)

            # Add a header label to each tab
            label = ttk.Label(tab, text=chapter, font=("Arial", 16))
            label.pack(pady=20)

        self.tab_control.pack(expand=1, fill="both")

if __name__ == "__main__":
    app = TabbedApplication()
    app.mainloop()