import sys
sys.path.append('/Users/ddeal/JWST-Treasure-Chest-2023/JWST-FGS-Spelunker-Repos/JWST-FGS-Spelunker/')
import spelunker.spelunker as spelunker

import tkinter as tk
from PIL import Image, ImageTk

window = tk.Tk()
window.iconbitmap('SPKLOGO.ico')
window.geometry("1200x800")

greeting = tk.Label(text="Hello, Tkinter")
greeting.pack()

image = Image.open('SPK.png')
test = ImageTk.PhotoImage(image)
label1 = tk.Label(image=test)
label1.image = test
label1.place(x=200, y=000)
window.mainloop()