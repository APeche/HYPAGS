# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 09:35:32 2022

@author: Peche.A
"""
import numpy as np
import matplotlib.pyplot as plt
from tkinter import *
from tkinter import messagebox
from PIL import ImageTk, Image
from pylab import text

def main():
            
    # GUI Body
    myGui = Tk()
    myGui.iconbitmap('HYPAGS_dependencies/HYPAGS.ico')
        
    # GUI title
    myGui.title("HYPAGS GUI (Version 1.2.2)")
    myGui.geometry("")

    # model out of bounds message  
    def messageBox01():
        messagebox.showwarning('WARNING:', 'Your input is out of model bounds. \n Model bounds:\n K [m/s]    - Min: 2.87*10^(-7), Max: 2.60*10^(-2)\n d_10 [m] - Min: 5.35*10^(-5), Max: 8.3*10^(-4) \n d_20 [m] - Min: 6.25*10^(-5), Max: 1.2*10^(-3)')
     
    
    # Calculation Functions
    # Input is K
    def Calculator(p, ip):
        # temproary workaround - assign global variables to alter use the in PlotSieveCurve()
        global Kt, d1, d2, d5, d6, net
        # constants:
        a1, a2, a3 = 1.004, 1.51e-4, 5.788e-3
        mu, g, rho = 1.1375e-3, 9.81, 999.7
        # set parameterization
        if p == 1:
            Pi, c1, c2 = .0009, 1.2, 1.13 # KC
        elif p == 2:
            Pi, c1, c2 = .0023, 1.5, 1.15 # BS
        elif p == 3:
            Pi, c1, c2 = .0019, 1.38, 1.17 # H
        if ip == 1:    
            Kt = InPar.get()
            # check for nonphysical input
            if Kt > 2.6e-2 or Kt < 2.87e-7:
                messageBox01()
            # calculate d1
            d1 = (Kt/Pi * (mu/(rho*g)))**.5
            # calculate d2
            d2 = c1 * d1
        elif ip ==2:
            d1 = InPar.get()
            # check for nonphysical input
            if d1 > 8.3e-4 or d1 < 5.35e-5:
                messageBox01()
            # calculate K 
            Kt = (Pi * rho * g * d1**2)/mu
            # calculate d2
            d2 = c1 * d1   
        elif ip == 3:
            d2 = InPar.get()
            # check for nonphysical input
            if d2 > 1.2e-3 or d2 < 6.25e-5:
                messageBox01()
            # calculate d1
            d1 = d2/c1
            # calculate K 
            Kt = (Pi * rho * g * d1**2)/mu
        # iteratively solve for d5,ne
        ne0 = a1 * Kt**a2
        d500 = a3 * d2
        net, d5 = iterator(ne0, d500, Kt)
        # d6
        d6 = c2 * d5
        # for output - fill output container
        K.set(['K','=', round(Kt, 9),'m/s'])
        d10.set(['d10','=', round(d1, 9),'m'])
        d20.set(['d20','=', round(d2, 9),'m'])
        d50.set(['d50','=', round(d5, 9),'m'])
        d60.set(['d60','=', round(d6, 9),'m'])
        ne.set(['ne','=', round(net, 3)])
        # Output to text file
        text_file = open("HYPAGS_result.csv", 'w')
        label = ('K[m/s];ne[-];d10[m];d20[m];d50[m];d60[m]'+'\n' + str(Kt) + ';' + str(net) + ';' + str(d1) + ';' + str(d2) + ';' + str(d5) + ';' + str(d6))
        text_file.write(label)
        text_file.close()
    
    
        
    def iterator(ne_i, d50_i, K):
      error = 10
      mu, g, rho = 1.1375e-3, 9.81, 999.7
      c = 180 * mu/(rho * g)
      n = np.zeros(100)
      d = np.zeros(100)
      d[0] = d50_i
      n[0] = ne_i
      j = 0
      while error > 1e-10:
          n[j+1] = (K/(d[j]**2) * c * (1 - n[j])**2)**(1/3)
          d[j+1] = (c * K * ( (1-n[j+1])**2 / (n[j+1]**3) ))**(1/2)
          error = np.abs(n[j+1]-n[j]+d[j+1]-d[j])
          j = j + 1
      n_r = n[j-1]
      d5 = d[j-1]
      return n_r, d5  
    
     
    def closeWin():
        myGui.destroy()  # Close Window Function
    
    
    def clearFunc():
        InPar.set("0")
        final.set("")
        K.set("")
        ne.set("")
        d10.set("")
        d20.set("")
        d50.set("")
        d60.set("")
        Kt = 0 
        d1 = 0 
        d2 = 0 
        d5 = 0 
        d6 = 0 
        net = 0
    

    # input variables
    final= StringVar()
    d10 = StringVar()
    K = StringVar()
    ne = StringVar()
    d20 = StringVar()
    d50 = StringVar()
    d60 = StringVar()
    
    def PlotSieveCurve(Kt, d1, d2, d5, d6, net):
        x = [10,20,50,60]
        fig = plt.figure("Sieve curve")
        ax = fig.add_subplot(1, 1, 1)
        plt.plot([d1,d2,d5,d6],x, 'k', linewidth=2.0, marker='o')
        ax.set_xscale('log')
        plt.ylim((0,100))
        plt.xlim((5e-7,4.3e-3))
        ax.axvspan(0, .002e-3, facecolor='gainsboro', alpha=0.5)
        ax.axvspan(.002e-3, 6e-5, facecolor='linen', alpha=0.5)
        ax.axvspan(6e-5, 2e-3, facecolor='peachpuff', alpha=0.5)
        ax.axvspan(2e-3, 9e-3, facecolor='lightsalmon', alpha=0.5)
        plt.text(8e-7,102,'Clay')
        plt.text(1.2e-5,102,'Silt')
        plt.text(2e-4,102,'Sand')
        plt.text(2.4e-3,102,'Gravel')
        posl, post = .02, .82
        label = ('K = ' + str(round(Kt,8)) + ' m/s' + '\n' + '$n_e$ = ' + str(round(net,3)) + '\n' '$d_{10}$ = ' + str(round(d1,6)) + ' m' + '\n' '$d_{20}$ = ' + str(round(d2,6)) + ' m' + '\n' '$d_{50}$ = ' + str(round(d5,6)) + ' m' + '\n' '$d_{60}$ = ' + str(round(d6,6)) + ' m') 
        text(posl, post, label, ha='left', va='center', transform=ax.transAxes, backgroundcolor='white') 
        plt.grid(True, which="both", ls="-")
        plt.ylabel('Mass/sieving fraction of sample [%]')
        plt.xlabel('Grain size [m]')
        plt.show()
    
    
    def about():
        global my_img
        top = Toplevel()
        citeas = 'The HYPAGS model is described in: Peche A., & Houben, G. (2022).\n Estimating Characteristic Grain Sizes\n and Effective Porosity from Hydraulic Conductivity Data.\n Groundwater, Wiley, DOI: 10.1111/gwat.13266.'
        info = 'HYPAGS GUI (Version 1.2.2) - Aaron Peche, Georg Houben, \n Federal Institute for Geosciences and Natural Resources, 2022.'
        info2 = 'HYPAGS (HYdraulic Parameters And Grain Sizes) model enables to retrieve \n virtual sieve curves (by means of representative grain sizes \n d10, d20, d50, d60 [m]) and hydraulically relevant parameters effective \n porosity ne [-] and hydraulic conductivity K [m/s] by single input values of K, d10, d20. \n PLEASE NOTE THAT THE MODEL IS VALID FOR COEFFICIENTS OF UNIFORMITY Cu <~ 3.'
        info3 = 'Results of the HYPAGS model are saved in the local directory in the file HYPAGS_result.csv. \n for further use with MS EXCEL etc.'
        newline = '\n \n'
        # label = info + newline + info2 + newline + info3 + newline + citeas
        Label(top, text=info).grid(column = 1 , row = 1)
        Label(top, text=info2).grid(column = 1 , row = 2)
        Label(top, text=info3).grid(column = 1 , row = 3)
        Label(top, text=citeas).grid(column = 1 , row = 4)
        my_img = ImageTk.PhotoImage(Image.open(r'HYPAGS_dependencies/Compi.png'))
        Label(top, image=my_img).grid(column = 1 , row = 6)
        top.title("About HYPAGS")
        
    def licence():
        paragraph='Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use,copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The ciation of the scientific work, above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.'
        messagebox.showinfo('Licence', paragraph)    
        
    # Radiobuttons for parameterization
    p = IntVar() 
    p.set(1) # set parameterization on Kozeny-Carman per default
    Label(myGui, text="Parameterization:").grid(column = 1 , row = 1)
    Radiobutton(myGui, text = "Kozeny-Carman-based", variable = p, value = 1).grid(column = 1, row = 2, sticky='w')
    Radiobutton(myGui, text = "Beyer-Schweiger-based", variable = p, value = 2).grid(column = 1, row = 3, sticky='w')
    Radiobutton(myGui, text = "Hazen-based", variable = p, value = 3).grid(column = 1, row = 4, sticky='w')
       
    # Radiobuttons for input parameter
    ip = IntVar() 
    ip.set(1) # set parameterization on Kozeny-Carman per default
    Label(myGui, text="Calculation based on:").grid(column = 1 , row = 6)
    Radiobutton(myGui, text = "K [m/s]", variable = ip, value = 1, command=lambda: clicked(ip.get())).grid(column = 1, row = 7, sticky='w')
    Radiobutton(myGui, text = "d10 [m]", variable = ip, value = 2, command=lambda: clicked(ip.get())).grid(column = 1, row = 8, sticky='w')
    Radiobutton(myGui, text = "d20 [m]", variable = ip, value = 3, command=lambda: clicked(ip.get())).grid(column = 1, row = 9, sticky='w')

    # Input entry widget    
    myLabel = Label(text='Input K [m/s]:', fg="black", justify='center').grid(column = 3, row = 1) 
    def clicked(ip):
        if ip == 2:
            myLabel = Label(text='Input d10 [m]:', fg="black", justify='center').grid(column = 3, row = 1) 
            clearFunc()
        elif ip == 3:
            myLabel = Label(text='Input d20 [m]:', fg="black", justify='center').grid(column = 3, row = 1) 
            clearFunc()
        elif ip == 1:
            myLabel = Label(text='Input K [m/s]:', fg="black", justify='center').grid(column = 3, row = 1)     
            clearFunc()
    InPar = DoubleVar()
    InParEntry = Entry(myGui, textvariable=InPar, justify='center').grid(column = 3, row = 2)
    # calculation button
    cp = IntVar() 
    cp.set(0) # set parameterization on Kozeny-Carman per default
    CalcButton = Button(text="Calculation", command=lambda: Calculator(p.get(), ip.get())).grid(column = 3, row = 3)
    # Output
    OutputLabel = Label(text ="Result :", textvariable=K, fg='blue', justify='center').grid(column = 3, row = 4)
    OutputLabel = Label(text ="Result :", textvariable=d10, fg='blue', justify='center').grid(column = 3, row = 5)
    OutputLabel = Label(text ="Result :", textvariable=d20, fg='blue', justify='center').grid(column = 3, row = 6)
    OutputLabel = Label(text ="Result :", textvariable=d50, fg='blue', justify='center').grid(column = 3, row = 7)
    OutputLabel = Label(text ="Result :", textvariable=d60, fg='blue', justify='center').grid(column = 3, row = 8)
    OutputLabel = Label(text ="Result :", textvariable=ne, fg='blue', justify='center').grid(column = 3, row = 9)
    
    
    # Buttons for 'Reset values', 'Program Info', and 'Quit program'
    reset = Button(text="Reset Input", command=clearFunc)
    reset.grid(column = 1, row = 16)
    
    abou = Button(text= "About", command=about)
    abou.grid(column = 1, row = 17)
    
    SievePlot = Button(text= "Plot sieve curve", command=lambda: PlotSieveCurve(Kt, d1, d2, d5, d6, net))
    SievePlot.grid(column = 3, row = 16)
    
    lic = Button(text= "Licence", command=licence)
    lic.grid(column = 3, row = 17)
    
    quit = Button(text="Quit", command=closeWin)
    quit.grid(column = 2, row = 18)
    
    myGui.mainloop()

main()

