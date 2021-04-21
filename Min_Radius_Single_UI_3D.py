# -*- coding: utf-8 -*-
'''
akayurin@gmail.com
Finding min curvature of 3D line  
'''
import pandas as pd
import numpy as np
import matplotlib.ticker
from matplotlib import pyplot as plt
import math
import re
from tkinter import *
from tkinter import filedialog as fd
from tkinter import messagebox as mb
import os

#%% Open file                     
def open_line_file():
    global line_filename
    global ol

    line_filename = fd.askopenfilename(title='Select line file')

    # Read file
    infile = open(line_filename,'r') # Format: Easting Northing Depth - any separator
    e_l = []
    n_l = []
    z_l = []
    k1_l = []
    k2_l = []
    k3_l = []
    k4_l = []
    
    # Segment lengths
    min_seg_l = 100000000
    max_seg_l = 0
    tot_seg_l = 0
    i = 0
    for line in infile:
        if line[0].isnumeric() or line[0] == '-':
            strn = re.sub(r'[;,:\s+]', ' ', line)
            strn = re.split(r'\s+', strn)
            east = float(strn[0])
            north = float(strn[1])
            z = float(strn[2])
            if i > 0:
                seg_l = ((east - e_l[i-1]) ** 2 + (north - n_l[i-1]) ** 2 + (z - z_l[i-1]) ** 2) ** 0.5
                tot_seg_l = tot_seg_l + seg_l
                if seg_l > max_seg_l:
                    max_seg_l = seg_l
                if seg_l < min_seg_l:
                    min_seg_l = seg_l
            i += 1
    # Coordinates of points of original line
            e_l.append(east)
            n_l.append(north)
            z_l.append(z)
    # Coeff's for circles equation        
            k1_l.append(2 * east) 
            k2_l.append(2 * north)
            k3_l.append(2 * z)
            k4_l.append(-(east ** 2 + north ** 2 + z ** 2))        
    infile.close()
    ol = pd.DataFrame({'E':e_l, 'N':n_l, 'Z':z_l, 'k1':k1_l, 'k2':k2_l, 'k3':k3_l, 'k4':k4_l}) # original line
        
    # Set sliders
    app.lb1var.set('No of points: ' + str((len(ol)))) # No of points in file
    app.lb3var.set(f'Mean segment 3D length: {(tot_seg_l / i):.2f}')
    app.lb4var.set(f'Min segment 3D length: {min_seg_l:.2f}')
    app.lb5var.set(f'Max segment 3D length: {max_seg_l:.2f}')    
    app.min_low_value.set(1)
    app.max_low_value.set(len(ol) - 3)
    app.min_high_value.set(3)
    app.max_high_value.set(len(ol))   
    app.dec_value.set(len(ol) / 2 - 1)
    app.tick_dec_interval.set(app.dec_value.get() / 5)
    app.scales_set()

#%% Main   
def minrad():
    dec = app.dec_scale.get()               
    # Plane coeff's
    ca = []
    cb = []
    cc = []
    # Circles
    ce = []
    cn = []
    cz = []
    r = []
    decimation = []
    p_no = []

    i = app.low_scale.get()  
    while i < app.high_scale.get(): # Triplets
        try:
            # Solve plane equation for triplets (find plane coeff's a*e + b*n + c = z)
            A = np.matrix([[ol['E'][i], ol['N'][i], 1],
                          [ol['E'][i + dec], ol['N'][i + dec], 1],
                          [ol['E'][i + 2 * dec], ol['N'][i + 2 * dec], 1]])
            B = np.matrix([[ol['Z'][i]], [ol['Z'][i + dec]], [ol['Z'][i + 2 * dec]]])
            X = np.linalg.solve(A, B)
            a = X[0, 0]
            b = X[1, 0]
            c = X[2, 0]
            ca.append(a)
            cb.append(b)
            cc.append(c)
            # Solve circle equation (find ce, cn, cz, r)
            C = np.matrix([[ol['k1'][i], ol['k2'][i], ol['k3'][i], 1],
                           [ol['k1'][i + dec], ol['k2'][i + dec], ol['k3'][i + dec], 1],
                           [ol['k1'][i + 2 * dec], ol['k2'][i + 2 * dec], ol['k3'][i + 2 * dec], 1],
                           [a, b, -1, 0]])
            D = np.matrix([[ol['k4'][i]], [ol['k4'][i + dec]], [ol['k4'][i + 2 * dec]], [c]])
            Y = np.linalg.solve(C, D)
            # Triplets circles
            ce.append(-Y[0, 0])
            cn.append(-Y[1, 0])
            cz.append(-Y[2, 0])
            r.append((Y[0, 0] ** 2 + Y[1, 0] ** 2 + Y[2, 0] ** 2 - Y[3, 0]) ** 0.5)
            decimation.append(dec)
            p_no.append(i + dec)

            i += 1
        except Exception: # if straight line
            i += 1
            continue
    
    # 1 - All circles for single decimation:
    all_circles_for_single_dec = pd.DataFrame({'E':ce, 'N':cn, 'Z':cz, 'R':r, 'a':ca, 'b':cb, 'c':cc, 'Dec':decimation, 'P_No':p_no})
    # 2- Min circle for single decimation:
    min_circle_in_single_dec = (all_circles_for_single_dec[all_circles_for_single_dec['R'] == min(all_circles_for_single_dec['R'])]) # Min circle      
    # 3 - Parameters of min circle 
    point_no = min_circle_in_single_dec['P_No'].item()
    ce_min = min_circle_in_single_dec['E'].item()
    cn_min = min_circle_in_single_dec['N'].item()
    cz_min = min_circle_in_single_dec['Z'].item()
    dec = min_circle_in_single_dec['Dec'].item()
    r_min = min_circle_in_single_dec['R'].item()
    a_min = min_circle_in_single_dec['a'].item()
    b_min = min_circle_in_single_dec['b'].item()
    c_min = min_circle_in_single_dec['c'].item()

    fname = os.path.basename(line_filename)
    no_extension = fname.split('.')[0]
    foldername = os.path.dirname(line_filename)
    circle_filename = foldername + '/' + no_extension + '_Dec_' + str(dec) + '_CIRCLE.dig'

    outfile = open(circle_filename,'w') # circle ENZ file
    outfile.write('#unit=m\n') # .dig file header
    # 3D circle min
    e_3d = []
    n_3d = []
    z_3d = []
    for i in range(0, 361):
        # r_min circle on EN plane
        e1 = ce_min + r_min * math.cos(math.radians(i)) 
        n1 = cn_min + r_min * math.sin(math.radians(i))
        z1 = (e1 - ce_min) * a_min + (n1 - cn_min) * b_min + cz_min
        # projection of r_min to min circle plane
        r2 = r_min * math.cos(math.atan((z1 - cz_min) / r_min)) 
        # r_min circle on min circle plane
        e2 = ce_min + r2 * math.cos(math.radians(i))
        n2 = cn_min + r2 * math.sin(math.radians(i))
        z2 = (e2 - ce_min) * a_min + (n2 - cn_min) * b_min + cz_min
        e_3d.append(e2)
        n_3d.append(n2)    
        z_3d.append(z2)
        # Write circle into ENZ file
        outstring = format(str(f'{e2:.5f}'))  + ' ' + str(f'{n2:.5f}') + ' ' + str(f'{z2:.5f}') + '\n'
        outfile.write(outstring)
    outfile.close()


    # Plot
    plt.figure(figsize=(7, 7), num='akayurin@gmail.com - Min radius')
    
    ax = plt.axes(xlim=(ce_min-5*r_min, ce_min+5*r_min), ylim=(cn_min-5*r_min, cn_min+5*r_min), aspect=1)
    ax.grid(which='major')
    formatter = matplotlib.ticker.FormatStrFormatter('%.0f')
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
    
    plt.plot(ol['E'], ol['N'], 'b.') # original line
    plt.plot(ol['E'], ol['N'], 'b--', lw=0.5)
    plt.plot(e_3d, n_3d, 'g') # 3D circle (projection to EN plane)
    circle1=plt.Circle((ce_min, cn_min), r_min, color='r', ls='--', fill=False) # 2D circle of r_min
    ax.add_patch(circle1)   
    plt.annotate(f'min R={r_min:.2f}m @ Point {point_no} \n      (Decimation={dec:.0f})', xy=(ol['E'][point_no], ol['N'][point_no]), xytext=(ol['E'][point_no], ol['N'][point_no]), arrowprops=dict(facecolor='r'))

    message = '3D circle file created:' + circle_filename
    mb.showinfo('Min radius', message)
    plt.show()


#%% Form
class HomePage:
    def __init__(self):
        self.root = Tk()
        self.root.title('Min radius 3D Single decimation akayurin@gmail.com')
        self.root.geometry('470x400+50+50')
        self.root.resizable(0,0)

        self.dec_value = IntVar()
        self.tick_dec_interval = IntVar()
        
        self.min_low_value = IntVar()
        self.max_low_value = IntVar()
        self.tick_low_interval = IntVar()
        
        self.min_high_value = IntVar()
        self.max_high_value = IntVar()  
        self.tick_high_interval = IntVar()
                
        self.lb1var = StringVar()
        self.lb1var.set('No file selected')
        self.lb2var = StringVar() 
        self.lb3var = StringVar() 
        self.lb4var = StringVar() 
        self.lb5var = StringVar() 

        self.lb_low = Label(self.root, text = 'From point:', width=21).place(relx=0.03, rely=0.3)
        self.lb_high = Label(self.root, text = 'To point:', width=21).place(relx=0.03, rely=0.45)
        self.lb_dec = Label(self.root, text = 'Decimation:', width=21).place(relx=0.03, rely=0.6)
       
        self.but1 = Button(self.root, width=16, height=3, text='Select line file', command=open_line_file).place(relx=0.1, rely=0.03)
        self.but2 = Button(self.root , width=16, height=3, text='Run', command=minrad).place(relx=0.1, rely=0.8)
        self.but3 = Button(self.root, width=16, height=3, text='Close', command=self.form_close).place(relx=0.53, rely=0.80)
        self.lb1 = Label(self.root, textvar = self.lb1var, anchor='w', width=25).place(relx=0.5, rely=0.02)
        self.lb2 = Label(self.root, textvar = self.lb2var, anchor='w', width=25).place(relx=0.5, rely=0.065)
        self.lb3 = Label(self.root, textvar = self.lb3var, anchor='w', width=25).place(relx=0.5, rely=0.11)        
        self.lb4 = Label(self.root, textvar = self.lb4var, anchor='w', width=25).place(relx=0.5, rely=0.155)  
        self.lb5 = Label(self.root, textvar = self.lb5var, anchor='w', width=25).place(relx=0.5, rely=0.2)  
        
        self.low_scale = Scale(self.root, from_=self.min_low_value.get(), to=self.max_low_value.get(), orient='horizontal', length=300, tickinterval= self.tick_low_interval.get(), command=self.low_scale_move)
        self.low_scale.place(relx=0.3, rely=0.25)        

        self.high_scale = Scale(self.root, from_=self.min_high_value.get(), to=self.max_high_value.get(), orient='horizontal', length=300, tickinterval=self.tick_high_interval.get(), command=self.high_scale_move)
        self.high_scale.place(relx=0.3, rely=0.4)        

        self.dec_scale = Scale(self.root, from_=1, to=self.dec_value.get(), orient='horizontal', length=300, tickinterval=self.tick_dec_interval.get(), command=self.dec_scale_move)
        self.dec_scale.place(relx=0.3, rely=0.55)
        
                
    def scales_set(self):                     
        self.tick_low_interval.set((self.max_low_value.get() - self.min_low_value.get()) / 5)
        self.low_scale.config(from_=self.min_low_value.get(), to=self.max_low_value.get(), tickinterval=self.tick_low_interval.get())
        self.tick_high_interval.set((self.max_high_value.get() -  self.min_high_value.get()) / 5)
        self.high_scale.config(from_=self.min_high_value.get(), to=self.max_high_value.get(), tickinterval=self.tick_high_interval.get())
        self.high_scale.set(self.max_high_value.get())       
        self.dec_scale.config(to=self.dec_value.get(), tickinterval=self.tick_dec_interval.get())                
        app.lb2var.set('No of iterations: ' + str(app.dec_scale.get() * (app.high_scale.get() - app.low_scale.get())))

    
    def high_scale_move(self, val):
        self.tick_low_interval.set((self.max_low_value.get() - self.min_low_value.get()) / 5)
        self.low_scale.config(to=(self.high_scale.get() - 3), tickinterval=self.tick_low_interval.get())        
        self.dec_value.set(self.high_scale.get() - self.low_scale.get())
        self.tick_dec_interval.set((self.high_scale.get() - self.low_scale.get()) / 5)
        self.dec_scale.config(to=((self.high_scale.get() - self.low_scale.get()) / 2 - 1), tickinterval=self.tick_dec_interval.get())
        app.lb2var.set('No of iterations: ' + str(app.dec_scale.get() * (app.high_scale.get() - app.low_scale.get())))

        
    def low_scale_move(self, val):
        self.tick_high_interval.set((self.max_high_value.get() -  self.min_high_value.get()) / 5)
        self.high_scale.config(from_=(self.low_scale.get() + 3), tickinterval=self.tick_high_interval.get())       
        self.dec_value.set(self.high_scale.get() - self.low_scale.get())
        self.tick_dec_interval.set((self.high_scale.get() - self.low_scale.get()) / 5)
        self.dec_scale.config(to=((self.high_scale.get() - self.low_scale.get()) / 2 - 1), tickinterval=self.tick_dec_interval.get())
        app.lb2var.set('No of iterations: ' + str(app.dec_scale.get() * (app.high_scale.get() - app.low_scale.get())))


    def dec_scale_move(self, val):
        app.lb2var.set('No of iterations: ' + str(app.dec_scale.get() * (app.high_scale.get() - app.low_scale.get())))
               
    
    def form_close(self):
        self.root.destroy()


if __name__ == '__main__':
    app = HomePage()
    app.root.mainloop()
