"""
Lump Model Simulation Script

This script utilizes the Tkinter library for GUI and matplotlib for visualization
to simulate a lumped-parameter model of particle temperature and liquid fraction.

Author: Azeddine Rachih
Date: August 10, 202
"""

from tkinter import *
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import matplotlib as mpl
from matplotlib import rc
mpl.rcParams['legend.numpoints'] = 1
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text')  # , usetex=True)

variables = {"rho_p": [8000., 'Particle density [kg/m3]'],
             "c_ps": [1000., 'Specific heat of the particle solid phase [J/K.kg]'],
             "c_pl": [1000., 'Specific heat of the particle liquid phase [J/K.kg]'],
             "d": [150.e-6, 'Particle diameter [m]'],
             "U0": [1., 'Relative particle velocity in the gas flow [m/s]'],
             "lambda_g": [0.0158, 'Gas thermal conductivity [W/m.K]'],
             "cp_g": [520.64, 'Gas specific heat [J/K.kg]'],
             "rho_g": [1.6228, 'Gas density [kg/m3]'],
             "mu_g": [2.125e-5, 'Gas dynamic viscosity [W/m.K]'],
             "T_L": [800., 'Liquidus [K]'],
             "T_S": [600., 'Solidus [K]'],
             "L": [200000., 'Enthalpy of fusion [J/kg]'],
             "T_inf": [300., 'Temperature of the surrounding gas [K]'],
             "T_0": [1000., 'Particle initial temperature [K]'],
             "t_0": [0., 'Initial instant of simulation [s]'],
             "t_f": [1., 'Final instant of simulation [s]'],
             "sigma": [5.67e-8, 'Stefan constant [W/m2.K4]'],
             "eps": [0.8, 'Particle thermal emissivity [-]']}


class MyWindow:
    def __init__(self, win):
        x0, xt0, y0, interline = 10, 300, 50, 40
        self.var = {}
        vari = {}
        for index, key in enumerate(variables):
            self.lbl = Label(win, text=variables[key][1])
            self.lbl.config(font=('Arial', 9))
            self.lbl.place(x=x0, y=y0 + index * interline)
            self.var[key] = getattr(self, key, Entry())
            self.var[key].place(x=xt0, y=y0 + index * interline)
            self.var[key].insert(END, str(variables[key][0]))
            vari[key] = float(self.var[key].get())

        self.radiation_active = IntVar()
        self.rad = Checkbutton(win, text="Activate radiation", variable=self.radiation_active)
        self.rad.config(font=('Arial', 9))
        self.rad.place(x=100, y=y0 + (index + 1) * interline)

        self.btn = Button(win, text='Compute')
        self.btn.bind('<Button-1>', self.plot)
        self.btn.place(x=200, y=y0 + (index + 2) * interline)

        self.figure = Figure(figsize=(5.5, 3), dpi=100)
        self.subplot1 = self.figure.add_subplot(211)
        self.subplot1.grid()
        self.subplot1.set_title('$Temperature\, \, evolution\, \,T_p(t)$', fontsize=14)
        self.subplot1.set_ylabel('$Particle\, \,Temperature\, \,(K)$', fontsize=11)
        self.subplot1.set_xlim(float(self.var['t_0'].get()), float(self.var['t_f'].get()))
        self.subplot1.set_ylim(float(self.var['T_inf'].get()), float(self.var['T_0'].get()))

        self.subplot2 = self.figure.add_subplot(212)
        self.subplot2.grid()
        self.subplot2.set_title('$Liquid\, \,fraction\, \,evolution\, \,\\alpha(t)$', fontsize=14)
        self.subplot2.set_xlabel('$Time(s)$', fontsize=11)
        self.subplot2.set_ylabel('$Liquid\, \,fraction$', fontsize=11)
        self.subplot2.set_xlim(float(self.var['t_0'].get()), float(self.var['t_f'].get()))
        self.subplot2.set_ylim(-1e-5, 1. + 1e-5)

        self.plots = FigureCanvasTkAgg(self.figure, win)
        self.plots.get_tk_widget().pack(side=RIGHT, fill=BOTH, expand=0)
        #------------ Nusselt number correlation ----------

    def nusselt_num(self, rho, mu, ur, d, cp, lamb):
        Re = rho * d * ur / mu
        Pr = mu * cp / lamb
        NuLam = 0.664 * Re**0.5 * Pr**(1. / 3.)
        NuTurb = 0.037 * Re**0.8 * Pr / \
            (1. + 2.443 * Re**(-0.1) * (Pr**(2. / 3.) - 1.))
        Nu = 2. + (NuLam * NuLam + NuTurb * NuTurb)**0.5
        return Nu
    #----------- Enthalpy in terms of temperature ------------

    def H(self, T):
        vari = {}
        for index, key in enumerate(variables):
            vari[key] = float(self.var[key].get())
        H_S = vari['c_ps'] * vari['T_S']
        c_pm = (vari['c_pl'] + vari['c_ps']) / 2.
        Leff = c_pm * (vari['T_L'] - vari['T_S']) + vari['L']
        H_L = H_S + Leff
        enthalpy = vari['c_ps'] * (T <= vari['T_S']) +\
            (c_pm * (vari['T_L'] - vari['T_S']) +
             H_S + ((T - vari['T_S']) / (vari['T_L'] - vari['T_S'])) * vari['L']) *\
            (T > vari['T_S']) * (T < vari['T_L']) +\
            (vari['c_pl'] * (T - vari['T_L']) + H_L) * (T >= vari['T_L'])
        return enthalpy
    #--------- Tempeature in terms of enthalpy -----------

    def T(self, H):
        vari = {}
        for index, key in enumerate(variables):
            vari[key] = float(self.var[key].get())
        H_S = vari['c_ps'] * vari['T_S']
        c_pm = (vari['c_pl'] + vari['c_ps']) / 2.
        Leff = c_pm * (vari['T_L'] - vari['T_S']) + vari['L']
        H_L = H_S + Leff
        temperature = (H / vari['c_ps']) * (H <= H_S) +\
            (vari['T_S'] + ((H - H_S) / (H_L - H_S)) * (vari['T_L'] - vari['T_S'])) *\
            (H > H_S) * (H < H_L) +\
            (vari['T_L'] + (H - H_L) / vari['c_pl']) * (H >= H_L)
        return temperature
    #-------------- Differential equation ---------------

    def lump(self, H, t, coeff):
        vari = {}
        for index, key in enumerate(variables):
            vari[key] = float(self.var[key].get())
        hp = vari['lambda_g'] * \
            self.nusselt_num(vari['rho_g'],
                             vari['mu_g'], vari['U0'], vari['d'],
                             vari['cp_g'], vari['lambda_g']) / vari['d']
        Volp = (1. / 6.) * np.pi * vari['d']**3
        Surfp = np.pi * vari['d']**2

        K_conv = (hp * Surfp) / (vari['rho_p'] * Volp)
        K_rad = vari['eps'] * vari['sigma'] * \
            (self.T(H) + vari['T_inf']) * (self.T(H)**2 + vari['T_inf']**2) *\
            Surfp / (vari['rho_p'] * Volp)
        K = K_conv + coeff * K_rad
        dHdt = -K * (self.T(H) - vari['T_inf'])
        return dHdt

    def resolution(self):
        vari = {}
        for index, key in enumerate(variables):
            vari[key] = float(self.var[key].get())
        #------------------- ODE resolution ------------------
        H_0 = self.H(vari['T_0'])
        t = np.linspace(vari['t_0'], vari['t_f'], 1000)
        H_p = odeint(self.lump, H_0, t, args=(self.radiation_active.get(),))
        T_p = self.T(H_p)
        return t, T_p

    def plot(self, event):
        t, T_p = self.resolution()

        vari = {}
        for index, key in enumerate(variables):
            vari[key] = float(self.var[key].get())

        alpha = 1. * (T_p >= vari['T_L']) + \
            ((T_p - vari['T_S']) / (vari['T_L'] - vari['T_S']))\
            * (T_p > vari['T_S']) * (T_p < vari['T_L']) + 0. * (T_p <= vari['T_S'])
        # self.plots.get_tk_widget().pack_forget()
        self.subplot1.clear()
        self.subplot2.clear()
        self.subplot1.grid()
        self.subplot2.grid()
        self.subplot1.set_title('$Temperature\, \, evolution\, \,T_p(t)$', fontsize=14)
        self.subplot1.set_ylabel('$Particle\, \,Temperature\, \,(K)$', fontsize=11)
        self.subplot2.set_title('$Liquid\, \,fraction\, \,evolution\, \,\\alpha(t)$', fontsize=14)
        self.subplot2.set_xlabel('$Time(s)$', fontsize=11)
        self.subplot2.set_ylabel('$Liquid\, \,fraction$', fontsize=11)
        self.subplot1.set_xlim(float(self.var['t_0'].get()), float(self.var['t_f'].get()))
        self.subplot1.set_ylim(float(self.var['T_inf'].get()), float(self.var['T_0'].get()))
        self.subplot2.set_xlim(float(self.var['t_0'].get()), float(self.var['t_f'].get()))
        self.subplot2.set_ylim(-1e-5, 1. + 1e-5)
        self.subplot1.plot(t, T_p, 'r', lw=2.5)
        self.subplot1.plot(t, vari['T_L'] * np.ones(np.shape(t)), 'k', lw=1.5)
        self.subplot1.plot(t, vari['T_S'] * np.ones(np.shape(t)), 'k', lw=1.5)
        self.subplot2.plot(t, alpha, 'b', lw=2.5)
        self.figure.savefig('plots.png', bbox_inches='tight')
        self.plots.draw()


if __name__ == "__main__":
    window = Tk()
    mywin = MyWindow(window)
    window.title('Lump model')
    window.geometry("1000x900+10+10")
    window.resizable(width=False, height=False)
    window.mainloop()
