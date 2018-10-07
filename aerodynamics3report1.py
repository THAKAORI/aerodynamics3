#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 13:19:51 2018

@author: takahiro
"""
import numpy as np
from numpy import pi, cos, sin
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D

class SIM3DOF:
    def __init__(self):
        self.tcos = 20971
        self.tsin = 18.7
        self.m = 7500
        self.dt = 0.1
        self.T = np.arange(0.0, 80.0 + self.dt, self.dt)
        self.N = self.T.shape[0]
        self.phi = np.zeros(self.N)
        for i in range(int(self.N / 4), self.N - int(self.N / 4)):
            self.phi[i] = 30 * pi / 180 * sin(pi / 20 * (self.T[i] - 20))
        self.X = np.zeros((self.N, 6))
        self.X[0, 2] = -5000
        self.X[0, 3] = 180
        self.runge()
    def runge(self):
        k1 = np.array(np.zeros(6))
        k2 = np.array(np.zeros(6))
        k3 = np.array(np.zeros(6))
        k4 = np.array(np.zeros(6))
        for i in range(self.N - 1):
                k1 = self._deriv(self.X[i,:], i)
                k2 = self._deriv(self.X[i,:] + 0.5 * self.dt * k1, i)
                k3 = self._deriv(self.X[i,:] + 0.5 * self.dt * k2, i)
                k4 = self._deriv(self.X[i,:] + self.dt * k3, i)
                self.X[i+1,:] = self.X[i,:] + self.dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    def _deriv(self, X, i):
        d_xe = X[3] * cos(X[4]) * cos(X[5])
        d_ye = X[3] * cos(X[4]) * sin(X[5])
        d_ze = -X[3] * sin(X[4])
        d_U = (self.tcos - self.getDrag(X[3])) / self.m - (9.8 * sin(X[4]))
        d_gamma = ((self.getLift(X[3]) + self.tsin) / self.m * cos(self.phi[i]) - (9.8 * cos(X[4]))) / X[3]
        d_psi = ((self.getLift(X[3]) + self.tsin) * sin(self.phi[i])) / (X[3] * self.m * cos(X[4]))
        temp = np.array([d_xe, d_ye, d_ze, d_U, d_gamma, d_psi])
        return temp
    def getLift(self, U):
        return (2.26827 * (U**2))
    def getDrag(self, U):
        return (0.6471 * (U**2))
    
def graph2d(X, T):
    fig,axes = plt.subplots(nrows=2,ncols=2,figsize=(10,8))
    pdf = PdfPages('aero2d.pdf')
    
    xe = X[:, 0]
    axes[0,0].plot(T, xe, linewidth=2)
    axes[0,0].set_title('xe')
    axes[0,0].set_xlabel('t')
    axes[0,0].set_ylabel('(m)')
    axes[0,0].grid(True)
    
    ye = X[:, 1]
    axes[0,1].plot(T, ye, linewidth=2)
    axes[0,1].set_title('ye')
    axes[0,1].set_xlabel('t')
    axes[0,1].set_ylabel('(m)')
    axes[0,1].grid(True)
    
    ze = X[:, 2]
    axes[1,0].plot(T, -ze, linewidth=2)
    axes[1,0].set_title('altitude')
    axes[1,0].set_xlabel('t')
    axes[1,0].set_ylabel('(m)')
    axes[1,0].grid(True)
    
    U = X[:, 3]
    axes[1,1].plot(T, U, linewidth=2)
    axes[1,1].set_title('velocity')
    axes[1,1].set_xlabel('t')
    axes[1,1].set_ylabel('(m)')
    axes[1,1].grid(True)

    pdf.savefig()
    pdf.close()

def graph3d(X):
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot(X[:, 0], -X[:, 1], -X[:, 2], "o-", color="#00aa00", ms=4, mew=0.5)
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    
    plt.show()
    fig.savefig("aero3d.pdf", bbox_inches='tight', pad_inches=0.3)
sim3dof=SIM3DOF()
graph2d(sim3dof.X, sim3dof.T)
graph3d(sim3dof.X)