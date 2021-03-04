import numpy as np
import matplotlib.pyplot as plt
from scinum import scinum as sn
from scipy.optimize import curve_fit
import pandas as pd
import scipy
plt.style.use('ggplot')

class dataa:
    def __init__(self, filename):
        self.filename = filename
        self.data = pd.read_csv(filename, delimiter= ',', header = 0)

    def show(self):
        return self.data

    def linear(self, dy = 0, dx = 0, title = "", ylabel = "", xlabel = "", yname = "", xname = ""):        
        Y = sn(self.data[yname], uncer = dv)        
        X = sn(self.data[xname], uncer = di)
        
        x = sn.f(self.up, X)
        y = sn.f(self.down, Y)
        z, cov = np.polyfit(x.n, y.n, 1, w= 1/y.u**2, cov = True)
        p = np.poly1d(z)
        
        plt.rcParams.update({'font.size': 11})
        fig = plt.figure(dpi = 300, figsize = (7,3))
        plt.title(title, fontsize = 16)
        plt.ylabel(ylabel, fontsize = 14)
        plt.errorbar(x.n, y.n, yerr = y.u, xerr = x.u, fmt = '.')
        plt.plot(x.n, p(x.n))
        plt.legend(fontsize = 10)
        ax1 = fig.add_axes((.124,-.1,.775,.2))
        ax1.errorbar(x.n, p(x.n) - y.n, yerr = y.u, xerr = x.u, label = 'Residuals', fmt = '.')
        ax1.set_xlabel(xlabel, fontsize = 14)
        plt.show()       

        print(z)
        fit_u = np.sqrt(np.diag(cov))
        print(fit_u)
        
        chi2 = self.linear_chisq(x.n, y.n, y.u, *z)
        
        print("Chi-Square: ", chi2)
        k = len(x.n) - len(z)
        print("degrees of freedom: ", k)
        print([k + (i+2)*np.sqrt(2*k) for i in range(3)])
        
    def nonlinear(self, f = None, dy = 0, dx = 0, title = "", ylabel = "", xlabel = "", yname = "", xname = "", p0 = None):
        
        Y = sn(self.data[yname], uncer = dv)        
        X = sn(self.data[xname], uncer = di)
        
        x = sn.f(self.up, X)
        y = sn.f(self.down, Y)
        
        best_params, fit_cov = curve_fit(f, x.n, y.n, p0)
        print(best_params)
        print(np.sqrt(np.diag(fit_cov)))

        print('chi sq: ', nonlinear_chisq(x.n, y.n, y.u, f, best_params))
        print('ratio: ', nonlinear_chisq(x.n, y.n, y.u, f, best_params)/(len(x.n) - len(best_params)))

        fit = [f(j, *best_params) for j in x.n]

        plt.rcParams.update({'font.size': 11})
        fig = plt.figure(dpi = 300, figsize = (7,3))
        plt.title(title, fontsize = 16)
        plt.ylabel(ylabel, fontsize = 14)
        plt.errorbar(x.n, y.n, yerr = y.u, xerr = x.u, fmt = '.', label = 'Data')
        plt.plot(x.n, fit, label = 'Fit')
        plt.legend(fontsize = 10)
        ax1 = fig.add_axes((.124,-.1,.775,.2))
        ax1.errorbar(x.n, fit - y.n, yerr = y.u, xerr = x.u, label = 'Residuals', fmt = '.')
        ax1.set_xlabel(xlabel, fontsize = 14)
        plt.show()

    def hist(self, label, bins = None):
        plt.hist(self.data[label], bins = bins)
        plt.xlabel(label)
        plt.ylabel('Count')
    
    def up(self, u):
        return u
    
    def down(self, d):
        return d
        
    def gaussian(self, x, mean, std):
        return 1/(std*np.sqrt(2*np.pi))*np.exp((-(x-mean)**2)/(2*std**2))
    
    def poisson(x, lamda, k):
        return lamda**k*np.exp(-lamda)/scipy.special.factorial(k)
    
    def get_a(self, x, y):
        return np.sum(x*y)/np.sum(x**2)

    def damped_trig(self, t, g, w, p):
        return np.exp(g*t)*np.cos(w*t+p)

    def nonlinear_chisq(self, x,y,ey,f, args):
        return np.sum((y-f(x, *args))**2/ey**2)
    
    def linear_chisq(self, x,y,ey, a, b):
        return np.sum((y-a*x-b)**2/ey**2)
    
    def get_chisq(self, x,y,ey,a, b, mode = "linear"):
        if mode == "linear":
            print("Chi-Square: ", self.linear_chisq(x, y, ey, a, b))
            k = len(x.n) - 2
            print("degrees of freedom: ", k)
            print([k + (i+2)*np.sqrt(2*k) for i in range(3)])
            print("CHISQ FACTOR: ", self.linear_chisq(x, y, ey, a, b)/k)
        else:
            print("Chi-Square: ", self.nonlinear_chisq(x, y, ey, a, b))
            k = len(x.n) - len(b)
            print("degrees of freedom: ", k)
            print([k + (i+2)*np.sqrt(2*k) for i in range(3)])
            print("CHISQ FACTOR: ", self.nonlinear_chisq(x, y, ey, a, b)/k)

    def weighted_x(self, x, u):
    	return np.sum(x/u**2)/(np.sum(1/u**2))

	def weighted_u(self, x, u):
	    return np.sqrt(1/(np.sum(1/u**2)))
