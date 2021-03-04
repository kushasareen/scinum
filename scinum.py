import numpy as np
from tabulate import tabulate
from scipy.misc import derivative
import sympy

class scinum:
    def __init__(self, arr, uncer = None): 
        if (uncer is not None):
            self.n = np.array(arr)
            self.len = len(np.array(arr))
            self.u = np.array([uncer for i in range(self.len)])
        else:    
            self.n = np.array(np.array(arr).T[0])
            self.u = np.array(np.array(arr).T[1])
            self.len = len(np.array(arr))
        
    def __mul__(self, data):
        if (isinstance(data, scinum)):
            if (self.len < data.len):
                self = scinum.expanddims(self, data.len)
            elif (self.len > data.len):
                data = scinum.expanddims(data, self.len)

            arr = []
            arr.append(self.n*data.n)
            arr.append(np.sqrt((self.u/self.n)**2 + (data.u/data.n)**2)*arr[0])

            return scinum(np.array(arr).T)
        else:
            arr = []
            arr.append(self.n*data)
            arr.append(self.u*data)

            return scinum(np.array(arr).T)
    
    def __truediv__(self, data):
        if (isinstance(data, scinum)):
            if (self.len < data.len):
                self = scinum.expanddims(self, data.len)
            elif (self.len > data.len):
                data = scinum.expanddims(data, self.len)

            arr = []
            arr.append(self.n/data.n)
            arr.append((self.u/self.n + data.u/data.n)*arr[0])

            return scinum(np.array(arr).T)
        else:
            arr = []
            arr.append(self.n/data)
            arr.append(self.u/data)

            return scinum(np.array(arr).T)
            
    
    def __add__(self, data):
        if (isinstance(data, scinum)):
            if (self.len < data.len):
                self = scinum.expanddims(self, data.len)
            elif (self.len > data.len):
                data = scinum.expanddims(data, self.len)

            arr = []
            arr.append(self.n+data.n)
            arr.append(np.sqrt(self.u**2 + data.u**2))

            return scinum(np.array(arr).T)
        else:
            arr = []
            arr.append(self.n + data)
            arr.append(self.u)

            return scinum(np.array(arr).T)
            
    
    def __sub__(self, data):
        if (isinstance(data, scinum)):
            if (self.len < data.len):
                self = scinum.expanddims(self, data.len)
            elif (self.len > data.len):
                data = scinum.expanddims(data, self.len)

            arr = []
            arr.append(self.n-data.n)
            arr.append(np.sqrt(self.u**2 + data.u**2))

            return scinum(np.array(arr).T)
        else:
            arr = []
            arr.append(self.n - data)
            arr.append(self.u)

            return scinum(np.array(arr).T)
        
    def __pow__(self, data):
        if (isinstance(data, scinum)):
            if (np.sum(data.u) == 0):
                arr = []
                arr.append(self.n**data.n)
                arr.append((data.n*self.u/self.n)*arr[0])
                return scinum(np.array(arr).T)
            else:
                print("scinum exponentiation not yet supported")
                return
        else:
            arr = []
            arr.append(self.n**data)
            arr.append((data*self.u/self.n)*arr[0])

            return scinum(np.array(arr).T)
    
    def __str__(self):
        return str(np.array([[self.n[i], self.u[i]] for i in range(self.len)]))
    
    def func(self, fn):
        arr = []
        arr.append(fn(self.n))
        arr.append(np.abs(derivative(fn, self.n, dx=1e-6))*self.u)

        return scinum(np.array(arr).T)
            
    def pr(self):
        print(np.array([[self.n[i], self.u[i]] for i in range(self.len)]))
        
    def present(self):
        x = []
        for i in range(self.len):
            if (-1*int(('%e' % self.u[i]).partition('e')[2]) > 0):
                x.append([round(self.n[i], -1*int(('%e' % self.u[i]).partition('e')[2])), 
                          round(self.u[i], -1*int(('%e' % self.u[i]).partition('e')[2]))])
            
            elif (-1*int(('%e' % self.u[i]).partition('e')[2]) <= 0):
                x.append([round(int(self.n[i]), -1*int(('%e' % self.u[i]).partition('e')[2])), 
                          round(int(self.u[i]), -1*int(('%e' % self.u[i]).partition('e')[2]))])
            else:
                x.append([round(self.n[i]), round(self.u[i])])

        return scinum(x)
    
    def checkzero(self, x):
        while (len(str(x[0]).split('.')[1]) != len(str(x[1]).split('.')[1])):
            x[0] = str(str(x[0])+"0")
            
        return x[0]
    
    def pres(self):
        x = []
        for i in range(self.len):
            if (-1*int(('%e' % self.u[i]).partition('e')[2]) > 0):
                x.append([self.checkzero([round(self.n[i], -1*int(('%e' % self.u[i]).partition('e')[2])), 
                                         round(self.u[i], -1*int(('%e' % self.u[i]).partition('e')[2]))]), 
                          round(self.u[i], -1*int(('%e' % self.u[i]).partition('e')[2]))])
            
            elif (-1*int(('%e' % self.u[i]).partition('e')[2]) <= 0):
                x.append([round(int(self.n[i]), -1*int(('%e' % self.u[i]).partition('e')[2])), 
                          round(int(self.u[i]), -1*int(('%e' % self.u[i]).partition('e')[2]))])
            else:
                x.append([round(self.n[i]), round(self.u[i])])

        return x              
        
    def gettable(self, columns = [], headers = None, fmt = "github"):
        if (headers is None):
            headers = ['' for i in range(len(columns) + 1)]
        
        x = self.pres()
        tab = np.array([[str(str(x[i][0]) + " (" + str(x[i][1]) + ")")] if str(x[i][1]) != "0" else [str(x[i][0])] for i in range(self.len)]) 
        
        for column in columns:
            x = column.pres()
            tab = np.hstack((tab, [[str(str(x[i][0]) + " (" + str(x[i][1]) + ")")] if str(x[i][1]) != "0" else [str(x[i][0])] for i in range(column.len)]))
            
#         print(tab)
        print(tabulate(tab, headers, tablefmt=fmt))
    
    @staticmethod
    def expanddims(num, dims):
        print("Warning, arrays are not the same size, expanding dims...")
        return scinum(np.array([[num.n[0], num.u[0]] for i in range(dims)]))
    
    @staticmethod
    def table(columns= [], headers = None, fmt = "github"):
        if (headers is None):
            headers = ['' for i in range(len(columns))]
            
        x = columns[0].pres()
        tab = np.array([[str(str(x[i][0]) + " (" + str(x[i][1]) + ")")] if str(x[i][1]) != "0" else [str(x[i][0])] for i in range(columns[0].len)]) 
        
        for column in columns[1:]:
            x = column.pres()
            tab = np.hstack((tab, [[str(str(x[i][0]) + " (" + str(x[i][1]) + ")")] if str(x[i][1]) != "0" else [str(x[i][0])] for i in range(column.len)]))
        
        print(tabulate(tab, headers, tablefmt=fmt))
        
    @staticmethod
    def f(fn, *snums):
        arr = []
        a = np.array([snums[i].n for i in range(len(snums))])
        arr.append(fn(*a))
        
        s = 0
        for i in range(len(snums)):
            s+= (scinum.df(fn, i, *a)**2)*(snums[i].u)**2
        
        arr.append(np.sqrt(s))
        return scinum(np.array(arr).T)
    
    @staticmethod
    def df(fn, n, *args):
        h = 1e-6
        var = list(args)
        var[n] = var[n] + h
        var = np.array(var)
        args = np.array(args)
        return (fn(*var) - fn(*args))/h
        