""" A little class for finding and propagating uncertainties in measured quantities, formatting them into LaTex tables
"""

import numpy as np
from tabulate import tabulate
from scipy.misc import derivative
import sympy

class scinum:

    """ A scinum object is defined by a set of values and respective undertainties and can be manipulated to propagate uncertainties through operations and produce data tables
    
    Attributes:
        len (int): Number of uncertain quantities
        n (list of float): list of uncertain numbers
        u (list of float): list of uncertainties
    """
    
    def __init__(self, arr, uncer = None): 
        """Initialize class variables
        
        Args:
            arr (list): Either 2D array shape (len, 2) containing numbers and uncertainty pairs or 1D array of numbers with u passed in 'uncer'
            uncer (None, optional): If all 'n' have same uncertainty, it is passed here
        """

        if (uncer is not None):
            self.n = np.array(arr)
            self.len = len(np.array(arr))
            self.u = np.array([uncer for i in range(self.len)])
        else:    
            self.n = np.array(np.array(arr).T[0])
            self.u = np.array(np.array(arr).T[1])
            self.len = len(np.array(arr))
        
    def __mul__(self, data):
        """Multiplication operator overloading
        
        Args:
            data (scium object or float/list): values to multiply with
        
        Returns:
            scinum object: result of multiplication
        """
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
        """Multiplication operator overloading
        
        Args:
            data (scium object or float/list): values to divide
        
        Returns:
            scinum object: result of division
        """
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
        """Addition operator overloading
        
        Args:
            data (scium object or float/list): values to add
        
        Returns:
            scinum object: result of addition
        """
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
        """Subtraction operator overloading
        
        Args:
            data (scium object or float/list): values to subtract
        
        Returns:
            scinum object: result of subtraction
        """
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
        """Exponentiation operator overloading
        
        Args:
            data (scium object or float/list): values to exponentiate
        
        Returns:
            scinum object: result of exponentiation
        """
        if (isinstance(data, scinum)):
            if (np.sum(data.u) == 0):
                arr = []
                arr.append(self.n**data.n)
                arr.append((data.n*self.u/self.n)*arr[0])
                return scinum(np.array(arr).T)
            else:
                def exp(a, b):
                    return a**b

                return scinum.f(exp, self, data)
        else:
            arr = []
            arr.append(self.n**data)
            arr.append((data*self.u/self.n)*arr[0])
            return scinum(np.array(arr).T)
    
    def __str__(self):
        """ Casts scinum object to string for printing
        
        Returns:
            str: list of self.n and self.u
        """
        return str(np.array([[self.n[i], self.u[i]] for i in range(self.len)]))
    
    def func(self, fn):
        """Apply a function depending only on 'self' to the existing object
        
        Args:
            fn (function): Function to apply
        
        Returns:
            scinum object: output of function
        """
        arr = []
        arr.append(fn(self.n))
        arr.append(np.abs(derivative(fn, self.n, dx=1e-6))*self.u)

        return scinum(np.array(arr).T)
            
    def pr(self):
        """Short function to print scinum object
        """
        print(np.array([[self.n[i], self.u[i]] for i in range(self.len)]))
        
    def present(self):
        """Appropriately rounds numbers to correct number of decimals for data presentation
        
        Returns:
            scinum object: scinum object containing rounded values
        """
        return scinum(self.pres())
    
    def checkzero(self, x):
        """Checks if value to be presented has a zero as the last digit
        
        Args:
            x (str): String of number
        
        Returns:
            str: number with 0 added
        """
        while (len(str(x[0]).split('.')[1]) != len(str(x[1]).split('.')[1])):
            x[0] = str(str(x[0])+"0")
            
        return x[0]
    
    def pres(self):
        """Appropriately rounds numbers to correct number of decimals for data presentation
        
        Returns:
            2D list of floats: list containing rounded values
        """
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
        
    def gettable(self, columns = [], headers = None, fmt = "latex"):
        """Wrapper for 'tabulate', returning a table containing 'self' values, in latex/markdown etc.
        
        Args:
            columns (list of scinum object, optional): Additonal scinum objects to add to table
            headers (list of string, optional): Headers for columns in the table
            fmt (str, optional): Table format, see 'tabulate' for options (ex. github/latex/markdown etc.)
        """
        if (headers is None):
            headers = ['' for i in range(len(columns) + 1)]
        
        x = self.pres()
        tab = np.array([[str(str(x[i][0]) + " (" + str(x[i][1]) + ")")] if str(x[i][1]) != "0" else [str(x[i][0])] for i in range(self.len)]) 
        
        for column in columns:
            x = column.pres()
            tab = np.hstack((tab, [[str(str(x[i][0]) + " (" + str(x[i][1]) + ")")] if str(x[i][1]) != "0" else [str(x[i][0])] for i in range(column.len)]))
            
        print(tabulate(tab, headers, tablefmt=fmt))
    
    @staticmethod
    def expanddims(num, dims):
        """Expanding a single valued scinum object for an operation with a larger scinum object
        
        Args:
            num (scinum object): Single value scinum object to be expanded
            dims (int): number of dimenions to expand to
        
        Returns:
            scinum object: expanded values
        """
        print("Warning, arrays are not the same size, expanding dims...")
        return scinum(np.array([[num.n[0], num.u[0]] for i in range(dims)]))
    
    @staticmethod
    def table(columns= [], headers = None, fmt = "github"):
        """Wrapper for 'tabulate', returning a table containing values for a list of scinum objects, in latex/markdown etc.
        
        Args:
            columns (list of scinum object, optional): Additonal scinum objects to add to table
            headers (list of string, optional): Headers for columns in the table
            fmt (str, optional): Table format, see 'tabulate' for options (ex. github/latex/markdown etc.)
        """
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
        """ Propagates uncertainty through ANY smooth function of scinum objects
        
        Args:
            fn (function): function to propagate uncertainties through
            *snums: list of scinum objects to be passed into the function
        
        Returns:
            scinum object: Output of function with uncertainty
        """
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
        """Computes the partial derivative of ANY smooth function with respect to a given variable
        
        Args:
            fn (function): function to find partial derivative
            n (int): index of variable of interest passed into function
            *args: function arguments
        
        Returns:
            float: numerical partial derivative
        """
        h = 1e-6
        var = list(args)
        var[n] = var[n] + h
        var = np.array(var)
        args = np.array(args)
        return (fn(*var) - fn(*args))/h
        