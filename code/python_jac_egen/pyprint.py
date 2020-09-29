# Plotting data from cpp-programs. Eigenvalue problems.
# Plotting numerical eigenvector vs analytic eigenvector for the
# lowest eigenvalue.

import numpy as np
import matplotlib.pyplot as plt

n = 100
N = int(n+1)
h = 1/N
d = 2/h**2
a = -1/h**2

eigval = [] # list containing eigenvalues.
eigvec = np.zeros((N-1,N-1)) # matrix containing N-1 eigenvectors with dim N-1.

Lambda = np.zeros(N) # array with analytical eigenvalues.

for j in range(N):
    Lambda[j] = d + 2*a*np.cos((j+1)*np.pi/N)

u_1 = np.zeros(N-1) # analytical eigenvector for lowest lambda, j = 1.
for p in range(N-1):
    u_1[p] = np.sin((p+1)*np.pi/N)
u_1 = u_1/np.linalg.norm(u_1)

infile = open("eigenpairs.txt","r") # file is produced from cpp-programs.
infile.readline()

i = 0 # number that corresponds to the line that we read.
for line in infile:
    values = line.split() # splitting line into words (numbers).
    infile.readline()
    eigval.append(eval(values[0]))
    for k in range(N-1):
        eigvec[i][k] = eval(values[k+1])
    i = i + 1 # next line.

infile.close()

plt.plot(u_1,label="Analytical")
plt.plot(eigvec[0][0:], label="numerical")
plt.title(r"Eigenvector for the lowest eigenvalue $\lambda_0={:}$".format(eigval[0]))
plt.legend()
plt.show()

error = np.zeros(N-1)
for i in range(N-1):
    error[i] = abs((eigvec[0][i]-u_1[i])/u_1[i])

plt.plot(error)
plt.title("Rel.Error between analytical and numerical eigvec")
plt.show()


"""
# Checking that the norm is 1:
norm_analytic = np.linalg.norm(u_1)
norm_numeric = np.linalg.norm(eigvec[0][0:])
print(norm_analytic)
print(norm_numeric)
"""
