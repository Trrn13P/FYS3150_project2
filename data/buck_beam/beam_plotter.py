import numpy as np
import matplotlib.pyplot as plt

path = "../buck_beam/"
file = "N100_jacobi.txt"

infile = open(path+file,"r")

first_line = infile.readline().split()
second_line = infile.readline().split()

eigenvalue = eval(second_line[1])

rho = []
u_a = []
u_n = []

for i in infile.readline():



print(eigenvalue)
plt.plot(rho,u_a)
plt.plot(rho,u_n)
plt.show()
