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
infile.readline()
for line in infile:
    line = line.split()
    rho.append(eval(line[0]))
    u_a.append(eval(line[1]))
    u_n.append(eval(line[2]))
infile.close()
print(eigenvalue)
plt.plot(rho,u_a)
plt.plot(rho,u_n)
plt.show()
