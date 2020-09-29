import numpy as np
import matplotlib.pyplot as plt

infile = open("../data/N500_jacobi.txt","r")
first_line = infile.readline().split()
tolerance = eval(first_line[0].split("=")[1])
iterations = eval(first_line[1].split("=")[1])
n = eval(first_line[2].split("=")[1])
infile.readline()
rho = infile.readline().split()
for i_,x_ in enumerate(rho):
    rho[i_] = eval(x_)

infile.readline()

Llambda = []
u = []

temp = []
for line in infile:
    line = line.split()
    Llambda.append(eval(line[0]))

    for u_val in line[1:]:
        temp.append(eval(u_val))
    u.append(temp)
    temp = []


for i in range(0,3):
    plt.plot(rho,u[i],label=r"$\lambda_%g=%.2f$"%(i+1,Llambda[i]))
    #plt.text(rho[15],u[i][15],r"$\lambda_%g=%.2f$"%(i+1,Llambda[i]),fontsize=14)
plt.xlabel(r"$\rho$",fontsize=12)
plt.ylabel("u",fontsize=12)
plt.legend()
plt.savefig("lambda_one_electron.png")
