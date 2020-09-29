import numpy as np
import matplotlib.pyplot as plt

path = "../data/"
from os import listdir
from os.path import isfile, join
filenames = [f for f in listdir(path) if isfile(join(path, f))]


for i_,x_ in enumerate(filenames):
    if (x_[-10:]!= "jacobi.txt"):
        filenames.pop(i_)
    filenames[i_] = eval(filenames[i_][1:-11])

filenames.sort()

for i in range(len(filenames)):
    filenames[i] = "N" + str(filenames[i]) + "_jacobi.txt"

n = []
iter = []
for file in filenames:
    infile = open(path+file,"r")
    first_line = infile.readline().split()
    n.append(eval(first_line[2][2:]))
    iter.append(eval(first_line[1][11:]))
    infile.close()
n = np.asarray(n)
iter = np.asarray(iter)


plt.figure(figsize=(9, 6))

# Remove the plot frame lines. They are unnecessary chartjunk.
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)


ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

range_ = range(n[0],n[-1])
for y in range(iter[2], iter[-1], int((iter[-1]-iter[0])/5)):
    plt.plot(range_, [y] * len(range_), "--", lw=1.0, color="black", alpha=0.3)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.tick_params(axis=u"both", which=u"both", bottom="off",length=0)


plt.xlabel("n",fontsize=14)
plt.ylabel("Iterations",fontsize=14)

clr_2 = (0,0.5,1)
plt.plot(n,iter,"o",color=clr_2)
plt.plot(n,iter,"--",lw=2,color=clr_2)
plt.text(n[-1]+10, iter[-1]-10, "Datapoints", fontsize=14, color=clr_2)

clr_2 = (0.5,0,0.5)
plt.plot(n,n**2,lw=2,color=clr_2)
plt.text(n[-1]+10, n[-1]**2-10, r"$n^2$", fontsize=14, color=clr_2)

plt.savefig("../figures/iteration_mot_n.png", bbox_inches="tight")
plt.clf()
