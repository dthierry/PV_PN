import matplotlib.pyplot as plt

with open("res.txt") as f:
    lines = f.readlines()

hour = []
pn = []
soc = []
sigma1 = []
sigma2 = []
sigma3 = []
u_value = []

for line in lines[1:]:
    s = line.split()
    sfloat = [float(i) for i in s]
    hour.append(int(sfloat[0]))
    pn.append(sfloat[1])
    soc.append(sfloat[2])
    sigma1.append(sfloat[3])
    sigma2.append(sfloat[4])
    sigma3.append(sfloat[5])
    u_value.append(sfloat[6])

d = [1 if i > 0 else 0 for i in sigma1]

fig, ax = plt.subplots()
lns0 = ax.plot(hour, pn, "b--", label="$P_N$")
ax1 = ax.twinx()

lns1 = ax1.plot(hour, soc, "r.-", label="$SOC$")

ax2 = ax.twinx()
lns2 = ax2.plot(hour, d, "cx", label="Diesel")

lns = lns0 + lns1 + lns2
labs = [l.get_label() for l in lns]

ax.legend(lns, labs, loc=5)


ax.set_xlabel("Hour of the year")
ax.set_ylabel("Net Power kW")
ax1.set_ylim([0, 1])
ax2.set_ylim([0, 1])
#ax2.get_yaxis().set_visible(False)
#ax2.get_yaxis().
ax1.set_ylabel("SOC")
ax2.set_ylabel("Diesel ON", loc="top")
ax.grid()
fig.savefig("my_res0.png")

ax.set_xlim([4000, 4730])
ax1.set_xlim([4000, 4730])
ax2.set_xlim([4000, 4730])
fig.savefig("my_res1.png")