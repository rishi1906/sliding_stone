import matplotlib.pyplot as plt
import csv

x = []
y = []

with open('output.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))
#print(x)
#print(y)
plt.figure(dpi=150)
plt.plot(x,y)
plt.xlabel('Time')
plt.ylabel('Distance')
plt.title('Block Move Dynamics')

plt.savefig('Block_Move.png')
plt.show()