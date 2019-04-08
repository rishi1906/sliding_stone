import matplotlib.pyplot as plt
import csv

x = []
y = []

with open('test_output.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))
#print(x)
#print(y)
plt.figure(dpi=150)
plt.plot(x,y)
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.title('Differentiation Matrix Validation')

plt.savefig('Differentiation_Validation.png')
plt.show()