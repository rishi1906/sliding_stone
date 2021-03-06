import matplotlib.pyplot as plt
import csv


x = []
y = []

with open('test_output_1.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))
#print(x)
#print(y)
# plt.figure(dpi=150)
# plt.plot(x,y,label = "Function Curve")

# plt.scatter(x, y,label= "CGL Nodes", color= "red", marker= "*", s=10) 
# plt.xlabel('Time')
# plt.ylabel('Amplitude')
# plt.title('Differentiation Matrix Validation')
# plt.legend()

# plt.savefig('Differentiation_Validation.png')
# plt.show()


a = []
b = []

with open('test_output_2.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        a.append(float(row[0]))
        b.append(float(row[1]))
#print(x)
#print(y)
plt.figure(dpi=150)
plt.plot(x,y,label = "Function Curve")
# plt.hold(True);
plt.plot(a,b,label = "Derivative Curve")
plt.scatter(a, b,label= "CGL Nodes", color= "green", marker= "*", s=10)
plt.scatter(x, y,label= "CGL Nodes", color= "red", marker= "*", s=10)  
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.title('Differentiation Matrix Validation')
plt.legend()

plt.savefig('Differentiation_Validation.png')
plt.show()