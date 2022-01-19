import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("data.csv")

N = data.iloc[:, 4]
F = data.iloc[:, 5]
V = data.iloc[:, 6]

X = [x/100 for x in range(len(N))]

plt.figure(1)
plt.plot(X, N)
plt.plot(X, F)
plt.xlabel('t, c')
plt.ylabel('F, H')
plt.legend(['N', 'F'])
plt.title('Силы реакции')
plt.show()

plt.figure(2)
plt.plot(X, V)
plt.xlabel('t, c')
plt.ylabel('V, m/c')
plt.title('Скорость')
plt.show()