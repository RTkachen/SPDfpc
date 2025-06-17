import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("sa_history.csv")
plt.plot(df["sampleIdx"], df["bestCost"], marker='o')
plt.xlabel("Numer prób temperaturowych")
plt.ylabel("Najlepszy C_max")
plt.title("Symulowane Wyżarzanie – zbieżność")
plt.grid(True)
plt.show()