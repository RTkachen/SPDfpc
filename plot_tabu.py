import pandas as pd
import matplotlib.pyplot as plt

# Wczytaj dane
df = pd.read_csv('tabu_history.csv')

# Rysuj wykres
plt.figure()
plt.plot(df['sampleIdx'], df['bestCost'])
plt.xlabel('Numer próbki')
plt.ylabel('Najlepszy Cmax')
plt.title('Postęp Tabu Search')
plt.grid(True)
plt.show()
