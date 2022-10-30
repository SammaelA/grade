import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('../opt_stat.csv')
print(df)
df["unit_id"] = df["unit_id"].map({0 : "0", 1 : "1", 2 : "2", 3 : "3", 4 : "4", 5 : "5", 6 : "6",
7 : "7", 8 : "8", 9 : "9", 10 : "10", 11 : "11", 12 : "12", 13 : "13", 14 : "14", 15 : "15", 16 : "16", 17 : "17", 18 : "18", 19 : "19", 20 : "20"})
fig = sns.relplot(data = df, kind="line", hue = "unit_id", x = "iteration", y = "loss")
plt.show()
