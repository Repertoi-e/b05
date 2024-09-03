import pandas as pd
import matplotlib.pyplot as plt

url = "https://docs.google.com/spreadsheets/d/1HbuwOGN-pKfY95QqNByAydbRJVicJ3oztc3d0fufQT8/export?format=csv&gid=0#gid=0"
df = pd.read_csv(url)

print(df)