
import pandas as pd

try:
    df = pd.read_excel("/Users/michaelwolfe/Documents/nirsdexata/data-raw/participantdata.xlsx")
    print(df[['id', 'time']].head(10))
except Exception as e:
    print(e)
