import pandas as pd

#Read excel file into a dataframe
df = pd.read_excel('Table.xlsx', 'Sheet.1', index_col=None)

#Write dataframe into csv
df.to_csv('table.tsv', sep='\t', encoding='utf-8',  index=False, line_terminator='\n')
