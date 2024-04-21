import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

counts_table_path = '/Users/stephaniecrilly/test/count_table.csv'

counts_table = pd.read_csv(counts_table_path, sep='\t')
print(counts_table['r1-hs-12-1'].sum())

#divide every column by the sum of the column
for col in counts_table.columns:
    print(col)
    print(counts_table[col].sum())
    if col == 'design':
        continue
    else:
        counts_table[col] = counts_table[col] / counts_table[col].sum() 

#select only control counts
list_of_ctrl_designs = ['bm01_ALFA_1c', 'bm01_ALFA_t2', 'bm01_ALFA_t3', 'bm01_ALFA_t6']
counts_table = counts_table[counts_table['design'].isin(list_of_ctrl_designs)]

#reorder the columns
cols = counts_table.columns.tolist()
#reorder list alphabetically
cols = sorted(cols)
counts_table = counts_table[cols]

#write the new table to a file
counts_table.to_csv('/Users/stephaniecrilly/test/normalized_counts_table_fraction_total_counts.tsv', sep='\t', index=False)

counts_table = pd.melt(counts_table, id_vars=['design'], value_vars=counts_table.columns[1:])
counts_table[['library_id', 'bin']] = counts_table['variable'].str.split('-hs-', expand=True)

sns.barplot(data=counts_table, x='bin', y='value', hue='design', palette=['cyan', 'red', 'orange', 'green'], order=['SD', '1-1', '1-2', '1-3', '1-4', '2-1', '2-2', '2-3', '2-4', '3-1', '3-2', '3-3', '3-4', '4-1', '4-2', '4-3', '4-4',
                                                                        '5-1', '5-2', '5-3', '5-4', '6-1', '6-2', '6-3', '6-4', '7-1', '7-2', '7-3', '7-4', '8-1', '8-2', '8-3', '8-4', '9-1', '9-2', '9-3', '10-1', '10-2', '10-3', '11-1', '11-2', '12-1'])
plt.show()
plt.clf()

fig, axes = plt.subplots(4, 1, figsize=(15, 5), sharex=True)
fig.suptitle('Normalized ngs counts control constructs')
fig.supxlabel('Bin')

sns.barplot(data=counts_table[counts_table['design'] == 'bm01_ALFA_1c'], ax=axes[0], x='bin', y='value', palette=['cyan'], order=['SD', '1-1', '1-2', '1-3', '1-4', '2-1', '2-2', '2-3', '2-4', '3-1', '3-2', '3-3', '3-4', '4-1', '4-2', '4-3', '4-4',
                        '5-1', '5-2', '5-3', '5-4', '6-1', '6-2', '6-3', '6-4', '7-1', '7-2', '7-3', '7-4', '8-1', '8-2', '8-3', '8-4', '9-1', '9-2', '9-3', '10-1', '10-2', '10-3', '11-1', '11-2', '12-1']).set(xlabel=None)
axes[0].set_title('bm01_ALFA_1c')

sns.barplot(data=counts_table[counts_table['design'] == 'bm01_ALFA_t2'], ax=axes[1], x='bin', y='value', palette=['red'], order=['SD', '1-1', '1-2', '1-3', '1-4', '2-1', '2-2', '2-3', '2-4', '3-1', '3-2', '3-3', '3-4', '4-1', '4-2', '4-3', '4-4',
                        '5-1', '5-2', '5-3', '5-4', '6-1', '6-2', '6-3', '6-4', '7-1', '7-2', '7-3', '7-4', '8-1', '8-2', '8-3', '8-4', '9-1', '9-2', '9-3', '10-1', '10-2', '10-3', '11-1', '11-2', '12-1']).set(xlabel=None)
axes[1].set_title('bm01_ALFA_t2')

sns.barplot(data=counts_table[counts_table['design'] == 'bm01_ALFA_t3'], ax=axes[2], x='bin', y='value', palette=['orange'], order=['SD', '1-1', '1-2', '1-3', '1-4', '2-1', '2-2', '2-3', '2-4', '3-1', '3-2', '3-3', '3-4', '4-1', '4-2', '4-3', '4-4',
                        '5-1', '5-2', '5-3', '5-4', '6-1', '6-2', '6-3', '6-4', '7-1', '7-2', '7-3', '7-4', '8-1', '8-2', '8-3', '8-4', '9-1', '9-2', '9-3', '10-1', '10-2', '10-3', '11-1', '11-2', '12-1']).set(xlabel=None)
axes[2].set_title('bm01_ALFA_t3')

sns.barplot(data=counts_table[counts_table['design'] == 'bm01_ALFA_t6'], ax=axes[3], x='bin', y='value', palette=['green'], order=['SD', '1-1', '1-2', '1-3', '1-4', '2-1', '2-2', '2-3', '2-4', '3-1', '3-2', '3-3', '3-4', '4-1', '4-2', '4-3', '4-4',
                        '5-1', '5-2', '5-3', '5-4', '6-1', '6-2', '6-3', '6-4', '7-1', '7-2', '7-3', '7-4', '8-1', '8-2', '8-3', '8-4', '9-1', '9-2', '9-3', '10-1', '10-2', '10-3', '11-1', '11-2', '12-1']).set(xlabel=None)
axes[3].set_title('bm01_ALFA_t6')

plt.show()




