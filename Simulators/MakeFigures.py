import numpy as np
import matplotlib as plt
import seaborn as sns
import pandas as pd

M1XV1_Cycles = [20562, 1041, 10448, 2048]
M1XV1_Memory = [310001, 62000, 9452, 22001]

M1XV2_Cycles = [10314,553, 5228, 2048]
M1XV2_Memory = [15001, 32000, 4762, 11501]

M9XV1_Cycles = [2070, 137, 1616, 2048]
M9XV1_Memory = [4001, 8000, 2555, 4001]

M1XM1_Cycles = [26572, 11287, 81933]
M1XM1_Memory = [31001,255714, 15000]

M1XM2_Cycles = [23166, 6027, 40848]
M1XM2_Memory = [31001, 157010, 1170]

M2XM3_Cycles = [13288, 6797, 329331]
M2XM3_Memory = [16001,126852,15090]

Mbeacxc_Cycles = [17565, 74224, 2309422]
Mbeacxc_Memory = [101533,12178040,541350]

gnutella_Cycles = [6210, 103419, 16241244]
gnutella_Memory = [361218, 1793376, 31179525]

pesa_Cycles = [193425, 55105, 3550849]
pesa_Memory = [250437, 1438816, 2554424]



df = pd.DataFrame([M1XV1_Cycles, M1XV2_Cycles, M9XV1_Cycles], index=["M1XV1", "M2XV2", "M4XV1"])

print(df)

sns.barplot(df)


import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

sns.set_theme(style="ticks")

diamonds = sns.load_dataset("diamonds")

f, ax = plt.subplots(figsize=(7, 5))

index=["M1XV1", "M2XV2", "M4XV1"]

ax.bar(index, np.log10(np.array(gnutella_Cycles)), width=1, edgecolor="white", linewidth=0.7)

ax.set(xlim=(-2, 3), xticks=[-1,0,1],
       yticks=(np.arange(0,10)), yticklabels=(10**np.arange(0,10)))

f.savefig("test")