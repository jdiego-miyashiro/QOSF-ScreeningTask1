# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 15:58:11 2020

@author: castromi
"""

from matplotlib import pyplot as plt
import pandas as pd


epsilon_data = []
epsilon_data2=[]

doc_name ="results/results_for_rzrx_gates/epsilon_results"
doc_name2 ="results/results_for_ryrx_gates/epsilon_results_ryrx_L"
lent=20
for i in range(1,lent):
    experiment_result=pd.read_csv(doc_name+str(i)+".csv")
    experiment_result2=pd.read_csv(doc_name2+str(i)+".csv")
    epsilon_data.append(experiment_result["epsilon"].iloc[-1])
    epsilon_data2.append(experiment_result2["epsilon"].iloc[-1])
    pass

print(epsilon_data)
fig=plt.figure()
ax1=fig.add_subplot(111)
plt.grid()
plt.title("Epsilon evolution by number of layers ")
plt.xticks(range(1,20))
plt.xlabel("Number of Layers")
plt.ylabel("Epsilon (%)")
plt.ylim(bottom=min(epsilon_data+epsilon_data2)*0.98,top=max(epsilon_data+epsilon_data2)*1.02)
ax1.plot(range(1,lent),epsilon_data,"-o")
ax1.plot(range(1,lent),epsilon_data2,"-s")
ax1.legend(["Circuit with Rz & Rx gates","Circuit with Ry & Rx gates"])
print(epsilon_data2)


plt.show()
