# for analyzing flow cytometry data
dependencies:
```
Python 3
numpy
pandas
seaborn
matplotlib
scipy 
FlowCytometryTools
```

Download the script into a folder, when you run make sure your target folder (the folder with fcs files) has a sampleName.csv with the following content:

for example:


| cells | MOI | Virus |
| :---- | :---- | :---- |
| HEK293T | 1000 | AAV2 |
| HEK293T | 5000 | AAV2 |

** make sure the all the characters are capitalized!


