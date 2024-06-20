'''
Data Extraction and Standardization
'''

import numpy as np

# This function excludes the first row of the nested list obtained with the readlines() method. 
# If there are gene names in the rows, it excludes them. If there are no gene names in the rows, 
# the else: statement comes into play, and it won't recalculate, but there won't be gene names in the table either. 
# We are using the version that contains gene names in the script later on.
def calculatorZ(exp_list): 
    if type(exp_list[0]) == str:
        z = [exp_list[0]]
        for i in range(1, len(exp_list)):
            x = (exp_list[i] - np.mean(exp_list[1:])) / np.std(exp_list[1:])
            z.append(x)
        return z
    else:
        z = []
        for i in range(len(exp_list)):
            x = (exp_list[i] - np.mean(exp_list)) / np.std(exp_list)
        return z

# I use "with" term to close the file automatically. It is a proper way of file processing.
with open("gene_expression.txt", "r") as file: 
    exp_values = file.readlines()

# Convert string type data obtained from the file into float type, reconstructing the list.
rows = []
for i in range(1, len(exp_values)):
    row = exp_values[i].replace("\n", "")
    row = row.split("\t")
    row[1:] = [float(j) for j in row[1:]]
    rows.append(row)

# Process the header row and convert it to a list format for later concatenation.
titles = exp_values[0].replace("\n", "")
titles = titles.split("\t")

# Using the calculatorZ function, assign the list containing gene names and z values to the variable table. 
# The first row of this variable contains the table headers.
table = [titles]
for i in rows:
    result = calculatorZ(i)
    table.append(result)

# Add necessary special characters to the data in list format to display the data in a table format in the .txt file.
table_with_special_characters = []
for i in table:
    for j in i:
        if j == i[-1]:
            x = str(j) + "\n"
        else:
            x = str(j) + "\t"
        table_with_special_characters.append(x)

# Convert the list-formatted data into a string type.
table_with_special_characters = "".join(table_with_special_characters)

# Write the data as tableZ.txt into a .txt file.
with open("tableZ.txt", "w") as file:
    file.write(table_with_special_characters)
