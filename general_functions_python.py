#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: miyang
"""

# control the number of digits
i = 0.6342
float('%.3g' % i)

# read txt file
pd.read_csv(path, sep="\t")


# create result folder
if os.path.exists(result_folder):
    os.chdir(result_folder)
else:
    os.makedirs(result_folder) ; os.chdir(result_folder)
    
    
# sample n elements from a list
random.sample([1,2,2,3,4,5], k=2)
    
# python equivalent of table in R
df.sample_group.value_counts()
    
# insert a first column
df.insert(0, "subjectID", subjectID_NP )
# drop a column
df = df.drop(["sample"], axis=1)
# remove rows containing NA for a specific column
df.dropna(subset = ["column2"], inplace=True)

# selecting rows based on multiple column values in pandas dataframe

# sort dataframe
df.sort_values(by='column_2', ascending=True)


# remove first few n character of a string
a_string = "abcde"
sliced = a_string[2:]
# take last 2 character of a list of string
[sub[ -2: ] for sub in list_string ]    
# split and take first/second element from list
[i.split('_', 1)[1] for i in temp]


# cbind 2 dataframes
pd.concat([df_a.reset_index(drop=True), df_b.reset_index(drop=True)], axis=1)
# rbind 2 dataframe
df3 = pd.concat( [ df1 , df2 ] ).reset_index(drop=True)
# append a list as a row to a DataFrame
df.loc[len(df)] = list
# cbind array
np.c_[x,y]



# aggregate single column
grouped_single = df.groupby('Team').agg({'Age': ['mean', 'min', 'max']}) #  pd.to_numeric if a column is not numeric
# aggregate dataframe based on 2 columns
grouped_multiple = df.groupby(["col_1","col_2"]).mean().reset_index()

# convert from coordinate table to matrix 
df.pivot_table(index='x', columns='y', values='z')
df.pivot(index="x", columns="y", values="z")

# drop a duplicate row, based on column name
df = df.drop_duplicates(subset='favorite_color', keep="first")

# Intersection of the lists in a list
common = set.intersection(*map(set,my_list))
# difference between 2 lists

# flatten list
regular_list = [[1, 2, 3, 4], [5, 6, 7], [8, 9]]
flat_list = [item for sublist in regular_list for item in sublist]

# select first 2 elment of a list of list
list( features_to_plot_list[i] for i in [0,1] )

# create a dataframe based on 1D lists
v1 = [1,2,3]
v2 = [1,2,3]
v3 = [1,2,3]
df = pd.DataFrame( list(zip(v1,v2,v3)), columns=['v1','v2','v3'] )

df = pd.DataFrame( {"v1":v1, "v2":v2 }  )



