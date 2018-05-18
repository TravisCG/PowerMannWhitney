# PowerMannWhitney

A large number of life scientist cannot understand statistics properly and they use
only non-parametric tests for eveything. Mann-Whitney test is one of them.

PowerMannWhitney is a highly optimized Mann-Whitney calculator for large data. It is using Bonferroni correction,
because I would like to avoid sorting the results.

## Usage
powermw onegroup|onevalue rowid path_to_groupfile path_to_valuefile outputfile

### onegroup

Choose one line from the group matrix and use this to perform Mann-Whitney on all the
rows in the value matrix.

### onevalue

Choose one line from the value matrix and use it to perform Mann-Whitney with every group
settings in the group matrix.

### rowid

The first column in both matrix is a rowid. This parameter helps to select one row from the
specified matrix.

Warning: The software does not check if there are
multiple rows with the same ID!

### path_to_groupfile / path_to_valuefile

Tab separated matrix files. The groupfile is a binary matrix. This specify the two groups
for the Mann-Whitney test.

## Alternatives

There is an alternative version (mutarget_core.c), which is used in www.mutarget.com. It has some additional
filtering options. I do not recommend the usage of this version in a standalone environment.

## Compiling

Simply type:

'''gcc powermw.c -lm -O3 -o powermw'''
