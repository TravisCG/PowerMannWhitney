# PowerMannWhitney
Highly optimized Mann-Whitney calculator for large data. It is using Bonferroni correction,
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

There is an alternative version, which is used in www.mutarget.com. It has some additional
filtering option.
