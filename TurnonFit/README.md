# Run the turnon fits 

`root -l -b -q turnon_chi2.C(${cat}\,${lower_bound}\,${higher_bound})`

# Check statistical info

`root -l -b -q getchi2.C(${cat}\,${lower_bound}\,${higher_bound}) # get raw chi2`
`root -l -b -q getprob.C(${cat}\,${lower_bound}\,${higher_bound})# convert to probability`