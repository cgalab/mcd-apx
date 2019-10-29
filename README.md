# MCD - Minimum Convex Decomposition by Martin

run using Knauer \& Spillner approach:
<code>mcd_apx --random [--seed 1] --obj --input FILE_IN --output FILE_OUT</code>

run using onion-based approximation:
<code>mcd_apx --onion --obj --input FILE_IN --output FILE_OUT</code>

| options       | description   |
| -------------:|:------------- |
|  --onion        | using onion-based approximation   |
|  --input &lt;file&gt;       | input file    |
|  --output &lt;file&gt;      | output file    |
|  --seed        | seed for random generator    |
|  --counter NUM        | set a specific counter of NUM    |
|  --timeout NUM       | set a timeout of NUM    |
|  --verbose       | verbose output    |
|  --random       | use a randomized approach    |
|  --index NUM       | set an index   |
|  --obj       | write object output format  |
|  --partition NUM       | partition into NUM blocks (use NUM=2^x, default: 1) |


## Branches

partition -- partition points into grids, apply acd_apx on each cell, then merge the resulting convex partition -- by gue
