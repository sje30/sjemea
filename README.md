# sjemea


Multielectrode array analysis in R

Please see the
[Gigascience](http://www.gigasciencejournal.com/content/3/1/3) paper
for recent details.

## Installing from github

## Citing
If you use this work in your publications, please cite it using the
output from the R function

```
citation("sjemea")
```

## TODO

### Portable C code

Update the C code tiling.c so that it avoids where possible any
R-specific headers/macros (e.g. R_NaN, Sfloat).  Then we can write a
mex function for calling it in matlab.

