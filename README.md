# sjemea


Multielectrode array analysis in R

Please see the
[Gigascience](http://www.gigasciencejournal.com/content/3/1/3) paper
for recent details.

# Installing this package

To install this package in R you need the devtools package; skip the
first line if you already have it installed:

	install.packages("devtools")
	devtools::install_github("sje30/sjemea")

Older versions of the package are available from my
[home page](http://www.damtp.cam.ac.uk/user/sje30/r).



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

