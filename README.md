# MAXIM-R
MAXIM-R is a collection of R programs intended to be used with the book *Modeling, Analysis, Design and Control of Stochastic Systems* by V G Kulkarni. The original collection is written in MATLAB, the translation is done by TJ Guo. This is a guide intended for STOR 445 students at UNC-CH on how to use MAXIM-R. 

# Prerequisite:
Download R and RStudio at the following links:

[Download R](https://cloud.r-project.org/)

[Download RStudio](https://www.rstudio.com/products/rstudio/download/#download)

After installing the programs, open RStudio and type in the following commands in the console to download devtools, a required package for installing this package.
```
install.packages("devtools")
library("devtools")
```
After loading devtools, type in the following command in the console to download this package:
```
install_github("cothurn/MAXIM-R")
```
Now you are ready to use this package!

# Using the package:
First, brush up on general R syntax:
[One of the many R tutortials out there](https://www.statmethods.net/r-tutorial/index.html)

To use the functions in the package, just type in the console/R file the function you want to call, followed by the arguments of the function in parenthesis. Example:
```R
> bincdf(10,0.5)
[1] 0.0009765625 0.0107421875 0.0546875000 0.1718750000 0.3769531250 0.6230468750 0.8281250000 0.9453125000
[9] 0.9892578125 0.9990234375 1.0000000000
```
You can use variables as arguments of the function/store the result of a function call as a variable if the result is not a plot. To see a complete list of available functions, please check the following section. You can also type in a question mark followed by the name of the function in RStudio console for a more detailed documentation and examples of its usage:

```R
?bincdf
```
## List of functions
### [probability](https://github.com/cothurn/MAXIM-R/blob/master/List%20of%20functions/Probability.md)
### [DTMC](https://github.com/cothurn/MAXIM-R/blob/master/List%20of%20functions/DTMC.md)
### [CTMC]() WIP
### [General models and queuing models]() WIP
