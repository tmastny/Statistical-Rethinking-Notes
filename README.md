# Statistical-Rethinking-Notes

This is my collection of notes on [*Statistical Rethinking: A Bayesian Course with Examples in R and Stan*](http://xcelab.net/rm/statistical-rethinking/). These are offered without assurances, besides the guarantee of mistakes and typos.

### How to use

There is a folder for each chapter and in each folder there are files for notes and exercises. There are two types of files, `.Rmd` and `.md`. To see the results in formatted html open the `.md` files on github. Some code or warnings may be suppressed. For the complete code, open the `.Rmd` files in RStudio.

### Other sets of notes

Self-study is difficult without feedback. I used the following sets of notes to check my understanding of the material.

https://github.com/cavaunpeu/statistical-rethinking

https://github.com/rlabbe/statistical-rethinking

http://scs.math.yorku.ca/index.php/SCS_Reads_2016-2017

### Why my notes are different

The book follows a traditional `r-base` style along with some added utility and Stan functions from the `rethinking` package. The other notes mostly conform to this style.  

As my notes progress I explore new packages while trying to verify the results using the expected methods. Following along with the book is a good way for me to build confidence and experience with new tools. 

Here are some of the packages that I incorporate:

1. [tidyverse](https://www.tidyverse.org/packages/): This is a collection of packages that attempt to make a consistent data science workflow. In particular, I explore

    a. [ggplot2](http://ggplot2.tidyverse.org/index.html) for visualization.

    b. [dplyr](http://dplyr.tidyverse.org/) for [tidying](https://cran.r-project.org/web/packages/tidyr/vignettes/tidy-data.html) data frames.

    c. [modelr](https://github.com/tidyverse/modelr) for `data_grid` to easily generate new data for models. 

2. [tidybayes](https://github.com/mjskay/tidybayes) to extend the tidy philosophy to Bayesian tools.

3. [brms](https://github.com/paul-buerkner/brms) as an alternative interface to Stan. This package allows you to construst full Bayesian models in Stan using the standard formula syntax, while offering some additional functionality to `rethinking::map2stan`. The package comes with a suite of helper functions that automate many of the common processes in the book, such as building the link function and calculating percentiles. For a comparsion to the rethinking package, check out the [brms overview](https://cran.r-project.org/web/packages/brms/vignettes/brms_overview.pdf). 


### Contact and Discussion

If you have any questions or suggestions you can reach me anytime at tim.mastny@gmail.com, or you can leave an issue on this repo. 

For more general *Rethinking* questions [Lotty Brand](https://twitter.com/LottyBrand22) recently started a *Rethinking* chat room on [slack](https://rethinkingstatistics.slack.com/). The author Richard McElreath has also joined the forum. 


