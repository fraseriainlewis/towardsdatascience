# Code for blog posts on towardsdatascience.com

## April 2021: "An alternative to Machine Learning: Bayesian diffusion modeling"
- See folder **BayesDiffOne**
- The code file is **blog_GompertzSDE.R** which when run should fit a model and generate a plot of the results.  
- The misc folder contains other files used to help in the blog but not essential

```bash
cd towardsdatascience/BayesDiffOne
R CMD BATCH --vanilla blog_GompertzSDE.R log.R # note as coded this takes ~ 1hr over 8 cpus
# or source blog_GompertzSDE.R in R Studio 
# after setting working directory to towardsdatascience/BayesDiffOne
```
