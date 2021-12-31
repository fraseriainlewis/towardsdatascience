# Code for blog posts on towardsdatascience.com

## Jan 2022: "Modeling biological systems in space and time"
## Folder **PDEOne**
 

## August 2021: "Variational Bayesian inference with normalizing flows: a simple example"
## Folder **BayesVIOne**
- This contains a Jupyter notebook **rats_lm_vi_github** which is a self-contained workbook contain all modeling code. 
- The notebook was originally run in Google Colab and is the recommended way to run the code. Save to google drive and open directly in colab. It should also work in Jupyter but this has not been tested
- A sub-folder **jags** has the various small files needed to run the jags MCMC analysis. These results - in CODAchain1.txt are hard coded into the Jupyter notebook above. For different ways to install jags see https://mcmc-jags.sourceforge.io/  



## April 2021: "An alternative to Machine Learning: Bayesian diffusion modeling"
## Folder **BayesDiffOne**
- The code file is **blog_GompertzSDE.R** which when run should fit a model and generate a plot of the results.  
- The misc folder contains other files used to help in the blog but not essential

```bash
cd towardsdatascience/BayesDiffOne
R CMD BATCH --vanilla blog_GompertzSDE.R log.R # note as coded this takes ~ 1hr over 8 cpus
# or source blog_GompertzSDE.R in R Studio 
# after setting working directory to towardsdatascience/BayesDiffOne
```
