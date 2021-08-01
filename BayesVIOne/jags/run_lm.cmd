model in "rats_lm.bug"
data in "rats-data_lm.R"
compile, nchains(2)
#inits in "rats-init.R"
initialize
update 10000 
monitor alpha, thin(10) 
monitor beta, thin(10) 
monitor std.c, thin(10)
update 10000
coda *
