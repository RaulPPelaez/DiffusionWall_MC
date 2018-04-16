


fit formula "y=a0*exp(-a1*x)"
fit with 2 parameters
fit prec 1e-6


a0 = 500
a0 constraints on
a0min = 0.0
a0max = 100

a1 = 0.0001
a1 constraints on
a1min = 0
a1max = 2.0

nonlfit(s0, 500)
nonlfit(s0, 500)
nonlfit(s0, 500)
nonlfit(s0, 500)

