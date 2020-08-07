##################################################################################################
# For use in python v. 3.7
#
# Originally created by Joshua I. Brown in August, 2020
#
# Description: Testing a Split with Migration model
#
#
# NOTE: This script depends on frequency spectrum file created by fs_from_nex_GBS.py
##################################################################################################



import dadi
import numpy
import sys
from numpy import array
import pylab
from dadi import Demographics2D


# Load the data
data = dadi.Spectrum.from_file('Dove.WWD.MD.spectra.projected.17.17.fs')
#Call number of samples (ns)
ns = data.sample_sizes
data.mask_corners()
# print number of samples to verify correct load
print (ns)
# Grid point settings will be used for extrapolation.
# Grid points need to be formated [n,n+10,n+20]. Needs to be bigger #than the number of samples you have (n>ns) and this will be a strong #determination as to how long your program will run
pts_l = [30,40,50]
# call particular model to run, the model choosen here is split.w.migration
func = dadi.Demographics2D.split_mig
#    params = (nu1,nu2,T,m)
#    ns = (n1,n2)
#
#    Split into two populations of specifed size, with migration.
#    nu1: Size of population 1 after split.
#    nu2: Size of population 2 after split.
#    T: Time in the past of split (in units of 2*Na generations) 
#    m: Migration rate between populations (2*Na*m)
#    n1,n2: Sample sizes of resulting Spectrum
#    pts: Number of grid points to use in integration.

# define your initial Params
params = array([2.64866114, 9.94528158, 5.49139634, 0.0381305])

#Set the upper and lower bounds to make sure that the bounderies are #there. Suggested time parameters: lower 0, upper 5, migration #parameters: lower 0, upper 10,size parameters: lower 1e-2, upper 100 
upper_bound = [15, 15, 20, 20]
lower_bound = [1e-3, 1e-3, 1e-3, 1e-3]


reps = 1
while reps > 0:
    func_ex = dadi.Numerics.make_extrap_func(func)


    p0 = dadi.Misc.perturb_params(params, fold=2, lower_bound=lower_bound,upper_bound=upper_bound)
    popt = dadi.Inference.optimize_log(p0, data, func_ex,
                                       pts_l,lower_bound=lower_bound,
                                       upper_bound=upper_bound,verbose= 1)
    model = func_ex(popt, ns, pts_l)
    theta = dadi.Inference.optimal_sfs_scaling(model, data)
    ll_opt = dadi.Inference.ll_multinom(model, data)

    outfile = "Dove.WWD.MD.Sp.M"
    pylab.figure(figsize = (8,6))
    dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=0.01, resid_range=1, 
                                        pop_ids =('WWD.SpM', 'MD.SpM'), show = False)
    pylab.savefig(outfile + str(reps) + ".pdf", dpi = 300, bbox_inches = "tight")
    pylab.savefig(outfile + str(reps) + ".png", dpi = 300, bbox_inches = "tight")

    print ('Optimized parameters', repr([theta,ll_opt, popt])) 
    file = open("Dove.WWD.MD.Sp.M.opt.theta.txt", "a")
    file.write("\nReplicate " + str(reps) + " " + repr([theta,ll_opt,popt]))
    file.close()
    reps = reps - 1
