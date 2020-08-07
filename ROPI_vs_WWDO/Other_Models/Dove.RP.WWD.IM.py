##################################################################################################
# For use in python v. 3.7
#
# Originally created by Joshua I. Brown in August, 2020
#
# Description: Testing a Isolation with Migration model
#
#
# NOTE: This script depends on frequency spectrum file created by fs_from_nex_GBS.py
##################################################################################################



import dadi
import numpy
import sys
from numpy import array
from dadi import Demographics2D


# Load the data
data = dadi.Spectrum.from_file('Dove.RP.WWD.spectra.projected.30.10.fs')
#Call number of samples (ns)
ns = data.sample_sizes
data.mask_corners()
# print number of samples to verify correct load
print (ns)
# Grid point settings will be used for extrapolation.
# Grid points need to be formated [n,n+10,n+20]. Needs to be bigger #than the number of samples you have (n>ns) and this will be a strong #determination as to how long your program will run
pts_l = [200,220,240]
# call particular model to run, the model choosen here is isolation #with migration (IM)
func = dadi.Demographics2D.IM
# Parameters:   ns = (n1,n2) params = (s,nu1,nu2,T,m12,m21)
#    s: Size of pop 1 after split. (Pop 2 has size 1-s.)
#    nu1: Final size of pop 1.
#    nu2: Final size of pop 2.
#    T: Time in the past of split (in units of 2*Na generations) 
#    m12: Migration from pop 2 to pop 1 (2*Na*m12)
#    m21: Migration from pop 1 to pop 2
#    n1,n2: Sample sizes of resulting Spectrum
#    pts: Number of grid points to use in integration.

# define your initial Params
params = array([1.00149291, 1.00623777, 1.99381135, 12.88572863, 1.55169363, 4.48731493])

#Set the upper and lower bounds to make sure that the bounderies are 
#there. Suggested time parameters: lower 0, upper 5, migration 
#parameters: lower 0, upper 10,size parameters: lower 1e-2, upper 100 
upper_bound = [5,5,5,20,20,20]
lower_bound = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]


reps = 5
while reps > 0:
    func_ex = dadi.Numerics.make_extrap_log_func(func)
    p0 = dadi.Misc.perturb_params(params, fold=2, lower_bound=lower_bound,upper_bound=upper_bound)
    popt = dadi.Inference.optimize_log(p0, data, func_ex,
    pts_l,lower_bound=lower_bound,upper_bound=upper_bound,
    verbose= 1)

    model = func_ex(params, ns, pts_l)

    theta = dadi.Inference.optimal_sfs_scaling(model, data)
    ll_opt = dadi.Inference.ll_multinom(model, data)
        


    outfile = "RP.WWD.IM." + str(reps) + ".pdf" 
    pylab.figure(figsize = (8,6))
    dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=0.01, resid_range=1,
                                        pop_ids =('RP.IM', 'WWD.IM'), show = False)
    pylab.savefig(outfile, dpi = 300, bbox_inches = "tight")
        
    print ('Optimized parameters', repr([theta,ll_opt,popt])) 
    file = open("Dove.RP.WWD.IM.txt", "a")
    file.write("\nReplicate " + str(reps) + " " + repr([theta, ll_opt, popt]))
    file.close()
    reps = reps - 1
