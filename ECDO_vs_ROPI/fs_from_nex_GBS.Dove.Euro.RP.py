##################################################################################################
# File: fs_from_nex_GBS.py
# For use in python v. 3.7
#
# Created by Joshua I. Brown in August, 2020
# Based on script originally writen by Michael G. Harvey on 21 September, 2012
#
# Description: generate multi-population allele frequency spectra from SNPs from Nexus dataset
#
# Usage: python fs_from_nex_GBS.py
#
# NOTE: This script depends on a modified version of the custom script "nex_mo.py" originally 
# written by R. Gutenkunst specifically for use with alignments in nexus format.
##################################################################################################


import dadi
import nex_mo



if __name__ == '__main__':
    pop_assignments = {
'DOV128a':'Euro',
'DOV128b':'Euro',
'DOV146a':'Euro',
'DOV146b':'Euro',
'DOV096a':'Euro',
'DOV096b':'Euro',
'DOV097a':'Euro',
'DOV097b':'Euro',
'DOV123a':'Euro',
'DOV123b':'Euro',
'DOV124a':'Euro',
'DOV124b':'Euro',
'DOV125a':'Euro',
'DOV125b':'Euro',
'DOV126a':'Euro',
'DOV126b':'Euro',
'DOV127a':'Euro',
'DOV127b':'Euro',
'DOV129a':'Euro',
'DOV129b':'Euro',
'DOV130a':'Euro',
'DOV130b':'Euro',
'DOV131a':'Euro',
'DOV131b':'Euro',
'DOV132a':'Euro',
'DOV132b':'Euro',
'DOV133a':'Euro',
'DOV133b':'Euro',
'DOV134a':'Euro',
'DOV134b':'Euro',
'DOV135a':'Euro',
'DOV135b':'Euro',
'DOV136a':'Euro',
'DOV136b':'Euro',
'DOV137a':'Euro',
'DOV137b':'Euro',
'DOV138a':'Euro',
'DOV138b':'Euro',
'DOV139a':'Euro',
'DOV139b':'Euro',
'DOV140a':'Euro',
'DOV140b':'Euro',
'DOV141a':'Euro',
'DOV141b':'Euro',
'DOV142a':'Euro',
'DOV142b':'Euro',
'DOV143a':'Euro',
'DOV143b':'Euro',
'DOV144a':'Euro',
'DOV144b':'Euro',
'DOV145a':'Euro',
'DOV145b':'Euro',
'DOV147a':'Euro',
'DOV147b':'Euro',
'DOV148a':'Euro',
'DOV148b':'Euro',
'DOV149a':'Euro',
'DOV149b':'Euro',
'DOV164a':'RP',
'DOV164b':'RP',
'DOV150a':'RP',
'DOV150b':'RP',
'DOV155a':'RP',
'DOV155b':'RP',
'DOV156a':'RP',
'DOV156b':'RP',
'DOV157a':'RP',
'DOV157b':'RP',
'DOV158a':'RP',
'DOV158b':'RP',
'DOV159a':'RP',
'DOV159b':'RP',
'DOV160a':'RP',
'DOV160b':'RP',
'DOV161a':'RP',
'DOV161b':'RP',
'DOV162a':'RP',
'DOV162b':'RP',
'DOV163a':'RP',
'DOV163b':'RP',
'DOV165a':'RP',
'DOV165b':'RP',
'DOV166a':'RP',
'DOV166b':'RP',
'DOV167a':'RP',
'DOV167b':'RP',
'DOV168a':'RP',
'DOV168b':'RP',
'DOV169a':'RP',
'DOV169b':'RP',
'DOV170a':'RP',
'DOV170b':'RP',
'DOV171a':'RP',
'DOV171b':'RP',
'DOV172a':'RP',
'DOV172b':'RP',
'DOV173a':'RP',
'DOV173b':'RP',
'DOV174a':'RP',
'DOV174b':'RP',
'DOV175a':'RP',
'DOV175b':'RP',
'DOV176a':'RP',
'DOV176b':'RP',
'DOV177a':'RP',
'DOV177b':'RP',
'DOV178a':'RP',
'DOV178b':'RP',
'DOV179a':'RP',
'DOV179b':'RP',
'DOV180a':'RP',
'DOV180b':'RP',
'DOV181a':'RP',
'DOV181b':'RP',
'DOV182a':'RP',
'DOV182b':'RP',
'DOV183a':'RP',
'DOV183b':'RP',
'DOV184a':'RP',
'DOV184b':'RP',
'DOV185a':'RP',
'DOV185b':'RP',
'DOV186a':'RP',
'DOV186b':'RP'
}


    # Generate data dictionary from NeXus alignment

    dd = nex_mo.data_dict_from_file('Euro.RP.nex', pop_assignments)

    # Generate frequency spectrum from data dictionary

    fs = dadi.Spectrum.from_data_dict(dd, ['Euro', 'RP'], [23,23], polarized = False)
    print (fs)
    # Print fs to file
    fs.tofile('Dove.Euro.RP.spectra.projection.xx.xx.fs')


    # Plot the fs
    import pylab
    pylab.figure(figsize = (8,6))
    dadi.Plotting.plot_single_2d_sfs(fs, vmin=0.01)
    pylab.savefig("Dove.Euro.RP.sepctra.projection.xx.xx.pdf", bbox_inches = 'tight')
    pylab.savefig("Dove.Euro.RP.sepctra.projection.xx.xx.png", dpi = 300, bbox_inches = 'tight')



