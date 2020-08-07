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
# NOTE: This script depends on a modified version of the custom script "NeXus.py" originally 
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
'DOV001a':'WWD',
'DOV001b':'WWD',
'DOV002a':'WWD',
'DOV002b':'WWD',
'DOV003a':'WWD',
'DOV003b':'WWD',
'DOV004a':'WWD',
'DOV004b':'WWD',
'DOV005a':'WWD',
'DOV005b':'WWD',
'DOV006a':'WWD',
'DOV006b':'WWD',
'DOV007a':'WWD',
'DOV007b':'WWD',
'DOV008a':'WWD',
'DOV008b':'WWD',
'DOV025a':'WWD',
'DOV025b':'WWD',
'DOV026a':'WWD',
'DOV026b':'WWD',
'DOV027a':'WWD',
'DOV027b':'WWD',
'DOV028a':'WWD',
'DOV028b':'WWD',
'DOV029a':'WWD',
'DOV029b':'WWD',
'DOV030a':'WWD',
'DOV030b':'WWD',
'DOV031a':'WWD',
'DOV031b':'WWD',
'DOV032a':'WWD',
'DOV032b':'WWD',
'DOV033a':'WWD',
'DOV033b':'WWD',
'DOV034a':'WWD',
'DOV034b':'WWD',
'DOV035a':'WWD',
'DOV035b':'WWD',
'DOV037a':'WWD',
'DOV037b':'WWD',
'DOV038a':'WWD',
'DOV038b':'WWD',
'DOV039a':'WWD',
'DOV039b':'WWD',
'DOV040a':'WWD',
'DOV040b':'WWD',
'DOV041a':'WWD',
'DOV041b':'WWD',
'DOV042a':'WWD',
'DOV042b':'WWD',
'DOV043a':'WWD',
'DOV043b':'WWD',
'DOV044a':'WWD',
'DOV044b':'WWD',
'DOV045a':'WWD',
'DOV045b':'WWD',
'DOV046a':'WWD',
'DOV046b':'WWD',
'DOV047a':'WWD',
'DOV047b':'WWD',
'DOV048a':'WWD',
'DOV048b':'WWD',
'DOV049a':'WWD',
'DOV049b':'WWD',
'DOV050a':'WWD',
'DOV050b':'WWD',
'DOV051a':'WWD',
'DOV051b':'WWD',
'DOV052a':'WWD',
'DOV052b':'WWD',
'DOV053a':'WWD',
'DOV053b':'WWD',
'DOV055a':'WWD',
'DOV055b':'WWD',
'DOV056a':'WWD',
'DOV056b':'WWD',
'DOV057a':'WWD',
'DOV057b':'WWD',
'DOV058a':'WWD',
'DOV058b':'WWD',
'DOV059a':'WWD',
'DOV059b':'WWD',
'DOV060a':'WWD',
'DOV060b':'WWD',
'DOV061a':'WWD',
'DOV061b':'WWD',
'DOV064a':'WWD',
'DOV064b':'WWD',
'DOV065a':'WWD',
'DOV065b':'WWD',
'DOV066a':'WWD',
'DOV066b':'WWD',
'DOV067a':'WWD',
'DOV067b':'WWD',
'DOV068a':'WWD',
'DOV068b':'WWD',
'DOV069a':'WWD',
'DOV069b':'WWD',
'DOV036a':'WWD',
'DOV036b':'WWD',
'DOV054a':'WWD',
'DOV054b':'WWD',
'DOV098a':'WWD',
'DOV098b':'WWD',
'DOV114a':'WWD',
'DOV114b':'WWD',
'DOV115a':'WWD',
'DOV115b':'WWD',
'DOV117a':'WWD',
'DOV117b':'WWD',
'DOV118a':'WWD',
'DOV118b':'WWD',
'DOV119a':'WWD',
'DOV119b':'WWD',
'DOV121a':'WWD',
'DOV121b':'WWD',
'DOV122a':'WWD',
'DOV122b':'WWD'
}


    # Generate data dictionary from NeXus alignment

    dd = nex_mo.data_dict_from_file('Euro.WWD.nex', pop_assignments)
    
    # Generate frequency spectrum from data dictionary

    fs = dadi.Spectrum.from_data_dict(dd, ['Euro', 'WWD'], [18,18], polarized=False)
    print (fs)
    # Print fs to file
    fs.tofile('Dove.Euro.WWD.spectra.projection.xx.xx.fs')


    # Plot the fs
    import pylab
    pylab.figure(figsize = (8,6))
    dadi.Plotting.plot_single_2d_sfs(fs, vmin=0.01)
    pylab.savefig("Dove.Euro.WWD.sepctra.projection.xx.xx.pdf", bbox_inches = 'tight')
    pylab.savefig("Dove.Euro.WWD.sepctra.projection.xx.xx.png", dpi = 300, bbox_inches = 'tight')



