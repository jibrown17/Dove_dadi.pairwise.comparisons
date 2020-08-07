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
'DOV009a':'MD',
'DOV009b':'MD',
'DOV010a':'MD',
'DOV010b':'MD',
'DOV011a':'MD',
'DOV011b':'MD',
'DOV012a':'MD',
'DOV012b':'MD',
'DOV013a':'MD',
'DOV013b':'MD',
'DOV014a':'MD',
'DOV014b':'MD',
'DOV015a':'MD',
'DOV015b':'MD',
'DOV016a':'MD',
'DOV016b':'MD',
'DOV017a':'MD',
'DOV017b':'MD',
'DOV019a':'MD',
'DOV019b':'MD',
'DOV020a':'MD',
'DOV020b':'MD',
'DOV021a':'MD',
'DOV021b':'MD',
'DOV022a':'MD',
'DOV022b':'MD',
'DOV023a':'MD',
'DOV023b':'MD',
'DOV024a':'MD',
'DOV024b':'MD',
'DOV070a':'MD',
'DOV070b':'MD',
'DOV071a':'MD',
'DOV071b':'MD',
'DOV018a':'MD',
'DOV018b':'MD',
'DOV072a':'MD',
'DOV072b':'MD',
'DOV090a':'MD',
'DOV090b':'MD',
'DOV108a':'MD',
'DOV108b':'MD',
'DOV073a':'MD',
'DOV073b':'MD',
'DOV074a':'MD',
'DOV074b':'MD',
'DOV075a':'MD',
'DOV075b':'MD',
'DOV076a':'MD',
'DOV076b':'MD',
'DOV077a':'MD',
'DOV077b':'MD',
'DOV078a':'MD',
'DOV078b':'MD',
'DOV079a':'MD',
'DOV079b':'MD',
'DOV080a':'MD',
'DOV080b':'MD',
'DOV081a':'MD',
'DOV081b':'MD',
'DOV082a':'MD',
'DOV082b':'MD',
'DOV083a':'MD',
'DOV083b':'MD',
'DOV084a':'MD',
'DOV084b':'MD',
'DOV085a':'MD',
'DOV085b':'MD',
'DOV086a':'MD',
'DOV086b':'MD',
'DOV087a':'MD',
'DOV087b':'MD',
'DOV088a':'MD',
'DOV088b':'MD',
'DOV089a':'MD',
'DOV089b':'MD',
'DOV091a':'MD',
'DOV091b':'MD',
'DOV092a':'MD',
'DOV092b':'MD',
'DOV093a':'MD',
'DOV093b':'MD',
'DOV094a':'MD',
'DOV094b':'MD',
'DOV095a':'MD',
'DOV095b':'MD',
'DOV099a':'MD',
'DOV099b':'MD',
'DOV100a':'MD',
'DOV100b':'MD',
'DOV101a':'MD',
'DOV101b':'MD',
'DOV102a':'MD',
'DOV102b':'MD',
'DOV103a':'MD',
'DOV103b':'MD',
'DOV104a':'MD',
'DOV104b':'MD',
'DOV105a':'MD',
'DOV105b':'MD',
'DOV106a':'MD',
'DOV106b':'MD',
'DOV107a':'MD',
'DOV107b':'MD',
'DOV109a':'MD',
'DOV109b':'MD',
'DOV110a':'MD',
'DOV110b':'MD',
'DOV111a':'MD',
'DOV111b':'MD',
'DOV112a':'MD',
'DOV112b':'MD',
'DOV113a':'MD',
'DOV113b':'MD',
'DOV151a':'MD',
'DOV151b':'MD',
'DOV152a':'MD',
'DOV152b':'MD',
'DOV153a':'MD',
'DOV153b':'MD',
'DOV154a':'MD',
'DOV154b':'MD'
}


    # Generate data dictionary from NeXus alignment

    dd = nex_mo.data_dict_from_file('Euro.MD.nex', pop_assignments)

    # Generate frequency spectrum from data dictionary

    fs = dadi.Spectrum.from_data_dict(dd, ['Euro', 'MD'], [20,15], polarized = False)
    print (fs)
    # Print fs to file
    fs.tofile('Dove.Euro.MD.spectra.projection.xx.xx.fs')


    # Plot the fs
    import pylab
    pylab.figure(figsize = (8,6))
    dadi.Plotting.plot_single_2d_sfs(fs, vmin=0.01)
    pylab.savefig("Dove.Euro.MD.sepctra.projection.xx.xx.pdf", bbox_inches = 'tight')
    pylab.savefig("Dove.Euro.MD.sepctra.projection.xx.xx.png", dpi = 300, bbox_inches = 'tight')



