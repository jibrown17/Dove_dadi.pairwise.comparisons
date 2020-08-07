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
'DOV122b':'WWD',
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

    dd = nex_mo.data_dict_from_file('WWD.MD.nex', pop_assignments)
    
    # Generate frequency spectrum from data dictionary

    fs = dadi.Spectrum.from_data_dict(dd, ['WWD', 'MD'], [17,17], polarized=False)
    print (fs)
    # Print fs to file
    fs.tofile('Dove.WWD.MD.spectra.projection.xx.xx.fs')


    # Plot the fs
    import pylab
    pylab.figure(figsize = (8,6))
    dadi.Plotting.plot_single_2d_sfs(fs, vmin=0.01)
    pylab.savefig("Dove.WWD.MD.sepctra.projection.xx.xx.pdf", bbox_inches = 'tight')
    pylab.savefig("Dove.WWD.MD.sepctra.projection.xx.xx.png", dpi = 300, bbox_inches = 'tight')



