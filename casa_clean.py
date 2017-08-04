import argparse
import subprocess as sp
import numpy as np

def clean(path, rms, mask):
         
    dirty_path = path + '.dirty'
    clean_path = path + '.clean'
    sp.call("rm -rf {{},{}}.*".format(dirty_path, clean_path), shell=True)
    
    print('Cleaning...')
    
    try: 
        rms = np.float(rms)
    except ValueError:
        #dirty clean to get rms
        sp.call("rm -rf {}.*".format(imagename), shell=True)
        tclean(vis=path + '.vis',
               imagename=dirty_path,
               imsize=512,
               cell='0.03arcsec',
               weighting='natural',
               niter=0)
        rms = imstat(imagename='{}.image'.format(dirty_path),
            region=rms, listit=False)['rms'][0]

    # if view:
    #     #Show dirty image, then clean up and delete all dirty clean files
    #     print 'dirty natural: {}'.format(dirty_natural_rms)
    #     viewer(infile=image + '.residual', displaytype='contour')
    #     raw_input('mask ready? ')
    sp.call("rm -rf {}.*".format(dirty_path), shell=True)

    # Clean with correct mask and rms
    tclean(vis=path + '.vis',
           imagename=clean_path,
           imsize=512,
           cell='0.03arcsec',
           weighting='natural',
           niter=100000000,
           threshold=rms / 2.,
           usemask='user',
           mask=mask,
           pbmask=None)

    # get rms
    clean_rms = imstat(imagename='{}.image'.format(clean_path),
        region=rms, listit=False)['rms'][0]
    print('Natural clean rms is {}'.format(natural_rms))
    
    # if view: viewer(infile=image + '.image')

    # Export to .fits
    exportfits(imagename='{}.image'.format(clean_path),
               fitsimage='{}.fits'.format(clean_path))
               
    return clean_rms

def concatenate(infiles, outfile):
    
    print('Concatenating measurement sets')
    print(np.array(infiles))
    print('into {}.ms'.format(outfile))
    
    sp.call("rm -rf {}.ms".format(outfile), shell=True)
    concat(vis=infiles, concatvis=outfile + '.ms', dirtol='2arcsec')
    return outfile
    


        

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Clean CASA measurement sets')
    parser.add_argument('vis', help='path(s) of visibilities to be cleaned')
    parser.add_argument('-concat', help='output path of concatenation')
    parser.add_argument('rms', help='either float or CASA region used to get rms with imstat')
    parser.add_argument('mask', help='clean mask file')
    args=parser.parse_args()
    
    vis_list = args.vis.split(',')
    if args.concat:
        concatenate(vis_list, args.concat)
        clean(args.conat, args.rms, args.mask)
    else:
        for vis in vis_list: clean(vis, rms, mask)
    
    # if len(vis_list) > 1:
    #     concat(vis_list,)
        
    # clean(vis_list, rms, mask)
