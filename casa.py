import subprocess as sp
import numpy as np

def pipe(commands):
    call_string = '\n'.join([command if type(command) is str else '\n'.join(command) for command in commands])
    
    print('Piping the following commands to CASA:\n')
    print(call_string)
    sp.call(['casa', '-c', call_string])
    
    # clean up .log files that casa poops out
    sp.call('rm -rf *.log', shell=True)
            
def concat(infiles, outfile):
    
    sp.call("rm -rf {}".format(outfile), shell=True)
    pipe(("concat(",
        "vis={},".format(infiles), 
        "concatvis='{}', dirtol='2arcsec')".format(outfile)))
    
def obs_clean(filename, rms, mask, 
    input_dir = '', output_dir = '', extension='.ms', 
    datacolumn='corrected', weighting='natural', uvtaper=None, clean_up=True):
         
    suffix = uvtaper[0] if uvtaper is not None else weighting
    dirty_filename = filename + '.' + suffix + '_dirty'
    clean_filename = filename + '.' + suffix + '_clean'
    
    commands = ()
    # get rms from region if necessary
    if type(rms) is str: 
        commands += \
        ("tclean(",
            "vis        = '{}',".format(input_dir + filename + extension),
            "datacolumn = '{}',".format(datacolumn),
            "imagename  = '{}',".format(output_dir + dirty_filename),
            "imsize     = 512,",
            "cell       = '0.03arcsec',",
            "weighting  = '{}',".format(weighting),
            "uvtaper    = {},".format(uvtaper),
            "niter      = 0)"),
        commands += \
        ("rms = imstat("
            "imagename='{}.image',".format(output_dir + dirty_filename),
            "region='{}', listit=False)['rms'][0]".format(rms))
        commands += \
        ("print('Dirty rms is {}'.format(rms))"),
    else:
        commands += ("rms = {}".format(rms),)
        
    # actually call tclean
    commands += \
    ("tclean(",
        "vis        = '{}',".format(input_dir + filename + extension),
        "datacolumn = '{}',".format(datacolumn),
        "imagename  = '{}',".format(output_dir + clean_filename),
        "imsize     = 512,",
        "cell       = '0.03arcsec',",
        "weighting  = '{}',".format(weighting),
        "uvtaper    = {},".format(uvtaper),
        "niter      = 100000000,",
        "threshold  = rms/2.,",
        "mask       = '{}')".format(mask)),
        
    # report clean rms if necessary
    if type(rms) is str: 
        commands += \
        ("clean_rms    = imstat(",
            "imagename = '{}.image',".format(output_dir + clean_filename),
            "region    = '{}', listit=False)['rms'][0]".format(rms)),
        commands += \
        ("print('Clean rms is {}'.format(clean_rms))"),
        
    # export to fits
    commands += \
    ("exportfits(",
        "imagename = '{}.image',".format(output_dir + clean_filename),
        "fitsimage = '{}.fits')".format(output_dir + clean_filename))
    
    
    sp.call('rm -rf {}*{{.fits,.image,.mask,.model,.pb,.psf,.residual,.sumwt}}'.format(output_dir + filename + '.' + suffix), shell=True)
    print('')
    print('=========================================================')
    print('=========================================================')
    print('Making {}...'.format(clean_filename))
    print('=========================================================')
    print('=========================================================')
    print('')
    pipe(commands)
        
    #clean up dirty files (hehe)
    if clean_up:
        sp.call('rm -rf {}*{{dirty.image,.mask,.model,.pb,.psf,.residual,.sumwt}}'.format(output_dir + filename + '.' + suffix), shell=True)
    
    # if view:
    #     #Show dirty image, then clean up and delete all dirty clean files
    #     print 'dirty natural: {}'.format(dirty_natural_rms)
    #     viewer(infile=image + '.residual', displaytype='contour')
    #     raw_input('mask ready? ')

    # # if view: viewer(infile=image + '.image')
    # 
    # sp.call("rm -rf {}.*".format(dirty_path), shell=True)
    # return clean_rms

def model_clean(path, rms, mask):
         
    dirty_path = path + '.dirty'
    clean_path = path + '.clean'
    sp.call('rm -rf {}*{{.image,.mask,.model,.pb,.psf,.residual,.sumwt}}'.format(path), shell=True)
    
    print('Cleaning...')
    pipe(
        ("tclean(",
            "vis='{}.ms',".format(path),
            "imagename='{}',".format(clean_path),
            "imsize=512,",
            "cell='0.03arcsec',",
            "weighting='natural',",
            "niter=100000000,",
            "threshold={}/2.,".format(rms),
            "mask='{}')".format(mask)),
        ("exportfits(",
            "imagename='{}.image',".format(clean_path),
            "fitsimage='{}.fits')".format(clean_path)))
    #clean up dirty files (hehe)
    sp.call("rm -rf {}.*".format(dirty_path), shell=True)

def to_fits(path):
    sp.call('rm -rf {}.fits'.format(path), shell=True)
    pipe("exportfits(",
        "imagename='{}.clean.image',".format(path),
        "fitsimage='{}.fits')".format(path))
               
def to_ms(paths):
    if type(paths) is str: paths = [paths]
    for path in paths: sp.call('rm -rf {}.ms'.format(path), shell=True)
    
    pipe(
    "for path in {}:".format(paths),
    ("    importuvfits(fitsfile = path+'.uvf', vis = path+'.ms')"))


        
# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description='Clean CASA measurement sets')
#     parser.add_argument('vis', help='path(s) of visibilities to be cleaned')
#     parser.add_argument('-concat', help='output path of concatenation')
#     parser.add_argument('rms', help='either float or CASA region used to get rms with imstat')
#     parser.add_argument('mask', help='clean mask file')
#     args=parser.parse_args()
#     
#     vis_list = args.vis.split(',')
#     if args.concat:
#         concatenate(vis_list, args.concat)
#         clean(args.concat, args.rms, args.mask)
#     else:
#         for vis in vis_list: clean(vis, rms, mask)
    
    # if len(vis_list) > 1:
    #     concat(vis_list,)
    #     
    # clean(vis_list, rms, mask)
