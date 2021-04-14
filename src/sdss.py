import typer
import pandas as pd
import aplpy
import os
import numpy as np
import warnings

from astroquery.sdss import SDSS
from astroquery.ned import Ned
import astropy.units as u
from astropy import coordinates as coords
from astropy.table import Table, Column
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

app = typer.Typer()

def fxn():
    warnings.warn("deprecated", DeprecationWarning)

def cleanup():
    '''
    A function to cleanup all FITS files downloaded and created in the process of generating images.
    '''
    os.remove('r.fits')
    os.remove('g.fits')
    os.remove('b.fits')
    os.remove('image.fits')
    os.remove('image_2d.fits')

@app.command()
def inputHandler(input_file: str = typer.Argument('input.csv', help='csv/txt file containing names and coordinates of targets (in degrees)'),
        output_dir:str = typer.Argument('output', help='Folder to write outputs to'),
        filters:str = typer.Argument('i,g,u', help='Choose three filters from z,i,r,g,u in the format "1,2,3". The order of the inputs decides the RGB bands of the output. Example: "i,g,u", where i is the r band, g is the g band, and u is the b band.'),
        output_format:str = typer.Argument('png', help='Choose your output format from any format supported by PIL'),
        pmin_r:float = typer.Option(0.5, help='pmin for r band'),
        pmax_r:float = typer.Option(99.5, help='pmax for r band'),
        pmin_g:float = typer.Option(0.5, help='pmin for g band)'),
        pmax_g:float = typer.Option(99.5, help='pmax for g band'),
        pmin_b:float = typer.Option(0.5, help='pmin for b band'),
        pmax_b:float = typer.Option(99.5, help='pmax for b band')):

    """
        A Python tool for plotting RGB images using data from SDSS. Ensure your input file is a csv/txt file in the format 'name, ra, dec'.
        The tool queries SDSS using the provided coordinates and names the output file according to the names provided.

        The wavelength of each available band is provided below (in angstroms):

            Ultraviolet (u):     3543

            Green (g):           4770

            Red (r):             6231

            Near Infrared (i):   7625

            Infrared (z):        9134

        You can choose any output format supported by PIL and set pmin and pmax for each band anywhere between 0 and 100. You can also save the FITS files for each image.

    """
    sdss(input_file, output_dir, filters, output_format, pmin_r,pmax_r, pmin_g, pmax_g, pmin_b, pmax_b)

def sdss(input_file,
        output_dir,
        filters,
        output_format,
        pmin_r,
        pmax_r,
        pmin_g,
        pmax_g,
        pmin_b,
        pmax_b):
    """
    Generates RGB images using SDSS data and APLPY. Information for parameters can be found in inputHandler()
    """
    
    # typer.secho(f'{exclude_bf}', fg=typer.colors.MAGENTA)

    warnings.filterwarnings("ignore")

    save_fits = False

    # get parameters and inputs into a suitable format
    df = pd.read_csv(input_file, header=None)
    output_format = output_format.lower()
    filters = filters.split(',')

    if len(filters) != 3:
        raise Exception('Error. Please provide 3 filters separated by commas')

    try:
        os.mkdir(output_dir)
    except:
        typer.echo(f'A dir with the name {typer.style(output_dir, bold=True, fg=typer.colors.RED)} already exists')

    for i, line in df.iterrows():
        try:

            name = line[0]
            ra = line[1]
            dec = line[2]

            styled_name = typer.style(name, fg=typer.colors.MAGENTA, bold=True)

            typer.echo(f'Plotting: ' + styled_name)

            pos = coords.SkyCoord(ra, dec, unit='deg')

            # query SDSS for each target and get the images
            try:
                xid = Table(SDSS.query_region(pos, spectro=False)[0])
            except:
                raise Exception(f'No images found on SDSS for target {name}')

            im = SDSS.get_images(matches=xid, band=filters)

            # raise exception if no images are found
            if len(im) == 0:
                raise Exception(f'No images found on SDSS for target {name}')

            # Obtain the PrimaryHDU from the HDUList for each band
            r, g, b = im[0][0], im[1][0], im[2][0]


            # save the fits files so they can be combined into a single rgb cube
            if not save_fits:
                r.writeto('r.fits', overwrite=True)
                g.writeto('g.fits', overwrite=True)
                b.writeto('b.fits', overwrite=True)
                aplpy.make_rgb_cube(['r.fits', 'g.fits', 'b.fits'], f'image.fits', north=True)
            else:
                pass
            
            aplpy.make_rgb_image(f'image.fits', fr'{output_dir}/{name}.{output_format}',  pmin_r=pmin_r, pmax_r=pmax_r, pmin_g=pmin_g, pmax_g=pmax_g, pmin_b=pmin_b, pmax_b=pmax_b)
            image = aplpy.FITSFigure(fr'{output_dir}/{name}.{output_format}')
            image.show_rgb()

            # add labels for filters used
            image.add_label(0.08, 0.15, "FILTERS:", relative=True, color='white')    

            filter_wavelengths = {'z': 913, 'i': 763, 'r': 623, 'g': 477, 'u': 354}
            bandColours = ['red', 'green', 'cyan']
            for i in range(3):
                w = filter_wavelengths[filters[i]]
                c = bandColours[i]
                image.add_label(0.08, 0.11 - 0.03*i, f"{filters[i]} ({w}nm)", relative=True, color=c)

            # add arrows pointing to object
            image.show_arrows(x=[ra + 0.01, ra - 0.01, ra + 0.01, ra - 0.01], 
                            y=[dec + 0.01, dec - 0.01, dec - 0.01, dec + 0.01], 
                            dx=[-0.007, 0.007, -0.007, 0.007], 
                            dy=[-0.007, 0.007, 0.007, -0.007], 
                            color='red')
            
            # add object name as label
            image.add_label(ra, dec + 0.02, name, color='red')

            # get redshift for scalebar
            resTable = Ned.query_region(pos, radius=0.01 * u.deg)
            resTable = resTable.to_pandas()
            redshift = 0
            for i, row in resTable.iterrows():
                if row['Redshift'] > 0:
                    redshift = row['Redshift']
                    break

            kpcArcmin = cosmo.kpc_proper_per_arcmin(float(redshift)).value
            length = round(kpcArcmin*2,1)  
            image.add_scalebar(0.03, label=f'120" | {length}kpc | r = {redshift}', color='red')

            image.save(fr'{output_dir}/{name}.{output_format}')
            
            typer.echo('Finished plotting: ' + styled_name)


        except Exception as e:
            typer.echo('Failed for: ' + styled_name)
            print(e)

    cleanup()

if __name__ == '__main__':
    app()