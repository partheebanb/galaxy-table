import typer
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import warnings
from pyvo.dal import tap
import csv
from PyAstronomy import pyasl
import datetime as dt
import glob
import matplotlib as mpl


from astropy.io import fits
from astropy.table import Table
from astropy.io import ascii
from astropy.utils.data import get_pkg_data_filename
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.modeling.fitting import LevMarLSQFitter

from pahfit.base import PAHFITBase

from astroquery.ned import Ned
from astroquery.alma import Alma
from astroquery.eso import Eso
from astroquery.vizier import Vizier
from astroquery.gemini import Observations

import plotSpectra
import visibility
import sdss
import pahfitHelper

app = typer.Typer()

def fxn():
    warnings.warn("deprecated", DeprecationWarning)

def getImagePath(dir, file):
    '''
    A helper to get image files from dir. Identifies the format present and returns an HTML img tag
    dir: directory to search for file in
    file: name of file to search for (without format)
    '''
    try:
        files = glob.glob(f'{dir}/{file}.*')
        file_format = files[0].split('.')[-1]
    except:
        file_format = ''
    return f'<img src="{dir}/{file}.{file_format}" alt="No Plot" height="300">'

def getVisibility(input_file):
    '''
    Plots visibility plots for each target in the input file at the given site during the given year. 
    Returns the directory containing the plots.

    input_file: is a csv/txt file in the format 'name,ra,dec', where ra and dec are ICRS J2000 coordinates in degrees.
    '''

    dataFile = open(input_file, 'r')
    dataFile = csv.reader(dataFile)

    typer.secho('\nProvide the following information to create the visibility plots\n', fg=typer.colors.GREEN)
    output_dir = typer.prompt('Enter the name of the directory to which the plots will be saved')
    pyasl.listObservatories()
    site = typer.prompt('Enter the code of the observatory from the list above')
    yearOrDay = typer.prompt('Is this plot for a whole year or a single day? (year/day default:year)').lower()

    if yearOrDay == 'day':
        date = typer.prompt('Enter the date (yyyy-mm-dd)')
    else:
        date = int(typer.prompt('Enter the year (yyyy)'))

    # convert the date parameter into Python datetime format 
    if yearOrDay == 'day':
        ymd = [int(i) for i in date.split('-')]
        date = dt.datetime(*ymd)

    # checks if output directory exists and creates it if it doesn't
    if os.path.isdir(output_dir):
        typer.echo('Directory named ' + typer.style(output_dir, fg=typer.colors.RED) + ' already exists')
    else:
        os.mkdir(output_dir)

    for line in dataFile:

        # prepare the data for visibility.py
        targets = []
        name = line[0]
        ra = line[1]
        dec = line[2]

        coord = SkyCoord(ra, dec, unit='deg')

        targets.append({'name': name,
                        'coord': coord})

        # get the plot
        if yearOrDay == 'day':
            fig = visibility.VisibilityPlot(date=date, targets=targets, observatory=site)
        else:
            fig = visibility.StarObsPlot(year=date, targets=targets, observatory=site, 
                        period=None, hover=False, sunless_hours=None,
                        remove_watermark=True)

        typer.echo('Plotting completed for '+ typer.style(name, fg=typer.colors.MAGENTA, bold=True))
        fig.savefig(fr'{output_dir}/{name.strip()}.png')
        plt.close(fig)
    return output_dir



def getSDSS_RGB(input_file):
    """
    Create RGB images using data from SDSS
    """

    typer.secho('\nProvide the requested information to generate RGB images using SDSS', fg=typer.colors.GREEN)

    output_dir = typer.prompt('Enter the name of the directory to which the images will be saved')
    filters = typer.prompt('Choose three filters from z,i,r,g,u in the format "1,2,3". The order of the inputs decides the RGB bands of the output. Example: "i,g,u", where i is the r band, g is the g band, and u is the b band.')
    output_format= typer.prompt('Choose either png or jpg as output format').lower()
    change_p = (typer.prompt('Would you like to change the pmin and pmax values? Defaults are pmin=0.5 and pmax=99.5 (y/n | Default: n)') == 'y')

    if change_p:
        pmin_r = float(typer.prompt('Enter pmin value for the r band'))
        pmax_r = float(typer.prompt('Enter pmax value for the r band'))
        pmin_g = float(typer.prompt('Enter pmin value for the g band'))
        pmax_g = float(typer.prompt('Enter pmax value for the g band'))
        pmin_b = float(typer.prompt('Enter pmin value for the b band'))
        pmax_b = float(typer.prompt('Enter pmax value for the b band'))
    else:
        pmin_r = pmin_g = pmin_b = 0.5
        pmax_r = pmax_g = pmax_b = 99.5 

    sdss.sdss(input_file, output_dir, filters, output_format, pmin_r, pmax_r, pmin_g, pmax_g, pmin_b, pmax_b)

    return output_dir

def getPahfit():

    typer.secho('\nProvide the requested information to perform PAHFits\n', fg=typer.colors.GREEN)

    input_dir = typer.prompt('Enter the name of the directory to load the Spizter spectra from')
    output_dir = typer.prompt('Enter the name of the directory to which the outputs will be saved')
    max_iter = int(typer.prompt('Enter the maximum number of iterations to be performed per PAHFit. This number should be an integer'))
    packfile = (typer.prompt('Enter the name of the packfile to be used for fitting'))
    error = float(typer.prompt('Enter the relative error to be used for early termination of PAHFit'))

    pahfitHelper.pahfit(input_dir, output_dir, max_iter, error, packfile, False)

    return output_dir, output_dir

def getSpectra():

    typer.secho('\nProvide the requested information to plot the Spizter spectra\n', fg=typer.colors.GREEN)

    input_dir = typer.prompt('Enter the name of the directory to load Spitzer spectra from')
    output_dir = typer.prompt('Enter the name of the directory to which the outputs will be saved')
    h = typer.prompt('Include vertical lines for Hydrogen emission spectra? (will be plotted as dotted red lines)'),
    f = typer.prompt('Include vertical lines forbidden emission spectra? (will be plotted as dotted red lines)'),
    pah = typer.prompt('Include vertical lines for PAH emission spectra? (will be plotted as dotted black lines)'),
    flags = typer.prompt('Include markers for data points with bit-flags? (will be plotted as green triangles)'),
    scale = typer.prompt('Scale the fluxes according to edge ratios of the observations specified in the header of IRS spectra file?')

    # 
    plotSpectra.plot(input_dir, output_dir, h, f, pah, flags, scale)

    return output_dir

    return

@app.command()
def buildTable():
    """
    Allows users to create a fully customizable HTML table containing information about a list of input galaxies.
    """
    input_file = typer.prompt('Enter the path to the csv/txt file containing the targets (format "target,ra,dec"), where ra, dec are the ICRS coordinates of target in degrees in J2000 equinox')
    output_name = typer.prompt('Enter the name of the output table')

    #  get the basic columns
    columns = ['NAME', 'RA', 'DEC', 'REDSHIFT', 'REFERENCE FOR REDSHIFT', 'NED LINK']
    
    # add the columns for images/plots/other files
    add_spectra = (typer.prompt('Add spectra plots to the table?(y/n)|(default: n)').lower() == 'y')
    create_spectra = False
    if add_spectra:
        create_spectra = (typer.prompt('Create spectra plots now? Type "n" if you already have them (y/n)|(default: n)').lower() == 'y')
        if not create_spectra:
            spectra_dir = typer.prompt('Enter directory with spectra plot files. File names must be in format (target).png/jpg')
        columns.append('SPECTRA')

    add_pahfitplot = (typer.prompt('Add pahfit plots and tables to the table?(y/n)|(default: n)').lower() == 'y')
    create_pahfit = False
    if add_pahfitplot:
        create_pahfit = (typer.prompt('Create pahfit plots and tables now? Type "n" if you already have them (y/n)|(default: n)').lower() == 'y')
        if not create_pahfit:
            pahfitplot_dir = typer.prompt('Enter directory with pahfit plot files. File names must be in format (target).png/jpg')
            pahfittable_dir = typer.prompt('Enter directory with pahfit output table files. File names must be in format (target)_output.ipac or (target).ipac')
        columns.extend(['PAHFIT PLOT', 'PAHFIT TABLE'])

    add_rgb = (typer.prompt('Add RGB images to the table(y/n)|(default: n)').lower() == 'y')
    create_rgb = False
    if add_rgb:
        create_rgb = (typer.prompt('Create rgb images now using SDSS? Type "n" if you already have the images(y/n)|(default: n)').lower() == 'y')
        if not create_rgb:
            rgb_dir = typer.prompt('Enter directory with the RGB images. File names musst be in format (target).png/jpg')
        columns.append('RGB IMAGE')

    add_visibility = (typer.prompt('Add visibility plots to the table?(y/n)|(default: n)').lower() == 'y')
    create_visibility = False
    if add_visibility:
        create_visibility = (typer.prompt('Create visibility plots now? Type "n" if you already have the plots(y/n)|(default: n)').lower() == 'y')
        if not create_visibility:
            visibility_dir = typer.prompt('Enter directory with the visibility plots. File names must be in format (target).png/jpg:')
        columns.append('VISIBILITY PLOT')

    # add other images and files the use might have
    other_image_dirs = []
    add_other_images = True
    while add_other_images:
        add_other_images = (typer.prompt('Add any other images/plots to the table?(y/n)|(default: n)').lower() == 'y')
        if add_other_images:
            other_image_dirs.append(typer.prompt('Enter directory with the other images. File names must be in format (target).png/jpg:'))
            columns.append(typer.prompt('Enter the name of the column for these images/plots').upper())

    other_files = []
    add_other_files = True
    while add_other_files:
        add_other_files = (typer.prompt('Add links to any other files?(y/n)|(default: n)').lower() == 'y')
        if add_other_files:
            format = typer.prompt('Enter the format of these files')
            dir = typer.prompt('Enter directory with the other images. File names must be in format (target).(format):')
            other_files.append({'format': format, 'dir': dir})

    # add all optional non-file columns

    add_eso = (typer.prompt('Add ESO links to the table?(y/n)|(default: n)').lower() == 'y')
    add_alma = (typer.prompt('Add ALMA links s to the table?(y/n)|(default: n)').lower() == 'y')
    add_manga = (typer.prompt('Add MANGA links to the table?(y/n)|(default: n)').lower() == 'y')
    add_gemini = (typer.prompt('Add GEMINI links to the table?(y/n)|(default: n)').lower() == 'y')
    add_guide_stars = (typer.prompt('Add guide star info to the table?(y/n)|(default: n)').lower() == 'y')


    if add_eso:
        columns.extend(['ESO INSTRUMENTS', 'ESO LINKS'])
    if add_alma:
        columns.append('ALMA LINK')
    if add_manga:
        columns.extend(['MaNGA ID', 'MaNGA PLATE IFUDesign', 'MaNGA LINK'])
    if add_gemini:
        columns.extend(['GEMINI INSTRUMENTS', 'GEMINI LINK'])
    if add_guide_stars:
        columns.extend(['HIGH STREHL', 'LOW STREHL'])

    # set up eso_tap_obs
    if add_eso:
        ESO_TAP_OBS = "http://archive.eso.org/tap_obs"
        tapobs = tap.TAPService(ESO_TAP_OBS)

    if add_manga:
        drpall = fits.open(r'C:\Users\bpart\OneDrive\Desktop\astro\MASSIVE\data\manga-drpall.fits')

    # set up Vizier so all the necessary columns can be obtained for the guide stars, including flags that signify galaxies
    if add_guide_stars:
        v = Vizier(columns=['RAJ2000', 'DEJ2000', 'pmRA', 'pmDE', 'UCAC4', 'LEDA', 'rmag', 'f_rmag', 'g', 'c', 'of'])



    # create spizter spectra plots
    if create_spectra:
        spectra_dir = getSpectra()


    # create pahfit plots and tables
    if create_pahfit:
        pahfitplot_dir, pahfittable_dir = getPahfit()

    # create rgb images
    if create_rgb:
        rgb_dir = getSDSS_RGB(input_file)

    # create visibility plots if user chose that option
    if create_visibility:
        visibility_dir = getVisibility(input_file)





    dataFile = open(input_file, 'r')
    dataFile = csv.reader(dataFile)

    output_name = 'index' if output_name == '' else output_name

    with open(f'{output_name}.csv', 'w') as outfile:
        csv_writer = csv.writer(outfile)
        csv_writer.writerow(columns)

        for line in dataFile:
            name = line[0]
            ra = line[1]
            dec = line[2]

            typer.secho('Building row for ' + typer.style(name, fg=typer.colors.MAGENTA))

            # get redshift from NED
            coord = SkyCoord(ra, dec, unit='deg')
            resTable = Ned.query_region(coord, radius=0.016 * u.deg)
            resTable = resTable.to_pandas()
            redshift = ''
            refCode, reflink = '', ''
            for i, row in resTable.iterrows():
                if row['Redshift'] > 0:
                    nName = row['Object Name']
                    redshift = row['Redshift']

                    # have to use name to query NED's redshift table to get reference
                    try:
                        redTable = Ned.get_table(nName, table='redshifts').to_pandas()

                        for j, r in redTable.iterrows():
                            if r['Published Redshift'] == redshift:
                                refCode = r['Refcode']
                                reflink = f'https://ui.adsabs.harvard.edu/abs/{refCode}/abstract'
                    except:
                        refCode = "No refcode found"

                    break;
            
            # get nedLink using coordinates
            nedLink = f'https://ned.ipac.caltech.edu/conesearch?search_type=Near%20Position%20Search&coordinates={ra}d%20{dec}d&radius=1' 

            # prepare row for required columns
            row = [name, ra, dec, redshift, f"<a href='{reflink}'>{refCode}</>", f"<a href='{nedLink}'>NED Link</>"]

            # prepare rows for other columns

            # Add spectra plot
            if add_spectra:
                row.append(getImagePath(spectra_dir, name))
                
            # Add PAHFIT plot and PAHFIT table, accounting for possible variations in naming conventions
            if add_pahfitplot:
                row.append(getImagePath(f'{pahfitplot_dir}', name))
                if os.path.exists(f'{pahfittable_dir}/{name}.ipac'):
                    row.append(f"<a href='{pahfittable_dir}/{name}.ipac'>PAHFIT Table</>")
                else:
                    row.append(f"<a href='{pahfittable_dir}/{name}_output.ipac'>PAHFIT Table</>")
            
            # Add visibility plot
            if add_visibility:
                row.append(getImagePath(visibility_dir, name))

            # Add RGB plots from SDSS
            if add_rgb:
                row.append(getImagePath(rgb_dir, name))

            # Add other plots/images
            for d in other_image_dirs:
                row.append(getImagePath(d, name))

            # ESO
            # tried using astroquery.eso but it would not work. I believe this is the best way to do it as of 04/2021
            if add_eso:
                
                query = f"""SELECT *
                FROM ivoa.ObsCore
                WHERE intersects(s_region, circle('', {ra}, {dec}, 0.02))=1""" 

                res = tapobs.search(query=query, maxrec=1000)

                # get all unique instruments
                esoInstruments = set(res['instrument_name'])
                esoInstruments = list(esoInstruments)
                esoInstruments = " ".join(i for i in esoInstruments).replace(",", "-")
                
                esoLink, esoText = '', ''

                # get the ESO link using coordinates
                if len(esoInstruments) > 1:
                    esoLink = f'https://archive.eso.org/scienceportal/home?pos={ra}%2C{dec}&r=0.02'
                    esoText = 'ESO Link'

                row.append(esoInstruments)
                row.append(f"<a href='{esoLink}'>{esoText}</>")
            
            # ALMA
            # there is a chance of astroquery.alma not working. Try installing the latest version of astroquery directly from the GitHub source in this case
            if add_alma:

                almaLink = ''
                almaText = ''
                try:
                    res = Alma.query_region(coord, 0.02*u.deg)
                    if len(res) > 0:
                        almaLink = f'https://almascience.eso.org/asax/?result_view=observation&raDec={ra}%20{dec}'
                        almaText = "ALMA"            
                except:
                    pass
                row.append(f"<a href='{almaLink}'>{almaText}</>")
            
            # MaNGA
            # the drpall file version 2_4_3 is required to obtain MaNGA data from release 16
            if add_manga:
                tbdata = drpall[1].data
                i = np.where(np.sqrt((tbdata['objra'] - float(ra))**2 + (tbdata['objdec'] - float(dec))**2) <= 0.02)
                mid, plate, ifudsgn, mlink, mtext = "", "", "", "", ""

                if len(i[0]) > 0:
                    # the mid, plate, and ifudsgn are used to locate manga data
                    # TODO provide more info on this
                    mid = tbdata['mangaid'][i][0]
                    plate = tbdata['plate'][i][0]
                    ifudsgn = tbdata['ifudsgn'][i][0]
                    mlink = f'https://data.sdss.org/sas/dr16/manga/spectro/redux/v2_4_3/{plate}/stack/'
                    mtext = "MaNGA Link"

                row.append(mid)
                row.append(ifudsgn)
                row.append(f'<a href="{mlink}">{mtext}</>')

            # GEMINI
            if add_gemini:
                gData = Observations.query_region(coordinates=SkyCoord(ra, dec, unit='deg'), radius=0.02*u.deg)
                gInstruments = set()
                gText, gLink = '', ''
                if len(gData) > 0:
                    gInstruments = set(gData['instrument'])

                    gText = "Gemini Link"
                    gLink = f'https://archive.gemini.edu/searchform/NotFail/ra={ra}/not_site_monitoring/notengineering/dec={dec}/cols=CTOWEQ'

                gInstruments = list(gInstruments)
                gInstruments = " ".join(i for i in gInstruments).replace(",", "-")

                row.append(gInstruments)
                row.append(f'<a href="{gLink}">{gText}</>')

            # Get guide stars using Vizier from the I/322A catalogue
            # Refer here for more information https://www.gemini.edu/instrumentation/altair
            if add_guide_stars:
                ucac4 = v.query_region(SkyCoord(ra, dec, unit='deg'), radius=25*u.arcsec, catalog='I/322A')
                highList, lowList = [], []
                if ucac4:
                    ucac4 = ucac4[0]
                    for star in ucac4:
                        dist = abs((np.float32(star['RAJ2000']) - np.float32(ra))**2 + abs(np.float32(star['DEJ2000']) - np.float32(dec))**2)**0.5 * 3600
                        magn = star['rmag']

                        # Objects with LEDA > 0 are galaxies and will be disregarded. 
                        if (3 < dist < 15 and 8.5 < magn < 15) and star['LEDA'] == star['g'] == star['c'] == 0:  
                            highList.append((star['UCAC4'], round(dist, 2), round(magn, 2), star['of'], star['f_rmag']))
                        elif (15 < magn < 18.5) and (3 < dist) and star['LEDA'] == star['g'] == star['c'] == 0:
                            lowList.append((star['UCAC4'], round(dist, 2), round(magn, 2), star['of'], star['f_rmag']))

                # prepare outputs for High and Low Strehl columns using data obtained from catalogue
                highList = " ".join((f'{str(h[0])} (rmag={str(h[2])}; f_rmag={str(h[4])}; distance={str(h[1])}"; of={str(h[3])})') for h in highList)
                lowList = " ".join((f'{str(h[0])} (rmag={str(h[2])}; f_rmag={str(h[4])}; distance={str(h[1])}"; of={str(h[3])})') for h in lowList)
                row.append(highList)
                row.append(lowList)
            
            # avoid adding empty strings to output as this can cause problems when converting to HTML 
            for i in range(len(row)):
                if row[i] == "":
                    row[i] = " "
            csv_writer.writerow(row)

    # create HTMl table from the csv table using pandas. The resulting table will look pretty ugly and as such it is recommended to add some basic styling.
    data = pd.read_csv(f'{output_name}.csv')
    file = open(f'{output_name}.html', 'w')

    file.writelines(data.to_html(escape=False, index=False))

if __name__ == '__main__':
    app()