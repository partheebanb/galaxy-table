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

from astropy.io import fits
from astropy.table import Table
from astropy.io import ascii
from astropy.utils.data import get_pkg_data_filename
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM

from astroquery.ned import Ned
from astroquery.alma import Alma
from astroquery.eso import Eso
from astroquery.vizier import Vizier
from astroquery.gemini import Observations

import visibility
import sdss

app = typer.Typer()

def fxn():
    warnings.warn("deprecated", DeprecationWarning)

def getImagePath(dir, file):
    '''A helper to get image files from dir'''
    if os.path.exists(f'{dir}/{file}.png'):
        return f'<img src="{dir}/{file}.png" alt="No Plot" height="300">'
    else:
        return f'<img src="{dir}/{file}.jpg" alt="No Plot" height="300">'

def getVisibility(input_file):
    '''
    Plots visibility plots for each target in the input file at the given site during the given year. Ensure input file is in the format 'name,ra,dec'.
    Returns the directory containing the plots.
    '''

    dataFile = open(input_file, 'r')
    dataFile = csv.reader(dataFile)

    typer.echo('Please provide the following information to create the visibility plots')
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
    output_dir = typer.prompt('Enter the name of the directory to which the images will be saved')
    filters = typer.prompt('Choose three filters from z,i,r,g,u in the format "1,2,3". The order of the inputs decides the RGB bands of the output. Example: "i,g,u", where i is the r band, g is the g band, and u is the b band.')
    output_format= typer.prompt('Choose either png or jpg as output format').lower()
    pmin_r = float(typer.prompt('Enter pmin value for the r band'))
    pmax_r = float(typer.prompt('Enter pmax value for the r band'))
    pmin_g = float(typer.prompt('Enter pmin value for the g band'))
    pmax_g = float(typer.prompt('Enter pmax value for the g band'))
    pmin_b = float(typer.prompt('Enter pmin value for the b band'))
    pmax_b = float(typer.prompt('Enter pmax value for the b band'))

    sdss.sdss(input_file, output_dir, filters, output_format, pmin_r, pmax_r, pmin_g, pmax_g, pmin_b, pmax_b)

    return output_dir



@app.command()
def buildTable():
    """
    Allows users to create a fully customizable HTML table containing information about a list of input galaxies.
    """
    input_file = typer.prompt('Enter the path to the csv file containing the targets (format "target,ra,dec")')
    output_name = typer.prompt('Enter the name of the output table')

    #  get the basic columns
    columns = ['NAME', 'RA', 'DEC', 'REDSHIFT', 'REFERENCE FOR REDSHIFT', 'NED LINK']
    
    # add the columns for images/plots/other files
    add_spectra = (typer.prompt('Add spectra plots to the table?(y/n)').lower() == 'y')
    if add_spectra:
        spectra_dir = typer.prompt('Enter directory with spectra plot files. File names must be in format (target).png/jpg')
        columns.append('SPECTRA')

    add_pahfitplot = (typer.prompt('Add pahfit plots to the table?(y/n)').lower() == 'y')
    if add_pahfitplot:
        pahfitplot_dir = typer.prompt('Enter directory with pahfit plot files. File names must be in format (target).png/jpg')
        columns.append('PAHFIT PLOT')

    add_pahfittable = (typer.prompt('Add pahfit tables to the table?(y/n)').lower() == 'y')
    if add_pahfittable:
        pahfittable_dir = typer.prompt('Enter directory with pahfit output table files. File names must be in format (target)_output.ipac or (target).ipac')
        columns.append('PAHFIT TABLE')

    add_rgb = (typer.prompt('Add RGB images to the table(y/n)').lower() == 'y')
    create_rgb = False
    if add_rgb:
        create_rgb = (typer.prompt('Create rgb images now using SDSS? Type "n" if you already have the images(y/n)').lower() == 'y')
        if not create_rgb:
            rgb_dir = typer.prompt('Enter directory with the RGB images. File names musst be in format (target).png/jpg')
        columns.append('RGB IMAGE')

    add_visibility = (typer.prompt('Add visibility plots to the table?(y/n)').lower() == 'y')
    create_visibility = False
    if add_visibility:
        create_visibility = (typer.prompt('Create visibility plots now? Type "n" if you already have the plots(y/n)').lower() == 'y')
        if not create_visibility:
            visibility_dir = typer.prompt('Enter directory with the visibility plots. File names must be in format (target).png/jpg:')
        columns.append('VISIBILITY PLOT')

    # add other images and files the use might have
    other_image_dirs = []
    add_other_images = True
    while add_other_images:
        add_other_images = (typer.prompt('Add any other images/plots to the table?(y/n)').lower() == 'y')
        if add_other_images:
            other_image_dirs.append(typer.prompt('Enter directory with the other images. File names must be in format (target).png/jpg:'))
            columns.append(typer.prompt('Enter the name of the column for these images/plots').upper())

    other_files = []
    add_other_files = True
    while add_other_files:
        add_other_files = (typer.prompt('Add links to any other files?(y/n)').lower() == 'y')
        if add_other_files:
            format = typer.prompt('Enter the format of these files')
            dir = typer.prompt('Enter directory with the other images. File names must be in format (target).(format):')
            other_files.append({'format': format, 'dir': dir})

    # add all optional non-file columns

    add_eso = (typer.prompt('Add ESO links to the table?(y/n)').lower() == 'y')
    add_alma = (typer.prompt('Add ALMA links s to the table?(y/n)').lower() == 'y')
    add_manga = (typer.prompt('Add MANGA links to the table?(y/n)').lower() == 'y')
    add_gemini = (typer.prompt('Add GEMINI links to the table?(y/n)').lower() == 'y')
    add_guide_stars = (typer.prompt('Add guide star info to the table?(y/n)').lower() == 'y')

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

    # set up Vizier so all the necessary columns can be obtained
    if add_guide_stars:
        v = Vizier(columns=['RAJ2000', 'DEJ2000', 'pmRA', 'pmDE', 'UCAC4', 'LEDA', 'rmag', 'f_rmag', 'g', 'c', 'of'])

    # create visibility plots if use chose that option
    if create_visibility:
        visibility_dir = getVisibility(input_file)

    # create rgb images
    if create_rgb:
        rgb_dir = getSDSS_RGB(input_file)


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
            
            
            nedLink = f'https://ned.ipac.caltech.edu/conesearch?search_type=Near%20Position%20Search&coordinates={ra}d%20{dec}d&radius=1' 

            row = [name, ra, dec, redshift, f"<a href='{reflink}'>{refCode}</>", f"<a href='{nedLink}'>NED Link</>"]

            # Add spectra plot
            if add_spectra:
                row.append(getImagePath(spectra_dir, name))
                
            # Add PAHFIT plot
            if add_pahfitplot:
                row.append(getImagePath(pahfitplot_dir, name))

            # Add PAHFIT table, accounting for possible variations in naming conventions
            if add_pahfittable:
                if os.path.exists(f'{pahfittable_dir}/{name}.ipac'):
                    row.append(f"<a href='{pahfittable_dir}/{name}.ipac'>PAHFIT Table</>")
                else:
                    row.append(f"<a href='{pahfittable_dir}/{name}_output.ipac'>PAHFIT Table</>")
            
            # Add visibility plot
            if add_visibility:
                row.append(getImagePath(visibility_dir, name))

            if add_rgb:
                row.append(getImagePath(rgb_dir, name))

            # Add other plots/images
            for d in other_image_dirs:
                row.append(getImagePath(d, name))

            # ESO
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

                # get the ESO link
                if len(esoInstruments) > 1:
                    esoLink = f'https://archive.eso.org/scienceportal/home?pos={ra}%2C{dec}&r=0.02'
                    esoText = 'ESO Link'

                row.append(esoInstruments)
                row.append(f"<a href='{esoLink}'>{esoText}</>")
            
            # ALMA
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
            if add_manga:
                tbdata = drpall[1].data
                i = np.where(np.sqrt((tbdata['objra'] - float(ra))**2 + (tbdata['objdec'] - float(dec))**2) <= 0.02)
                mid, plate, ifudsgn, mlink, mtext = "", "", "", "", ""

                if len(i[0]) > 0:
                    # the mid, plate, and ifudsgn are used to locate manga data
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
                gData = Observations.query_region(coordinates=SkyCoord(ra, dec, unit='deg'), radius=0.016*u.deg)
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
            if add_guide_stars:
                ucac4 = v.query_region(SkyCoord(ra, dec, unit='deg'), radius=25*u.arcsec, catalog='I/322A')
                highList, lowList = [], []
                if ucac4:
                    ucac4 = ucac4[0]
                    for star in ucac4:
                        dist = abs((np.float32(star['RAJ2000']) - np.float32(ra))**2 + abs(np.float32(star['DEJ2000']) - np.float32(dec))**2)**0.5 * 3600
                        magn = star['rmag']

                        # Objects with LEDA > 0 are galaxies and will disregarded
                        if (3 < dist < 15 and 8.5 < magn < 15) and star['LEDA'] == star['g'] == star['c'] == 0:  
                            highList.append((star['UCAC4'], round(dist, 2), round(magn, 2), star['of'], star['f_rmag']))
                        elif (15 < magn < 18.5) and (3 < dist) and star['LEDA'] == star['g'] == star['c'] == 0:
                            lowList.append((star['UCAC4'], round(dist, 2), round(magn, 2), star['of'], star['f_rmag']))

                highList = " ".join((f'{str(h[0])} (rmag={str(h[2])}; f_rmag={str(h[4])}; distance={str(h[1])}"; of={str(h[3])})') for h in highList)
                lowList = " ".join((f'{str(h[0])} (rmag={str(h[2])}; f_rmag={str(h[4])}; distance={str(h[1])}"; of={str(h[3])})') for h in lowList)
                row.append(highList)
                row.append(lowList)
            
            for i in range(len(row)):
                if row[i] == "":
                    row[i] = " "
            csv_writer.writerow(row)

        # create HTMl table
    data = pd.read_csv(f'{output_name}.csv')
    file = open(f'{output_name}.html', 'w')

    file.writelines(data.to_html(escape=False, index=False))

if __name__ == '__main__':
    app()