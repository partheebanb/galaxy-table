import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import warnings
import typer
import glob

from astropy.io import ascii

app = typer.Typer()

def fxn():
    warnings.warn("deprecated", DeprecationWarning)

# set size and style for plt
plt.style.use('seaborn')
plt.rcParams['figure.figsize'] = [10.0, 4.0]
plt.rcParams['figure.dpi'] = 140

# get pah features
pah = open(r'data/pah.txt', 'r')
pahFeatures = []
for line in pah:
    pahFeatures.append(float(line))

# get Hydrogen and forbidden features
features = ascii.read(r'data/Smith2007_Table2_H2_nebular_wavelengths.txt')
features = features.to_pandas()

featureLabels = features['line']

def applyScale(fluxes, wavelengths, scale, lowerBound, upperBound):
    '''
    Scale fluxes using the edge ratios found in the IRS spectra files' headers. The scale is applied to fluxes with wavelengths between lowerBound and upperBound.
    This is done to help eliminate problems such as massive jumps and drops in flux at the wavelengths where the data is stitched together.

    fluxes: the observed fluxes
    wavelengths: the observed wavelengths. The wavelengths correspond to the observed fluxes
    scale: the scale to be applied to the fluxes of observations with wavlengths between lowerBound and upperBound
    lowerBound: the shortest wavelength to apply the scaling to
    upperBound: the longest wavelength to apply the scaling to

    returns: scaled fluxes
    '''
    for i in range(len(fluxes)):
        
        if wavelengths[i] >= lowerBound and wavelengths[i] <= upperBound:
            fluxes[i] = fluxes[i] * scale
        
    return fluxes

@app.command()
def plot(
    input_dir: str = typer.Argument('input', help='Folder to read Spitzer data from.'),
    output_dir:str = typer.Argument('output', help='Directory to write outputs to'),
    h:bool = typer.Option(True, help='Include vertical lines for Hydrogen emission spectra? (will be plotted as dotted red lines)'),
    f:bool = typer.Option(True, help='Include vertical lines forbidden emission spectra? (will be plotted as dotted red lines)'),
    pah:bool = typer.Option(True, help='Include vertical lines for PAH emission spectra? (will be plotted as dotted black lines)'),
    flags:bool = typer.Option(True, help='Include markers for data points with bit-flags? (will be plotted as green triangles)'),
    scale:bool = typer.Option(True, help='Scale the fluxes according to edge ratios of the observations specified in the header of IRS spectra file?')
    ):
    '''
    Plots Spitzer IRS spectra and annotates emission lines and bit-flags. Also scales segments based on their edge ratios if required
    '''
    # TODO add redshifts?
    plot(input_dir, output_dir, h, f, pah, flags, scale)


def plot(input_dir, output_dir, h, f, pah, flags, scale):
    try:
        os.mkdir(output_dir)
    except:
        typer.echo(f'A dir with the name {typer.style(output_dir, bold=True, fg=typer.colors.RED)} already exists')

    for subdir, dirs, files in os.walk(fr'./{input_dir}'):
        for filename in files:
            # try: 
                filepath = subdir + os.sep + filename
                redshift = 0

                # spizterMap = open(r'C:\Users\bpart\OneDrive\Desktop\astro\MASSIVE\data\data-maps\spitzer.txt', 'r')

                data = ascii.read(filepath, format='ipac')
                npData = data.to_pandas().to_numpy().transpose()

                for name, keyword in data.meta['keywords'].items():
                    if name == "OBJECT":
                        objName = keyword['value'].split("'")[1]
                        print(objName)
                    elif name == "RA":
                        ra = keyword['value'].split(" ")[0]
                    elif name == "DEC":
                        dec = keyword['value'].split(" ")[0]
                    elif name == "RATIO_SL1_LL2":
                        ll2 = float(keyword['value'].split(" ")[0])
                    elif name == "RATIO_LL2_LL1":
                        ll1 = float(keyword['value'].split(" ")[0])
                    elif name == "RATIO_SL2_SL1":
                        sl1 = float(keyword['value'].split(" ")[0])

                map = False
                if map:
                    # TODO
                    spizterMap = open('map.txt', 'r')
                    for line in spizterMap:
                        if objName in line:
                            publishedName = line.split(",")[0]
                            break
                else:
                    publishedName = filename.split('/')[-1].split('.')[0]

                # TODO get redshift

                npFeatures = features.to_numpy().transpose()
                redShiftedFeatures = npFeatures[1] + npFeatures[1]*redshift

                redShiftedPahFeatures = []
                for i in range(len(pahFeatures)):
                    redShiftedPahFeatures.append(pahFeatures[i] + pahFeatures[i] * redshift)


                wavelengths = npData[1]
                fluxes = npData[2]
                errors = npData[3]
                flags = npData[4]

                if sl1 > 0:
                    fluxes = applyScale(fluxes, wavelengths, sl1, 7.56, 14.28)
                    ll2 *= sl1

                if ll2 > 0:
                    fluxes = applyScale(fluxes,  wavelengths, ll2, 14.28, 20.66)
                    ll1 *= ll2

                if ll1 > 0:
                    fluxes = applyScale(fluxes,  wavelengths, ll1, 20.68, 100)

                

                plt.step(wavelengths, fluxes, color="#1F77B4")
                plt.errorbar(wavelengths, fluxes, yerr=errors, fmt="none", ecolor="#1F77B4", elinewidth=0.5, alpha=0.7)
                
                plt.xlabel('Wavelength (microns)')
                plt.ylabel('Flux density (Jy)')
                plt.title(f'{publishedName}, z = {redshift}')

                # if redshift > 0:
                for i in range(len(redShiftedFeatures)):
                    if redShiftedFeatures[i] > max(wavelengths):
                        break
                    elif redShiftedFeatures[i] >= min(wavelengths):
                        plt.axvline(x=redShiftedFeatures[i], ymax=1, ymin=0, color='red', linestyle='dotted', alpha=0.3, label=1)
                        x_bounds = plt.xlim()

                        # code to place Ar II and Fe II labels on bottom to avoid overlap
                        if featureLabels[i] == '[Ar II]' or featureLabels[i] == '[Fe II]':
                            plt.annotate(text=featureLabels[i], xy =(((redShiftedFeatures[i]-x_bounds[0])/(x_bounds[1]-x_bounds[0]))+ 0.015,0.2), xycoords='axes fraction', verticalalignment='baseline', horizontalalignment='right' , rotation = 90, fontsize='x-small')
                        else:
                            plt.annotate(text=featureLabels[i], xy =(((redShiftedFeatures[i]-x_bounds[0])/(x_bounds[1]-x_bounds[0])),0.8), xycoords='axes fraction', verticalalignment='baseline', horizontalalignment='right' , rotation = 90, fontsize='x-small')

                
                for p in redShiftedPahFeatures:
                    if p > max(wavelengths):
                        break
                    elif p > min(wavelengths):
                        plt.axvline(x=p, ymax=0.2, ymin=0.05, color='black',linestyle='dotted', alpha=0.3, label=1)
                

                for i in range(len(flags)):
                    if flags[i] > 0:
                        plt.plot(wavelengths[i], fluxes[i] + max(fluxes) /20, color='green', marker=7)


                plt.annotate(text=f'SL1 (7.56-14.28): {sl1 if sl1 > 0 else 0}\nLL2 (14.28-20.66): {round(ll2, 4) if ll2 > 0 else 0} \nLL1 (14.28-): {round(ll1, 4) if ll1 > 0 else 0}', xy=(0.97, 1), xycoords='axes fraction', fontsize='x-small')

                plotName = fr"{output_dir}/{publishedName}.png"
                plt.savefig(plotName)
                plt.close()
            # except:
            #     pass
    
if __name__ == '__main__':
    app()