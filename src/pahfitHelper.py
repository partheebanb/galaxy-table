import matplotlib.pyplot as plt
import matplotlib as mpl

import pandas as pd
import numpy as np
import os
import os.path
from os import path

from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
from astropy.modeling.fitting import LevMarLSQFitter

from pahfit.base import PAHFITBase

import typer

app = typer.Typer()

# chiSquaredFile = open('data/chi_squared.txt', 'w')

def calc_reduced_chi_square(fit, x, y, yerr, N, n_free):
    '''
    fit (array) values for the fit
    x,y,yerr (arrays) data
    N total number of points
    n_free number of parameters we are fitting
    '''
    return 1.0/(N - n_free)*sum(((fit - y)/yerr)**2)


@app.command()
def input_handler(inp: str = typer.Argument('input', help='Folder to read Spitzer data from.'), 
        output:str = typer.Argument('output', help='Folder to write outputs to'),
        max_iter: int = typer.Argument(300, help='Number of iterations to run the fitter'),
        error:float = typer.Argument(1e-10, help='Relative error to use for early termination'),
        packfile:str = typer.Argument('packfile.ipac', help='Name of the packfile used for fitting'),
        exclude_bf:bool = typer.Option(False, help='Exclude data points with bit flags from PAHFIT')):
    """
        PAHFit works by fitting Spitzer IRS spectra using models for blackbody, hydrogen, PAH and forbidden emissions.
        More info here http://tir.astro.utoledo.edu/jdsmith/research/pahfit.php and here https://pahfit.readthedocs.io/en/latest/
    """
    
    pahfit(inp, output, max_iter, error, packfile, exclude_bf)



def pahfit(input_dir, output, max_iter, error, packfile, exclude_bf):

    try:
        os.mkdir(output)
    except:
        pass

    for subdir, dirs, files in os.walk(fr'./{input_dir}'):
        for specfile in files:
            try:
                if specfile.split('.')[-1]  == 'txt' or specfile.split('.')[-1]  == 'ipac':

                    filepath = subdir + os.sep + specfile

                    outputname = specfile.split('.')[0]
                    typer.echo(f'Performing PAHFit for {typer.style(outputname, bold=True, fg=typer.colors.MAGENTA)}')
                        
                    obs_spectrum = Table.read(filepath, format='ipac')
                    obs_x = obs_spectrum['wavelength'].to(u.micron,
                                                        equivalencies=u.spectral())
                    obs_y = obs_spectrum['flux'].to(u.Jy,
                                                    equivalencies=u.spectral_density(obs_x))
                    obs_unc = obs_spectrum['sigma'].to(u.Jy,
                                                    equivalencies=u.spectral_density(obs_x))

                    obs_x = obs_x.value
                    obs_y = obs_y.value
                    weights = 1. / obs_unc.value

                    # packfile = fr'C:\Users\bpart\OneDrive\Desktop\astro\MASSIVE\data\packFile.ipac'

                    pmodel = PAHFITBase(obs_x, obs_y, filename=packfile)

                    # pick the fitter
                    fit = LevMarLSQFitter()

                    # fit
                    obs_fit = fit(pmodel.model, obs_x, obs_y, weights=weights,
                                maxiter=max_iter, epsilon=1e-10, acc=error)

                    pmodel.save(obs_fit, fr'{output}\{outputname}','ipac')

                    # plot result
                    fontsize = 18
                    font = {'size': fontsize}
                    mpl.rc('font', **font)
                    mpl.rc('lines', linewidth=2)
                    mpl.rc('axes', linewidth=2)
                    mpl.rc('xtick.major', width=2)
                    mpl.rc('ytick.major', width=2)

                    fig, ax = plt.subplots(figsize=(15, 10))

                    # modelfig, modelax = plt.subplots(figsize=(15, 10))

                    # modelax.plot(obs_x, obs_fit(obs_x) /obs_x, "g-")

                    # modelfig.tight_layout()

                    modelcsv = open(fr'{output}\{outputname}_model.csv', 'w')

                    modelcsv.write('observed wavelength (microns), observed flux (Jy), model flux (Jy)\n')
                    ys = obs_fit(obs_x)
                    for i in range(len(obs_x)):

                        modelcsv.write(f'{str(obs_x[i])}, {str(obs_y[i])}, {str(ys[i])}\n')
                    modelcsv.close()

                    pmodel.plot(ax, obs_x, obs_y, obs_fit)
                    ax.plot(obs_x, obs_y/obs_x, "ks", fillstyle="full")

                    ax.set_yscale('linear')
                    ax.set_xscale('log')

                    # use the whitespace better
                    fig.tight_layout()

                    fig.savefig(fr'{output}\{outputname}.png')

                    # chiSquared = calc_reduced_chi_square(obs_fit(obs_x), obs_x, obs_y, obs_unc, len(obs_x), 139)

                    print("SUCCESS", outputname)
            except:
                print('PAHFIT FAILED FOR', specfile)


if __name__ == '__main__':
    app()