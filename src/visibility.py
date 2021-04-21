# -*- coding: utf-8 -*-
# code from https://ia-observationtools.readthedocs.io/en/latest/visibility.html | https://github.com/iastro-pt/ObservationTools

# used under MIT License

# The MIT License (MIT)

# Copyright (c) 2016 Daniel Thaagaard Andreasen

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
from __future__ import print_function
import sys
import numpy as np
import datetime as dt
from dateutil import tz
import pickle
from random import choice
from PyAstronomy import pyasl
from astropy.coordinates import SkyCoord
from astropy.coordinates import name_resolve
import ephem
import argparse
import calendar
import typer
import os

try:
  from tqdm import tqdm
except ImportError:
  tqdm = lambda x: x

import io
import matplotlib.pyplot as plt
import matplotlib
replace_figure = True
try:
    from PySide.QtGui import QApplication, QImage
except ImportError:
    try:
        from PyQt4.QtGui import QApplication, QImage
    except ImportError:
        try:
            from PyQt5.QtWidgets import QApplication
            from PyQt5.QtGui import QImage
        except ImportError:
            replace_figure = False


app = typer.Typer()

import warnings

def fxn():
    warnings.warn("deprecated", DeprecationWarning)


def add_clipboard_to_figures():
    # replace the original plt.figure() function with one that supports 
    # clipboard-copying
    oldfig = plt.figure

    def newfig(*args, **kwargs):
        fig = oldfig(*args, **kwargs)
        def clipboard_handler(event):
            if event.key == 'ctrl+c':
                # store the image in a buffer using savefig(), this has the
                # advantage of applying all the default savefig parameters
                # such as background color; those would be ignored if you simply
                # grab the canvas using Qt
                buf = io.BytesIO()
                fig.savefig(buf)
                QApplication.clipboard().setImage(QImage.fromData(buf.getvalue()))
                buf.close()
                print('Ctrl+C pressed: image is now in the clipboard')

        fig.canvas.mpl_connect('key_press_event', clipboard_handler)
        return fig

    plt.figure = newfig

if replace_figure: add_clipboard_to_figures()


def decdeg2dms(dd):
    """ Convert decimal degrees to deg,min,sec """
    is_positive = dd >= 0
    dd = abs(dd)
    minutes,seconds = divmod(dd*3600,60)
    degrees,minutes = divmod(minutes,60)
    degrees = degrees if is_positive else -degrees
    return (degrees,minutes,seconds)

class CacheSkyCoord(SkyCoord):
  @classmethod
  def from_name(cls, name, frame='icrs'):
    try:
      cached = pickle.load(open('CachedSkyCoords.pickle', 'rb'))
    except FileNotFoundError:
      cached = {}
    
    if name in cached:
      return cached[name]
    else:
      original = super(CacheSkyCoord, cls).from_name(name, frame)
      # keep the cached dict manageable
      n = len(cached)
      if n>100:
        # remove a random cached target
        cached.pop(choice(list(cached.keys())))
      cached.update({name: original})
      pickle.dump(cached, open('CachedSkyCoords.pickle', 'wb'))
      return original


ESO_periods = {
  104 : [ (2019, 10, 1),  (2020, 3, 31)],
  103 : [ (2019, 4,  1),  (2019, 9, 30)],
  102 : [ (2018, 10, 1),  (2019, 3, 31)],
  101 : [ (2018, 4,  1),  (2018, 9, 30)],
  100 : [ (2017, 10, 1),  (2018, 3, 31)],
  99  : [ (2017, 4,  1),  (2017, 9, 30)],
  98  : [ (2016, 10, 1),  (2017, 3, 31)],
  97  : [ (2016, 4,  1),  (2016, 9, 30)],
  96  : [ (2015, 10, 1),  (2016, 3, 31)],
  95  : [ (2015, 4,  1),  (2015, 9, 30)],
  94  : [ (2014, 10, 1),  (2015, 3, 31)],
  93  : [ (2014, 4,  1),  (2014, 9, 30)],
  92  : [ (2013, 10, 1),  (2014, 3, 31)],
}

def get_ESO_period(period):
  """ Return the JD of start and end of ESO period """
  assert isinstance(period, str) or isinstance(period, int)
  P = int(period)  

  getjd = lambda y,m,d: pyasl.jdcnv(dt.datetime(y, m, d))
  jd_start, jd_end = [getjd(*d) for d in ESO_periods[P]] 

  return jd_start, jd_end


def StarObsPlot(year=None, targets=None, observatory=None, period=None, 
                hover=False, sunless_hours=None, remove_watermark=False):
  """
    Plot the visibility of target.

    Parameters
    ----------
    year: int
        The year for which to calculate the visibility.
    targets: list
        List of targets.
        Each target should be a dictionary with keys 'name' and 'coord'.
        The key 'name' is a string, 'coord' is a SkyCoord object.
    observatory: string
        Name of the observatory that pyasl.observatory can resolve.
        Basically, any of pyasl.listObservatories().keys()
    period: string, optional
        ESO period for which to calculate the visibility. Overrides `year`.
    hover: boolean, optional
        If True, color visibility lines when mouse over.
    sunless_hours: float, optional
        If not None, plot sunless hours above this airmass
  """
  date = year
  from mpl_toolkits.axes_grid1 import host_subplot
  from matplotlib.ticker import MultipleLocator
  from matplotlib.font_manager import FontProperties
  from matplotlib import rcParams
  rcParams['xtick.major.pad'] = 12
  font0 = FontProperties()
  font1 = font0.copy()
  font0.set_family('sans-serif')
  font0.set_weight('light')
  font1.set_family('sans-serif')
  font1.set_weight('medium')

  # set the observatory
  if isinstance(observatory, dict):
    obs = observatory
  else:
    obs = pyasl.observatory(observatory)


  fig = plt.figure(figsize=(15,10))
  fig.subplots_adjust(left=0.07, right=0.8, bottom=0.15, top=0.88)
  
  # watermak
  if not remove_watermark:
    fig.text(0.99, 0.99, 'Created with\ngithub.com/iastro-pt/ObservationTools',
            fontsize=10, color='gray',
            ha='right', va='top', alpha=0.5)

  # plotting sunless hours?
  shmode = False
  if sunless_hours is not None:
    shmode = True
    # limit in airmass (assumed plane-parallel atm)
    shairmass = sunless_hours
    # correspoing limit in altitude
    from scipy.optimize import bisect
    shalt = 90 - bisect(lambda alt: pyasl.airmassPP(alt) - shairmass, 0, 89)

  if shmode:
    fig.subplots_adjust(hspace=0.35)
    ax = host_subplot(211)
    axsh = host_subplot(212)
    plt.text(0.5, 0.47,
          "- sunless hours above airmass {:.1f} - \n".format(shairmass),
          transform=fig.transFigure, ha='center', va='bottom', fontsize=12)
    plt.text(0.5, 0.465,
          "the thick line above the curves represents the total sunless hours "\
          "for each day of the year",
          transform=fig.transFigure, ha='center', va='bottom', fontsize=10)

  else:
    ax = host_subplot(111)

  for n, target in enumerate(targets):

    target_coord = target['coord']
    target_ra = target_coord.ra.deg
    target_dec = target_coord.dec.deg

    
    if period is not None:
      jd_start, jd_end = get_ESO_period(period)
    else:
      jd_start = pyasl.jdcnv(dt.datetime(year, 1, 1))
      jd_end = pyasl.jdcnv(dt.datetime(year, 12, 31))
    
    jdbinsize = 1 # every day
    each_day = np.arange(jd_start, jd_end, jdbinsize)
    jds = []

    ## calculate the mid-dark times
    sun = ephem.Sun()
    for day in each_day:
      date_formatted = '/'.join([str(i) for i in pyasl.daycnv(day)[:-1]])
      s = ephem.Observer();
      s.date = date_formatted;
      s.lat = ':'.join([str(i) for i in decdeg2dms(obs['latitude'])])
      s.lon = ':'.join([str(i) for i in decdeg2dms(obs['longitude'])])
      jds.append(ephem.julian_date(s.next_antitransit(sun)))


    jds = np.array(jds)

    # Get JD floating point
    jdsub = jds - np.floor(jds[0])

    # Get alt/az of object
    altaz = pyasl.eq2hor(jds, np.ones_like(jds)*target_ra, np.ones_like(jds)*target_dec, \
                        lon=obs['longitude'], lat=obs['latitude'], alt=obs['altitude'])
    ax.plot( jdsub, altaz[0], '-', color='k')

    # label for each target
    plabel = "[{0:2d}]  {1!s}".format(n+1, target['name'])

    # number of target at the top of the curve
    ind_label = np.argmax(altaz[0])
    # or at the bottom if the top is too close to the corners
    # if jdsub[ind_label] < 5 or jdsub[ind_label] > jdsub.max()-5:
    #   ind_label = np.argmin(altaz[0])
    ax.text( jdsub[ind_label], altaz[0][ind_label], str(n+1), color="b", fontsize=14, \
             fontproperties=font1, va="bottom", ha="center")


    if n+1 == 29:
      # too many?
      ax.text(1.1, 1.0-float(n+1)*0.04, "too many targets", ha="left", va="top", transform=ax.transAxes, \
              fontsize=10, fontproperties=font0, color="r")
    else:
      ax.text(1.1, 1.0-float(n+1)*0.04, plabel, ha="left", va="top", transform=ax.transAxes, \
              fontsize=12, fontproperties=font0, color="b")

  if shmode:
    sunless_hours = []
    for day in each_day:
      date_formatted = '/'.join([str(i) for i in pyasl.daycnv(day)[:-1]])
      s = ephem.Observer();
      s.date = date_formatted;
      s.lat = ':'.join([str(i) for i in decdeg2dms(obs['latitude'])])
      s.lon = ':'.join([str(i) for i in decdeg2dms(obs['longitude'])])
      # hours from sunrise to sunset
      td = pyasl.daycnv(ephem.julian_date(s.next_setting(sun)), mode='dt') \
            - pyasl.daycnv(ephem.julian_date(s.next_rising(sun)), mode='dt')
      sunless_hours.append(24 - td.total_seconds() / 3600)

    days = each_day - np.floor(each_day[0])
    axsh.plot(days, sunless_hours, '-', color='k', lw=2)
    axsh.set(ylim=(0, 15), yticks=range(1,15), ylabel='Useful hours',
             yticklabels=[r'${}^{{\rm h}}$'.format(n) for n in range(1,15)])


  ax.text(1.1, 1.03, "List of targets", ha="left", va="top", transform=ax.transAxes, \
          fontsize=12, fontproperties=font0, color="b")


  axrange = ax.get_xlim()
  
  
  if period is None:
    months = range(1, 13)
    ndays = [0] + [calendar.monthrange(date, m)[1] for m in months]
    ax.set_xlim([0, 366])
    ax.set_xticks(np.cumsum(ndays)[:-1] + (np.array(ndays)/2.)[1:])
    ax.set_xticklabels(map(calendar.month_abbr.__getitem__, months), fontsize=10)
    if shmode:
      axsh.set_xlim([0, 366])
      axsh.set_xticks(np.cumsum(ndays)[:-1] + (np.array(ndays)/2.)[1:])
      axsh.set_xticklabels(map(calendar.month_abbr.__getitem__, months), fontsize=10)
  else:
    if int(period) % 2 == 0:
      # even ESO period, Oct -> Mar
      months = [10, 11, 12, 1, 2, 3]
      ndays = [0] + [calendar.monthrange(date, m)[1] for m in months]
      ax.set_xlim([0, 181])
      ax.set_xticks(np.cumsum(ndays)[:-1] + (np.array(ndays)/2.)[1:])
      ax.set_xticklabels(map(calendar.month_abbr.__getitem__, months), fontsize=10)
      if shmode:
        axsh.set_xlim([0, 181])
        axsh.set_xticks(np.cumsum(ndays)[:-1] + (np.array(ndays)/2.)[1:])
        axsh.set_xticklabels(map(calendar.month_abbr.__getitem__, months), fontsize=10)
    else:
      # odd ESO period, Apr -> Sep
      months = range(4, 10)
      ndays = [0] + [calendar.monthrange(date, m)[1] for m in months]
      ax.set_xlim([0, 182])
      ax.set_xticks(np.cumsum(ndays)[:-1] + (np.array(ndays)/2.)[1:])
      ax.set_xticklabels(map(calendar.month_abbr.__getitem__, months), fontsize=10)
      if shmode:
        axsh.set_xlim([0, 182])
        axsh.set_xticks(np.cumsum(ndays)[:-1] + (np.array(ndays)/2.)[1:])
        axsh.set_xticklabels(map(calendar.month_abbr.__getitem__, months), fontsize=10)

  if axrange[1]-axrange[0] <= 1.0:
    jdhours = np.arange(0,3,1.0/24.)
    utchours = (np.arange(0,72,dtype=int)+12)%24
  else:
    jdhours = np.arange(0,3,1.0/12.)
    utchours = (np.arange(0,72, 2, dtype=int)+12)%24


  # Make ax2 responsible for "top" axis and "right" axis
  ax2 = ax.twin()
  # Set upper x ticks
  ax2.set_xticks(np.cumsum(ndays))
  ax2.set_xlabel("Day")

  # plane-parallel airmass
  airmass_ang = np.arange(10, 81, 5)
  geo_airmass = pyasl.airmass.airmassPP(airmass_ang)[::-1]
  ax2.set_yticks(airmass_ang)
  airmassformat = []
  for t in range(geo_airmass.size):
    airmassformat.append("{0:2.2f}".format(geo_airmass[t]))
  ax2.set_yticklabels(airmassformat)#, rotation=90)
  ax2.set_ylabel("Relative airmass", labelpad=32)
  ax2.tick_params(axis="y", pad=6, labelsize=8)
  plt.text(1.02,-0.04, "Plane-parallel", transform=ax.transAxes, ha='left', \
           va='top', fontsize=10, rotation=90)

  ax22 = ax.twin()
  ax22.set_xticklabels([])
  ax22.set_frame_on(True)
  ax22.patch.set_visible(False)
  ax22.yaxis.set_ticks_position('right')
  ax22.yaxis.set_label_position('right')
  ax22.spines['right'].set_position(('outward', 30))
  ax22.spines['right'].set_color('k')
  ax22.spines['right'].set_visible(True)
  airmass2 = list(map(lambda ang: pyasl.airmass.airmassSpherical(90. - ang, obs['altitude']), airmass_ang))
  ax22.set_yticks(airmass_ang)
  airmassformat = []
  for t in range(len(airmass2)): 
    airmassformat.append(" {0:2.2f}".format(airmass2[t]))
  ax22.set_yticklabels(airmassformat, rotation=90)
  ax22.tick_params(axis="y", pad=8, labelsize=8)
  plt.text(1.05,-0.04, "Spherical+Alt", transform=ax.transAxes, ha='left', va='top', \
           fontsize=10, rotation=90)


  ax.set_ylim([0, 91])
  ax.yaxis.set_major_locator(MultipleLocator(15))
  ax.yaxis.set_minor_locator(MultipleLocator(5))
  yticks = ax.get_yticks()
  ytickformat = []
  for t in range(yticks.size):
    ytickformat.append(str(int(yticks[t]))+r"$^\circ$")
  ax.set_yticklabels(ytickformat, fontsize=16)
  ax.set_ylabel("Altitude", fontsize=18)
  yticksminor = ax.get_yticks(minor=True)
  ymind = np.where( np.array(yticksminor) % 15. != 0. )[0]
  yticksminor = [yticksminor[i] for i in ymind]
  yticksminor = np.array(yticksminor)
  ax.set_yticks(yticksminor, minor=True)
  m_ytickformat = []
  for t in range(yticksminor.size):
    m_ytickformat.append(str(int(yticksminor[t]))+r"$^\circ$")
  ax.set_yticklabels(m_ytickformat, minor=True)
  ax.set_ylim([0, 91])

  ax.yaxis.grid(color='gray', linestyle='dashed')
  ax.yaxis.grid(color='gray', which="minor", linestyle='dotted')
  ax2.xaxis.grid(color='gray', linestyle='dotted')

  if period is not None:
    plt.text(0.5, 0.95,
             "Visibility over P{0!s}\n - altitudes at mid-dark time -".format(period),
             transform=fig.transFigure, ha='center', va='bottom', fontsize=12)
  else:
    plt.text(0.5, 0.95,
             "Visibility over {0!s}\n - altitudes at mid-dark time -".format(date),
             transform=fig.transFigure, ha='center', va='bottom', fontsize=12)


  obsco = "Obs coord.: {0:8.4f}$^\circ$, {1:8.4f}$^\circ$, {2:4f} m".format(obs['longitude'], obs['latitude'], obs['altitude'])

  plt.text(0.01,0.97, obsco, transform=fig.transFigure, ha='left', va='center', fontsize=10)
  plt.text(0.01,0.95, obs['name'], transform=fig.transFigure, ha='left', va='center', fontsize=10)

  # interactive!
  if hover:
    main_axis = fig.axes[0]
    all_lines = set(main_axis.get_lines())
    def on_plot_hover(event):
      for line in main_axis.get_lines():
          if line.contains(event)[0]:
            line.set_color('red') # make this line red
            # and all others black
            all_other_lines = all_lines - set([line])
            for other_line in all_other_lines:
              other_line.set_color('black')
            fig.canvas.draw_idle()
    fig.canvas.mpl_connect('motion_notify_event', on_plot_hover)

  return fig



def VisibilityPlot(date=None, targets=None, observatory=None, plotLegend=True, 
                   showMoonDist=True, print2file=False, remove_watermark=False):
  """
    Plot the visibility of target.

    Parameters
    ----------
    date: datetime
        The date for which to calculate the visibility.
    targets: list
        List of targets.
        Each target should be a dictionary with keys 'name' and 'coord'.
        The key 'name' is aa string, 'coord' is a SkyCoord object.
    observatory: string
        Name of the observatory that pyasl.observatory can resolve.
        Basically, any of pyasl.listObservatories().keys()
    plotLegend: boolean, optional
        If True (default), show a legend.
    showMoonDist : boolean, optional
        If True (default), the Moon distance will be shown.
  """

  from mpl_toolkits.axes_grid1 import host_subplot
  from matplotlib.ticker import MultipleLocator
  from matplotlib.font_manager import FontProperties
  from matplotlib import rcParams
  rcParams['xtick.major.pad'] = 12


  if isinstance(observatory, dict):
    obs = observatory
  else:
    obs = pyasl.observatory(observatory)

  fig = plt.figure(figsize=(15,10))
  fig.subplots_adjust(left=0.07, right=0.8, bottom=0.15, top=0.88)

  # watermak
  if not remove_watermark:
    fig.text(0.99, 0.99, 'Created with\ngithub.com/iastro-pt/ObservationTools',
            fontsize=10, color='gray',
            ha='right', va='top', alpha=0.5)

  ax = host_subplot(111)

  font0 = FontProperties()
  font1 = font0.copy()
  font0.set_family('sans-serif')
  font0.set_weight('light')
  font1.set_family('sans-serif')
  font1.set_weight('medium')


  for n, target in enumerate(targets):

    target_coord = target['coord']
    target_ra = target_coord.ra.deg
    target_dec = target_coord.dec.deg

    # JD array
    jdbinsize = 1.0/24./20.
    # jds = np.arange(allData[n]["Obs jd"][0], allData[n]["Obs jd"][2], jdbinsize)
    jd = pyasl.jdcnv(date)
    jd_start = pyasl.jdcnv(date)-0.5
    jd_end = pyasl.jdcnv(date)+0.5
    jds = np.arange(jd_start, jd_end, jdbinsize)
    # Get JD floating point
    jdsub = jds - np.floor(jds[0])
    # Get alt/az of object
    altaz = pyasl.eq2hor(jds, np.ones(jds.size)*target_ra, np.ones(jds.size)*target_dec, \
                        lon=obs['longitude'], lat=obs['latitude'], alt=obs['altitude'])
    # Get alt/az of Sun
    sun_position = pyasl.sunpos(jd)
    sun_ra, sun_dec = sun_position[1], sun_position[2]
    sunpos_altaz = pyasl.eq2hor(jds, np.ones(jds.size)*sun_ra, np.ones(jds.size)*sun_dec, \
                                lon=obs['longitude'], lat=obs['latitude'], alt=obs['altitude'])

    # Define plot label
    plabel = "[{0:2d}]  {1!s}".format(n+1, target['name'])

    # Find periods of: day, twilight, and night
    day = np.where( sunpos_altaz[0] >= 0. )[0]
    twi = np.where( np.logical_and(sunpos_altaz[0] > -18., sunpos_altaz[0] < 0.) )[0]
    night = np.where( sunpos_altaz[0] <= -18. )[0]

    if (len(day) == 0) and (len(twi) == 0) and (len(night) == 0):
      print
      print("VisibilityPlot - no points to draw")
      print

    mpos = pyasl.moonpos(jds)
    # mpha = pyasl.moonphase(jds)
    # mpos_altaz = pyasl.eq2hor(jds, mpos[0], mpos[1],
    #                            lon=obs['longitude'], lat=obs['latitude'], alt=obs['altitude'])
    # moonind = np.where( mpos_altaz[0] > 0. )[0]

    if showMoonDist:
      mdist = pyasl.getAngDist(mpos[0], mpos[1], np.ones(jds.size)*target_ra, \
                              np.ones(jds.size)*target_dec)
      bindist = int((2.0/24.)/jdbinsize)
      firstbin = np.random.randint(0,bindist)
      for mp in range(0, int(len(jds)/bindist)):
        bind = firstbin+mp*bindist
        if altaz[0][bind]-1. < 5.: continue
        ax.text(jdsub[bind], altaz[0][bind]-1., str(int(mdist[bind]))+r"$^\circ$", ha="center", va="top", \
                fontsize=8, stretch='ultra-condensed', fontproperties=font0, alpha=1.)


    if len(twi) > 1:
      # There are points in twilight
      linebreak = np.where( (jdsub[twi][1:]-jdsub[twi][:-1]) > 2.0*jdbinsize)[0]
      if len(linebreak) > 0:
        plotrjd = np.insert(jdsub[twi], linebreak+1, np.nan)
        plotdat = np.insert(altaz[0][twi], linebreak+1, np.nan)
        ax.plot( plotrjd, plotdat, "-", color='#BEBEBE', linewidth=1.5)
      else:
        ax.plot( jdsub[twi], altaz[0][twi], "-", color='#BEBEBE', linewidth=1.5)

    ax.plot( jdsub[night], altaz[0][night], '.k', label=plabel)
    ax.plot( jdsub[day], altaz[0][day], '.', color='#FDB813')

    altmax = np.argmax(altaz[0])
    ax.text( jdsub[altmax], altaz[0][altmax], str(n+1), color="b", fontsize=14, \
             fontproperties=font1, va="bottom", ha="center")

    if n+1 == 29:
      ax.text( 1.1, 1.0-float(n+1)*0.04, "too many targets", ha="left", va="top", transform=ax.transAxes, \
              fontsize=10, fontproperties=font0, color="r")
    else:
      ax.text( 1.1, 1.0-float(n+1)*0.04, plabel, ha="left", va="top", transform=ax.transAxes, \
              fontsize=12, fontproperties=font0, color="b")

  ax.text( 1.1, 1.03, "List of targets", ha="left", va="top", transform=ax.transAxes, \
          fontsize=12, fontproperties=font0, color="b")

  axrange = ax.get_xlim()
  ax.set_xlabel("UT [hours]")

  if axrange[1]-axrange[0] <= 1.0:
    jdhours = np.arange(0,3,1.0/24.)
    utchours = (np.arange(0,72,dtype=int)+12)%24
  else:
    jdhours = np.arange(0,3,1.0/12.)
    utchours = (np.arange(0,72, 2, dtype=int)+12)%24
  ax.set_xticks(jdhours)
  ax.set_xlim(axrange)
  ax.set_xticklabels(utchours, fontsize=18)

  # Make ax2 responsible for "top" axis and "right" axis
  ax2 = ax.twin()
  # Set upper x ticks
  ax2.set_xticks(jdhours)
  ax2.set_xticklabels(utchours, fontsize=18)
  ax2.set_xlabel("UT [hours]")

  # Horizon angle for airmass
  airmass_ang = np.arange(5.,90.,5.)
  geo_airmass = pyasl.airmass.airmassPP(90.-airmass_ang)
  ax2.set_yticks(airmass_ang)
  airmassformat = []
  for t in range(geo_airmass.size):
    airmassformat.append("{0:2.2f}".format(geo_airmass[t]))
  ax2.set_yticklabels(airmassformat, rotation=90)
  ax2.set_ylabel("Relative airmass", labelpad=32)
  ax2.tick_params(axis="y", pad=10, labelsize=10)
  plt.text(1.015,-0.04, "Plane-parallel", transform=ax.transAxes, ha='left', \
           va='top', fontsize=10, rotation=90)

  ax22 = ax.twin()
  ax22.set_xticklabels([])
  ax22.set_frame_on(True)
  ax22.patch.set_visible(False)
  ax22.yaxis.set_ticks_position('right')
  ax22.yaxis.set_label_position('right')
  ax22.spines['right'].set_position(('outward', 25))
  ax22.spines['right'].set_color('k')
  ax22.spines['right'].set_visible(True)
  airmass2 = list(map(lambda ang: pyasl.airmass.airmassSpherical(90. - ang, obs['altitude']), airmass_ang))
  ax22.set_yticks(airmass_ang)
  airmassformat = []
  for t in airmass2:
      airmassformat.append("{0:2.2f}".format(t))
  ax22.set_yticklabels(airmassformat, rotation=90)
  ax22.tick_params(axis="y", pad=10, labelsize=10)
  plt.text(1.045,-0.04, "Spherical+Alt", transform=ax.transAxes, ha='left', va='top', \
           fontsize=10, rotation=90)

  ax3 = ax.twiny()
  ax3.set_frame_on(True)
  ax3.patch.set_visible(False)
  ax3.xaxis.set_ticks_position('bottom')
  ax3.xaxis.set_label_position('bottom')
  ax3.spines['bottom'].set_position(('outward', 50))
  ax3.spines['bottom'].set_color('k')
  ax3.spines['bottom'].set_visible(True)

  ltime, ldiff = pyasl.localtime.localTime(utchours, np.repeat(obs['longitude'], len(utchours)))
  jdltime = jdhours - ldiff/24.
  ax3.set_xticks(jdltime)
  ax3.set_xticklabels(utchours)
  ax3.set_xlim([axrange[0],axrange[1]])
  ax3.set_xlabel("Local time [hours]")

  ax.set_ylim([0, 91])
  ax.yaxis.set_major_locator(MultipleLocator(15))
  ax.yaxis.set_minor_locator(MultipleLocator(5))
  yticks = ax.get_yticks()
  ytickformat = []
  for t in range(yticks.size): ytickformat.append(str(int(yticks[t]))+r"$^\circ$")
  ax.set_yticklabels(ytickformat, fontsize=16)
  ax.set_ylabel("Altitude", fontsize=18)
  yticksminor = ax.get_yticks(minor=True)
  ymind = np.where( np.array(yticksminor) % 15. != 0. )[0]
  yticksminor = [yticksminor[i] for i in ymind]
  yticksminor = np.array(yticksminor)
  ax.set_yticks(yticksminor, minor=True)
  m_ytickformat = []
  for t in range(yticksminor.size): m_ytickformat.append(str(int(yticksminor[t]))+r"$^\circ$")
  ax.set_yticklabels(m_ytickformat, minor=True)
  ax.set_ylim([0, 91])

  ax.yaxis.grid(color='gray', linestyle='dashed')
  ax.yaxis.grid(color='gray', which="minor", linestyle='dotted')
  ax2.xaxis.grid(color='gray', linestyle='dotted')

  plt.text(0.5,0.95,"Visibility on {0!s}".format(date.date()), \
           transform=fig.transFigure, ha='center', va='bottom', fontsize=20)

  if plotLegend:
    line1 = matplotlib.lines.Line2D((0,0),(1,1), color='#FDB813', linestyle="-", linewidth=2)
    line2 = matplotlib.lines.Line2D((0,0),(1,1), color='#BEBEBE', linestyle="-", linewidth=2)
    line3 = matplotlib.lines.Line2D((0,0),(1,1), color='k', linestyle="-", linewidth=2)

    lgd2 = plt.legend((line1,line2,line3),("day","twilight","night",), \
                        bbox_to_anchor=(0.88, 0.13), loc='best', borderaxespad=0.,prop={'size':12}, fancybox=True)
    lgd2.get_frame().set_alpha(.5)

  obsco = "Obs coord.: {0:8.4f}$^\circ$, {1:8.4f}$^\circ$, {2:4f} m".format(obs['longitude'], obs['latitude'], obs['altitude'])

  plt.text(0.01,0.97, obsco, transform=fig.transFigure, ha='left', va='center', fontsize=10)
  plt.text(0.01,0.95, obs['name'], transform=fig.transFigure, ha='left', va='center', fontsize=10)

  return fig


import astropy.units as u
import csv



@app.command()
def list_obs():
  ''' Lists all available observatories
  '''
  pyasl.listObservatories()


@app.command()
def vis_year(
    input_file:str = typer.Argument('input.csv', help='Path to the csv/txt file containing the targets (format: "target,ra,dec"), where ra, dec are the ICRS coordinates of target in degrees in J2000 equinox'),
    output_dir:str = typer.Argument('output', help='Folder to write outputs to'),
    site:str = typer.Argument('cfht', help='Observatory. Refer to list_obs for a complete list'),
    year:int = typer.Argument(2021, help='Year (format: yyyy)')
  ):

  '''
  Plots visibility plots for each target in the input file at the given site during the given year. Ensure input file is in the format 'name,ra,dec'
  '''

  dataFile = open(input_file)
  dataFile = csv.reader(dataFile)

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
      fig = StarObsPlot(year=year, targets=targets, observatory=site, 
                        period=None, hover=False, sunless_hours=None,
                        remove_watermark=True)

      typer.echo('Plotting completed for '+ typer.style(name, fg=typer.colors.MAGENTA, bold=True))
      fig.savefig(fr'{output_dir}/{name.strip()}.png')
      plt.close(fig)


@app.command()
def vis_day(
    input_file:str = typer.Argument('input.csv', help='Path to the csv/txt file containing the targets (format: "target,ra,dec"), where ra, dec are the ICRS coordinates of target in degrees in J2000 equinox'),
    output_dir:str = typer.Argument('output', help='Folder to write outputs to'),
    site:str = typer.Argument('cfht', help='Observatory'),
    date:str = typer.Argument(None, help='Date (format: yyyy-mm-dd)')
  ):

  '''
  Plots visibility plots for each target in the input file at the given site during the given day. 

  '''

  # convert the date parameter into Python datetime format
  ymd = [int(i) for i in date.split('-')]
  date = dt.datetime(*ymd)

  dataFile = open(input_file)
  dataFile = csv.reader(dataFile)

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
      fig = VisibilityPlot(date=date, targets=targets, observatory=site)

      typer.echo('Plotting completed for '+ typer.style(name, fg=typer.colors.MAGENTA, bold=True))
      fig.savefig(fr'{output_dir}/{name.strip()}.png')
      plt.close(fig)

if __name__ == '__main__':
  '''
  This module contains a wrapper on top of visibility.py (https://ia-observationtools.readthedocs.io/en/latest/visibility.html | https://github.com/iastro-pt/ObservationTools) to allow users to 
  obtain visibility plots for large numbers of objects on separate plots without having to write any code themselves. Whereas the original plots the visibilities for all targets
  on the same plot, this version plots them separately. The original also uses object names to get the coordinates, but this version directly uses the 
  ICRS J2000 coordinates  of the targets to reduce the chances of objects being confused with other objects.
  '''
  app()