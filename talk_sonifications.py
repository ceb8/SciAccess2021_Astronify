
# Imports
import numpy as np
from astropy.table import Table, Column
from astropy.modeling import models
import astropy.units as u
import matplotlib.pyplot as plt

from astronify.series import SoniSeries
from astronify.simulator import simulated_lc


# Sonification definition plot
# Slide 7

test_data = simulated_lc("sine", lc_length=300, sine_period=300)

test_obj = SoniSeries(test_data)
test_obj.sonify()

test_obj.write("sonifications/sine_ex.wav")


# The basic usage example
# Slide 8-9

# Data
test_data = Table({"time":[0,1,2,3,4,5,6,7,8,9,10],
                   "flux":[0,1,2,3,4,5,6,7,8,9,10]})

# Sonification
test_obj = SoniSeries(test_data)
test_obj.note_spacing = 0.2  # sec
test_obj.note_duration = 0.3  # sec

test_obj.sonify()
test_obj.write("sonifications/line.wav")

# Visualization
fig, ax = plt.subplots(figsize=(12,8), facecolor="black", dpi=150)

ax.set_facecolor('black')
ax.spines['bottom'].set_color('white')
ax.spines['top'].set_color('white') 
ax.spines['right'].set_color('white')
ax.spines['left'].set_color('white')
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel("Time", fontsize=24, color="white")
ax.set_ylabel("Flux", fontsize=24, color="white")
    
ax.plot(test_data["time"],test_data["flux"], color="white",  linewidth=2)

fig.savefig("plots/line.jpg")


# Sonifying a GALEX flare
# Slide 12-13

# Data
galex_lc = Table.read("data/6371232513623328240_LC.csv")

# Sonification
galex_obj = SoniSeries(galex_lc[:153], time_col="t_mean", val_col="flux_bgsub")
galex_obj.note_duration = 0.8
galex_obj.note_spacing = 0.04
galex_obj.sonify()

galex_obj.write("sonifications/6371232513623328240.wav")

# Visualization
times = galex_lc['t_mean'][:153]
fluxes = galex_lc['flux_bgsub'][:153]
kid = 10263691

f, ax = plt.subplots(figsize=(9, 8))   
ax.tick_params(axis='both', which='major', labelsize=14)

# Turn times into relative time and into minutes
stInt = times[0]
times = [(x-stInt)/60 for x in times]

ax.plot(times, fluxes, color='#460099', label=f"Kepler ID: {kid}") 
ax.axhline(np.median(fluxes),ls="--", color='#717171',linewidth=2)    
    
ax.set_xlabel('Time (Min)',fontsize=16)
ax.set_ylabel(r'Flux ($erg\cdot sec^{-1} cm^{-2} \AA^{-1}$)',fontsize=16)
ax.legend(prop={'size': 14})

f.savefig("plots/6371232513623328240.jpg")


# Waveband comparison
# Slide 14-15


# Data
scale = 1
bb_9000 = models.BlackBody(9000*u.K, scale=scale * u.erg / (u.cm ** 2 * u.AA * u.s * u.sr))
xvals = np.linspace(1000, 10000, 100)
yvals =  bb_9000(xvals*u.AA)*np.pi*u.sr

bb_table = Table(names=["lambda","flux"], data=[xvals, yvals])

# Adding column to show where telesopes observe
bb_table.add_column(Column(name="telescopes", length=len(bb_table), data=np.nan))

bb_table["telescopes"][0] = 0 # this is for timing, will have to manually remove
bb_table["telescopes"][(bb_table["lambda"] > 1700) & (bb_table["lambda"] < 3000)] = 800
bb_table["telescopes"][(bb_table["lambda"] > 4300) & (bb_table["lambda"] < 8900)] = 400

# Sonification

# Black body curve
bb = SoniSeries(bb_table, time_col="lambda")
bb.note_spacing = 0.05
bb.sonify()

bb.write("sonifications/blackbody.wav")

# Telescope passbands

bb = SoniSeries(bb_table, time_col="lambda", val_col="telescopes")
bb.note_spacing = 0.05
bb.sonify()

bb.write("sonifications/telescopes.wav")

# I then used Audacity to combine the two wave files and silence the dummy tone
# that makes sure they are in alignment
# The resulting sonification is blackbody_telescopes.wav

# Visualization
f, ax = plt.subplots(figsize=(7, 5))
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)

ax.plot(xvals, yvals, label="9000K Blackbody", color="black", ls=":", alpha=1) 
   
# GALEX
ax.axvspan(1700, 3000, alpha=0.2, color='#1581e1')
ax.text(x=2350, y=8*10**7, s="GALEX\nNUV", va='center', ha="center", color="#1581e1", fontsize=14)

# Kepler
ax.axvspan((430*u.nm).to(u.angstrom).value, (890*u.nm).to(u.angstrom).value, alpha=0.2, color='#df0040')
ax.text(x=6600, y=8*10**7, s="<------------ Kepler ------------>", va='center', ha="center", color="#df0040", fontsize=14)


ax.set_xlabel("Wavelength ($\AA$)", fontsize=14)
ax.set_ylabel("Flux density ($erg\; s^{-1} cm^{-2} \AA^{-1}$)", fontsize=14)
ax.set_ylim(0, 9*10**7)
ax.legend(fontsize=14, loc="lower right")

f.tight_layout()
f.savefig("plots/multiwavelength.jpg")


# Galex flare with Kepler data
# Slide 16-17

# Data

g_sub_lc = Table.read("data/galex_lc.ecsv")
k_sub_lc = Table.read("data/kepler_lc.ecsv")

# Dictionary with extra into about the Kepler light curve
kep_info = {"k_ind": 8,
            "k_qui_m": -50.35789232862335,
            "k_qui_b": 123652643.70482102,
            "k_qui_std": 6.1929104641125345,
            "flux_min": 4664.331,  # e-/s
            "flux_max": 5206.7856}  # e-/s

# Sonification

# Kepler light curve
kepler = SoniSeries(k_sub_lc)

kepler.note_spacing = 0.2
kepler.note_duration = 0.9
kepler.pitch_mapper.pitch_map_args["minmax_value"] = [kep_info["flux_min"], kep_info["flux_max"]]

kepler.sonify()

kepler.write("sonifications/Kepler.wav")

# Kepler quiescent flux min/max lines
qui_flux = kep_info["k_qui_m"]*k_sub_lc['time'].jd+kep_info["k_qui_b"]

min_qui = Table(names=["time","flux"], data=[k_sub_lc['time'].jd, qui_flux-kep_info["k_qui_std"]])

qui_l = SoniSeries(min_qui)
qui_l.note_spacing = 0.2
qui_l.note_duration = 0.9
qui_l.pitch_mapper.pitch_map_args["minmax_value"] = [kep_info["flux_min"], kep_info["flux_max"]]
qui_l.sonify()

qui_l.write("sonifications/qui_min.wav")

max_qui = Table(names=["time","flux"], data=[k_sub_lc['time'].jd, qui_flux+kep_info["k_qui_std"]])

qui_h = SoniSeries(max_qui)
qui_h.note_spacing = 0.2
qui_h.note_duration = 0.9
qui_h.pitch_mapper.pitch_map_args["minmax_value"] = [kep_info["flux_min"], kep_info["flux_max"]]
qui_h.sonify()

qui_h.write("sonifications/qui_max.wav")

# GALEX light curve

# Calculating note spacing to match the Kepler sonification
delt_k = np.mean(k_sub_lc["time"].jd[1:]-k_sub_lc["time"].jd[:-1])
delt_g = np.median(g_sub_lc["time"].jd[1:]-g_sub_lc["time"].jd[:-1])
g_x = 0.2*delt_g/delt_k

# Add a dummy point at the start of the kepler time frame for alignment
# Must be silenced during post-processing
galex_tab = Table.copy(g_sub_lc[20:39]["time","flux_mcatbgsub_acorr"])
galex_tab.add_row([k_sub_lc["time"][0],0])
galex_tab.sort("time")

galex = SoniSeries(galex_tab, val_col='flux_mcatbgsub_acorr')

galex.note_spacing = 0.001
galex.note_duration = 0.2
galex.pitch_mapper.pitch_map_args["minmax_value"] = [g_sub_lc["flux_mcatbgsub_acorr"].min(),
                                                     g_sub_lc["flux_mcatbgsub_acorr"].max()]
galex.pitch_mapper.pitch_map_args["pitch_range"] = [100,800]
galex.sonify()

galex.write("sonifications/galex.wav")

# I then used Audacity to combine all of these sonifications into one, adjusting volumes
# to make the GALEX flare blip sufficiently audible (and silencing the dummy point
# added at the start). The resulting file is keple_galex.wav.

# I also made another sonification withincluding only the Kepler data, so we can
# hear how there is no evidence of the flare in the Kepler light curve
# This file is full_kepler.wav.

# Visualization
f, axs = plt.subplots(nrows=2, sharex=True, figsize=(6, 4))

axs[0].plot(g_sub_lc['time'].jd, g_sub_lc['flux_mcatbgsub_acorr'], color="#a3003e")
axs[0].fill_between(g_sub_lc['time'].jd, 
                    g_sub_lc['flux_mcatbgsub_acorr']-g_sub_lc['flux_mcatbgsub_err'], 
                    g_sub_lc['flux_mcatbgsub_acorr']+g_sub_lc['flux_mcatbgsub_err'],
                    color="#a3003e", alpha=0.5)

axs[1].plot(k_sub_lc['time'].jd, k_sub_lc['flux'], color="#0089aa")
k_ind = kep_info["k_ind"]
axs[1].errorbar(k_sub_lc['time'][k_ind].jd, k_sub_lc['flux'][k_ind], yerr=k_sub_lc['flux_err'][k_ind], 
                color='black', capsize=3)

qui_flux = kep_info["k_qui_m"]*k_sub_lc['time'].jd+kep_info["k_qui_b"]
axs[1].fill_between(k_sub_lc['time'].jd, qui_flux-kep_info["k_qui_std"], qui_flux+kep_info["k_qui_std"],
                    color="#0089aa", alpha=0.5)

axs[0].axes.xaxis.set_visible(False)
axs[0].axes.yaxis.set_visible(False)

axs[1].axes.xaxis.set_visible(False)
axs[1].axes.yaxis.set_visible(False)

f.tight_layout()
f.savefig("plots/galex_kepler_flare.jpg")


# GALEX flare inset
# Slide 16

# Displaying the GALEX light curve on its own so we can see the flare shape

# Sonification
galex = SoniSeries(g_sub_lc[20:39], val_col='flux_mcatbgsub_acorr')

galex.note_spacing = 0.2
galex.note_duration = 0.6
galex.pitch_mapper.pitch_map_args["minmax_value"] = [g_sub_lc["flux_mcatbgsub_acorr"].min(),
                                                     g_sub_lc["flux_mcatbgsub_acorr"].max()]
galex.pitch_mapper.pitch_map_args["pitch_range"] = [100,800]

galex.sonify()

galex.write("sonifications/galex_zoom.wav")

# Visualization
f, ax = plt.subplots(figsize=(6, 4))

ax.plot(g_sub_lc['time'][20:39].jd, g_sub_lc['flux_mcatbgsub_acorr'][20:39], color="#a3003e", linewidth=3)
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)

f.tight_layout()
f.savefig("plots/galex_flare_zoom.jpg")
