import numpy as np
import astropy.units as u
import pickle
from scipy.interpolate import griddata
import os
from .LMCtoGalCoord import LMCtoGal
from astropy.coordinates import SkyCoord


directory = os.path.dirname(__file__)
def get_egb_extra_term(egb_template_filename = None):
    """
    returns phi from extragalactic background to be subtracted 

    Parameters
    ----------
    egb_template_filename:
        filename of EGB ionizing flux data
    """

    # get data
    if egb_template_filename == None:
        egb_template_filename = os.path.join(directory, "data/egb_Fox2014.sed")

    nu, fnu = np.loadtxt(egb_template_filename, unpack = True)
    nu*= u.Hz
    fnu*= u.erg * u.cm**-2 * u.s**-1 * u.Hz**-1
    nu_e = nu.to(u.Ry, u.spectral())
    mask = nu_e > 1*u.Ry
    phi_nu = fnu / nu.to(u.erg, u.spectral()) * u.photon
    return np.trapz(phi_nu[mask], nu[mask])

def get_MW_IonizingFlux(c, data_filename = None):
    """
    Returns ionizing photon flux at provided 3D coordinate from Milky Way
    
    Parameters
    ----------
    c: `astropy.coordinates.SkyCoord`
        Coordinate must have 3D info
    data_filename:
        filename of MW ionizing flux data
    """
    # get data
    if data_filename == None:
        # load from default location
        data_filename = os.path.join(directory,"data/mw_ionizing_field.txt")
        
    data = np.loadtxt(data_filename, unpack = True)
    
    # transform coordinates
    c = c.transform_to("galactocentric")
    r = np.sqrt(c.x**2 + c.y**2 + c.z**2).to(u.kpc).value
    rho = np.sqrt(c.x**2 + c.y**2).to(u.kpc).value
    th = np.arctan2(np.abs(c.z.to(u.kpc).value), rho)
    
    PHI = griddata(data[:2,:].T, data[2,:], 
                   (r, th), 
                   method = 'cubic', 
                   rescale = True)
    
    res =  10**PHI * u.photon * u.cm**-2 * u.s**-1

    # egb = get_egb_extra_term()
    # out = res - egb
    # mask = out<0* u.photon * u.cm**-2 *u.s**-1
    # out[mask] = 0.1* u.photon * u.cm**-2 *u.s**-1

    return res


    

def get_LMC_IonizingFlux(c, lmc_par = None, LMC_coord = None):
    """
    Returns ionizing photon flux at provided 3D coordinate from LMC

    Parameters
    ----------
    c: `astropy.coordinates.SkyCoord`
        Coordinate must have 3D info
    lmc_par: `number`, optional, must be keyword
        parameter for inverse squared function, default to result from DK
    LMC_coord: `astropy.coordinates.SkyCoord`, optional, must be keyword
        coordinate of LMC in 3D
    """
    def inverse_r2(r, lmc_par):
        return lmc_par/r**2 
    # set lmc_par
    if lmc_par == None:
        lmc_par = 17735370.35471095


    # transform_coordinates
    if LMC_coord is None:
        if hasattr(c, "LMC_coord"):
            if c.LMC_coord is not None:
                LMC_coord = c.LMC_coord
            else:
                LMC_coord = SkyCoord(ra = 80.89416667 * u.deg, 
                                 dec = -69.75611111 * u.deg,
                                 distance = 50.0 * u.kpc, 
                                 frame = "icrs")
        else:
            LMC_coord = SkyCoord(ra = 80.89416667 * u.deg, 
                                 dec = -69.75611111 * u.deg,
                                 distance = 50.0 * u.kpc, 
                                 frame = "icrs")


    LMC_distance = c.separation_3d(LMC_coord).to(u.kpc).value


    res =  inverse_r2(LMC_distance, lmc_par) * u.photon * u.cm**-2 *u.s**-1

    # egb = get_egb_extra_term()
    # out = res - egb
    # mask = out<0* u.photon * u.cm**-2 *u.s**-1
    # out[mask] = 0.1* u.photon * u.cm**-2 *u.s**-1

    return res
    


def get_SMC_IonizingFlux(c, smc_par = None, SMC_coord = None):
    """
    Returns ionizing photon flux at provided 3D coordinate from LMC

    Parameters
    ----------
    c: `astropy.coordinates.SkyCoord`
        Coordinate must have 3D info
    smc_par: `number`, optional, must be keyword
        parameter for inverse squared function, default to result from DK
    SMC_coord: `astropy.coordinates.SkyCoord`, optional, must be keyword
        coordinate of SMC in 3D

    """

    def inverse_r2(r, smc_par):
        return smc_par/r**2 
    # set lmc_par
    if smc_par == None:
        smc_par = 10398694.06801058

    # SMC coordinate
    if SMC_coord is None:
        SMC_coord = SkyCoord(ra = 13.15833333*u.deg, dec = -72.80027778*u.deg, 
                             distance = 60*u.kpc, frame = "icrs")
    

    # transform coordinates
    SMC_distance = c.separation_3d(SMC_coord).to(u.kpc).value

    res =  inverse_r2(SMC_distance, smc_par) * u.photon * u.cm**-2 *u.s**-1

    # egb = get_egb_extra_term()
    # out = res - egb
    # mask = out<0* u.photon * u.cm**-2 *u.s**-1
    # out[mask] = 0.1* u.photon * u.cm**-2 *u.s**-1

    return res

def MW_LMC_SMC_IonizingFlux(c, data_filename = None, lmc_par = None, smc_par =None, 
                            SMC_coord = None, LMC_coord = None, as_dict = False):
    """
    Returns ionizing photon flux at provided 3D coordinate from combination of MW, LMC, and SMC, 
    as a dictionary of each and a total

    Parameters
    ----------
    c: `astropy.coordinates.SkyCoord`
        Coordinate must have 3D info
    data_filename:
        filename of MW ionizing flux data
    smc_par: `number`, optional, must be keyword
        SMC parameter for inverse squared function, default to result from DK
    lmc_par: `number`, optional, must be keyword
        LMC parameter for inverse squared function, default to result from DK
    SMC_coord: `astropy.coordinates.SkyCoord`, optional, must be keyword
        coordinate of SMC in 3D
    LMC_coord: `astropy.coordinates.SkyCoord`, optional, must be keyword
        coordinate of LMC in 3D
    as_dict: `bool`, optional, must be keyword
        if True, returns dictionary of individual and total contributions
        if False, returns total contribution as single value
    """

    res = {}

    res["MW"] = get_MW_IonizingFlux(c, data_filename = data_filename)
    res["LMC"] = get_LMC_IonizingFlux(c, lmc_par = lmc_par, LMC_coord = LMC_coord)
    res["SMC"] = get_SMC_IonizingFlux(c, smc_par = smc_par, SMC_coord = SMC_coord)

    res["TOTAL"] = res["MW"] + res["LMC"] + res["SMC"]

    if as_dict:
        return res

    else:
        return res["TOTAL"]

def get_input_spectra(c, spectra_template_filename = None, egb_template_filename = None, **kwargs):
    """
    Returns input spectra, as individual components, and combined total for MW, LMC, SMC,
    as a dictionary 

    Parameters
    ----------
    c: `astropy.coordinates.SkyCoord`
        Coordinate must have 3D info
    spectra_template_filename: `str`, optional, must be keyword
        template tabulated spectrum, defaults to that from Fox et al. 2005 (Figure 8)
    data_filename:
        filename of MW ionizing flux data
    smc_par: `number`, optional, must be keyword
        SMC parameter for inverse squared function, default to result from DK
    lmc_par: `number`, optional, must be keyword
        LMC parameter for inverse squared function, default to result from DK
    SMC_coord: `astropy.coordinates.SkyCoord`, optional, must be keyword
        coordinate of SMC in 3D
    LMC_coord: `astropy.coordinates.SkyCoord`, optional, must be keyword
        coordinate of LMC in 3D
    as_dict: `bool`, optional, must be keyword
        if True, returns dictionary of individual and total contributions
        if False, returns total contribution as single value
    """

    res = {}

    
    kwargs["as_dict"] = True

    norms = MW_LMC_SMC_IonizingFlux(c, **kwargs)

    # get spectral shape
    if spectra_template_filename == None:
        spectra_template_filename = os.path.join(directory,"data/MW.0kpc.dat")

    freq, F_nu = np.loadtxt(spectra_template_filename, unpack = True)
    freq *= u.Hz
    F_nu *= u.erg * u.cm**-2 * u.s**-1 * u.Hz**-1

    # flip if needed
    if np.argmin(freq) != 0:
        freq = np.flip(freq)
        F_nu = np.flip(F_nu)

    # convert to photon flux
    phi_nu = F_nu / freq.to(u.erg, u.spectral()) * u.photon


    # get min frequency to integrate over
    f_0 = 1*u.Ry
    f_0 = f_0.to(u.Hz, u.spectral())
    mask = freq > f_0

    # Get default PHI
    PHI_0 = np.trapz(phi_nu[mask], freq[mask])

    # get normalization correction_factors
    correction_factors = {}
    for key in ["MW","SMC","LMC","TOTAL"]:
        correction_factors[key] = norms[key]/PHI_0
        res[key] = F_nu[:,None] * correction_factors[key]
        res[key] = res[key].T

    res["nu"] = freq

    return res, norms



    






    