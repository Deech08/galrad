import astropy.coordinates as coord
from astropy.coordinates import frame_transform_graph, SkyCoord

import astropy.units as u

  
import numpy as np

# Coordinate Frame for connecting LMC and Milky Way, with Southern Galactic Pole
# as used in Bland-Hawthorn et al. (2019)

class LMCtoGal(coord.BaseCoordinateFrame):
    """
    A cartesian coordinate system that has the X-Z plane such that it crosses
    through Galactic Center, the LMC, and the Southern Galactic Pole. 
    The z-axis is perpendicular to the Galactic plane, with positive Z pointed
    towards the Southern Galactic Pole.
    Requires knowing the coordinates and distance to the LMC, as well as a distance
    to galactic center of the Sun's position in the Galaxy
    Will default to using the same values used in Bland-Hawthorn et al. (2019)

    Attributes
    ----------
    LMC_coord: `astropy.coordinates.SkyCoord`, optional, must be keyword
        Coordinate of LMC, including distance from sun
    galcen_distance: `astropy.units.Quantity`, optional, must be keyword
        distance of Sun to Galactic Center

    Parameters
    ----------
    x: `astropy.units.Quantity`, optional, must be keyword
        The x coordinate in the LMCtoGal coordinate system
    y: `astropy.units.Quantity`, optional, must be keyword
        The y coordinate in the LMCtoGal coordinate system
    z: `astropy.units.Quantity`, optional, must be keyword
        The z coordinate in the LMCtoGal coordinate system
    """
    default_representation = coord.CartesianRepresentation

    # Specify frame attributes required to fully specify the frame
    galcen_distance = coord.QuantityAttribute(default = 8.122 * u.kpc, unit = u.kpc)

    LMC_coord = coord.CoordinateAttribute(coord.Galactic, 
                                          default = SkyCoord(ra = 80.89416667 * u.deg, 
                                                             dec = -69.75611111 * u.deg,
                                                             distance = 50.0 * u.kpc, 
                                                             frame = "icrs"))

    



def get_transformation_matrix(LMCtoGal_frame, inverse = False):
    """
    Create coordinate transformation matrix for converting between the LMCtoGal frame and Galactocentric frame
    Parameters
    ----------
    LMCtoGal_frame: LMCtoGal class Coordinate frame
    inverse: 'bool', optional, must be keyword
        if True, return the transposed matrix for converting from the Galactocentric frame to the LMCtoGal frame 
    """
    galcen_distance = LMCtoGal_frame.galcen_distance.value
    LMC_coord = LMCtoGal_frame.LMC_coord.transform_to(coord.Galactocentric(galcen_distance = galcen_distance*u.kpc))

    #theta - angle between galcen x axis and LMCtoGal x axis
    
    # get lmc_calcen coords

    theta = np.arctan2(LMC_coord.y, LMC_coord.x).value

    # beta - inclination angle
    beta = np.pi # inclination remains unchanged - z axis is perpendicular to galactic plane but flips sign

    # alpha - tilt angle
    alpha = 0. # no tilt along this axis

    # Generate rotation matrix for coordinate transformation into coord.Galactocentric
    R_matrix = np.array([np.cos(beta)*np.cos(theta), np.cos(beta)*np.sin(theta), -np.sin(beta), 
                        -np.cos(theta)*np.sin(alpha)*-np.sin(beta) - np.cos(alpha)*np.sin(theta), 
                        np.cos(alpha)*np.cos(theta) + np.sin(alpha)*np.sin(beta)*np.sin(theta), 
                        np.cos(beta)*np.sin(alpha), 
                        np.cos(alpha)*np.cos(theta)*np.sin(beta) + np.sin(alpha)*np.sin(theta), 
                        -np.cos(theta)*np.sin(alpha) + np.cos(alpha)*np.sin(beta)*np.sin(theta), 
                        np.cos(alpha)*np.cos(beta)]).reshape(3,3)
    if inverse:
        return R_matrix.transpose()
    else:
        return R_matrix

@frame_transform_graph.transform(coord.DynamicMatrixTransform, LMCtoGal, coord.Galactocentric)
def td_to_galactocentric(LMCtoGal_coord, galactocentric_frame):
    """ Compute the transformation matrix from the Tilted Disk 
        coordinates to Galactocentric coordinates.
    """
    return get_transformation_matrix(LMCtoGal_coord)
    
@frame_transform_graph.transform(coord.DynamicMatrixTransform, coord.Galactocentric, LMCtoGal)
def galactocentric_to_td(galactocentric_coord, LMCtoGal_frame):
    """ Compute the transformation matrix from Galactocentric coordinates to
        Tilted Disk coordinates.
    """
    return get_transformation_matrix(LMCtoGal_frame, inverse = True)


