import sys
sys.path.append('/global/homes/r/rhliu/projects/repos/ThumbStack')
from headers import *
from catalog import Catalog

def addIntegratedY(cat, U, nu=150.e9):
    """Integrated tSZ signal: int d^2theta n_e sigma_T (k_B T_e / m_e c^2)
    in [sr].
    To get dT in muK*sr, multiply by Tcmb * f(nu).
    Simple power-law fit to Greco et al 2014, fig4.
    """
    print("- add integrated y")

    # in arcmin^2
    yCcyltilda = (cat.Mstellar/1.e11)**3.2 * 1.e-6

    # in arcmin^2
    yCcyl = yCcyltilda * (U.hubble(cat.Z) / U.hubble(0.))**(2./3.)
    yCcyl /= (U.bg.comoving_distance(cat.Z) / (500.*U.bg.h))**2
    # in sr
    yCcyl *= (np.pi/180./60.)**2

    # cat['integratedY'] = yCcyl
    cat.integratedY = yCcyl


def make_Catalog(U, MassConversion, df, name, 
                 out_dir='/pscratch/sd/r/rhliu/projects/Weak_lensing/desi/spec_Y3/',
                 fig_dir='/pscratch/sd/r/rhliu/projects/Weak_lensing/figs/', 
                 cat_fn="/catalog.txt"):
    """Hard coded method for making a ThumbStack Catalog from a dict or DataFrame.

    Args:
        U (_type_): _description_
        MassConversion (_type_): _description_
        df (_type_): _description_
        name (_type_): _description_

    Returns:
        _type_: _description_
    """

    cat = Catalog(U, MassConversion, name=name, save=False, load=False, 
                  out_dir=out_dir, fig_dir=fig_dir, cat_fn=cat_fn)

    cat.nObj = len(df['RA'])
    cat.RA = df['RA'].values # [deg]
    cat.DEC = df['DEC'].values    # [deg]
    cat.Z = df['Z'].values 

    cat.coordX = np.zeros_like(cat.RA)   # [Mpc/h]
    cat.coordY = np.zeros_like(cat.RA)  # [Mpc/h]
    cat.coordZ = np.zeros_like(cat.RA)   # [Mpc/h]
    #
    # displacement from difference,
    # not including the Kaiser displacement,
    # from differences of the observed and reconstructed fields
    cat.dX = np.zeros_like(cat.RA)   # [Mpc/h]
    cat.dY = np.zeros_like(cat.RA)   # [Mpc/h]
    cat.dZ = np.zeros_like(cat.RA)   # [Mpc/h]
    #
    # Kaiser-only displacement
    # originally from differences of the observed and reconstructed fields
    cat.dXKaiser = np.zeros_like(cat.RA)   # [Mpc/h] from cartesian catalog difference
    cat.dYKaiser = np.zeros_like(cat.RA)   # [Mpc/h]
    cat.dZKaiser = np.zeros_like(cat.RA)   # [Mpc/h]
    #
    # velocity in cartesian coordinates
    cat.vX = np.zeros_like(cat.RA)   #[km/s]
    cat.vY = np.zeros_like(cat.RA)   #[km/s]
    cat.vZ = np.zeros_like(cat.RA)   #[km/s]
    #
    # velocity in spherical coordinates,
    # from catalog of spherical displacements
    # cat.vR = np.zeros_like(cat.RA)  # [km/s]   from spherical catalog, >0 away from us
    try:
        cat.vR = df['vR'].values  # [km/s]
    except:
        print('no velocity column for catalogue. Replaced with zeros')
        cat.vR = np.zeros_like(cat.RA)  # [km/s]
    cat.vTheta = np.zeros_like(cat.RA)   # [km/s]
    cat.vPhi = np.zeros_like(cat.RA)  # [km/s]
    #
    # Stellar masses
    cat.Mstellar = np.ones_like(cat.RA) * MassConversion.fmVirTomStar(2e13)  # [M_sun], placeholder
    cat.Mvir = np.ones_like(cat.RA) * 2e13   # [M_sun]
    # cat.Mstellar = df['Mstellar'].values   # [M_sun], from Maraston et al
    # cat.Mvir = df['Mvir'].values
    
    # cat.integratedKSZ = np.zeros_like(cat.RA)
    
    # addIntegratedY(cat, U)
    cat.addHaloMass()
    cat.addIntegratedTau()
    cat.addIntegratedKSZ()
    cat.addIntegratedY()

    return cat