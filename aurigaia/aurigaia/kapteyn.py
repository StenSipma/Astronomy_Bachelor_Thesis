from . import read_aurigaia as ra
from . import coords as cc
available_props = ["AccretedFlag", "Age",
                   "EffectiveTemperature", "EffectiveTemperatureObs", "EffectiveTemperatureError",
                   "Extinction31",
                   "GBmagnitude", "GBmagnitudeError", "GBmagnitudeObs",
                   "GRmagnitude", "GRmagnitudeError", "GRmagnitudeObs",
                   "Gmagnitude", "GmagnitudeError", "GmagnitudeObs",
                   "GravPotential",
                   "HCoordinates", "HCoordinatesObs", "HCoordinateErrors",  # ra, dec, parallax
                   "HVelocities", "HVelocitiesObs", "HVelocityErrors",  # pm_ra_cosdec, pm_dec, vr
                   "Pos", "PosObs",  # Note conversion, added!
                   "Vel", "VelObs",
                   "Mass",
                   "Metallicity",
                   "ParticleID",
                   "SurfaceGravity", "SurfaceGravityError", "SurfaceGravityObs",
                   "VabsMagnitude", "IabsMagnitude",
                   # Sub grouped under "Magnitude"
                   "Umagnitude", "Bmagnitude", "Rmagnitude", "Jmagnitude", 
                   "Hmagnitude", "Kmagnitude", "Vmagnitude", "Imagnitude", 
                   ]
available_halos = [6, 16, 21, 23, 24, 27]
available_angles = [30, 120, 210, 300]
available_ext = [True, False]

import copy
from .read_aurigaia import RangeFilter, DifferenceFilter


def kapteyn_aurigaia(halo_n, datasets, filters=None, angle=30, ext=True, fsample=None, verbose=False, units=True):
    data_dir = "/net/gaia2/data/users/callingham/data/aurigaia/data/ICC/v3/chabrierIMF/level3_MHD"
    halo_dir = f"/halo_{halo_n}/mockdir_angle{angle:03d}"
    base_dir = data_dir + halo_dir
    if ext:
        base_name = f"mock_{angle:03d}"
    else:
        base_name = f"mock_noex_{angle:03d}"

    print(f"Reading data from mock {halo_dir}, {base_name}")

    if fsample is not None:
        print(f"Randomly sampling {fsample} fraction")
        if filters is None:
            filters = [ra.RandomSampleFilter(fsample)]
        else:
            filters.append(ra.RandomSampleFilter(fsample))

    checked_datasets = check_posvel_datasets(datasets)

    star_data = ra.read_aurigaia(
        base_dir, base_name, datasets=checked_datasets, filters=filters, verbose=verbose)

    star_data = add_posvel_data(star_data, datasets)

    if not units:
        for p in list(star_data.keys()):
            try:
                star_data[p] = star_data[p].value
            except Exception:
                pass

    # "Header":star_data["Header"],
    data = {"Parameters": star_data["Parameters"], "Header":star_data["Header"]}
    del star_data["Header"]
    del star_data["Parameters"]
    data["stars"] = star_data

    return data


def check_posvel_datasets(datasets):
    checked_datasets = copy.deepcopy(datasets)
    for p in ["", "Obs"]:
        if "Pos" + p in checked_datasets:
            if "HCoordinates" + p not in checked_datasets:
                checked_datasets.append("HCoordinates" + p)
            checked_datasets.remove("Pos" + p)
        if "Vel" + p in checked_datasets:
            if "HVelocities" + p not in checked_datasets:
                checked_datasets.append("HVelocities" + p)
            checked_datasets.remove("Vel" + p)
    return checked_datasets


def add_posvel_data(data, datasets):
    checked_datasets = list(data.keys())
    for p in ["", "Obs"]:
        if "Pos" + p in datasets and "Vel" + p in datasets:
            print("Applying Cartesian_Conversion")
            data["Pos" + p], data["Vel" + p] = cc.galactic_cartesian_coordinates(data["HCoordinates" + p],
                                                                                 data["HVelocities" + p])
        elif "Pos" + p in datasets:
            print("Applying Cartesian_Conversion")
            data["Pos" +
                 p] = cc.galactic_cartesian_coordinates(data["HCoordinates" + p])
        for x in ["HVelocities", "HCoordinates"]:
            if x + p not in datasets and x + p in checked_datasets:
                del data[x + p]
    return data
