import numpy as np
from netCDF4 import Dataset

# ----------------------------------------------------------------------------------------------------------------------

#                                             file creation

# ----------------------------------------------------------------------------------------------------------------------
Vario_ncdf = Dataset("Vario_uncor_ncdf.nc", "w", format="NETCDF4")

# ----------------------------------------------------------------------------------------------------------------------

#                                             group creation

# ----------------------------------------------------------------------------------------------------------------------

uncor_data_grp = Vario_ncdf.createGroup("Uncor_data")

bins_incor_grp = Vario_ncdf.createGroup("Uncor_data/bins_uncor data")
expe_data_uncor_grp = Vario_ncdf.createGroup("Uncor_data/experimental_data_uncor_group")
res_data_uncor_grp = Vario_ncdf.createGroup("Uncor_data/residual_uncor_data_group")
param_uncor_grp = Vario_ncdf.createGroup("Uncor_data/parameters_uncor")


cor_data_grp = Vario_ncdf.createGroup("cor_data")

bins_cor_grp = Vario_ncdf.createGroup("cor_data/bins_cor data")
expe_data_cor_grp = Vario_ncdf.createGroup("cor_data/experimental_data_cor_group")
res_data_cor_grp = Vario_ncdf.createGroup("cor_data/residual_cor_data_group")
param_cor_grp = Vario_ncdf.createGroup("cor_data/parameters_cor")

# ----------------------------------------------------------------------------------------------------------------------

#                                             dimension creation

# ----------------------------------------------------------------------------------------------------------------------

bins_cor = bins_incor_grp.createDimension()

