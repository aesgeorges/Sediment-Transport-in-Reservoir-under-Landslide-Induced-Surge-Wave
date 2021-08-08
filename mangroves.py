import os
import gc
import earthpy as et
import rasterio as rio
import rioxarray as rxr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pylab as pl
import matplotlib.pyplot as plt

gc.collect()

mangrove_path = 'CMS_Global_Map_Mangrove_Canopy_1665/data/Mangrove_agb_Haiti.tif'
mangrove_canopy = rxr.open_rasterio(mangrove_path)
bounds = mangrove_canopy.rio.bounds()

lon_min = bounds[0]
lon_max = bounds[2]
lat_min = bounds[1]
lat_max = bounds[3]

f, ax = plt.subplots(figsize=(16,8))
cmap = pl.cm.viridis

ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
#ax.add_feature(cfeature.LAND)
#ax.add_feature(cfeature.OCEAN)
ax.set_extent((lon_min, lon_max, lat_min, lat_max))

mangrove_canopy.plot(cmap='Oranges')
plt.show()