![](../docs/banner.png)

# Chapter 2: Visualizing and modelling spatial data

**Tomas Beuzen, 2021**

![](img/space.png)

## Chapter Outline
<hr>

<div class="toc"><ul class="toc-item"><li><span><a href="#Chapter-Learning-Objectives" data-toc-modified-id="Chapter-Learning-Objectives-2">Chapter Learning Objectives</a></span></li><li><span><a href="#Imports" data-toc-modified-id="Imports-3">Imports</a></span></li><li><span><a href="#1.-Spatial-visualization" data-toc-modified-id="1.-Spatial-visualization-4">1. Spatial visualization</a></span></li><li><span><a href="#2.-Spatial-modelling" data-toc-modified-id="2.-Spatial-modelling-5">2. Spatial modelling</a></span></li></ul></div>

## Chapter Learning Objectives
<hr>

- Make informed choices about how to plot your spatial data, e.g., scattered, polygons, 3D, etc..
- Plot spatial data using libraries such as `geopandas`, `plotly`, and `keplergl`.
- Interpolate unobserved spatial data using deterministic methods such as nearest-neighbour interpolation.
- Interpolate data from one set of polygons to a different set of polygons using areal interpolation.

## Imports
<hr>

import warnings
import keplergl
import numpy as np
import osmnx as ox
import pandas as pd
import geopandas as gpd
import plotly.express as px
from skgstat import Variogram
import matplotlib.pyplot as plt
from shapely.geometry import Point
from pykrige.ok import OrdinaryKriging
from scipy.interpolate import NearestNDInterpolator
from tobler.area_weighted import area_interpolate
# Custom functions
from scripts.utils import pixel2poly
# Plotting defaults
plt.style.use('ggplot')
px.defaults.height = 400; px.defaults.width = 620
plt.rcParams.update({'font.size': 16, 'axes.labelweight': 'bold', 'figure.figsize': (6, 6), 'axes.grid': False})

## 1. Spatial visualization
<hr>

### 1.1. Geopandas

We saw last chapter how to easily plot geospatial data using the `geopandas` method `.plot()`. This workflow is useful for making quick plots, exploring your data, and easily layering geometries. Let's import some data of UBC buildings using `osmnx` (our Python API for accessing OpenStreetMap data) and make a quick plot:

ubc = (ox.geometries_from_place("University of British Columbia, Canada", tags={'building':True})
         .loc[:, ["geometry"]]                 # just keep the geometry column for now
         .query("geometry.type == 'Polygon'")  # only what polygons (buidling footprints)
         .assign(Label="Building Footprints")  # assign a label for later use
         .reset_index(drop=True)               # reset to 0 integer indexing
      )
ubc.head()

Recall that we can make a plot using the `.plot()` method on a `GeoDataFrame`:

ax = ubc.plot(figsize=(8, 8), column="Label", legend=True,
              edgecolor="0.2", markersize=200, cmap="rainbow")
plt.title("UBC");

Say I know the "point" location of my office but I want to locate the building footprint (a "polygon"). That's easily done with `geopandas`!

First, I'll use `shapely` (the Python geometry library `geopandas` is built on) to make my office point, but you could also use the `geopandas` function `gpd.points_from_xy()` like we did last chapter:

point_office = Point(-123.2522145, 49.2629555)
point_office

Now, I can use the `.contains()` method to find out which building footprint my office resides in:

ubc[ubc.contains(point_office)]

Looks like it's index 48! I'm going to change the label of that one to "Tom's Office":

ubc.loc[48, "Label"] = "Tom's Office"

Now let's make a plot!

ax = ubc.plot(figsize=(8, 8), column="Label", legend=True,
              edgecolor="0.2", markersize=200, cmap="rainbow")
plt.title("UBC");

We can add more detail to this map by including a background map. For this, we need to install the [contextily](https://github.com/geopandas/contextily) package. Note that most web providers use the Web Mercator projection, "EPSG:3857" ([interesting article on that here](https://www.esri.com/arcgis-blog/products/arcgis-pro/mapping/mercator-its-not-hip-to-be-square/)) so I'll convert to that before plotting:

import contextily as ctx

ax = (ubc.to_crs("EPSG:3857")
         .plot(figsize=(10, 8), column="Label", legend=True,
               edgecolor="0.2", markersize=200, cmap="rainbow")
     )
ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik)  # I'm using OSM as the source. See all provides with ctx.providers
plt.axis("off")
plt.title("UBC");

### 1.2. Plotly

The above map is nice, but it sure would be helpful in some cases to have interactive functionality for our map don't you think? Like the ability to zoom and pan? Well there are several packages out there that can help us with that, but `plotly` is one of the best for this kind of mapping. `plotly` supports plotting of maps backed by [MapBox](https://www.mapbox.com/) (a mapping and location data cloud platform).

Here's an example using the `plotly.express` function `px.choropleth_mapbox()`:

fig = px.choropleth_mapbox(ubc, geojson=ubc.geometry, locations=ubc.index, color="Label",
                           center={"lat": 49.261, "lon": -123.246}, zoom=12.5,
                           mapbox_style="open-street-map")
fig.update_layout(margin=dict(l=0, r=0, t=30, b=10))

You can pan and zoom above as you desire! How about we add a terrain map instead:

fig = px.choropleth_mapbox(ubc, geojson=ubc.geometry, locations=ubc.index, color="Label",
                           center={"lat": 49.261, "lon": -123.246}, zoom=12.5,
                           mapbox_style="stamen-terrain")
fig.update_layout(margin=dict(l=0, r=0, t=30, b=10))

We are just colouring our geometries based on a label at the moment, but we could of course do some cooler things. Let's colour by building area:

# Calculate area
ubc["Area"] = ubc.to_crs(epsg=3347).area  # note I'm projecting to EPSG:3347 (the projected system in meters that Statistics Canada useshttps://epsg.io/3347)

# Make plot
fig = px.choropleth_mapbox(ubc, geojson=ubc.geometry, locations=ubc.index, color="Area",
                           center={"lat": 49.261, "lon": -123.246}, zoom=12.5,
                           mapbox_style="carto-positron")
fig.update_layout(margin=dict(l=0, r=0, t=30, b=10))

Check out the `plotly` [documentation](https://plotly.com/python/maps/) for more - there are many plotting options and examples to learn from! Other popular map plotting options include [altair](https://altair-viz.github.io/gallery/index.html#maps) ([doesn't support interactivity yet](https://github.com/vega/vega-lite/issues/3306)), [folium](https://python-visualization.github.io/folium/), and [bokeh](https://docs.bokeh.org/en/latest/docs/user_guide/geo.html).

### 1.3. Kepler.gl

The above mapping was pretty cool, but are you ready for more power?

![](img/thanos.gif)

Time to introduce [kepler.gl](https://docs.kepler.gl/)! `keplergl` is a web-based tool for visualing spatial data. Luckily, it has a nice Python API and Jupyter extension for us to use (see the [install instructions](https://docs.kepler.gl/docs/keplergl-jupyter#install)). The basic way `keplergl` works is:
1. We create an instance of a map with `keplergl.KeplerGl()`
2. We add as much data to the map as we like with the `.add_data()` method
3. We customize and configure the map in any way we like using the GUI (graphical user interface)

ubc_map = keplergl.KeplerGl(height=500)
ubc_map.add_data(data=ubc.copy(), name="Building heights")
ubc_map

%%html
<iframe src="../_images/ubc-2d.html" width="80%" height="500"></iframe>

I'll do you one better than that! Let's add a 3D element to our plot! I'm going to load in some data of UBC building heights I downloaded from the [City of Vancouver Open Data Portal](https://opendata.vancouver.ca/explore/dataset/building-footprints-2009/information/):

ubc_bldg_heights = gpd.read_file("data-spatial/ubc-building-footprints-2009")
ubc_bldg_heights.head()

I'm going to combine this with our `ubc` data and do a bit of clean-up and wrangling. You'll do this workflow in your lab too so I won't spend too much time here:

ubc_bldg_heights = (gpd.sjoin(ubc, ubc_bldg_heights[["hgt_agl", "geometry"]], how="inner")
                       .drop(columns="index_right")
                       .rename(columns={"hgt_agl": "Height"})
                       .reset_index()
                       .dissolve(by="index", aggfunc="mean")  # dissolve is like "groupby" in pandas. We use it because it retains geometry information
                   )

Now I'll make a new map and configure it to be in 3D using the GUI!

ubc_height_map = keplergl.KeplerGl(height=500)
ubc_height_map.add_data(data=ubc_bldg_heights.copy(), name="Building heights")
ubc_height_map

%%html
<iframe src="../_images/ubc-3d.html" width="80%" height="500"></iframe>

You can save your configuration and customization for re-use later ([see the docs here](https://docs.kepler.gl/docs/keplergl-jupyter#5-save-and-load-config))

## 2. Spatial modelling
<hr>

There’s typically two main ways we might want to "model" spatial data:
1. Spatial interpolation: use a set of observations in space to estimate the value of a spatial field
2. Areal interpolation: project data from one set of polygons to another set of polygons

Both are based on the fundamental premise: “*everything is related to everything else, but near things are more related than distant things.*” ([Tobler's first law of geography](https://en.wikipedia.org/wiki/Tobler%27s_first_law_of_geography)). To demonstrate some of these methods, we'll look at the annual average air pollution (PM 2.5) recorded at stations across BC during 2020 (which I downloaded from [DataBC](https://catalogue.data.gov.bc.ca/dataset/air-quality-monitoring-verified-hourly-data/resource/54492296-b999-4c53-b42c-21927d025e5a)):

pm25 = pd.read_csv("data/bc-pm25.csv")
pm25.head()

fig = px.scatter_mapbox(pm25, lat="Lat", lon="Lon", color="PM_25", size="PM_25",
                        color_continuous_scale="RdYlGn_r",
                        center={"lat": 52.261, "lon": -123.246}, zoom=3.5,
                        mapbox_style="carto-positron", hover_name="Station Name")
fig.update_layout(margin=dict(l=0, r=0, t=30, b=10))
fig.show()

The goal is to interpolate these point measurements to estimate the air pollution for all of BC.

### 2.1. Deterministic spatial interpolation

Create surfaces directly from measured points using a mathematical function. Common techniques (all available in the `scipy` module `interpolate`) are:
- Inverse distance weighted interpolation
- Nearest neighbour interpolation
- Polynomial interpolation
- Radial basis function (RBF) interpolation

Let's try nearest neighbour interpolation now using `scipy.interpolate.NearestNDInterpolator`. As angular coordinates (lat/lon) are not good for measuring distances, I'm going to first convert my data to the linear, meter-based Lambert projection recommend by Statistics Canada and extract the `x` and `y` locations as columns in my `GeoDataFrame` ("Easting" and "Northing" respectively):

gpm25 = (gpd.GeoDataFrame(pm25, crs="EPSG:4326", geometry=gpd.points_from_xy(pm25["Lon"], pm25["Lat"]))
            .to_crs("EPSG:3347")
        )
gpm25["Easting"], gpm25["Northing"] = gpm25.geometry.x, gpm25.geometry.y
gpm25.head()

Now let's create a grid of values (a raster) to interpolate over. I'm just going to make a square grid of fixed resolution that spans the bounds of my observed data points (we'll plot this shortly so you can see what it looks like):

resolution = 25_000  # cell size in meters
gridx = np.arange(gpm25.bounds.minx.min(), gpm25.bounds.maxx.max(), resolution)
gridy = np.arange(gpm25.bounds.miny.min(), gpm25.bounds.maxy.max(), resolution)

So now let's interpolate using a nearest neighbour method (the code below is straight from the [scipy docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.NearestNDInterpolator.html)):

model = NearestNDInterpolator(x = list(zip(gpm25["Easting"], gpm25["Northing"])),
                              y = gpm25["PM_25"])
z = model(*np.meshgrid(gridx, gridy))
plt.imshow(z);

Okay so it looks like we successfully interpolated our made-up grid, but let's re-project it back to our original map. To do this I need to convert each cell in my raster to a small polygon using a simple function I wrote called `pixel2poly()` which we imported at the beginning of the notebook:

polygons, values = pixel2poly(gridx, gridy, z, resolution)

Now we can convert that to a `GeoDataFrame` and plot using `plotly`:

pm25_model = (gpd.GeoDataFrame({"PM_25_modelled": values}, geometry=polygons, crs="EPSG:3347")
                 .to_crs("EPSG:4326")
             )

fig = px.choropleth_mapbox(pm25_model, geojson=pm25_model.geometry, locations=pm25_model.index,
                           color="PM_25_modelled", color_continuous_scale="RdYlGn_r", opacity=0.5,
                           center={"lat": 52.261, "lon": -123.246}, zoom=3.5,
                           mapbox_style="carto-positron")
fig.update_layout(margin=dict(l=0, r=0, t=30, b=10))
fig.update_traces(marker_line_width=0)

Very neat! Also, notice how the grid distorts a bit once we project it back onto an angular (degrees-based) projection from our linear (meter-based) projection.

### 2.2. Probabilistic (geostatistical) spatial interpolation

Geostatistical interpolation (called "Kriging") differs to deterministic interpolation in that we interpolate using statistical models that include estimates of spatial autocorrelation. There's a lot to read about [kriging](https://desktop.arcgis.com/en/arcmap/10.3/tools/3d-analyst-toolbox/how-kriging-works.htm), I just want you to be aware of the concept in case you need to do something like this in the future. Usually I'd do this kind of interpolation in GIS software like ArcGIS or QGIS, but it's possible to do in Python too.

The basic idea is that if we have a set of observations $Z(s)$ at locations $s$, we estimate the value of an unobserved location ($s_0$) as a weighted sum:
$$\hat{Z}(s_0)=\sum_{i=0}^{N}\lambda_iZ(s_i)$$

Where $N$ is the size of $s$ (number of observed samples) and $\lambda$ is an array of weights. The key is deciding which weights $\lambda$ to use. Kriging uses the spatial autocorrelation in the data to determine the weights. Spatial autocorrelation is calculated by looking at the squared difference (the variance) between points at similar distances apart, let's see what that means:

warnings.filterwarnings("ignore")  # Silence some warnings
vario = Variogram(coordinates=gpm25[["Easting", "Northing"]],
                  values=gpm25["PM_25"],
                  n_lags=20)
vario.distance_difference_plot();

The idea is to then fit a model to this data that describes how variance (spatial autocorrelation) changes with distance ("lag") between locations. We look at the average variance in bins of the above distances/pairs and fit a line through them. This model is called a "variogram" and it's analogous to the autocorrelation function for time series. It defines the variance (autocorrelation structure) as a function of distance:

vario.plot(hist=False);

The above plot basically shows that:
- at small distances (points are close together), variance is reduced because points are correlated
- but at a distance around 400,000 m the variance flattens out indicating points are two far away to have any impactful spatial autocorrelation. This location is called the "range" while the variance at the "range" is called the "sill" (like a ceiling). We can extract the exact range:

vario.describe()["effective_range"]

By the way, we call it the semi-variance because there is a factor of 0.5 in the equation to account for the fact that variance is calculated twice for each pair of points ([read more here](https://www.aspexit.com/en/variogram-and-spatial-autocorrelation/)).

Remember, our variogram defines the spatial autocorrelation of the data (i.e., how the locations in our region affect one another). Once we have a variogram model, we can use it to estimate the weights in our kriging model. I won't go into detail on how this is done, but there is a neat walkthrough in the [scikit-gstat docs here](https://scikit-gstat.readthedocs.io/en/latest/userguide/kriging.html).

Anyway, I'll briefly use the [pykrige](https://github.com/GeoStat-Framework/PyKrige) library to do some kriging so you can get an idea of what it looks like:

krig = OrdinaryKriging(x=gpm25["Easting"], y=gpm25["Northing"], z=gpm25["PM_25"], variogram_model="spherical")
z, ss = krig.execute("grid", gridx, gridy)
plt.imshow(z);

Now let's convert our raster back to polygons so we can map it. I'm also going to load in a polygon of BC using `osmnx` to clip my data so it fits nicely on my map this time:

polygons, values = pixel2poly(gridx, gridy, z, resolution)
pm25_model = (gpd.GeoDataFrame({"PM_25_modelled": values}, geometry=polygons, crs="EPSG:3347")
                 .to_crs("EPSG:4326")
                 )
bc = ox.geocode_to_gdf("British Columbia, Canada")
pm25_model = gpd.clip(pm25_model, bc)

fig = px.choropleth_mapbox(pm25_model, geojson=pm25_model.geometry, locations=pm25_model.index,
                           color="PM_25_modelled", color_continuous_scale="RdYlGn_r",
                           center={"lat": 52.261, "lon": -123.246}, zoom=3.5,
                           mapbox_style="carto-positron")
fig.update_layout(margin=dict(l=0, r=0, t=30, b=10))
fig.update_traces(marker_line_width=0)

I used an "ordinary kriging" interpolation above which is the simplest implementation of kriging. The are many other forms of kriging too that can account for underlying trends in the data ("universal kriging"), or even use a regression or classification model to make use of additional explanatory variables. `pykrige` [supports most variations](https://geostat-framework.readthedocs.io/projects/pykrige/en/stable/examples/index.html). In particular for the latter, `pykrige` can accept `sklearn` models which is useful!

### 2.3. Areal interpolation

Areal interpolation is concerned with mapping data from one polygonal representation to another. Imagine I want to map the air pollution polygons I just made to FSA polygons (recall FSA is "forward sortation area", which are groups of postcodes). The most intuitive way to do this is to distribute values based on area proportions, hence "areal interpolation".

I'll use the [tobler](https://github.com/pysal/tobler) library for this. First, load in the FSA polygons:

van_fsa = gpd.read_file("data-spatial/van-fsa")
ax = van_fsa.plot(edgecolor="0.2")
plt.title("Vancouver FSA");

Now I'm just going to made a higher resolution interpolation using kriging so we can see some of the details on an FSA scale:

resolution = 10_000  # cell size in meters
gridx = np.arange(gpm25.bounds.minx.min(), gpm25.bounds.maxx.max(), resolution)
gridy = np.arange(gpm25.bounds.miny.min(), gpm25.bounds.maxy.max(), resolution)
krig = OrdinaryKriging(x=gpm25["Easting"], y=gpm25["Northing"], z=gpm25["PM_25"], variogram_model="spherical")
z, ss = krig.execute("grid", gridx, gridy)
polygons, values = pixel2poly(gridx, gridy, z, resolution)
pm25_model = (gpd.GeoDataFrame({"PM_25_modelled": values}, geometry=polygons, crs="EPSG:3347")
                 .to_crs("EPSG:4326")
                 )

Now we can easily do the areal interpolation using the function `area_interpolate()`:

areal_interp = area_interpolate(pm25_model.to_crs("EPSG:3347"),
                                van_fsa.to_crs("EPSG:3347"),
                                intensive_variables=["PM_25_modelled"]).to_crs("EPSG:4326")
areal_interp.plot(column="PM_25_modelled", figsize=(8, 8),
                  edgecolor="0.2", cmap="RdBu", legend=True)
plt.title("FSA Air Pollution");

There are other methods you can use for areal interpolation too, that include additional variables or use more advanced interpolation algorithms. The [tobbler documentation](https://pysal.org/tobler/) describes some of these.

![](img/bye.png)