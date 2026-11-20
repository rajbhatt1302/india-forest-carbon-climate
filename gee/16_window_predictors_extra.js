// Cropland NDVI (rabi, dry season), canopy height, soil moisture, wind and population per window.
// Output: grid05_predictors_extra.csv
var A = 'projects/YOUR_PROJECT/assets/';
var india = ee.FeatureCollection(A + 'india_pov_boundary').geometry();
var lc = ee.Image(A + 'india_lc_154m');
var crop = lc.select('fCrop').gte(50), nat = lc.select('fNat').gte(50), plant = lc.select('fPlant').gte(50);

var mod = ee.ImageCollection('MODIS/061/MOD13Q1').filterDate('2019-01-01', '2024-01-01').select('NDVI');
var rabi = mod.filter(ee.Filter.calendarRange(1, 3, 'month')).mean().multiply(0.0001).updateMask(crop).rename('NDVI_rabi');
var dry = mod.filter(ee.Filter.calendarRange(4, 5, 'month')).mean().multiply(0.0001).updateMask(crop).rename('NDVI_dry');
var ch = ee.Image('users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1').select('b1');
var chN = ch.updateMask(nat).rename('CH_nat'), chP = ch.updateMask(plant).rename('CH_plant');
var era = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR').filterDate('2019-01-01', '2024-01-01');
var sm = era.select('volumetric_soil_water_layer_1').mean().rename('SM');
var wind = era.map(function(i){ return i.select('u_component_of_wind_10m').hypot(i.select('v_component_of_wind_10m')); })
  .mean().rename('WIND');
var pop = ee.ImageCollection('JRC/GHSL/P2023A/GHS_POP').filterDate('2020-01-01', '2020-12-31').first().rename('POP');

var xs = ee.List.sequence(68.0, 97.0, 0.5), ys = ee.List.sequence(6.5, 37.0, 0.5);
var grid = ee.FeatureCollection(xs.map(function(x){ return ys.map(function(y){
  x = ee.Number(x); y = ee.Number(y);
  var id = x.multiply(100).round().format('%d').cat('_').cat(y.multiply(100).round().format('%d'));
  return ee.Feature(ee.Geometry.Rectangle([x, y, x.add(0.5), y.add(0.5)], null, false), {gid: id});
});}).flatten()).filterBounds(india);
var inIndia = ee.Image.constant(1).clip(india);

var m1 = ee.Image.cat([rabi, dry, chN, chP]).updateMask(inIndia).reduceRegions({collection: grid, reducer: ee.Reducer.mean(),
  crs: 'EPSG:4326', crsTransform: [1/720, 0, 68.0, 0, -1/720, 37.5], tileScale: 8});
var m2 = ee.Image.cat([sm, wind]).updateMask(inIndia).reduceRegions({collection: m1, reducer: ee.Reducer.mean(), scale: 5000,
  tileScale: 4});
var m3 = pop.updateMask(inIndia).reduceRegions({collection: m2, reducer: ee.Reducer.sum(), scale: 100, tileScale: 16})
  .map(function(f){ return f.set('POP', f.get('sum')).setGeometry(null); });
Export.table.toDrive({collection: m3, description: 'grid05_predictors_extra', folder: 'india_forest_climate', fileFormat: 'CSV',
  selectors: ['gid', 'NDVI_rabi', 'NDVI_dry', 'CH_nat', 'CH_plant', 'SM', 'WIND', 'POP']});
