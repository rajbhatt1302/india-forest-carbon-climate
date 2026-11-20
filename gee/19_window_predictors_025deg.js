// Predictors of 15_window_predictors.js on the 0.25 deg grid. Output: grid025_predictors.csv
var A = 'projects/YOUR_PROJECT/assets/';
var india = ee.FeatureCollection(A + 'india_pov_boundary').geometry();

// TerraClimate 1991-2020 climatology
var tc = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE').filterDate('1991-01-01', '2021-01-01');
var P = tc.select('pr').sum().divide(30).rename('P');                      // mm/yr
var PET = tc.select('pet').sum().multiply(0.1).divide(30).rename('PET');   // mm/yr
var AI = P.divide(PET).rename('AI');
var Tmax = tc.select('tmmx').mean().multiply(0.1).rename('Tmax');          // deg C
var Tmin = tc.select('tmmn').mean().multiply(0.1).rename('Tmin');
var VPD = tc.select('vpd').mean().multiply(0.01).rename('VPD');            // kPa
// share of annual rain in JJAS
var jjas = tc.filter(ee.Filter.calendarRange(6, 9, 'month')).select('pr').sum().divide(30);
var Pjjas = jjas.divide(P.multiply(1)).rename('Pjjas');

// ERA5-Land downward shortwave, 2019-2023 mean (W/m2)
var era = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR').filterDate('2019-01-01', '2024-01-01')
  .select('surface_solar_radiation_downwards_sum');
var SW = era.sum().divide(5 * 365.25 * 86400).rename('SWdown');

// terrain
var dem = ee.Image('USGS/SRTMGL1_003');
var elev = dem.rename('elev');
var slope = ee.Terrain.slope(dem).rename('slope');

// background surface state
var clim = ee.Image(A + 'india_1km_climate_v2').select(['LSTd', 'LSTn', 'ALB', 'ET', 'AGB'],
  ['LSTd_bg', 'LSTn_bg', 'ALB_bg', 'ET_bg', 'AGB_bg']);

var img = ee.Image.cat([P, PET, AI, Tmax, Tmin, VPD, Pjjas, SW, elev, slope, clim]).updateMask(
  ee.Image.constant(1).clip(india));

var G = 0.25;
var xs = ee.List.sequence(68.0, 97.5 - G, G), ys = ee.List.sequence(6.5, 37.5 - G, G);
var grid = ee.FeatureCollection(xs.map(function(x){ return ys.map(function(y){
  x = ee.Number(x); y = ee.Number(y);
  var id = x.multiply(100).round().format('%d').cat('_').cat(y.multiply(100).round().format('%d'));
  return ee.Feature(ee.Geometry.Rectangle([x, y, x.add(G), y.add(G)], null, false), {gid: id});
});}).flatten()).filterBounds(india);

var out = img.reduceRegions({collection: grid, reducer: ee.Reducer.mean(), scale: 2000, tileScale: 4})
  .map(function(f){ return f.setGeometry(null); });
var elevSD = elev.rename('elevSD').reduceRegions({collection: grid, reducer: ee.Reducer.stdDev(), scale: 1000, tileScale: 4})
  .map(function(f){ return ee.Feature(null, {gid: f.get('gid'), elevSD: f.get('stdDev')}); });
var joined = ee.Join.inner().apply(out, elevSD, ee.Filter.equals({leftField: 'gid', rightField: 'gid'}))
  .map(function(p){ return ee.Feature(p.get('primary')).set('elevSD', ee.Feature(p.get('secondary')).get('elevSD')); });

Export.table.toDrive({collection: joined, description: 'grid025_predictors', folder: 'india_forest_climate', fileFormat: 'CSV',
  selectors: ['gid', 'P', 'PET', 'AI', 'Tmax', 'Tmin', 'VPD', 'Pjjas', 'SWdown', 'elev', 'elevSD', 'slope',
              'LSTd_bg', 'LSTn_bg', 'ALB_bg', 'ET_bg', 'AGB_bg']});
