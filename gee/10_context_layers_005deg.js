// Aridity index (TerraClimate 1991-2020) and elevation (GMTED2010) at 0.05 deg.
var india = ee.FeatureCollection('projects/YOUR_PROJECT/assets/india_pov_boundary').geometry();
var tc = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE').filterDate('1991-01-01', '2021-01-01');
var P = tc.select('pr').sum().divide(30);                // mm/yr
var PET = tc.select('pet').sum().multiply(0.1).divide(30);   // pet scale factor 0.1 -> mm/yr
var AI = P.divide(PET).multiply(1000).rename('AI_x1000');
var elev = ee.Image('USGS/GMTED2010_FULL').select('mea').rename('elev_m');
var img = AI.addBands(elev.reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 1024}))
  .round().toInt16().clip(india);
Export.image.toDrive({image: img, description: 'india_context_005deg', folder: 'india_forest_climate',
  fileNamePrefix: 'india_context_005deg', crs: 'EPSG:4326', crsTransform: [0.05, 0, 68.0, 0, -0.05, 37.5],
  region: ee.Geometry.Rectangle([68.0, 6.5, 97.5, 37.5], null, false), maxPixels: 1e9, priority: 500});
