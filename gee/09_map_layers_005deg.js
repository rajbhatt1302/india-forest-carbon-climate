// Land cover, forest change, LST, albedo, ET and biomass at 0.05 deg.
var A = 'projects/YOUR_PROJECT/assets/';
var india = ee.FeatureCollection(A + 'india_pov_boundary').geometry();
var CT = [0.05, 0, 68.0, 0, -0.05, 37.5];
var region = ee.Geometry.Rectangle([68.0, 6.5, 97.5, 37.5], null, false);

var lc = ee.Image(A + 'india_lc_154m')
  .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 2048})
  .round().toUint8().clip(india);
Export.image.toDrive({image: lc, description: 'india_lc_005deg', folder: 'india_forest_climate', fileNamePrefix: 'india_lc_005deg',
  crs: 'EPSG:4326', crsTransform: CT, region: region, maxPixels: 1e9});

var cl = ee.Image(A + 'india_1km_climate_v2');
var bio = ee.Image.cat([
  cl.select(['LSTd','LSTn','LSTd_MAM','LSTn_MAM','LSTd_JJAS','LSTd_ONDJF']).multiply(100),
  cl.select('ALB').multiply(10000), cl.select('ET'), cl.select('AGB').multiply(10)
]).reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 1024, bestEffort: true})
  .round().toInt16().unmask(-32768).clip(india);
Export.image.toDrive({image: bio, description: 'india_bio_005deg', folder: 'india_forest_climate', fileNamePrefix: 'india_bio_005deg',
  crs: 'EPSG:4326', crsTransform: CT, region: region, maxPixels: 1e9});
