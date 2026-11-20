// WorldCover 2021 class fractions at 0.05 deg, trees split as in 02_landcover_154m.js.
var A = 'projects/YOUR_PROJECT/assets/';
var india = ee.FeatureCollection(A + 'india_pov_boundary').geometry();
var wc = ee.ImageCollection('ESA/WorldCover/v200').first();
var nf = ee.ImageCollection('projects/nature-trace/assets/forest_typology/natural_forest_2020_v1_0_collection')
  .mosaic().unmask(0);                       // probability 0-250; 125 = 0.5
var tree = wc.eq(10);
var ind = ee.Image.cat([tree.and(nf.gte(125)), tree.and(nf.lt(125)), wc.eq(20), wc.eq(30), wc.eq(40), wc.eq(50), wc.eq(60),
  wc.eq(70), wc.eq(80), wc.eq(90), wc.eq(95), wc.eq(100)])
  .rename(['nat', 'plant', 'shrub', 'grass', 'crop', 'built', 'bare', 'snow', 'water', 'wetland', 'mangrove', 'moss']).toFloat();
var P150 = ee.Projection('EPSG:4326').atScale(150);
var P005 = ee.Projection('EPSG:4326', [0.05, 0, 68.0, 0, -0.05, 37.5]);
var f = ind.reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 512}).reproject(P150)
  .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 2048}).reproject(P005)
  .multiply(100).round().toUint8().clip(india);
Export.image.toDrive({image: f, description: 'india_worldcover_classes_005deg', folder: 'india_forest_climate',
  fileNamePrefix: 'india_worldcover_classes_005deg', crs: 'EPSG:4326', crsTransform: [0.05, 0, 68.0, 0, -0.05, 37.5],
  region: ee.Geometry.Rectangle([68.0, 6.5, 97.5, 37.5], null, false), maxPixels: 1e9, shardSize: 64});
