// Forest change at 0.05 deg and summed per 0.5 deg window (km2).
var A = 'projects/YOUR_PROJECT/assets/';
var india = ee.FeatureCollection(A + 'india_pov_boundary').geometry();
var hz = ee.Image(A + 'india_hansen_154m_v113').divide(100);

// 0.05 deg
var m = hz.reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 2048}).multiply(100).round().toUint8().clip(india);
Export.image.toDrive({image: m, description: 'india_hansen_v113_005deg', folder: 'india_forest_climate', fileNamePrefix: 'india_hansen_v113_005deg',
  crs: 'EPSG:4326', crsTransform: [0.05, 0, 68.0, 0, -0.05, 37.5],
  region: ee.Geometry.Rectangle([68.0, 6.5, 97.5, 37.5], null, false), maxPixels: 1e9});

// window sums (km2)
var xs = ee.List.sequence(68.0, 97.0, 0.5), ys = ee.List.sequence(6.5, 37.0, 0.5);
var grid = ee.FeatureCollection(xs.map(function(x){ return ys.map(function(y){
  x = ee.Number(x); y = ee.Number(y);
  var id = x.multiply(100).round().format('%d').cat('_').cat(y.multiply(100).round().format('%d'));
  return ee.Feature(ee.Geometry.Rectangle([x, y, x.add(0.5), y.add(0.5)], null, false), {gid: id});
});}).flatten()).filterBounds(india);
var km2 = hz.multiply(ee.Image.pixelArea().divide(1e6));
var sums = km2.reduceRegions({collection: grid, reducer: ee.Reducer.sum(), crs: 'EPSG:4326',
  crsTransform: [1/720, 0, 68.0, 0, -1/720, 37.5], tileScale: 4}).map(function(f){ return f.setGeometry(null); });
Export.table.toDrive({collection: sums, description: 'grid05_hansen_v113', folder: 'india_forest_climate', fileFormat: 'CSV',
  selectors: ['gid', 'hF2000', 'hLoss0125', 'hLoss0110', 'hLoss1120', 'hLossNat2125', 'hLossOth2125', 'hGain', 'hGainNonNat']});
