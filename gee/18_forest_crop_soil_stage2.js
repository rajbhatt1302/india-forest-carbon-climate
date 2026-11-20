// Stage 2: regional (100 and 250 km) and window values of the stage 1 predictors.
// Output: grid05_predictors_type.csv, grid025_predictors_type.csv
var A = 'projects/YOUR_PROJECT/assets/';
var india = ee.FeatureCollection(A + 'india_pov_boundary').geometry();
var P005 = ee.Projection('EPSG:4326', [0.05, 0, 68.0, 0, -0.05, 37.5]);
var src = ee.Image.cat([ee.Image(A + 'india_type_EVG_005deg'), ee.Image(A + 'india_type_FSEAS_005deg'),
  ee.Image(A + 'india_type_CROP_005deg'), ee.Image(A + 'india_type_NDVIK_005deg'), ee.Image(A + 'india_type_soil_005deg')]);
var inIndia = src.select('CLAY').mask();
var names = ['EVG', 'FSEAS', 'CNC', 'CGSL', 'CAMP', 'NDVI_kharif'];
var bands = [src.select(['CLAY', 'SAND'])];
names.forEach(function(nm){
  var s = src.select(nm + '_s').unmask(0), n = src.select(nm + '_n').unmask(0);
  [[18, '_R100'], [45, '_R250']].forEach(function(k){
    var ker = ee.Kernel.circle(k[0], 'pixels', false);
    var ns = n.convolve(ker);
    bands.push(s.convolve(ker).divide(ns).updateMask(ns.gt(0.5)).updateMask(inIndia).reproject(P005).rename(nm + k[1]));
  });
  bands.push(src.select([nm + '_s', nm + '_n']));
});
var img = ee.Image.cat(bands);
function grid(G){
  var xs = ee.List.sequence(68.0, 97.5 - G, G), ys = ee.List.sequence(6.5, 37.5 - G, G);
  return ee.FeatureCollection(xs.map(function(x){ return ys.map(function(y){
    x = ee.Number(x); y = ee.Number(y);
    var id = x.multiply(100).round().format('%d').cat('_').cat(y.multiply(100).round().format('%d'));
    return ee.Feature(ee.Geometry.Rectangle([x, y, x.add(G), y.add(G)], null, false), {gid: id});
  });}).flatten()).filterBounds(india);
}
var cols = ['gid'];
names.forEach(function(nm){ cols.push(nm + '_s', nm + '_n', nm + '_R100', nm + '_R250'); });
cols.push('CLAY', 'SAND');
[[0.5, 'grid05_predictors_type'], [0.25, 'grid025_predictors_type']].forEach(function(g){
  var t = img.reduceRegions({collection: grid(g[0]), reducer: ee.Reducer.mean(), crs: 'EPSG:4326',
    crsTransform: [0.05, 0, 68.0, 0, -0.05, 37.5], tileScale: 4}).map(function(f){ return f.setGeometry(null); });
  Export.table.toDrive({collection: t, description: g[1], folder: 'india_forest_climate', fileFormat: 'CSV', selectors: cols});
});
