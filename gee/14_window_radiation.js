// ERA5-Land downward shortwave per window, 2019-2023. Output: grid05_swdown.csv
var A = 'projects/YOUR_PROJECT/assets/';
var india = ee.FeatureCollection(A + 'india_pov_boundary').geometry();
var sw = ee.Image(A + 'era5land_swdown_2019_2023');
var xs = ee.List.sequence(68.0, 97.0, 0.5), ys = ee.List.sequence(6.5, 37.0, 0.5);
var grid = ee.FeatureCollection(xs.map(function(x){ return ys.map(function(y){
  x = ee.Number(x); y = ee.Number(y);
  var id = x.multiply(100).round().format('%d').cat('_').cat(y.multiply(100).round().format('%d'));
  return ee.Feature(ee.Geometry.Rectangle([x, y, x.add(0.5), y.add(0.5)], null, false), {gid: id});
});}).flatten()).filterBounds(india);
var out = sw.rename('SWdown').reduceRegions({collection: grid, reducer: ee.Reducer.mean(), scale: 11132})
  .map(function(f){ return f.setGeometry(null); });
Export.table.toDrive({collection: out, description: 'grid05_swdown', folder: 'india_forest_climate', fileFormat: 'CSV',
  selectors: ['gid', 'SWdown']});
