// CMIP6 annual temperature and aridity per 0.5 deg window (area-weighted mean of the 0.25 deg cells).
var A = 'projects/YOUR_PROJECT/assets/';
var india = ee.FeatureCollection(A + 'india_pov_boundary').geometry();
var xs = ee.List.sequence(68.0, 97.0, 0.5), ys = ee.List.sequence(6.5, 37.0, 0.5);
var grid = ee.FeatureCollection(xs.map(function(x){ return ys.map(function(y){
  x = ee.Number(x); y = ee.Number(y);
  var id = x.multiply(100).round().format('%d').cat('_').cat(y.multiply(100).round().format('%d'));
  return ee.Feature(ee.Geometry.Rectangle([x, y, x.add(0.5), y.add(0.5)], null, false), {gid: id});
});}).flatten()).filterBounds(india);
// ACCESS_CM2 and EC_Earth3: one asset per model (05_cmip6_annual.js)
var SINGLE = ['ACCESS_CM2', 'EC_Earth3'];
var MODELS = SINGLE.concat(['GFDL_ESM4', 'INM_CM5_0', 'IPSL_CM6A_LR', 'MIROC6', 'MPI_ESM1_2_HR', 'MRI_ESM2_0',
                            'NorESM2_MM', 'CMCC_ESM2']);
MODELS.forEach(function(m){
  ['hist', 'ssp245', 'ssp585'].forEach(function(s){
    var asset = SINGLE.indexOf(m) >= 0 ? 'cmip6_annual_' + m : 'cmip6_annual_' + m + '_' + s;
    var img = ee.Image(A + asset).select(['T_' + s + '_.*', 'AI_' + s + '_.*']).clip(india);
    var t = img.reduceRegions({collection: grid, reducer: ee.Reducer.mean(), crs: 'EPSG:4326',
      crsTransform: [0.25, 0, 68.0, 0, -0.25, 37.5], tileScale: 2}).map(function(f){ return f.setGeometry(null); });
    Export.table.toDrive({collection: t, description: m + '_' + s, folder: 'india_forest_climate', fileFormat: 'CSV'});
  });
});
