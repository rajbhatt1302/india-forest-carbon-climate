// Window regressions on the 0.25 deg grid. Output: grid025_window_results.csv
var A = 'projects/YOUR_PROJECT/assets/';
var india = ee.FeatureCollection(A + 'india_pov_boundary').geometry();
// 154 m fractions -> 1/120 deg
var lc = ee.Image(A + 'india_lc_154m').divide(100)
  .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 64})
  .addBands(ee.Image('USGS/SRTMGL1_003').rename('elev').toFloat()
    .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 65535, bestEffort: true}));
var cl = ee.Image(A + 'india_1km_climate_v2');
var soc = ee.Image('projects/soilgrids-isric/ocs_mean').select('ocs_0-30cm_mean').rename('SOC').toFloat()
  .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 64});     // 250 m -> 1/120 deg
var CRS = 'EPSG:4326', CT = [0.008333333333333333, 0, 68.0, 0, -0.008333333333333333, 37.5];

// grid
var G = 0.25;
var xs = ee.List.sequence(68.0, 97.5 - G, G), ys = ee.List.sequence(6.5, 37.5 - G, G);
var grid = ee.FeatureCollection(xs.map(function(x){ return ys.map(function(y){
  x = ee.Number(x); y = ee.Number(y);
  var id = x.multiply(100).round().format('%d').cat('_').cat(y.multiply(100).round().format('%d'));
  return ee.Feature(ee.Geometry.Rectangle([x, y, x.add(G), y.add(G)], null, false), {lon: x.add(G / 2), lat: y.add(G / 2), gid: id});
});}).flatten()).filterBounds(india).map(function(f){
  return f.set('fIndia', f.geometry().intersection(india, 1000).area(1000).divide(f.geometry().area(1000)));
});

// masks
var inIndia = ee.Image.constant(1).clip(india).mask();
var fNat = lc.select('fNat'), fPlant = lc.select('fPlant'), fCrop = lc.select('fCrop');
var clean = fNat.add(fPlant).add(fCrop).gte(0.85).and(lc.select('fBuilt').lt(0.05)).and(lc.select('fWater').lt(0.05)).and(inIndia);
var X = ee.Image.constant(1).rename('c').addBands(fNat).addBands(fPlant).addBands(lc.select('elev').divide(1000).rename('elevkm'));

function regress(fc, Yimg, yNames, tag){
  var k = yNames.length;
  var img = X.addBands(Yimg.select(yNames)).updateMask(clean);
  var valid = img.mask().reduce(ee.Reducer.min());
  img = img.updateMask(valid);
  var out = img.reduceRegions({collection: fc, reducer: ee.Reducer.linearRegression(4, k), crs: CRS, crsTransform: CT, tileScale: 4});
  var cov = img.select(['fNat','fPlant','elevkm']).toArray().reduceRegions({collection: out,
    reducer: ee.Reducer.covariance().unweighted(), crs: CRS, crsTransform: CT, tileScale: 4});
  var cnt = img.select('c').rename('nreg').reduceRegions({collection: cov, reducer: ee.Reducer.count(), crs: CRS, crsTransform: CT, tileScale: 4});
  return cnt.map(function(f){
    var has = f.get('coefficients');
    var coef = ee.Array(ee.Algorithms.If(has, f.get('coefficients'), ee.Array(ee.List.repeat(ee.List.repeat(-9999, k), 4))));
    var res = ee.Array(ee.Algorithms.If(has, f.get('residuals'), ee.Array(ee.List.repeat(-9999, k))));
    var n = ee.Number(f.get('count'));
    var C3 = ee.Array(ee.Algorithms.If(f.get('covariance'), f.get('covariance'), ee.Array([[1,0,0],[0,1,0],[0,0,1]])));
    var det = C3.matrixDeterminant();
    // slope standard errors
    var inv = ee.Array(ee.Algorithms.If(n.gt(10).and(det.abs().gt(1e-12)),
      C3.multiply(n.subtract(1)).matrixInverse(), ee.Array([[1e9,0,0],[0,1e9,0],[0,0,1e9]])));
    var d = {}; d['n_' + tag] = n;
    yNames.forEach(function(nm, j){
      var rms = res.get([j]);
      d['bN_' + nm] = coef.get([1, j]);
      d['bP_' + nm] = coef.get([2, j]);
      d['sN_' + nm] = ee.Number(rms).multiply(ee.Number(inv.get([0, 0])).sqrt());
      d['sP_' + nm] = ee.Number(rms).multiply(ee.Number(inv.get([1, 1])).sqrt());
      // standard error of bP - bN
      d['sD_' + nm] = ee.Number(rms).multiply(ee.Number(inv.get([0, 0])).add(ee.Number(inv.get([1, 1])))
        .subtract(ee.Number(inv.get([0, 1])).multiply(2)).max(0).sqrt());
    });
    return ee.Feature(f.geometry(), d).copyProperties(f, ['gid','lon','lat','fIndia']);
  });
}

var r1 = regress(grid, cl, ['LSTd','LSTn','ALB','ET'], 'r1');
var r2 = regress(grid, cl, ['LSTd_MAM','LSTn_MAM','LSTd_JJAS','LSTd_ONDJF'], 'r2');
var r3 = regress(grid, cl, ['AGB'], 'r3');
var r4 = regress(grid, soc, ['SOC'], 'r4');
function carry(src, dst){
  var j = ee.Join.saveFirst('m').apply(dst, src, ee.Filter.equals({leftField: 'gid', rightField: 'gid'}));
  return ee.FeatureCollection(j).map(function(f){ return ee.Feature(f).copyProperties(ee.Feature(f.get('m'))).set('m', null); });
}
var allR = carry(r1, carry(r2, carry(r3, r4)));
var lcm = ee.Image.cat([fNat, fPlant, fCrop]).updateMask(inIndia);
var cnt = ee.Image.cat([fNat.gt(0.5).rename('nNat'), fPlant.gt(0.5).rename('nPlant'), fCrop.gt(0.5).rename('nCrop'),
                        clean.rename('nClean')]).updateMask(clean);
var w1 = lcm.reduceRegions({collection: allR, reducer: ee.Reducer.mean(), crs: CRS, crsTransform: CT, tileScale: 4});
var w2 = cnt.reduceRegions({collection: w1, reducer: ee.Reducer.sum(), crs: CRS, crsTransform: CT, tileScale: 4});
Export.table.toDrive({collection: w2.map(function(f){ return f.setGeometry(null); }), description: 'grid025_window_results',
  folder: 'india_forest_climate', fileFormat: 'CSV'});
