// As 05_cmip6_annual.js, one asset per model and scenario.
var C = ee.ImageCollection('NASA/GDDP-CMIP6');
var MODELS = ['ACCESS-CM2','EC-Earth3','GFDL-ESM4','INM-CM5-0','IPSL-CM6A-LR','MIROC6','MPI-ESM1-2-HR','MRI-ESM2-0','NorESM2-MM','CMCC-ESM2'];
var lat = ee.Image.pixelLonLat().select('latitude').multiply(Math.PI / 180);
var region = ee.Geometry.Rectangle([68.0, 6.5, 97.5, 37.5], null, false);
var CT = [0.25, 0, 68.0, 0, -0.25, 37.5];
var RUNS = [['historical', 'hist', 1995, 2014], ['ssp245', 'ssp245', 2015, 2099], ['ssp585', 'ssp585', 2015, 2099]];

function daily(i){
  var doy = ee.Date(i.get('system:time_start')).getRelative('day', 'year').add(1);
  var ang = doy.multiply(2 * Math.PI / 365);
  var dr = ang.cos().multiply(0.033).add(1), dec = ang.subtract(1.39).sin().multiply(0.409);
  var ws = lat.tan().multiply(dec.tan().multiply(-1)).clamp(-1, 1).acos();
  var ra = ws.multiply(lat.sin()).multiply(dec.sin()).add(lat.cos().multiply(dec.cos()).multiply(ws.sin()))
    .multiply(dr.multiply(24 * 60 / Math.PI * 0.0820 * 0.408));          // extraterrestrial radiation, mm/day (FAO-56)
  var tx = i.select('tasmax').subtract(273.15), tn = i.select('tasmin').subtract(273.15), tm = tx.add(tn).divide(2);
  var pet = ra.multiply(0.0023).multiply(tm.add(17.8)).multiply(tx.subtract(tn).max(0).sqrt());
  return i.select('pr').multiply(86400).rename('P').addBands(pet.rename('PET')).addBands(tm.rename('T'));
}

function runStack(model, r){
  var scen = r[0], tag = r[1];
  var cm = C.filter(ee.Filter.eq('model', model)).filter(ee.Filter.eq('scenario', scen)).select(['pr','tasmax','tasmin']);
  return ee.ImageCollection(ee.List.sequence(r[2], r[3]).map(function(y){
    y = ee.Number(y);
    var c = cm.filterDate(ee.Date.fromYMD(y, 1, 1), ee.Date.fromYMD(y.add(1), 1, 1)).map(daily);
    var P = c.select('P').sum(), PET = c.select('PET').sum(), T = c.select('T').mean();
    var ys = ee.String(tag).cat('_').cat(y.format('%d'));
    return P.addBands(PET).addBands(T).addBands(P.divide(PET))
      .rename(ee.List(['P_','PET_','T_','AI_']).map(function(b){ return ee.String(b).cat(ys); })).toFloat();
  })).toBands().regexpRename('^[0-9]+_', '');
}
var TODO = ['GFDL-ESM4','INM-CM5-0','IPSL-CM6A-LR','MIROC6','MPI-ESM1-2-HR','MRI-ESM2-0','NorESM2-MM','CMCC-ESM2'];
TODO.forEach(function(m){
  RUNS.forEach(function(r){
    var name = 'cmip6_annual_' + m.replace(/-/g, '_') + '_' + r[1];
    Export.image.toAsset({image: runStack(m, r), description: name,
      assetId: 'projects/YOUR_PROJECT/assets/' + name,
      crs: 'EPSG:4326', crsTransform: CT, region: region, maxPixels: 1e9, pyramidingPolicy: {'.default': 'mean'}});
  });
});
print('SUBMITTED 05b_cmip6_annual_by_scenario.js');
