// Forest type, crop phenology and soil predictors, stage 1: class-weighted sums at 0.05 deg.
// EVG evergreen share, FSEAS leaf seasonality, CNC/CGSL/CAMP crop cycles, NDVI_kharif, CLAY/SAND.
var A = 'projects/YOUR_PROJECT/assets/';
var india = ee.FeatureCollection(A + 'india_pov_boundary').geometry();
var lc = ee.Image(A + 'india_lc_154m');
var nat = lc.select('fNat').gte(50), crop = lc.select('fCrop').gte(50);
var P005 = ee.Projection('EPSG:4326', [0.05, 0, 68.0, 0, -0.05, 37.5]);
var inIndia = ee.Image.constant(1).clip(india).reproject(P005).mask();

var q1c = ee.ImageCollection('MODIS/061/MCD12Q1').filterDate('2019-01-01', '2024-01-01').select('LC_Type1');
var q2c = ee.ImageCollection('MODIS/061/MCD12Q2').filterDate('2019-01-01', '2024-01-01');
var m13 = ee.ImageCollection('MODIS/061/MOD13Q1').filterDate('2019-01-01', '2024-01-01').select('NDVI');
var P500 = ee.Image(q1c.first()).projection(), P250 = ee.Image(m13.first()).projection();

// [mean(V*M), mean(M)] at 0.05 deg
function agg(V, M, proj, name){
  var ok = V.mask().unmask(0).and(M.unmask(0)).setDefaultProjection(proj);
  var fVM = V.updateMask(ok).unmask(0).setDefaultProjection(proj).reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 1024})
    .reproject(P005).rename(name + '_s');
  var fM = ok.toFloat().reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 1024}).reproject(P005).rename(name + '_n');
  return [fVM.updateMask(inIndia), fM.updateMask(inIndia)];
}

// forest type and leaf seasonality
var lc1 = q1c.mode();
var modF = lc1.gte(1).and(lc1.lte(5));
var evg = lc1.eq(1).or(lc1.eq(2)).toFloat();
var natF = nat.and(modF);
var ndW = m13.filter(ee.Filter.calendarRange(9, 11, 'month')).median().multiply(0.0001);
var ndD = m13.filter(ee.Filter.calendarRange(3, 4, 'month')).median().multiply(0.0001);
var fseas = ee.Image(1).subtract(ndD.divide(ndW)).updateMask(ndW.gt(0.2));
// crop phenology
var q2 = q2c.map(function(i){
  var good = i.select('QA_Overall_1').lte(2);
  return i.select('NumCycles').toFloat().rename('CNC')
    .addBands(i.select('Dormancy_1').subtract(i.select('Greenup_1')).toFloat().updateMask(good).rename('CGSL'))
    .addBands(i.select('EVI_Amplitude_1').multiply(0.0001).toFloat().updateMask(good).rename('CAMP'));
}).mean();
var ndK = m13.filter(ee.Filter.calendarRange(8, 10, 'month')).mean().multiply(0.0001);

// SoilGrids clay and sand (%), 0-30 cm
function soil(prop){
  var sg = ee.Image('projects/soilgrids-isric/' + prop + '_mean');
  return sg.select(prop + '_0-5cm_mean').multiply(5).add(sg.select(prop + '_5-15cm_mean').multiply(10))
    .add(sg.select(prop + '_15-30cm_mean').multiply(15)).divide(30 * 10);
}
var Psg = ee.Image('projects/soilgrids-isric/clay_mean').projection();
var soilImg = ee.Image.cat([soil('clay').rename('CLAY'), soil('sand').rename('SAND')]).setDefaultProjection(Psg)
  .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 1024}).reproject(P005).updateMask(inIndia);

function cat(parts){ return ee.Image.cat(parts).toFloat(); }
var groups = {
  EVG: cat(agg(evg, natF, P500, 'EVG')),
  FSEAS: cat(agg(fseas, nat, P250, 'FSEAS')),
  CROP: cat(agg(q2.select('CNC'), crop, P500, 'CNC').concat(agg(q2.select('CGSL'), crop, P500, 'CGSL'))
            .concat(agg(q2.select('CAMP'), crop, P500, 'CAMP'))),
  NDVIK: cat(agg(ndK, crop, P250, 'NDVI_kharif')),
  soil: soilImg.toFloat()
};
var region = ee.Geometry.Rectangle([68.0, 6.5, 97.5, 37.5], null, false);
Object.keys(groups).forEach(function(g){
  Export.image.toAsset({image: groups[g], description: 'india_type_' + g + '_005deg',
    assetId: A + 'india_type_' + g + '_005deg', region: region, crs: 'EPSG:4326',
    crsTransform: [0.05, 0, 68.0, 0, -0.05, 37.5], maxPixels: 1e10, shardSize: 64});
});
