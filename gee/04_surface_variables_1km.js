// LST (annual and seasonal), albedo, ET and GEDI biomass on the MODIS 1 km grid, 2019-2023.
var y0 = '2019-01-01', y1 = '2024-01-01';
function lstColl(band, qc){
  return ee.ImageCollection('MODIS/061/MYD11A2').filterDate(y0, y1).map(function(i){
    // QA: good or other quality, LST error <= 2 K
    var q = i.select(qc).bitwiseAnd(3).lte(1).and(i.select(qc).rightShift(6).bitwiseAnd(3).lte(1));
    return i.select(band).multiply(0.02).subtract(273.15).updateMask(q).copyProperties(i, ['system:time_start']);
  });
}
var lstD = lstColl('LST_Day_1km','QC_Day'), lstN = lstColl('LST_Night_1km','QC_Night');   // Aqua ~13:30 / ~01:30
function seas(c, months){ return c.filter(ee.Filter.calendarRange(months[0], months[1], 'month')).mean(); }
var winF = ee.Filter.or(ee.Filter.calendarRange(10,12,'month'), ee.Filter.calendarRange(1,2,'month'));
var lstDann = lstD.mean().rename('LSTd'), lstNann = lstN.mean().rename('LSTn');
var lstDmam = seas(lstD,[3,5]).rename('LSTd_MAM'), lstNmam = seas(lstN,[3,5]).rename('LSTn_MAM');
var lstDjjas = seas(lstD,[6,9]).rename('LSTd_JJAS'), lstNjjas = seas(lstN,[6,9]).rename('LSTn_JJAS');
var lstDwin = lstD.filter(winF).mean().rename('LSTd_ONDJF'), lstNwin = lstN.filter(winF).mean().rename('LSTn_ONDJF');
var alb = ee.ImageCollection('MODIS/061/MCD43A3').filterDate(y0, y1).map(function(i){
  return i.select('Albedo_WSA_shortwave').multiply(0.001)
    .updateMask(i.select('BRDF_Albedo_Band_Mandatory_Quality_shortwave').eq(0));   // full BRDF inversions only
}).mean().rename('ALB');
var et = ee.ImageCollection('MODIS/061/MOD16A2GF').filterDate(y0, y1).select('ET').mean()
  .multiply(0.1 * 365 / 8).rename('ET');                              // mm per year
var agb = ee.Image('LARSE/GEDI/GEDI04_B_002').select('MU').rename('AGB');   // Mg/ha

var region = ee.Geometry.Rectangle([68.0, 6.5, 97.5, 37.5], null, false);
var stackCL = ee.Image.cat([lstDann, lstNann, lstDmam, lstNmam, lstDjjas, lstNjjas, lstDwin, lstNwin, alb, et, agb]).toFloat();
Export.image.toAsset({image: stackCL, description: 'india_1km_climate_v2',
  assetId: 'projects/YOUR_PROJECT/assets/india_1km_climate_v2', crs: 'SR-ORG:6974', scale: 926.625433055833,
  region: region, maxPixels: 1e10, pyramidingPolicy: {'.default': 'mean'}});
