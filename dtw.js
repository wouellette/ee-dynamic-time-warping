// ****************************************************************************************************************** //
// ********** Module to implement Time-Weighted or Time-Constrained Dynamic Time Warping (TW/TC-DTW) **************** //
// ****************************************************************************************************************** //

/**
 * Compute Dynamic Time Warping dissimilarity for each pixel in the multi-dimensional image array using a list of patterns/signatures.
 * For a deeper understanding of how arrays in GEE work, check out: https://medium.com/google-earth/runs-with-arrays-400de937510a
 * The time-weighted approach is taken from: Maus, V., Câmara, G., Cartaxo, R., Sanchez, A., Ramos, F. M., & De Queiroz, G. R. (2016).
 *                                           A time-weighted dynamic time warping method for land-use and land-cover mapping.
 *                                           IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 9(8), 3729-3739.
 * The time-constrained approach is taken from: Csillik, O., Belgiu, M., Asner, G. P., & Kelly, M. (2019).
 *                                              Object-based time-constrained dynamic time warping classification of crops using Sentinel-2.
 *                                              Remote sensing, 11(10), 1257.
 * The angular distance approach is taken from: Teke, Mustafa, and Yasemin Y. Çetin.
 *                                             "Multi-year vector dynamic time warping-based crop mapping."
 *                                             Journal of Applied Remote Sensing 15.1 (2021): 016517.
 * @param {Array} patterns_arr: An array of dimension [k, n, t],
 *                              with k the number of patterns, n number of bands,
 *                              and t the number of timestamps in the time series.
 * @param {Array} stack: An image collection with t number of images. Each image has n bands.
 * @param {Dictionary} options: The options consist of the following parameters:
 *                              - @param {Number} patterns_no: Number of patterns to iterate over.
 *                                If not specified, the patterns number is computed from the patterns array
 *                                (as long as function is not used inside of a mapping routine, will fail if not provided).
 *                              - @param {Number} band_no: Number of bands (excluding the Day of Year band) to iterate over.
 *                                If not specified, the bands number computed from the patterns array
 *                                (as long as function is not used inside of a mapping routine, will fail if not provided).
 *                              - @param {Number} timeseries_len: The length of the image time series.
 *                                Also computed from the image array if not provided
 *                                (as long as function is not used inside of a mapping routine, will fail if not provided).
 *                              - @param {Number} patterns_len: The length of the reference pattern time series.
                                  Also computed from the patterns attributes if not provided
                                  (as long as function is not used inside of a mapping routine, will fail if not provided).
 *                              - @param {String} constraint_type: The type of time constraint to use,
 *                                whether 'time-weighted' or 'time-constrained'.
 *                              - @param {String} weight_type: The type of weight to apply for the 'time-weighted' approach.
 *                                Defaults to 'logistic' as it represents natural and phenological cycles better.
 *                                Ignored if 'time-constrained' is chosen as constraint type.
 *  *                           - @param {String} distance_type: The type of distance to apply for the dissimilarity calculation.
 *                                Defaults to 'euclidean', but 'angular' can also be provided to compute the angular distance
 *                                as specified by the Spectral Angle Mapper (SAM).
 *                              - @param {Number} beta: The Beta parameter of the time-weighted approach,
 *                                                      defining the tolerance of the function.
 *                                                      If 'time-constrained' is defined, beta is the constraint period.
 *                                                      Defaults to 50 (days).
 *                              - @param {Number} alpha: The Alpha parameter of the time-weighted approach,
 *                                                       defining the steepness of the logistic function used.
 *                                                       Defaults to 0.1.
 * @returns {Image}
 * @ignore
 */
exports.DTWDist = function(patterns_arr, stack, options){
  stack=ee.Image(stack);
  
  var bandNames=stack.bandNames();
  patterns_arr = ee.Array(patterns_arr);
  var patterns_no = options.patterns_no || patterns_arr.length().get([0]);
  var band_no = options.band_no || patterns_arr.length().get([1]).subtract(1);
  var timeseries_len = options.timeseries_len;
  var patterns_len = options.patterns_len || patterns_arr.length().get([2]);
  var constraint_type = options.constraint_type || 'time-weighted';
  var weight_type = options.weight_type || 'logistic';
  var distance_type = options.distance_type || 'euclidean';
  var beta = options.beta || 50;
  var alpha = options.alpha || 0.1;

  // This step calculates distance between NDVI_series and mean_NDVI of croppoints using DTW.
  var dis_arr;
  var time_arr;
  var dis_mat;
  var dis_list;
  var cost_weight;
  var x1;
  var y1;
  var t1;
  var x2;
  var y2;
  var t2;
  
  var matrix=ee.List.sequence(1, ee.Number(timeseries_len).subtract(1)).map(function(i){
    var matrix_tmp=ee.List.sequence(1, ee.Number(patterns_len).subtract(1)).map(function(j){
      return ee.List([i,j]);
    });
    return matrix_tmp;
  });
  matrix=ee.List(matrix.iterate(function(x,previous){
    return ee.List(previous).cat(ee.List(x));
  },ee.List([])));
  
          
  var dtw_image_list = ee.List.sequence(1, patterns_no).map(function(k){
    
    if (distance_type === 'euclidean') {
      if (constraint_type === 'time-constrained') {
    
      var dt_list=ee.List.sequence(1, timeseries_len).map(function (i) {
        i = ee.Number(i);
        var dt_list = ee.List.sequence(1, patterns_len).map(function (j) {
          var time_fc=ee.List([]);
          j = ee.Number(j);
          //ASSUMPTION: timeseries_len IS EQUAL TO patterns_len AND THE SAME TIME STEPS ARE USED
          var t1 = patterns_arr.get(ee.List([ee.Number(k).subtract(1),-1,i.subtract(1)]));
          var t2 = patterns_arr.get(ee.List([ee.Number(k).subtract(1),-1,j.subtract(1)]));
          var time_arr= t1.subtract(t2).abs();
          time_fc = time_fc.add(ee.Feature(null, {
            'dt': time_arr,
            'i': i,
            'j': j
          }));
          return time_fc;
        });
        return dt_list;
      });
      dt_list=dt_list.flatten();
      
      dis_mat=ee.List(ee.List.sequence(1, timeseries_len).iterate(function(i,previous0){
        i=ee.Number(i);
        previous0=ee.List(previous0);
        //filter the patterns within dt<=beta
        //iterate over all time steps dt<=beta
        var patterns_tmp=ee.FeatureCollection(dt_list).filter(ee.Filter.and(ee.Filter.eq('i',i),ee.Filter.lte('dt',beta))).aggregate_array('j');
        var dis_list0=ee.List(patterns_tmp.iterate(function(j,previous){
          j=ee.Number(j);
          previous=ee.List(previous);

          var dis_sum=ee.Image(ee.List.sequence(1, ee.Number(band_no)).iterate(function(n,previous2){
            n=ee.Number(n);
            previous2=ee.Image(previous2);
            var x1=stack.select(ee.String(bandNames.get(n.subtract(1).multiply(ee.Number(timeseries_len)).add(i).subtract(1))));
            x1=x1.toArray().toInt16();
            var y1 = patterns_arr.get(ee.List([ee.Number(k).subtract(1), n.subtract(1), j.subtract(1)]));
            var dis_arr = x1.subtract(y1).pow(2);
            return dis_arr.add(previous2);
          },ee.Image(0)));
          
          previous = previous.add(dis_sum.sqrt().set('j',j));
          
          return previous;
        },ee.List([])));
        //iterate over all time steps dt>beta
        patterns_tmp=ee.FeatureCollection(dt_list).filter(ee.Filter.and(ee.Filter.eq('i',i),ee.Filter.gt('dt',beta))).aggregate_array('j');
        var dis_list=ee.List(patterns_tmp.iterate(function(j,previous){
          j=ee.Number(j);
          previous=ee.List(previous);
          previous = previous.add(ee.Image(1e6).set('j',j));
          return previous;
        },dis_list0));
        //now sort in chronological order
        dis_list=ee.ImageCollection(dis_list).sort('j').toList(dis_list.length());
      
        return previous0.add(dis_list);
      
      },ee.List([])));

    } else if (constraint_type === 'time-weighted') {
      
      dis_mat=ee.List(ee.List.sequence(1, timeseries_len).iterate(function(i,previous0){
        i=ee.Number(i);
        previous0=ee.List(previous0);      
        var dis_list=ee.List(ee.List.sequence(1, patterns_len).iterate(function(j,previous){
          j=ee.Number(j);
          previous=ee.List(previous);

          var dis_sum=ee.Image(ee.List.sequence(1, ee.Number(band_no)).iterate(function(n,previous2){
            n=ee.Number(n);
            previous2=ee.Image(previous2);
            var x1=stack.select(ee.String(bandNames.get(n.subtract(1).multiply(ee.Number(timeseries_len)).add(i).subtract(1))));
            x1=x1.toArray().toInt16();
            var y1 = patterns_arr.get(ee.List([ee.Number(k).subtract(1), n.subtract(1), j.subtract(1)]));
            var dis_arr = x1.subtract(y1).pow(2);
            return dis_arr.add(previous2);
          },ee.Image(0)));
          
          //the doy bands come last in the stack
          var t1=stack.select(ee.String(bandNames.get(ee.Number(band_no).multiply(ee.Number(timeseries_len)).add(i).subtract(1))));
          t1=t1.toArray().toInt16();
          var t2 = patterns_arr.get(ee.List([ee.Number(k).subtract(1), -1, j.subtract(1)]));
          var time_arr = t1.subtract(t2).abs();
  
          var cost_weight;
          if (weight_type === 'logistic') {
            cost_weight = ee.Image(1).divide(ee.Image(1)
                          .add(ee.Image(alpha).multiply(time_arr.subtract(beta))).exp());
          } else if (weight_type === 'linear') {
            cost_weight = ee.Image(alpha).multiply(time_arr).add(beta);
          }
          previous = previous.add(dis_sum.sqrt().add(cost_weight));
        
          return previous;
        },ee.List([])));
      
        return previous0.add(dis_list);
      
      },ee.List([])));
    }} 
    //distance_type 'angular' not yet implemented as a server-side only function
    /*
    dis_mat = ee.List([]);
    for(var i=1;i<=timeseries_len;i++){
      dis_list = ee.List([]);
      for(var j=1;j<=patterns_len;j++){
        var dis_sum = ee.Image(0);
        for(var n=1; n<=band_no;n++){

          x1 = image_arr.arraySlice(0, i-1, i).arraySlice(1, n-1, n);
          y1 = patterns_arr.get(ee.List([ee.Number(k).subtract(1), n-1, j-1]));

          if (distance_type === 'angular' && i > 1) {
            x2 = image_arr.arraySlice(0, i-2, i-1).arraySlice(1, n-1, n);
            y2 = patterns_arr.get(ee.List([ee.Number(k).subtract(1), n-1, j-2]));
            dis_arr = x1.multiply(y1).add(x2.multiply(y2))
                      .divide(x1.pow(2).add(x2.pow(2)).sqrt().multiply(y1.pow(2).add(y2.pow(2)).sqrt()))
                      .acos();
          } else if (distance_type === 'euclidean') {
            dis_arr = x1.subtract(y1).pow(2);
          }

          dis_sum = (distance_type === 'angular' && i === 1) ? dis_sum : dis_arr.add(dis_sum);
        }

        t1 = image_arr.arraySlice(0, i-1, i).arraySlice(1, -1);
        t2 = patterns_arr.get(ee.List([ee.Number(k).subtract(1), -1, j-1]));
        time_arr = t1.subtract(t2).abs();

        if (constraint_type === 'time-weighted') {
          if (weight_type === 'logistic') {
            cost_weight = ee.Image(1).divide(ee.Image(1)
                          .add(ee.Image(alpha).multiply(time_arr.subtract(beta))).exp());
          } else if (weight_type === 'linear') {
            cost_weight = ee.Image(alpha).multiply(time_arr).add(beta);
          }
          dis_list = dis_list.add(dis_sum.sqrt().add(cost_weight));
        } else if (constraint_type === 'time-constrained') {
          dis_list = dis_list.add(ee.Algorithms.If(time_arr.lte(beta),
                                                   dis_sum.sqrt(),
                                                   dis_sum.add(1e6)));
        } else {
          dis_list = dis_list.add(dis_sum.sqrt());
        }
      }
      dis_mat = dis_mat.add(dis_list);
    }
    */
    var D_mat = ee.List([]);
    var dis = ee.List(dis_mat.get(0)).get(0);
    var d_mat = D_mat.add(dis);

    d_mat=ee.List(ee.List.sequence(1, ee.Number(patterns_len).subtract(1)).iterate(function(j,previous){
      j=ee.Number(j);
      previous=ee.List(previous);
      var dis = ee.Image(previous.get(j.subtract(1))).add(ee.List(dis_mat.get(0)).get(j));
      return previous.add(dis);
    },d_mat));
    D_mat = D_mat.add(d_mat);
    
    D_mat=ee.List(ee.List.sequence(1, ee.Number(timeseries_len).subtract(1)).iterate(function(i,previous){
      i=ee.Number(i);
      previous=ee.List(previous);
      var dis = ee.Image(ee.List(previous.get(i.subtract(1))).get(0)).add(ee.List(dis_mat.get(i)).get(0));
      return previous.add(ee.List(previous.get(i.subtract(1))).set(0, dis));
    },D_mat));
    
    D_mat=ee.List(ee.List.sequence(0, matrix.length().subtract(1)).iterate(function(x,previous){
      var i=ee.Number(ee.List(matrix.get(x)).get(0));
      var j=ee.Number(ee.List(matrix.get(x)).get(1));
      previous=ee.List(previous);
      var dis = ee.Image(ee.List(previous.get(i.subtract(1))).get(j)).min(ee.List(previous.get(i)).get(j.subtract(1)))
        .min(ee.List(previous.get(i.subtract(1))).get(j.subtract(1))).add(ee.List(dis_mat.get(i)).get(j));
      return previous.set(i, ee.List(previous.get(i)).set(j, dis));
    },D_mat));

    return ee.Image(ee.List(D_mat.get(-1)).get(-1))
                                          .arrayProject([0])
                                          .arrayFlatten([['DTW']]);
  });

  return ee.ImageCollection(dtw_image_list).min().toUint16();
};


/**
 * A utility that converts the feature collection containing the signatures/patterns to an dtw-ready array.
 * @param {FeatureCollection} signatures: A feature collection containing the signatures/patterns to be used as input to DTW.
 * @param {String} class_name: The property name of the label class containing the class values as integers.
 *
 * @param {Number} class_no: Integer of the label class to retrieve from the feature collection.
 * @param {Number} band_no: Number of bands (excluding the Day of Year band).
 * @param {Number} patterns_len:  The length of the reference pattern time series.
 * @param {List} band_names: The list of band names containing the pattern/signature values to retrieve from the feature collection.
 *                           The Day of Year band (doy) must be placed last; for instance for ndvi, VV and day of year (doy) bands:
 *                           [ndvi, ndvi_1, ndvi_2, ndvi_n, evi, evi_1, evi_2, evi_n, doy, doy_1, doy_2, doy_n].
 * @returns {Array}
 * @ignore
 */
exports.prepareSignatures = function(signatures, class_name, class_no, band_no, patterns_len, band_names){
  var train_points = signatures.filter(ee.Filter.eq(class_name, class_no)).select(band_names);
  var feature_list = train_points.toList(train_points.size());
  var table_points_list = feature_list.map(function (feat){
    return ee.List(band_names).map(function (band){return ee.Feature(feat).get(band)});
  });

  return ee.Array(table_points_list).reshape([-1, band_no+1, patterns_len]).toInt16();
};

/**
 * A utility that converts a multi-band image containing the time series' timestamps to a dtw-ready array.
 * @param {Image} image: The multi-band image containing the time series' timestamps.
 * @param {Number} band_no: Number of bands (excluding the Day of Year band) to iterate over.
 * @param {Number} timeseries_len: The length of the image time series.
 * @param {List} band_names: The list of band names containing the pattern/signature values to retrieve from the feature collection.
 *                           The Day of Year band (doy) must be placed last; for instance for ndvi, VV and day of year (doy) bands:
 *                           [ndvi, ndvi_1, ndvi_2, ndvi_n, evi, evi_1, evi_2, evi_n, doy, doy_1, doy_2, doy_n].
 * @returns {Array}
 * @ignore
 */
exports.prepareBands = function(image, band_no, timeseries_len, band_names){
  var arr_list = ee.List([]);
  var band_image_arr = ee.List.sequence({start: timeseries_len, step: timeseries_len, count: band_no})
                       .iterate(function(i, img){
    var image_arr = image.select(band_names.slice(i, ee.Number(i).add(timeseries_len))).toArray().toArray(1).toInt16();
    return ee.Image(img).arrayCat(image_arr, 1)},
    image.select(band_names.slice(0, timeseries_len)).toArray().toArray(1).toInt16());

  return band_image_arr;
};
