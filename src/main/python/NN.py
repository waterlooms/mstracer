#!/usr/bin/python

import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()

from tensorflow import keras
import numpy as np
import csv
import math
import os

tf.logging.set_verbosity(tf.logging.ERROR)

class NN:
  def __init__(self):
    self.ID = 0
    self.mz = 1
    self.rt = 2
    self.z = 3
    self.isonum = 4
    self.int_shape = 5
    self.iso_distr = 6
    self.intensity_area_percentage = 7
    self.scan_num = 8
    self.quantification_peaks_sum = 9
    self.quantification_peaks_area = 10
    self.svr_score = 11
    self.quality_score = 12

  def run(self, feature_file):
    filepath = os.path.splitext(feature_file)[0]
    with open(filepath + "_featureOneZ.tsv", 'r') as file:
      reader = csv.reader(file, delimiter='\t')
      test_data = list(reader)
    path = "/mstracer/src/main/python"
    X_nn = test_data[1: len(test_data)]
    features = [self.isonum, self.int_shape, self.iso_distr,
                self.intensity_area_percentage]
    for i in range(0, len(X_nn)):
      for j in range(0, len(X_nn[i])):
        X_nn[i][j] = float(X_nn[i][j])

    for i in range(0, len(X_nn)):
      arr = [0] * len(features)
      for j in range(0, len(features)):
        arr[j] = X_nn[i][features[j]]
      X_nn[i] = arr
    X_nn = np.array(X_nn)
    predict_nn = NN.get_prediction(X_nn, path + "/model/NN_feature_rank")
    with open(filepath + "_NNScore.tsv", "w+") as outfile:
      outfile.write("id\tmz\trt\tz\tisotope_num\tintensity_shape_score\tisotope_distribution_score\tintensity_area_percentage\tscan_num\tquantification_peaks_sum\tquantification_peaks_area\tsvr_score\tquality_score\tmzs,rts,ints\n")
      for i in range(1, len(X_nn) + 1):
        scannum = int(test_data[i][self.scan_num])
        mzsStr = str(test_data[i][self.quality_score + 1])
        rtsStr =  str(test_data[i][self.quality_score + 1 + scannum])
        intsStr =  str(test_data[i][self.quality_score + 1 + scannum * 2])
        for j in range(1, scannum):
          mzsStr += "\t%s" % test_data[i][self.quality_score + 1 + j]
          rtsStr += "\t%s" % test_data[i][self.quality_score + 1 + scannum + j]
          intsStr += "\t%s" % test_data[i][self.quality_score + 1 + scannum * 2 + j]

        outfile.write("%s\t" % test_data[i][self.ID])
        outfile.write("%s\t" % test_data[i][self.mz])
        outfile.write("%s\t" % test_data[i][self.rt])
        outfile.write("%s\t" % test_data[i][self.z])
        outfile.write("%s\t" % test_data[i][self.isonum])
        outfile.write("%s\t" % test_data[i][self.int_shape])
        outfile.write("%s\t" % test_data[i][self.iso_distr])
        outfile.write("%s\t" % test_data[i][self.intensity_area_percentage])
        outfile.write("%s\t" % scannum)
        outfile.write("%s\t" % test_data[i][self.quantification_peaks_sum])
        outfile.write("%s\t" % test_data[i][self.quantification_peaks_area])
        outfile.write("%s\t" % test_data[i][self.svr_score])
        outfile.write("%f\t" % predict_nn[i - 1])
        outfile.write(mzsStr + "\t")
        outfile.write(rtsStr + "\t")
        outfile.write(intsStr + "\n")
    print("Assaigning NN Score Finished!")

  def get_prediction(X_test, fname):
    model = keras.models.load_model(fname)
    probability_model = tf.keras.Sequential([model, tf.keras.layers.Softmax()])
    predictions = probability_model.predict(X_test)
    prediction = []
    for i in range(len(predictions)):
      if predictions[i][0] == 0:
        log_score = 100000
      else:
        log_score = math.log(predictions[i][1] / predictions[i][0])
      prediction.append(log_score)
    return prediction


# Read in arguments
import argparse  # add to beginning of file

parser = argparse.ArgumentParser(
  description='Train on features extracted from mzXML.')
parser.add_argument("-features",
                    help="Features computed from FeatureDetect."
                    )
parser.add_argument("-parameters",
                    default="/data/featuredetect.params",
                    help="FeatureDetect parameter file."
                    )
args = parser.parse_args()

# Parse the .env like parameter file
# from dotenv import dotenv_values

# parsed = dotenv_values(stream=args.parameters)

# Invoke the training
training = NN()
# TO-DO: change to training.run(args.features, parsed)
training.run(args.features)
