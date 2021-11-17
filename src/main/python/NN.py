#!/usr/bin/python

import tensorflow.compat.v1 as tf

tf.disable_v2_behavior()
from tensorflow import keras
import numpy as np
import csv
import math


class NN:
  def __init__(self):
    # TODO: Setup the following variables to read from the parameter_dict
    self.ID = 0
    self.mz = 1
    self.rt = 2
    self.z = 3
    self.isonum = 4
    self.int_shape = 5
    self.iso_distr = 6
    self.intensity_window_avg = 7
    self.intensity_area_percentage = 8
    self.rt_start = 9
    self.rt_end = 10
    self.scan_num = 11
    self.intensity_sum = 12
    self.svr = 13

  def run(self, feature_file):
    filepath = feature_file
    with open(filepath + "_feature_one_z", 'r') as file:
      reader = csv.reader(file, delimiter='\t')
      test_data = list(reader)
    path = "/waterlooms/src/main/python"
    X_nn = test_data[1: len(test_data)]
    #intensity_lst = np.empty(len(X_nn))
    features = [self.isonum, self.int_shape, self.iso_distr,
                self.intensity_area_percentage]
    for i in range(0, len(X_nn)):
      for j in range(0, len(X_nn[i])):
        X_nn[i][j] = float(X_nn[i][j])
      #intensity_lst[i] = X_nn[i][self.intensity_sum]
    #np.sort(intensity_lst)
    #scale_int = intensity_lst[99]
    for i in range(0, len(X_nn)):
      #X_nn[i][self.intensity_sum] /= scale_int
      arr = [0] * len(features)
      for j in range(0, len(features)):
        arr[j] = X_nn[i][features[j]]
      X_nn[i] = arr
    X_nn = np.array(X_nn)
    predict_nn = NN.get_prediction(X_nn, path + "/model/NN_feature_rank")
    with open(filepath + "_nn_score", "w+") as outfile:
      outfile.write(
        "id\tmz\trt\tz\tisotope_num\tintensity_shape_score\tisotope_distribution_score\tintensity_window_avg\tintensity_intensity_area_percentage\trt_start\trt_end\tscan_num\tintensity_sum\tSVRscore\tquality_score\n")
      for i in range(len(X_nn)):
        outfile.write("%s\t" % test_data[i + 1][self.ID])
        outfile.write("%s\t" % test_data[i + 1][self.mz])
        outfile.write("%s\t" % test_data[i + 1][self.rt])
        outfile.write("%s\t" % test_data[i + 1][self.z])
        outfile.write("%s\t" % test_data[i + 1][self.isonum])
        outfile.write("%s\t" % test_data[i + 1][self.int_shape])
        outfile.write("%s\t" % test_data[i + 1][self.iso_distr])
        #test_data[i + 1][self.intensity_window_avg] *= scale_int
        outfile.write("%s\t" % test_data[i + 1][self.intensity_window_avg])
        outfile.write("%s\t" % test_data[i + 1][self.intensity_area_percentage])
        outfile.write("%s\t" % test_data[i + 1][self.rt_start])
        outfile.write("%s\t" % test_data[i + 1][self.rt_end])
        outfile.write("%s\t" % test_data[i + 1][self.scan_num])
        outfile.write("%s\t" % test_data[i + 1][self.intensity_sum])
        outfile.write("%s\t" % test_data[i + 1][self.svr])
        outfile.write("%f\n" % predict_nn[i])
    print("Done!")

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