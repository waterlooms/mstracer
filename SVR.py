#!/usr/bin/python

import csv
import pickle
import os

class SVR:
    def __init__(self):
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

    def run(self, feature_file):
        filepath = feature_file + "_feature_all_z"
        features = [self.isonum, self.int_shape, self.iso_distr]
        with open(filepath, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            test_data = list(reader)
        X_svr = test_data[1: len(test_data)]
        for i in range(0, len(X_svr)):
            arr = [0] * len(features)
            for j in range(0, len(features)):
                arr[j] = X_svr[i][features[j]]
            X_svr[i] = arr
        for i in range(0, len(X_svr)):
            for j in range(0, len(X_svr[i])):
                X_svr[i][j] = float(X_svr[i][j])
        clf = pickle.load(open(os.getcwd() + "/model/SVR_z_selection", 'rb'))
        # Prediction
        predict_svr = clf.predict(X_svr)
        # Write file
        with open(feature_file + "_svr_score", "w+") as outfile:
            outfile.write("id\tmz\trt\tz\tisotope_num\tintensity_shape_score\tisotope_distribution_score\tintensity_window_avg\tintensity_area_percentage\trt_start\trt_end\tscan_num\tintensity_sum\tSVRscore\n")
            for i in range(len(X_svr)):
                outfile.write("%s\t" % test_data[i + 1][self.ID])
                outfile.write("%s\t" % test_data[i + 1][self.mz])
                outfile.write("%s\t" % test_data[i + 1][self.rt])
                outfile.write("%s\t" % test_data[i + 1][self.z])
                outfile.write("%s\t" % test_data[i + 1][self.isonum])
                outfile.write("%s\t" % test_data[i + 1][self.int_shape])
                outfile.write("%s\t" % test_data[i + 1][self.iso_distr])
                outfile.write("%s\t" % test_data[i + 1][self.intensity_window_avg])
                outfile.write("%s\t" % test_data[i + 1][self.intensity_area_percentage])
                outfile.write("%s\t" % test_data[i + 1][self.rt_start])
                outfile.write("%s\t" % test_data[i + 1][self.rt_end])
                outfile.write("%s\t" % test_data[i + 1][self.scan_num])
                outfile.write("%s\t" % test_data[i + 1][self.intensity_sum])
                outfile.write("%f\n" % predict_svr[i])
        print("SVR Done!")



# Read in arguments
import argparse # add to beginning of file
parser = argparse.ArgumentParser(description='Train on features extracted from mzXML.')
parser.add_argument("-features",
                    help="Features computed from FeatureDetect."
                    )
parser.add_argument("-parameters",
                    default="/src/src/main/resources/featuredetect.params",
                    help="FeatureDetect parameter file."
                    )
args = parser.parse_args()

# Invoke the training
training = SVR()
training.run(args.features)
