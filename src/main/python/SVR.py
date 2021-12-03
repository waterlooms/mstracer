#!/usr/bin/python

# TODO: Suppress warnings but not errors in production

import csv
import pickle
import sys
import os
sys.path.append("/mstracer/src/main/python") #TODO Possibly remove this

import argparse

class SVR:
    def __init__(self):
        # TODO: Setup the following variables to read from the parameter_dict
        self.ID = 0
        self.mz = 1
        self.rt = 2
        self.z = 3
        self.isonum = 4
        self.int_shape = 5
        self.iso_distr = 6
        self.intensity_area_percentage = 7
        self.rt_start = 8
        self.rt_end = 9
        self.scan_num = 10
        self.quantification_peaks_sum = 11
        self.quantification_peaks_area = 12

    # TO-DO: change to run(self, feature_file, parameter_dict)
    def run(self, feature_file):
        start = self.isonum
        end = self.iso_distr
        filepath = os.path.splitext(feature_file)[0]
        with open(filepath + "_featureAllZ.tsv", 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            test_data = list(reader)
        X_svr = test_data[1: len(test_data)]  # jump header
        for i in range(0, len(X_svr)):
            X_svr[i] = X_svr[i][start:end+1]

        for i in range(0, len(X_svr)):
            for j in range(0, len(X_svr[i])):
                X_svr[i][j] = float(X_svr[i][j])

        path = "/mstracer/src/main/python" #TODO Finalize the pathing
        clf = pickle.load(open(path + "/model/SVR_z_selection", 'rb'))
        predict_svr = clf.predict(X_svr)

        with open(filepath + "_SVRScore.tsv", "w+") as outfile:
            outfile.write("id\tmz\trt\tz\tisotope_num\tintensity_shape_score\tisotope_distribution_score\tintensity_area_percentage\trt_start\trt_end\tscan_num\tquantification_peaks_sum\tquantification_peaks_area\tSVRscore\n")
            for i in range(len(X_svr)):
                outfile.write("%s\t" % test_data[i + 1][self.ID])
                outfile.write("%s\t" % test_data[i + 1][self.mz])
                outfile.write("%s\t" % test_data[i + 1][self.rt])
                outfile.write("%s\t" % test_data[i + 1][self.z])
                outfile.write("%s\t" % test_data[i + 1][self.isonum])
                outfile.write("%s\t" % test_data[i + 1][self.int_shape])
                outfile.write("%s\t" % test_data[i + 1][self.iso_distr])
                outfile.write("%s\t" % test_data[i + 1][self.intensity_area_percentage])
                outfile.write("%s\t" % test_data[i + 1][self.rt_start])
                outfile.write("%s\t" % test_data[i + 1][self.rt_end])
                outfile.write("%s\t" % test_data[i + 1][self.scan_num])
                outfile.write("%s\t" % test_data[i + 1][self.quantification_peaks_sum])
                outfile.write("%s\t" % test_data[i + 1][self.quantification_peaks_area])
                outfile.write("%f\n" % predict_svr[i])
        print("Done!")



# Read in arguments
import argparse # add to beginning of file
parser = argparse.ArgumentParser(description='Train on features extracted from mzXML.')
parser.add_argument("-features",
    help="Features computed from FeatureDetect."
    )
parser.add_argument("-parameters",
    default="/data/featuredetect.params",
    help="FeatureDetect parameter file."
    )
args = parser.parse_args()



# Parse the .env like parameter file
# TO-DO uncomment this
# from dotenv import dotenv_values
# parsed = dotenv_values(stream=args.parameters)


# Invoke the training
training = SVR()
# TO-DO: change to training.run(args.features, parsed)
training.run(args.features)