#!/usr/bin/python

# TODO: Suppress warnings but not errors in production

import csv
import pickle
import sys
sys.path.append("/waterlooms/src/main/python") #TODO Possibly remove this

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
        self.intensity_window_avg = 7
        self.intensity_area_percentage = 8
        self.rt_start = 9
        self.rt_end = 10
        self.scan_num = 11
        self.intensity_sum = 12

    # TO-DO: change to run(self, feature_file, parameter_dict)
    def run(self, feature_file):
        start = self.isonum
        end = self.iso_distr
        filepath = feature_file
        with open(filepath + "_feature_all_z", 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            test_data = list(reader)
        X_svr = test_data[1: len(test_data)]  # jump header
        for i in range(0, len(X_svr)):
            X_svr[i] = X_svr[i][start:end+1]

        for i in range(0, len(X_svr)):
            for j in range(0, len(X_svr[i])):
                X_svr[i][j] = float(X_svr[i][j])

        path = "/waterlooms/src/main/python" #TODO Finalize the pathing
        clf = pickle.load(open(path + "/model/SVR_z_selection", 'rb'))
        predict_svr = clf.predict(X_svr)

        with open(filepath + "_svr_score", "w+") as outfile:
            outfile.write("id\tmz\trt\tz\tisotope_num\tintensity_shape_score\tisotope_distribution_score\tintensity_intensity_window_avg\tintensity_intensity_area_percentage\trt_start\trt_end\tscan_num\tintensity_sum\tSVRscore\n")
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