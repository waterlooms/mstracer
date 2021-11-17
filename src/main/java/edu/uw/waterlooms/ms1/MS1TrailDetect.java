package edu.uw.waterlooms.ms1;

import edu.uw.waterlooms.entity.XIC;
import edu.uw.waterlooms.service.ParameterService;

import java.io.*;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class MS1TrailDetect {
    private int id_index;
    private int mz_index;
    private int rt_index;
    private int z_index;
    private int isotope_num_index;
    private int invalidVal;

    public MS1TrailDetect() {
        setParameters();
    }

    public void matchFeatureToTrails(String oldFilePath, String workingDirectory, String rawFileName, ArrayList<XIC> xics) throws IOException {
        List<Double[]> featuresMonoIsotope = readFile(oldFilePath + ".precursors");
        featuresMonoIsotope.sort(Comparator.comparing(l -> l[rt_index]));
        featuresMonoIsotope.sort(Comparator.comparing(l -> l[mz_index]));
        List<Double[]> featuresIsotope1 = getIsotopeFromIndex(featuresMonoIsotope, 1);
        List<Double[]> featuresIsotope2 = getIsotopeFromIndex(featuresMonoIsotope, 2);

        xics.sort(XIC::compareRtXIC);
        xics.sort(XIC::compareMzXIC);

        List<String[]> trailMonoisotope = matchFeaturesForIsotope(featuresMonoIsotope, xics);
        List<String[]> trailIsotope1 = matchFeaturesForIsotope(featuresIsotope1, xics);
        List<String[]> trailIsotope2 = matchFeaturesForIsotope(featuresIsotope2, xics);

        ArrayList<List<String[]>> trailInfoComponent = new ArrayList<>();
        trailInfoComponent.add(trailMonoisotope);
        trailInfoComponent.add(trailIsotope1);
        trailInfoComponent.add(trailIsotope2);

//        writeMS1TrailsData(oldFilePath + "_ms1_trails.tsv", xics);
        writeFile(oldFilePath + "_precursors_trails.tsv", featuresMonoIsotope, trailInfoComponent);
    }

    private List<String[]> matchFeaturesForIsotope (List<Double[]> featureIsotope, ArrayList<XIC> xics) {
        List<String[]> trailInfoStr = new ArrayList<>();
        int count = 0;
        double mzPPMTol = 8 * 10e-6;
        double rtTol = 0.1;
        String[] emptyStr = new String[]{"", "", "", "", ""};
        int j = 0;
        for (Double[] feature: featureIsotope) {
            double mz1 = feature[mz_index];
            double rt1 = feature[rt_index];
            if (mz1 == invalidVal) {
                trailInfoStr.add(emptyStr);
                continue;
            }
            // Debug
//            if (mz1 > 312.84 && mz1 < 312.85 && rt1 > 0.156 && rt1 < 0.157) {
//                int iii=0;
//            }
            if (j == xics.size()) {
                break;
            }
            while (j > 0 && mz1 * (1 - mzPPMTol) < xics.get(j).getMassChargeRatios().get(0)) {
                j--;
            }
            for (; j < xics.size(); j++) {
                ArrayList<Double> mzs = xics.get(j).getMassChargeRatios();
                ArrayList<Double> rts = xics.get(j).getRetentionTimes();
                ArrayList<Double> ints = xics.get(j).getIntensities();
                boolean if_found = false;
                boolean if_no_match = false;
                for (int m = 0; m < mzs.size(); m++) {
                    double mz2 = mzs.get(m);
                    double rt2 = rts.get(m);
                    if (m == mzs.size() - 1 && mz2 > mz1 * (1 + mzPPMTol)) {
                        trailInfoStr.add(emptyStr);
                        if_no_match = true;
                        break;
                    } else if (mz2 <=  mz1 * (1 + mzPPMTol)
                            && mz2 >= mz1 * (1 - mzPPMTol)
                            && rt2 <= rt1 + rtTol
                            && rt2 >= rt1 - rtTol) {
                        Double[] calculate_quantification = calculate_quantification(rts, ints);
                        String mzsStr = mzs.get(0).toString();
                        String rtsStr = rts.get(0).toString();
                        String intsStr = ints.get(0).toString();
                        String peaksSumStr = calculate_quantification[0].toString();
                        String peaksAreaStr = calculate_quantification[1].toString();
                        for (int k = 1; k < mzs.size(); k++) {
                            mzsStr +=  "," + mzs.get(k).toString();
                            rtsStr +=  "," + rts.get(k).toString();
                            intsStr += "," + ints.get(k).toString();
                        }
                        trailInfoStr.add(new String[]{mzsStr, rtsStr, intsStr, peaksSumStr, peaksAreaStr});
                        count++;
                        if_found = true;
                        break;
                    }
                }
                if (if_found == true || if_no_match == true) {
                    break;
                }
            }
        }
        System.out.println("Number of matched MS1 trails to ms-tracer results: " + count);
        return trailInfoStr;
    }

    private void setParameters() {
        id_index = ParameterService.getIdIndex();
        mz_index = ParameterService.getMzIndex();
        rt_index = ParameterService.getRtIndex();
        z_index = ParameterService.getzIndex();
        isotope_num_index = ParameterService.getIsonumIndex();
        invalidVal = ParameterService.getInvalidVal();
    }

    /**
     * Quantify trails w.r.t peak sum and peak area
     * @param rts
     * @param ints
     * @return [peaks_sum, peaks_area]
     */
    private Double[] calculate_quantification(ArrayList<Double> rts, ArrayList<Double> ints) {
        Double peaks_sum = ints.get(0);
        Double peaks_area = 0.0;
        Double prev_rt = rts.get(0);
        Double prev_int = ints.get(0);
        for (int i = 1; i < rts.size(); i++) {
            Double cur_rt = rts.get(i);
            Double cur_int = ints.get(i);
            peaks_sum += cur_int;
            peaks_area += (cur_rt - prev_rt) * (cur_int + prev_int) / 2;
        }
        return new Double[] {peaks_sum, peaks_area};
    }

    private List<Double[]> getIsotopeFromIndex(List<Double[]> isoList, int isotopeIndex){
        List<Double[]> result = new ArrayList<>();
        for (int i = 0; i < isoList.size(); i++) {
            Double[] line =  isoList.get(i);
            double cur_mz = line[mz_index]+ isotopeIndex * 1.003355/line[z_index];
            Double[] new_line = line.clone();
            if (isotopeIndex >= line[isotope_num_index]) {
                new_line[mz_index] = (double)invalidVal;
            } else {
                new_line[mz_index] = cur_mz;
            }
            result.add(new_line);
        }
        return result;
    }

    static List<Double[]> readFile(String filename) throws IOException {
        List<Double[]> result = new ArrayList<>();
        FileReader fileReader = new FileReader(filename);
        BufferedReader bufferedReader = new BufferedReader(fileReader);
        String line;
        boolean IsfirstLine = true;
        while ((line = bufferedReader.readLine()) != null) {
            // process first line
            if (IsfirstLine == true) {
                IsfirstLine = false;
            } else {
                // process line
                String[] words = line.split("\\t");
                Double[] nums = new Double[words.length];
                for (int i = 0; i < words.length; i++) {
                    nums[i] = Double.valueOf(words[i]);
                }
                result.add(nums);
            }
        }
        return result;
    }

    void writeFile(String filename, List<Double[]> oldList, ArrayList<List<String[]>> trailInfoComponent) throws IOException {
        File file = new File(filename);
        try {
            // create FileWriter object with file as parameter
            FileWriter outputfile = new FileWriter(file);
            // PrintWriter
            PrintWriter printWriter = new PrintWriter(outputfile);
            String header =
                    "id"
                            + '\t'
                            + "mz"
                            + '\t'
                            + "rt"
                            + '\t'
                            + "z"
                            + '\t'
                            + "isotope_num"
                            + '\t'
                            + "intensity_shape_score"
                            + '\t'
                            + "isotope_distribution_score"
                            + '\t'
                            + "intensity_window_evg"
                            + '\t'
                            + "intensity_area_percentage"
                            + '\t'
                            + "rt_start"
                            + '\t'
                            + "rt_end"
                            + '\t'
                            + "scan_num"
                            + '\t'
                            + "intensity_sum"
                            + '\t'
                            + "svr_score"
                            + '\t'
                            + "quality_score"
                            + '\t'
                            + "mzs"
                            + '\t'
                            + "rts"
                            + '\t'
                            + "ints"
                            + '\t'
                            + "peaks_sum"
                            + '\t'
                            + "peaks_area"
                            + '\t'
                            + "mzs_isotope+1"
                            + '\t'
                            + "rts_isotope+1"
                            + '\t'
                            + "ints_isotope+1"
                            + '\t'
                            + "peaks_sum_isotope+1"
                            + '\t'
                            + "peaks_area_isotope+1"
                            + '\t'
                            + "mzs_isotope+2"
                            + '\t'
                            + "rts_isotope+2"
                            + '\t'
                            + "ints_isotope+2"
                            + '\t'
                            + "peaks_sum_isotope+2"
                            + '\t'
                            + "peaks_area_isotope+2"
                            + '\n';
            printWriter.print(header);

            int length = oldList.get(0).length;
            int trailInfoLength = trailInfoComponent.get(0).get(0).length;
            for (int i = 0; i < oldList.size(); i++) {
                StringBuilder data = new StringBuilder();

                data.append(oldList.get(i)[id_index].intValue());
                for (int j = 1; j < length; j++) {
                    data.append('\t').append(oldList.get(i)[j].toString());
                }
                for (int j = 0; j < trailInfoComponent.size(); j++) {
                    for (int k = 0; k < trailInfoLength; k++) {
                        data.append('\t').append(trailInfoComponent.get(j).get(i)[k]);
                    }
                }

                data.append('\n');
                printWriter.print(data);
            }
            // closing writer connection
            printWriter.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
