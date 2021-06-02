package edu.uw.waterlooms.msutil;

import org.apache.commons.math3.util.Pair;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.fileio.exceptions.FileParsingException;

import java.util.ArrayList;
import java.util.Map;
import java.util.Set;

// *** this class holds extracted mzML information ***
// *** including:
//        Isolation window ranges
//        MS2 scans divided according to isolation windows

/**
 * Class holds extracted mzXML information.
 * Isolation Window ranges.
 * MS_x scans divided according to Isolation Window ranges.
 * x is the MS level, one of [1,2].
 */
public class MZXMLReader {
  protected OpenMzxml of;
  protected Scans msx_scans;
  protected int msLevel;
  // isolation window info
  ArrayList<Pair<Double, Double>> isolation_windows = new ArrayList<>();

  // MS Info
  protected Set<Map.Entry<Integer, IScan>> scans;
  ArrayList<ArrayList<ArrayList<ScanEntry>>> spectrumPoints = new ArrayList<>(); // spectrum points
  ArrayList<ArrayList<IScan>> fullScans = new ArrayList<>();

  /**
   * MZXML Reader Constructor
   * @param mzXML string containing ABSOLUTE path to mzXML file.
   * @param msLevel integer representing the msLevel to parse at. [1, 2]
   */
  public MZXMLReader(String mzXML, int msLevel) {
    // read mzml
    this.of = new OpenMzxml(mzXML);
    this.msx_scans = new Scans(of, msLevel);
    this.msLevel = msLevel;
    this.scans = msx_scans.scanEntries;
    read_ms2();
  }

  public ArrayList<Pair<Double, Double>> getIsolation_windows() {
    return this.isolation_windows;
  }

  public ArrayList<ArrayList<ArrayList<ScanEntry>>> getSpectrumPoints() {
    return this.spectrumPoints;
  }

  private void read_ms2() {
    for (Map.Entry<Integer, IScan> scan : scans) {
      // get window size
      double mz_start = scan.getValue().getPrecursor().getMzRangeStart();
      double mz_end = scan.getValue().getPrecursor().getMzRangeEnd();
      Pair win_range = new Pair<>(mz_start, mz_end);
      // get scan info
      IScan i_scan = scan.getValue();
      ISpectrum spectrum = null;
      try {
        spectrum = i_scan.fetchSpectrum();
      } catch (FileParsingException e) {
        e.printStackTrace();
      }
      double[] MZs = spectrum.getMZs();
      double[] intensities = spectrum.getIntensities();
      double rt = i_scan.getRt();
      ArrayList<ScanEntry> spec = new ArrayList<>(MZs.length); // Attempt to pre-define size of ArrayList to avoid incurring heavy copy costs when doing .add operation(s)
      for (int j = 0; j < MZs.length; ++j) {
        ScanEntry se = new ScanEntry(MZs[j], rt, intensities[j]);
        spec.add(se);
      }
      // find isolation window info
      if (isolation_windows.contains(win_range)) {
        int window_index = isolation_windows.indexOf(win_range);
        spectrumPoints.get(window_index).add(spec);
        fullScans.get(window_index).add(i_scan);
      } else {
        ArrayList<ArrayList<ScanEntry>> new_window = new ArrayList<>();
        new_window.add(spec);
        isolation_windows.add(win_range);
        spectrumPoints.add(new_window);
        ArrayList<IScan> new_ms2s = new ArrayList<>();
        new_ms2s.add(i_scan);
        fullScans.add(new_ms2s);
      }
    }
  }
}
