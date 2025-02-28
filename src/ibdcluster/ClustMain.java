/*
 * Copyright 2023-2025 Brian L. Browning
 *
 * This file is part of the ibd-cluster program.
 *
 * Licensed under the Apache License, Version 2.0 (the License);
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an AS IS BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package ibdcluster;

import blbutil.Const;
import blbutil.FileUtil;
import blbutil.Utilities;
import java.io.File;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Locale;

/**
 * <p>Class {@code ClustMain} contains the main() method for the ibd-cluster
 * program.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ClustMain {

    private static final String EXECUTABLE = "ibd-cluster.jar";
    private static final String PROGRAM = EXECUTABLE + "  [ version 0.2.0, 27Feb25.e01 ]";
    private static final String COPYRIGHT = "Copyright (C) 2023 Brian L. Browning";
    static final String COMMAND = "java -jar " + EXECUTABLE;

    private static final String HELP_MESSAGE = "Enter \"" + COMMAND
            + "\" to print a list of command line arguments";

    private ClustMain() {
        // private constructor to prevent instantiation
    }

    /**
     * Entry point to the ibd-cluster program.  See the {@code ClustPar}
     * class for details of the program arguments.
     *
     * @param args the command line arguments
     * @throws NullPointerException if {@code args = null}
     */
    public static void main(String[] args) {
	Locale.setDefault(Locale.US);
        ClustPar par = getPar(args);
        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism",
                String.valueOf(par.nthreads()));

        try (PrintWriter log = FileUtil.printWriter(new File(par.out() + ".log"))) {
            long t0 = System.nanoTime();
            Utilities.duoPrintln(log, startInfo(par));
            ClustStats stats = ClustAnalysis.run(par);
            Utilities.duoPrintln(log, statistics(par, stats));
            Utilities.duoPrintln(log, endInfo(t0));
        }
    }

    private static ClustPar getPar(String[] args) {
        if (args.length==0 || args[0].toLowerCase().startsWith("help")) {
            System.out.println(PROGRAM);
            System.out.println();
            System.out.println(ClustPar.usage());
            System.exit(0);
        }
        ClustPar par = new ClustPar(args);
        checkOutputPrefix(par);
        checkOutputFilename(par, par.out() + ".ibdclust.gz");
        checkOutputFilename(par, par.out() + ".log");
        return par;
    }

    private static void checkOutputPrefix(ClustPar par) {
        File outPrefix = new File(par.out());
        if (outPrefix.isDirectory()) {
            String err = "The out parameter cannot be a directory";
            String info = Const.nl + "Error      :  " + err
                    + Const.nl     + "Parameter  :  " + "out=" + par.out()
                    + Const.nl     + "Directory  :  " + outPrefix;
            Utilities.exit(new Throwable(err), info);
        }
    }

    private static void checkOutputFilename(ClustPar par, String filename) {
        File file = new File(filename);
        if (file.equals(par.gt()) || file.equals(par.map())
                || file.equals(par.excludesamples()) || file.equals(par.excludemarkers())) {
            String err = "An output file has the same name as an input file";
            String info = Const.nl + "Error      :  " + err
                    + Const.nl     + "Filename   :  " + filename;
            Utilities.exit(new Throwable(err), info);
        }
    }

    private static String startInfo(ClustPar par) {
        StringBuilder sb = new StringBuilder(300);
        sb.append(COPYRIGHT);
        sb.append(Const.nl);
        sb.append(HELP_MESSAGE);
        sb.append(Const.nl);
        sb.append(Const.nl);
        sb.append("Program             :  ");
        sb.append(PROGRAM);
        sb.append(Const.nl);
        sb.append("Start Time          :  ");
        sb.append(Utilities.timeStamp());
        sb.append(Const.nl);
        sb.append("Max Memory          :  ");
        long maxMemory = Runtime.getRuntime().maxMemory();
        if (maxMemory != Long.MAX_VALUE) {
            sb.append(maxMemory >> 30);
            sb.append(" GiB");
        }
        else {
            sb.append("[no limit])");
        }
        sb.append(Const.nl);
        sb.append(Const.nl);
        par.appendDocumentedParameters(sb);
        par.appendNonDefaultUndocumentedParameters(sb);
        return sb.toString();
    }

    private static String statistics(ClustPar par, ClustStats stats) {
        int nSamples = stats.nSamples();
        long nMarkers = stats.nMarkers();
        long nFilteredMarkers = stats.nFilteredMarkers();
        long nOutputPositions = stats.nOutputPositions();
        double filteredPercent = (100.0*nFilteredMarkers) / nMarkers;
        double ibdSetsPerLocus = (double) stats.nIbdSets() / nOutputPositions;
        DecimalFormat DF1 = new DecimalFormat("0.0");
        DecimalFormat DFE2 = new DecimalFormat("0.00E0");
        StringBuilder sb = new StringBuilder(300);
        sb.append(Const.nl);
        sb.append("Analysis summary");
        sb.append(Const.nl);
        sb.append("  samples           :  ");
        sb.append(stats.nSamples());
        sb.append(Const.nl);
        sb.append("  haplotypes        :  ");
        sb.append(nSamples << 1);
        sb.append(Const.nl);
        sb.append("  input VCF records :  ");
        sb.append(nMarkers);
        sb.append(Const.nl);
        sb.append("  filtered records  :  ");
        sb.append(nFilteredMarkers);
        sb.append("  (");
        sb.append( DF1.format(filteredPercent) );
        sb.append("% of records)");
        sb.append(Const.nl);
        sb.append("  output positions  :  ");
        sb.append(nOutputPositions);
        sb.append(Const.nl);
        sb.append("  clusters/position :  ");
        sb.append( (int) Math.rint(ibdSetsPerLocus));
        sb.append(Const.nl);
        double discordRate = stats.discordRate();
        sb.append("  discordance rate  :  ");
        sb.append(DFE2.format(discordRate));
        sb.append(Const.nl);
        return sb.toString();
    }

    private static String endInfo(long startNanoTime) {
        StringBuilder sb = new StringBuilder(300);
        long wallclockNanos = System.nanoTime() - startNanoTime;
        sb.append(Const.nl);
        sb.append("Wallclock Time:     :  ");
        sb.append(Utilities.elapsedNanos(wallclockNanos));
        sb.append(Const.nl);
        sb.append("End Time            :  ");
        sb.append(Utilities.timeStamp());
        return sb.toString();
    }
}
