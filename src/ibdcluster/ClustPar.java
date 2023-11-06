/*
 * Copyright 2023 Brian L. Browning
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

import beagleutil.ChromInterval;
import blbutil.Const;
import blbutil.Utilities;
import blbutil.Validate;
import java.io.File;
import java.util.Map;

/**
 * <p>Class {@code ClustPar} represents the parameters for an ibd-cluster
 * analysis.</p>
 *
 * <p>Instances of class {@code ClustPar} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class ClustPar {

    private final String[] args;

    private final File gt;
    private final File map;
    private final String out;
    private final String chrom;
    private final ChromInterval chromInt;
    private final File excludesamples;
    private final File excludemarkers;
    private final float min_maf;

    private final float length;
    private final float trim;
    private final float aggregate;  // max cM in an aggregate marker
    private final int nthreads;

    private static final float DEF_MIN_MAF = 0.1f;
    private static final float DEF_LENGTH = 2.0f;
    private static final float DEF_TRIM = 1.0f;
    private static final float DEF_AGGREGATE = 0.005f;
    private static final int   DEF_NTHREADS = Runtime.getRuntime().availableProcessors();

    /**
     * Constructs a new {@code ClustPar} instance from the specified
     * command line arguments.  The JVM will exit with an error message if the
     * {@code length} parameter is not at least {@code 1.9999f} times
     * the {@code trim} parameter.
     * @param args the command line arguments
     * @throws IllegalArgumentException if a command line argument
     * is incorrectly specified
     * @throws NumberFormatException if a numeric value for a parameter
     * is incorrectly specified
     * @throws NullPointerException if {@code args == null} or if there
     * exists {@code j} such that {@code (0 <= j && j < args.length)} and
     * {@code (args[j] == null)}
     */
    public ClustPar(String[] args) {
        int IMAX = Integer.MAX_VALUE;
        float FMIN = Float.MIN_VALUE;
        float FMAX = Float.MAX_VALUE;

        this.args = args.clone();
        Map<String, String> argsMap = Validate.argsToMap(args, '=');

        gt = Validate.getFile(
                Validate.stringArg("gt", argsMap, true, null, null));
        map = Validate.getFile(Validate.stringArg("map", argsMap, true, null, null));
        out = Validate.stringArg("out", argsMap, true, null, null);
        chrom = Validate.stringArg("chrom", argsMap, false, null, null);
        chromInt = parseChromInt(chrom);
        excludesamples = Validate.getFile(
                Validate.stringArg("excludesamples", argsMap, false, null, null));
        excludemarkers = Validate.getFile(
                Validate.stringArg("excludemarkers", argsMap, false, null, null));
        min_maf = Validate.floatArg("min-maf", argsMap, false, DEF_MIN_MAF, -FMAX, Math.nextDown(0.5f));

        length = Validate.floatArg("length", argsMap, false, DEF_LENGTH, FMIN, FMAX);
        trim = Validate.floatArg("trim", argsMap, false, DEF_TRIM, FMIN, FMAX);
        aggregate = Validate.floatArg("aggregate", argsMap, false, DEF_AGGREGATE, 0.0f, FMAX);
        nthreads = Validate.intArg("nthreads", argsMap, false, DEF_NTHREADS,
                1, IMAX);

        if (length < 1.99999f*trim) {
            System.out.println(usage());
            System.out.println();
            String s = "ERROR: the \"length\" parameter must be at least two times the \"trim\" parameter";
            Utilities.exit(s);
        }
        Validate.confirmEmptyMap(argsMap);
    }

    private static ChromInterval parseChromInt(String chrom) {
        ChromInterval ci = ChromInterval.parse(chrom);
        if (chrom!=null && ci==null) {
            throw new IllegalArgumentException("Invalid chrom parameter: " + chrom);
        }
        return ci;
    }

    /**
     * Returns the command line arguments.
     * @return the command line arguments
     */
    public String[] args() {
        return args.clone();
    }

    /**
     * Returns the ibd-cluster usage instructions.
     * @return the ibd-cluster usage instructions
     */
    public static String usage() {
        String nl = Const.nl;
        return  "Syntax: " + ClustMain.COMMAND + " [arguments in format: parameter=value]" + nl
                + nl
                + "  gt=<VCF file with phased genotypes>                (required)" + nl
                + "  map=<PLINK map file with cM units>                 (required)" + nl
                + "  out=<output file prefix>                           (required)" + nl
                + nl
                + "  chrom=< [chrom] or [chrom]:[start]-[end] >         (optional)" + nl
                + "  excludesamples=<file with 1 sample ID per line>    (optional)" + nl
                + "  excludemarkers=<file with 1 marker ID per line>    (optional)" + nl
                + "  min-maf=<min frequency of each non-major allele>   (default=" + DEF_MIN_MAF + ")" + nl
                + nl
                + "  length=<minimum IBS cM length>                     (default=" + DEF_LENGTH + ")" + nl
                + "  trim=<cM trimmed from each IBS segment end>        (default=" + DEF_TRIM + ")" + nl
                + "  aggregate=<maximum cM in an aggregate marker>      (default=" + DEF_AGGREGATE + ")" + nl
                + "  nthreads=<number of threads>                       (default: all CPU cores)" + nl
                + nl;
    }

    /**
     * Returns the gt parameter.
     * @return the gt parameter
     */
    public File gt() {
        return gt;
    }

    /**
     * Returns the map parameter or {@code null} if no map parameter was
     * specified.
     * @return the map parameter or {@code null} if no map parameter was
     * specified
     */
    public File map() {
        return map;
    }

    /**
     * Returns the out parameter.
     * @return the out parameter
     */
    public String out() {
        return out;
    }

    /**
     * Returns the chrom parameter or {@code null}
     * if no chrom parameter was specified.
     *
     * @return the chrom parameter or {@code null}
     * if no chrom parameter was specified
     */
    public String chrom() {
        return chrom;
    }

    /**
     * Returns the chromosome interval or {@code null} if no chrom
     * parameter was specified.
     *
     * @return the chromosome interval or {@code null} if no chrom
     * parameter was specified
     */
    public ChromInterval chromInt() {
        return chromInt;
    }

    /**
     * Returns the excludesamples parameter or {@code null}
     * if no excludesamples parameter was specified.
     *
     * @return the excludesamples parameter or {@code null}
     * if no excludesamples parameter was specified
     */
    public File excludesamples() {
        return excludesamples;
    }

    /**
     * Returns the excludemarkers parameter or {@code null}
     * if no excludemarkers parameter was specified.
     *
     * @return the excludemarkers parameter or {@code null}
     * if no excludemarkers parameter was specified
     */
    public File excludemarkers() {
        return excludemarkers;
    }

    /**
     * Returns the min-maf parameter.
     * @return the min-maf parameter
     */
    public float min_maf() {
        return min_maf;
    }

    /**
     * Returns the length parameter.
     * @return the length parameter
     */
    public float length() {
        return length;
    }

    /**
     * Returns the trim parameter.
     * @return the trim parameter
     */
    public float trim() {
        return trim;
    }

    /**
     * Returns the aggregate parameter.
     * @return the aggregate parameter
     */
    public float aggregate() {
        return aggregate;
    }

    /**
     * Returns the nthreads parameter.
     * @return the nthreads parameter
     */
    public int nthreads() {
        return nthreads;
    }
}
