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

import beagleutil.ChromInterval;
import blbutil.Const;
import blbutil.Validate;
import java.io.File;
import java.util.Map;

/**
 * <p>Class {@code ClustPar} represents the analysis parameters
 * for an ibd-ends analysis.</p>
 *
 * <p>Class {@code ClustPar} is immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class ClustPar {

    private final String[] args;

    // required parameters
    private final File gt;
    private final File map;
    private final String out;

    // optional parameters
    private final String chrom;
    private final ChromInterval chromInt;
    private final File excludesamples;
    private final File excludemarkers;

    private final float min_maf;
    private final float min_ibs_cm;
    private final float min_ibd_cm;
    private final int pbwt;
    private final float trim;
    private final float discord;
    private final float out_cm;
    private final int nthreads;

    private static final float DEF_MIN_MAF = 0.1f;
    private static final float DEF_MIN_IBS_CM = 1.0f;
    private static final float DEF_MIN_IBD_CM = 1.0f;
    private static final int DEF_PBWT = 4;
    private static final float DEF_TRIM = 0.5f;
    private static final float DEF_DISCORD = 0.0005f;
    private static final float DEF_OUT_CM = 0.02f;
    private static final int DEF_NTHREADS = Runtime.getRuntime().availableProcessors();

    // undocumented parameters
    private final float ne;              // effective population size
    private final float quantile;        // quantile of endpoint distribution
    private final int gc_bases;          // max bases in a gene conversion tract
    private final float gc_discord;      // allele discord probability for gc_bases following discordant alleles
    private final int local_segments;    // haplotype pairs used used to estimate local IBS length destribution
    private final float local_max_cdf;   // max local IBS length CDF probability
    private final int global_loci;       // loci used to estimated one-sided global IBS length distribution
    private final int global_segments;   // haplotype pairs per locus used to estimate one-sided global IBS length distribution
    private final float global_quantile; // global one-sided IBS length destribution quantile used for filtering loci
    private final float global_multiple; // max permitted multiple of mean (over loci) global_quantile value
    private final float min_cdf_ratio;   // min CDF ratio for two markers with discordant alleles.
    private final int max_its;           // max number of iterative updates to an IBD segment's end points
    private final float end_morgans;     // morgans between first/last marker and nominal preceding/next discordance
    private final boolean fix_focus;     // if true, do not iteratively update focus
    private final float prefocus_quantile; // IBD segment endpoint distribution quantile used for prefocus length
    private final float max_rel_change;  // max permitted relative change in pre-focus segment length
    private final int out_window_size;   // number of output positions per window
    private final long seed;             // seed for random number generation

    private static final float DEF_NE = 10000f;
    private static final float DEF_QUANTILE = 0.5f;
    private static final int DEF_GC_BASES = 1000;
    private static final float DEF_GC_DISCORD = 0.1f;
    private static final int DEF_LOCAL_SEGMENTS = 10000;
    private static final int MAX_LOCAL_SEGMENTS = 45000; // MAX_LOCAL_SEGMENTS*(MAX_LOCAL_SEGMENTS-1) <= Integer.MAX_VALUE
    private static final float DEF_LOCAL_MAX_CDF = 0.999f;
    private static final int DEF_GLOBAL_LOCI = 1000;
    private static final int DEF_GLOBAL_SEGMENTS = 2000;
    private static final float DEF_GLOBAL_QUANTILE = 0.9f;
    private static final float DEF_GLOBAL_MULTIPLE = 3.0f;
    private static final float DEF_MIN_CDF_RATIO = 0.001f;
    private static final int DEF_MAX_ITS = 2;
    private static final float DEF_END_MORGANS = 1.0f;
    private static final boolean DEF_FIX_FOCUS = false;
    private static final float DEF_PREFOCUS_QUANTILE = 0.05f;
    private static final float DEF_MAX_REL_CHANGE = 0.1f;
    private static final int DEF_OUT_WINDOW_SIZE = 500;
    private static final int DEF_SEED = -99999;

    /**
     * Constructs an {@code ClustPar} object that represents the
     * analysis parameters for an ibd-cluster analysis.  See the
     * {@code usage()} method for a description of the command line
     * parameters.
     *
     * @param args the command line arguments
     * @throws IllegalArgumentException if the command line arguments
     * are incorrectly specified
    */
    public ClustPar(String[] args) {
        float MIN_PROP = Float.MIN_VALUE;
        float MAX_PROP = Math.nextDown(1.0f);
        int IMAX = Integer.MAX_VALUE;
        float FMIN = Float.MIN_VALUE;
        float FMAX = Float.MAX_VALUE;
        long LMIN = Long.MIN_VALUE;
        long LMAX = Long.MAX_VALUE;

        Map<String, String> argsMap = Validate.argsToMap(args, '=');
        this.args = args.clone();

        // data parameters
        gt = Validate.getFile(Validate.stringArg("gt", argsMap, true, null, null));
        map = Validate.getFile(Validate.stringArg("map", argsMap, true, null, null));
        out = Validate.stringArg("out", argsMap, true, null, null);

        chrom = Validate.stringArg("chrom", argsMap, false, null, null);
        chromInt = parseChromInt(chrom);
        excludesamples = Validate.getFile(Validate.stringArg("excludesamples",
                argsMap, false, null, null));
        excludemarkers = Validate.getFile(
                Validate.stringArg("excludemarkers", argsMap, false, null, null));

        // algorithm parameters
        min_maf = Validate.floatArg("min-maf", argsMap, false, DEF_MIN_MAF, FMIN, 0.5f);
        min_ibs_cm = Validate.floatArg("min-ibs-cm", argsMap, false, DEF_MIN_IBS_CM, FMIN, FMAX);
        min_ibd_cm = Validate.floatArg("min-ibd-cm", argsMap, false, DEF_MIN_IBD_CM, FMIN, FMAX);
        pbwt = Validate.intArg("pbwt", argsMap, false, DEF_PBWT, 1, IMAX);
        trim = Validate.floatArg("trim", argsMap, false, DEF_TRIM, 0.0f, FMAX);
        discord = Validate.floatArg("discord", argsMap, false, DEF_DISCORD, MIN_PROP, MAX_PROP);
        out_cm = Validate.floatArg("out-cm", argsMap, false, DEF_OUT_CM, FMIN, FMAX);
        nthreads = Validate.intArg("nthreads", argsMap, false, DEF_NTHREADS, 1, IMAX);

        // undocumented parameters
        ne = Validate.floatArg("ne", argsMap, false, DEF_NE, 1, LMAX);
        quantile = Validate.floatArg("quantile", argsMap, false, DEF_QUANTILE, MIN_PROP, MAX_PROP);
        gc_bases = Validate.intArg("gc-bases", argsMap, false, DEF_GC_BASES, 0, IMAX);
        gc_discord = Validate.floatArg("gc-discord", argsMap, false, DEF_GC_DISCORD, MIN_PROP, MAX_PROP);
        local_segments = Validate.intArg("local-segments", argsMap, false, DEF_LOCAL_SEGMENTS, 1, MAX_LOCAL_SEGMENTS);
        local_max_cdf = Validate.floatArg("local-max-cdf", argsMap, false,
                DEF_LOCAL_MAX_CDF, MIN_PROP, MAX_PROP);
        global_loci = Validate.intArg("global-loci", argsMap, false, DEF_GLOBAL_LOCI, 1, IMAX);
        global_segments = Validate.intArg("global-segments", argsMap, false, DEF_GLOBAL_SEGMENTS, 1, IMAX);
        global_quantile = Validate.floatArg("global-quantile", argsMap, false,
                DEF_GLOBAL_QUANTILE, MIN_PROP, MAX_PROP);
        global_multiple = Validate.floatArg("global-multiple", argsMap, false,
                DEF_GLOBAL_MULTIPLE, FMIN, FMAX);
        min_cdf_ratio = Validate.floatArg("min-cdf-ratio", argsMap, false,
                DEF_MIN_CDF_RATIO, FMIN, FMAX);
        max_its = Validate.intArg("max-its", argsMap, false, DEF_MAX_ITS, 1, IMAX);
        end_morgans = Validate.floatArg("end-morgans", argsMap, false,
                DEF_END_MORGANS, 0.0f, Float.MAX_VALUE);

        fix_focus = Validate.booleanArg("fix-focus", argsMap, false, DEF_FIX_FOCUS);
        prefocus_quantile = Validate.floatArg("prefocus-quantile", argsMap, false,
                DEF_PREFOCUS_QUANTILE, MIN_PROP, MAX_PROP);
        max_rel_change = Validate.floatArg("max-rel-diff", argsMap, false, DEF_MAX_REL_CHANGE,
                MIN_PROP, MAX_PROP);
        out_window_size = Validate.intArg("out-window-size", argsMap, false, DEF_OUT_WINDOW_SIZE, 1, IMAX);
        seed = Validate.longArg("seed", argsMap, false, DEF_SEED, LMIN, LMAX);

        Validate.confirmEmptyMap(argsMap);
    }

    private static ChromInterval parseChromInt(String str) {
        ChromInterval chromInt = ChromInterval.parse(str);
        if (str!=null && str.length()>0 && chromInt==null) {
            throw new IllegalArgumentException("Invalid chrom parameter: " + str);
        }
        return chromInt;
    }

    /**
     * Returns the command line arguments.
     * @return the command line arguments
     */
    public String[] args() {
        return args.clone();
    }

    /**
     * Returns a string describing the command line arguments.
     * The format of the returned string is unspecified and subject to change.
     *
     * @return a string describing the command line arguments.
     */
    public static String usage() {
        String nl = Const.nl;
        return "Usage: " + ClustMain.COMMAND + " [parameters]" + nl
                + nl
                + "Required parameters: " + nl
                + nl
                + "  gt=<VCF file with phased genotypes>                   (required)" + nl
                + "  map=<PLINK map file with cM units>                    (required)" + nl
                + "  out=<output file prefix>                              (required)" + nl
                + nl
                + nl
                + "Optional parameters: " + nl
                + "  chrom=<[chrom] or [chrom]:[start]-[end]>              (optional)" + nl
                + "  excludesamples=<file with 1 sample ID per line>       (optional)" + nl
                + "  excludemarkers=<file with 1 marker ID per line>       (optional)" + nl
                + nl
                + "  min-maf=<min frequency of each non-major allele>      (default=" + DEF_MIN_MAF + ")" + nl
                + "  min-ibs-cm=<minimum IBS segment cM length>            (default="+ DEF_MIN_IBS_CM + ")" + nl
                + "  min-ibd-cm=<minimum IBD segment cM length>            (default="+ DEF_MIN_IBD_CM + ")" + nl
                + "  pbwt=<number of interleaved PBWT analyses>            (default="+ DEF_PBWT + ")" + nl
                + "  trim=<cM trimmed from each IBD segment end>           (default=" + DEF_TRIM + ")" + nl
                + "  discord=<allele discord probability in IBD segment>   (default: " + DEF_DISCORD + ")" + nl
                + "  out-cm=<cM between output positions>                  (default: " + DEF_OUT_CM + ")" + nl
                + "  nthreads=<number of computational threads>            (default: all CPU cores)" + nl
                + nl;
    }

    // data input/output parameters

    /**
     * Returns the gt parameter.
     * @return the gt parameter
     */
    public File gt() {
        return gt;
    }

    /**
     * Returns the map parameter.
     * @return the map parameter
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
     * if no excludemarkeres parameter was specified
     */
    public File excludemarkers() {
        return excludemarkers;
    }

    // algorithm parameters

    /**
     * Returns the min-maf parameter.
     * @return the min-maf parameter
     */
    public float min_maf() {
        return min_maf;
    }

    /**
     * Returns the min-ibs-cm parameter.
     * @return the min-ibs-cm parameter
     */
    public float min_ibs_cm() {
        return min_ibs_cm;
    }

    /**
     * Returns the min-ibd-cm parameter.
     * @return the min-ibd-cm parameter
     */
    public float min_ibd_cm() {
        return min_ibd_cm;
    }

    /**
     * Returns the pbwt parameter.
     * @return the pbwt parameter
     */
    public int pbwt() {
        return pbwt;
    }


    /**
     * Returns the trim parameter
     * @return the trim parameter
     */
    public float trim() {
        return trim;
    }

    /**
     * Returns the discord parameter.
     * @return the discord parameter
     */
    public float discord() {
        return discord;
    }

    /**
     * Returns the out-cm parameter.
     * @return the out-cm parameter
     */
    public float out_cm() {
        return out_cm;
    }

    /**
     * Returns the nthreads parameter.
     * @return the nthreads parameter
     */
    public int nthreads() {
        return nthreads;
    }

    // undocumented parameters

    /**
     * Returns the ne parameter.
     * @return the ne parameter
     */
    public float ne() {
        return ne;
    }

    /**
     * Returns the quantile parameter
     * @return the quantile parameter
     */
    public float quantile() {
        return quantile;
    }
    /**
     * Returns the gc-bases parameter.
     * @return the gc-bases parameter
     */
    public int gc_bases() {
        return gc_bases;
    }

    /**
     * Returns the gc-discord parameter.
     * @return the gc-discord parameter
     */
    public float gc_discord() {
        return gc_discord;
    }

    /**
     * Returns the local-segments parameter.
     * @return the local-segments parameter
     */
    public int local_segments() {
        return local_segments;
    }

    /**
     * Returns the maximum permitted number of local segments.
     *
     * @return the maximum permitted number of local segments
     */
    public static int max_local_segments() {
        return MAX_LOCAL_SEGMENTS;
    }

    /**
     * Returns the local-max-cdf parameter.
     * @return the local-max-cdf parameter
     */
    public float local_max_cdf() {
        return local_max_cdf;
    }

    /**
     * Returns the global-loci parameter.
     * @return the global-loci parameter
     */
    public int global_loci() {
        return global_loci;
    }

    /**
     * Returns the global-segments parameter.
     * @return the global-segments parameter
     */
    public int global_segments() {
        return global_segments;
    }

    /**
     * Returns the global-quantile parameter.
     * @return the global-quantile parameter
     */
    public float global_quantile() {
        return global_quantile;
    }

    /**
     * Returns the global-multiple parameter.
     * @return the global-multiple parameter
     */
    public float global_multiple() {
        return global_multiple;
    }

    /**
     * Returns the min-cdf-ratio parameter.
     * @return the min-cdf-ratio parameter
     */
    public float min_cdf_ratio() {
        return min_cdf_ratio;
    }

    /**
     * Returns the max-its parameter.
     * @return the max-its parameter
     */
    public int max_its() {
        return max_its;
    }

    /**
     * Returns the end-morgans parameter.
     * @return the end-morgans parameter
     */
    public float end_morgans() {
        return end_morgans;
    }

    /**
     * Returns the fix-focus parameter.
     * @return the fix-focus parameter
     */
    public boolean fix_focus() {
        return fix_focus;
    }

    /**
     * Returns the prefocus-quantile parameter.
     * @return the prefocus-quantile parameter
     */
    public float prefocus_quantile() {
        return prefocus_quantile;
    }

    /**
     * Returns the max-rel-change parameter.
     * @return the max-rel-change parameter
     */
    public float max_rel_change() {
        return max_rel_change;
    }

    /**
     * Returns the out-window-size parameter
     * @return the out-window-size parameter
     */
    public int out_window_size() {
        return out_window_size;
    }

    /**
     * Returns the seed parameter.
     * @return the seed parameter
     */
    public long seed() {
        return seed;
    }

    void appendDocumentedParameters(StringBuilder sb) {
        sb.append("Parameters");
        sb.append(Const.nl);
        sb.append("  gt                :  ");
        sb.append(gt);
        sb.append(Const.nl);
        sb.append("  map               :  ");
        sb.append(map!=null ? map : " [1 cM = 1 Mb]");
        sb.append(Const.nl);
        sb.append("  out               :  ");
        sb.append(out);
        if (chrom != null) {
            sb.append(Const.nl);
            sb.append("  chrom             :  ");
            sb.append(chrom);
        }
        if (excludesamples!=null) {
            sb.append(Const.nl);
            sb.append("  excludesamples    :  ");
            sb.append(excludesamples);
        }
        if (excludemarkers!=null) {
            sb.append(Const.nl);
            sb.append("  excludemarkers    :  ");
            sb.append(excludemarkers);
        }
        sb.append(Const.nl);
        sb.append("  min-maf           :  ");
        sb.append(min_maf);
        sb.append(Const.nl);
        sb.append("  min-ibs-cm        :  ");
        sb.append(min_ibs_cm);
        sb.append(Const.nl);
        sb.append("  min-ibd-cm        :  ");
        sb.append(min_ibd_cm);
        sb.append(Const.nl);
        sb.append("  pbwt              :  ");
        sb.append(pbwt);
        sb.append(Const.nl);
        sb.append("  trim              :  ");
        sb.append(trim);
        sb.append(Const.nl);
        sb.append("  discord           :  ");
        sb.append(discord);
        sb.append(Const.nl);
        sb.append("  out-cm            :  ");
        sb.append(out_cm);
        sb.append(Const.nl);
        sb.append("  nthreads          :  ");
        sb.append(nthreads);
    }

    void appendNonDefaultUndocumentedParameters(StringBuilder sb) {
        if (ne != DEF_NE) {
            sb.append(Const.nl);
            sb.append("  ne                :  ");
            sb.append(ne);
        }
        if (quantile != DEF_QUANTILE) {
            sb.append("  quantile          :  ");
            sb.append(quantile);
            sb.append(Const.nl);
        }
        if (gc_bases != DEF_GC_BASES) {
            sb.append("  gc-bases          :  ");
            sb.append(gc_bases);
            sb.append(Const.nl);
        }
        if (gc_discord != DEF_GC_DISCORD) {
            sb.append("  gc-discord        :  ");
            sb.append(gc_discord);
            sb.append(Const.nl);
        }
        if (local_segments != DEF_LOCAL_SEGMENTS) {
            sb.append(Const.nl);
            sb.append("  local-haps        :  ");
            sb.append(local_segments);
        }
        if (local_max_cdf != DEF_LOCAL_MAX_CDF) {
            sb.append(Const.nl);
            sb.append("  local-max-cdf     :  ");
            sb.append(local_max_cdf);
        }
        if (global_loci != DEF_GLOBAL_LOCI) {
            sb.append(Const.nl);
            sb.append("  global-pos        :  ");
            sb.append(global_loci);
        }
        if (global_segments != DEF_GLOBAL_SEGMENTS) {
            sb.append(Const.nl);
            sb.append("  global-segments   :  ");
            sb.append(global_segments);
        }
        if (global_quantile != DEF_GLOBAL_QUANTILE) {
            sb.append(Const.nl);
            sb.append("  global-quantile   :  ");
            sb.append(global_quantile);
        }
        if (global_multiple != DEF_GLOBAL_MULTIPLE) {
            sb.append(Const.nl);
            sb.append("  global-multiple   :  ");
            sb.append(global_multiple);
        }
        if (min_cdf_ratio != DEF_MIN_CDF_RATIO) {
            sb.append(Const.nl);
            sb.append("  min-cdf-ratio     :  ");
            sb.append(min_cdf_ratio);
        }
        if (max_its != DEF_MAX_ITS) {
            sb.append(Const.nl);
            sb.append("  max-its           :  ");
            sb.append(max_its);
        }
        if (end_morgans != DEF_END_MORGANS) {
            sb.append(Const.nl);
            sb.append("  end-morgans       :  ");
            sb.append(end_morgans);
        }
        if (fix_focus != DEF_FIX_FOCUS) {
            sb.append(Const.nl);
            sb.append("  fix-focus         :  ");
            sb.append(fix_focus);
        }
        if (prefocus_quantile != DEF_PREFOCUS_QUANTILE) {
            sb.append(Const.nl);
            sb.append("  prefocus-quantile :  ");
            sb.append(prefocus_quantile);
        }
        if (max_rel_change != DEF_MAX_REL_CHANGE) {
            sb.append(Const.nl);
            sb.append("  max-rel-change    :  ");
            sb.append(max_rel_change);
        }
        if (out_window_size != DEF_OUT_WINDOW_SIZE) {
            sb.append(Const.nl);
            sb.append("  out-window-ize    :  ");
            sb.append(out_window_size);
        }
        if (seed != DEF_SEED) {
            sb.append(Const.nl);
            sb.append("  seed              :  ");
            sb.append(seed);
        }
    }
}