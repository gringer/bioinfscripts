package affytools;

import java.util.zip.*;
import java.util.Vector;
import java.util.HashMap;
import java.io.*;

/**
 * @author David Eccles (gringer)
 * @version 0.2
 *
 * Created 2008-Sep-23
 *
 * This program processes files, converting the affymetrix data into a simplegt
 * format. That is, convert from a n*m-line file with 4 columns (marker,
 * individual, genotype, value) to an m line file with n+1 columns (marker,
 * genotype(s)).
 *
 * Changed 2009-Sep-15
 * - added dots to report progress
 * - fixed bug checking array size when setting genotypes of individuals
 * - individuals not genotyped for a particular marker will now get
     '--' as their genotype (previously caused a null pointer
     exception)
 */
public class Affy2SimpleGT {
    /**
     * @param args
     *            List of files to process
     */
    static int     POS_GT     = 2;

    static int     POS_IND    = 1;

    static int     POS_MARKER = 0;

    static boolean DEBUG      = false;

    public static void main(String[] args) {
        long startTime = System.currentTimeMillis();
        long linesread = 0;
        boolean firstTime = true;
        if (Affy2SimpleGT.DEBUG) {
            System.out.println("Reading in files");
        }
        StreamTokenizer inputScanner = null;
        /* Assume 17 (or fewer) possible genotypes initially */
        Vector<String> gtVals = new Vector<String>(17);
        /*
         * Assume 10,000 (or fewer) individuals initially. HashMap allows quick
         * checking
         */
        HashMap<String, Integer> hIndVals = new HashMap<String, Integer>(10000);
        /* Reverse Map for easier/faster iteration over individuals */
        Vector<String> rIndVals = new Vector<String>(10000);
        /*
         * Assume 500,000 (or fewer) markers initially. HashMap allows quick
         * checking
         */
        HashMap<String, Integer> markVals = new HashMap<String, Integer>(500000);
        /* Reverse Map for easier/faster iteration over markers */
        Vector<String> rMarkVals = new Vector<String>(500000);
        /*
         * markGenotypes is a vector(markers) of vectors(individuals) of
         * genotypes
         */
        Vector<Vector<Integer>> markGenotypes = new Vector<Vector<Integer>>(
                500000);
        for (int argPos = 0; argPos < args.length; argPos++) {
            int tokenCounter = 0;
            System.err.println("Reading " + args[argPos] + " (one '.' per 1,000,000 lines)");
            startTime = System.currentTimeMillis();
            try {
                inputScanner = new StreamTokenizer(new InputStreamReader(
                        new GZIPInputStream(new FileInputStream(args[argPos]))));
            } catch (FileNotFoundException e) {
                System.out.println("Error: file not found '" + args[argPos]
                        + "'");
                System.exit(1);
            } catch (IOException e1) {
                if (e1.getMessage().contains("Not in GZIP format")) {
                    try {
                        inputScanner = new StreamTokenizer(
                                new InputStreamReader(new FileInputStream(
                                        args[argPos])));
                    } catch (FileNotFoundException e) {
                        System.out.println("Error: file not found '"
                                + args[argPos] + "'");
                        System.exit(1);
                    }
                } else {
                    e1.printStackTrace();
                }
            }
            /* Consider end-of-line to be a special token */
            inputScanner.eolIsSignificant(true);
            /*
             * Consider non-numeric printable ASCII characters to be word/number
             * characters
             */
            inputScanner.wordChars(33, 126);
            if (inputScanner != null) {
                try {
                    int tokenPos = 0;
                    int curMarker = -1;
                    int curIndiv = -1;
                    int tTokenID = StreamTokenizer.TT_EOL;
                    while (tTokenID != StreamTokenizer.TT_EOF) {
                        /*
                         * The general process is to read the input files, line
                         * by line, and store the individual/marker/genotype
                         * combinations into an array for layer output. Data are
                         * stored in a vector of vectors (top level: by marker,
                         * next level: by individual) with vector location being
                         * based on the order in which they have appeared in the
                         * input file(s). i.e. individuals will be unsorted, so
                         * it is necessary to generate a header indicating this
                         * ordering (although it is still a good idea to
                         * generate a header even when individuals are sorted).
                         */
                        try {
                            tTokenID = inputScanner.nextToken();
                        } catch (IOException e) {
                            /*
                             * workaround a Java bug caused by GZIP files being
                             * >2GB
                             * http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=4262583
                             * ... Of course, there's a problem here if the GZIP
                             * file actually <em>does</em> have a corrupt
                             * trailer.
                             */
                            if (e.getMessage().contains("Corrupt GZIP trailer")) {
                                tTokenID = StreamTokenizer.TT_EOF;
                            } else {
                                throw e;
                            }
                        }
                        /*
                         * Number: probably QC/reliability value or similar (not
                         * currently used)
                         */
                        if (tTokenID == StreamTokenizer.TT_NUMBER) {
                            if (Affy2SimpleGT.DEBUG) {
                                System.out.println("Read number [" + tokenPos
                                        + "]: " + inputScanner.nval);
                            }
                            tokenCounter++;
                            tokenPos++;
                        }
                        /*
                         * Word: Some combination of alphanumeric characters,
                         * hyphens, and other odd symbols. This is treated as
                         * either a marker name, individual name, or a genotype,
                         * depending on the field/column number. The only way to
                         * change this (currently) is to modify the POS_X static
                         * integers at the top of this code. A command-line
                         * method to alter these numbers could be easily
                         * implemented if the file format changes in the future.
                         */
                        if (tTokenID == StreamTokenizer.TT_WORD) {
                            if (Affy2SimpleGT.DEBUG) {
                                System.err.println("Read word [" + tokenPos
                                        + "]: " + inputScanner.sval);
                            }
                            /*
                             * A list of markers is maintained, with each marker
                             * (even across different input files) having its
                             * own unique numerical ID for referencing into the
                             * vector.
                             */
                            if (tokenPos == Affy2SimpleGT.POS_MARKER) {
                                Integer markObject = markVals
                                        .get(inputScanner.sval);
                                /*
                                 * If the marker has been seen before, use that
                                 * ID, otherwise generate a new one.
                                 */
                                int markID = (markObject == null) ? -1
                                        : markObject;
                                if (markID == -1) {
                                    rMarkVals.add(inputScanner.sval);
                                    markGenotypes.add(new Vector<Integer>());
                                    markID = rMarkVals.size() - 1;
                                    if (Affy2SimpleGT.DEBUG) {
                                        System.err.println("New Marker: "
                                                + inputScanner.sval + "["
                                                + markID + "]");
                                    }
                                    markVals.put(inputScanner.sval, markID);
                                } else {
                                    if (Affy2SimpleGT.DEBUG) {
                                        System.err.println("Found Marker: "
                                                + markVals.get(markID));
                                    }
                                }
                                curMarker = markID;
                            }
                            /*
                             * A list of individuals is also maintained, with
                             * each individual having its own unique numerical
                             * ID for referencing into the vector. It is
                             * expected that the individuals in the first file
                             * compose the entire set of individuals that have
                             * genotypes.
                             */
                            if (tokenPos == Affy2SimpleGT.POS_IND) {
                                Integer indObject = hIndVals
                                        .get(inputScanner.sval);
                                /*
                                 * If the individual has been seen before, use
                                 * that ID, otherwise generate a new one
                                 */
                                int indID = (indObject == null) ? -1
                                        : indObject;
                                if (indID == -1) {
                                    if (!firstTime) {
                                        /*
                                         * This is not a good situation to be
                                         * in: there is more than one input
                                         * file, and an individual has been seen
                                         * that was not observed in the first
                                         * input file. The procedure will
                                         * continue, with a warning, and
                                         * evidence (via an increase in the
                                         * number of genotypes on a line) that
                                         * this has happened.
                                         */
                                        System.err
                                                .println("Warning: new individual ("
                                                        + inputScanner.sval
                                                        + ") added after first file");
                                    }
                                    rIndVals.add(inputScanner.sval);
                                    indID = rIndVals.size() - 1;
                                    if (Affy2SimpleGT.DEBUG) {
                                        System.err.println("New Individual: "
                                                + inputScanner.sval + "["
                                                + indID + "]");
                                    }
                                    hIndVals.put(inputScanner.sval, indID);
                                }
                                curIndiv = indID;
                            }
                            if (tokenPos == Affy2SimpleGT.POS_GT) {
                                int gtID = gtVals.indexOf(inputScanner.sval);
                                if (gtID == -1) {
                                    gtID = gtVals.size();
                                    if (Affy2SimpleGT.DEBUG) {
                                        System.out.println("New Genotype: "
                                                + inputScanner.sval + "["
                                                + gtID + "]");
                                    }
                                    gtVals.add(inputScanner.sval);
                                } else {
                                    if (Affy2SimpleGT.DEBUG) {
                                        System.out.println("Found Genotype: "
                                                + gtVals.get(gtID));
                                    }
                                }
                                /*
                                 * The genotype is expected to be the last
                                 * "word" on the line. The individual/marker
                                 * ordering can change, however. This seems to
                                 * be typical for most genotype input files of
                                 * this type. A completely order-agnostic
                                 * implementation would store the genotype once
                                 * it had all three values, rather than at the
                                 * time of receipt of the genotype as is the
                                 * case here.
                                 */
                                if ((curIndiv == -1) || (curMarker == -1)) {
                                    System.err
                                            .println("Error: genotype before individual/marker on line "
                                                    + inputScanner.lineno()
                                                    + " of '"
                                                    + args[argPos]
                                                    + "'");
                                } else {
                                    Vector<Integer> tMarkGenotype = markGenotypes
                                            .get(curMarker);
                                    if (tMarkGenotype.size() <= curIndiv) {
                                        tMarkGenotype.setSize(curIndiv + 1);
                                    }
                                    if (tMarkGenotype.size() <= curIndiv) {
                                        System.err
                                            .println("Error: size is still less than individual number on line "
                                                     + inputScanner.lineno()
                                                     + " of '"
                                                     + args[argPos]
                                                     + "'");
                                    }
                                    tMarkGenotype.set(curIndiv, gtID);
                                }
                            }
                            tokenCounter++;
                            tokenPos++;
                        }
                        if (tTokenID == StreamTokenizer.TT_EOL) {
                            if((inputScanner.lineno() % 1000000) == 0){
                                System.err.print(".");
                            }
                            if (Affy2SimpleGT.DEBUG) {
                                System.err.println("[New Line]");
                            }
                            tokenPos = 0;
                            curMarker = -1;
                            curIndiv = -1;
                        }
                    }
                    /* Subtract 1 because the last line in a text file is (or should be) empty */
                    System.err.println((inputScanner.lineno() - 1) + " lines read");
                    /*
                     * A header containing identifiers for the individuals may
                     * be provided on the first line. Each population is
                     * identified by '## <Individual/Column IDs: ', followed by
                     * a space-separated list of individual labels, followed by ' >
                     * ##'. The header for multiple populations can be
                     * concatenated on the first line, possibly excluding the
                     * initial '##' (this means the GNU 'join' program can be
                     * used for combining two separate simplegt files.
                     */
                    if (firstTime) {
                        /*
                         * The header is only output once, at the end of the
                         * reading of the first input file. If there are more
                         * individuals in subsequent files that were not present
                         * in the first file, they will be nameless in the
                         * output.
                         */
                        System.out.print("## <Individual/Column IDs: ");
                        for (int indPos = 0; indPos < rIndVals.size(); indPos++) {
                            System.out.print(rIndVals.get(indPos) + " ");
                        }
                        System.out.println(" > ##");
                        firstTime = false;
                    }
                    /*
                     * This short loop generates the "simplegt" formatted file.
                     * Each line in the file is composed of a marker, followed
                     * by a space-separated list of genotypes.
                     */
                    System.err.print("Writing genotypes to output file (one '.' per 1000 markers)");
                    for (int markerPos = 0; markerPos < markGenotypes.size(); markerPos++) {
                        Vector<Integer> tMarkGenotype = markGenotypes.get(markerPos);
                        if(tMarkGenotype != null){
                            if(((markerPos+1) % 1000) == 0){
                                System.err.print(".");
                            }
                            /* rMarkVals is a vector of marker name
                             * strings, markGenotypes is where the
                             * actual genotypes are stored (as a
                             * vector of individual genotypes inside a
                             * vector of markers) */
                            String tMarker = rMarkVals.get(markerPos);
                            System.out.print(tMarker + " ");
                            for (int gtPos = 0; gtPos < tMarkGenotype.size(); gtPos++) {
                                Integer tGtVal = tMarkGenotype.get(gtPos);
                                if(tGtVal == null){
                                    /* Each un-genotyped individual
                                     * for a known marker gets '--' as
                                     * their genotype */
                                    System.err.println("Warning: no genotype for individual #" + gtPos + ", marker " + tMarker);
                                    System.out.print("--" + " ");
                                } else {
                                    System.out.print(gtVals.get(tMarkGenotype
                                                                .get(gtPos))
                                                     + " ");
                                }
                            }
                            System.out.println();
                            /* free memory for markers that have been output */
                            markGenotypes.setElementAt(null, markerPos);
                        }
                    }
                    System.err
                            .println("Completed in "
                                    + ((System.currentTimeMillis() - startTime) / 1000.0)
                                    + " seconds");
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }
}
