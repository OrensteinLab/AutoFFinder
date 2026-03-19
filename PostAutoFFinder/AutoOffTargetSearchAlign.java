package PostAutoFFinder;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.nio.file.Path;
import java.nio.file.Files;
import java.io.File;
import java.io.IOException;
import java.io.FileWriter;
import java.io.BufferedReader;
import java.io.FileReader;

import java.time.Instant;
import java.time.Duration;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;


public class AutoOffTargetSearchAlign {
    private static final String CHROMOSOME_STR = "Chromosome";
    private static final String END_POS_STR = "EndPosition";
    private static final String SEQ_STR = "SiteSeqPlusMaxEditsBefore";
    private static final String TARGET_STR = "Target";
    private static final String EDIT_NUM_STR = "#Edit";
    private static final String STRAND_STR = "Strand";
    private static final String ALIGNED_TARGET_STR = "AlignedTarget";
    private static final String ALIGNED_TEXT_STR = "AlignedText";
    private static final String MISMATCHES_NUM_STR = "#Mismatches";
    private static final String BUGLES_NUM_STR = "#Bulges";
    private static String OUTPUT_PATH = "./output";
    private static int NUM_THREADS = 8;
    private static int SITE_WINDOW_SIZE = 50;
    private static int MAX_EDITS = 7;
    private static int MAX_MISMATCHES_WITHOUT_BULGES = MAX_EDITS;
    private static int MAX_MISMATCHES_WITH_BULGES = 6;
    private static int MAX_BULGES = 1;
    private static int EFFECTIVE_MAX_EDIT_WITH_BULGE = Math.min(MAX_EDITS, MAX_MISMATCHES_WITH_BULGES + MAX_BULGES);
    private static boolean ALLOW_PAM_EDITS = false;
    // shared DP matrix for all threads to reduce memory allocation, since we only compute one alignment at a time, we can reuse the same matrix
    private static final ThreadLocal<int[][]> sharedDpMatrix = ThreadLocal.withInitial(() -> new int[60][60]);
    

    public static void setOutputPath(String value) {
        OUTPUT_PATH = value;
    }

    public static void setNumThreads(int value) {
        NUM_THREADS = value;
    }

    public static void setSiteWindowSize(int value) {
        SITE_WINDOW_SIZE = value;
    }

    public static void setMaxEdits(int value) {
        MAX_EDITS = value;
        EFFECTIVE_MAX_EDIT_WITH_BULGE = Math.min(MAX_EDITS, MAX_MISMATCHES_WITH_BULGES + MAX_BULGES);
    }

    public static void setMaxMismatchesWithoutBulges(int value) {
        MAX_MISMATCHES_WITHOUT_BULGES = value;
    }

    public static void setMaxMismatchesWithBulges(int value) {
        MAX_MISMATCHES_WITH_BULGES = value;
        EFFECTIVE_MAX_EDIT_WITH_BULGE = Math.min(MAX_EDITS, MAX_MISMATCHES_WITH_BULGES + MAX_BULGES);
    }

    public static void setMaxBulges(int value) {
        MAX_BULGES = value;
        EFFECTIVE_MAX_EDIT_WITH_BULGE = Math.min(MAX_EDITS, MAX_MISMATCHES_WITH_BULGES + MAX_BULGES);
    }

    public static void setAllowPamEdits(boolean value) {
        ALLOW_PAM_EDITS = value;
    }

    /**
     * Creates a FileWriter to write the output file for the given target.
     * The output file will be a CSV file with a header row containing specific
     * column names.
     *
     * @param target The target string used to name the output file.
     * @return A FileWriter object for the specified output file.
     * @throws IOException If an I/O error occurs while creating the FileWriter.
     */
    public static FileWriter writeOutputFile() throws IOException {
        FileWriter writer = new FileWriter(OUTPUT_PATH + ".csv");
        String[] cols = new String[] { CHROMOSOME_STR, STRAND_STR, END_POS_STR, TARGET_STR, SEQ_STR, EDIT_NUM_STR,
                ALIGNED_TARGET_STR, ALIGNED_TEXT_STR, MISMATCHES_NUM_STR, BUGLES_NUM_STR };
        for (int colIndex = 0; colIndex < cols.length; colIndex++) {
            writer.append(cols[colIndex]);
            if (colIndex < (cols.length - 1)) {
                writer.append(",");
            }
        }
        writer.append("\n");

        return writer;
    }

    private static boolean charactersMatch(char target, char text, boolean allowNsInText) {
        if (allowNsInText) {
            return target == text || target == 'N' || text == 'N';
        } else {
            return target == text || target == 'N';
        }
    }
    
    private static boolean pamRegionsMatch(String target, String text, String pam, boolean allowNsInText) {
        int pamSize = pam.length();
        if (target.length() < pamSize || text.length() < pamSize) {
            return false;
        }
        
        String targetPam = target.substring(target.length() - pamSize);
        String textPam = text.substring(text.length() - pamSize);
        
        for (int i = 0; i < pamSize; i++) {
            if (!charactersMatch(targetPam.charAt(i), textPam.charAt(i), allowNsInText)) {
                return false;
            }
        }
        return true;
    }
    
    public static Alignment find_alignment(
        String pam, Boolean allowNsInText, Boolean allowPamEdits, String target, String text) {
        // Handle PAM edits constraint
        String workingTarget = target;
        String workingText = text;
        String targetPamSuffix = "";
        String textPamSuffix = "";
        
        if (!allowPamEdits) {
            // Check if PAM regions match
            if (!pamRegionsMatch(target, text, pam, allowNsInText)) {
                return null;
            }
            
            // Extract PAM suffix and work with truncated sequences
            int pamSize = pam.length();
            targetPamSuffix = target.substring(target.length() - pamSize);
            textPamSuffix = text.substring(text.length() - pamSize);
            workingTarget = target.substring(0, target.length() - pamSize);
            workingText = text.substring(0, text.length() - pamSize);
        }
        
        return computeAlignment(workingTarget, workingText, allowNsInText, targetPamSuffix, textPamSuffix);
    }

    private static void smithWatermanLastCellLogic(
        int[][] M, String target,String text, int row, int col, Boolean allowNsInText) {
        int matchCell = M[row - 1][col - 1];
        int targetBulgeCell = M[row - 1][col];
        int textBulgeCell = M[row][col - 1];
        if (allowNsInText) {
            M[row][col] = Math.min(Math.min(
                (target.charAt(row - 1) != text.charAt(col - 1) &&
                 target.charAt(row - 1) != 'N' &&
                 text.charAt(col - 1) != 'N') ? matchCell + 1 : matchCell,
                 targetBulgeCell + 1), textBulgeCell + 1
                );
        }
        else {
            M[row][col] = Math.min(Math.min(
                (target.charAt(row - 1) != text.charAt(col - 1) &&
                 target.charAt(row - 1) != 'N') ? matchCell + 1 : matchCell,
                 targetBulgeCell + 1), textBulgeCell + 1
                );
        }
    }

    private static void smithWatermanRowsFill(int[][] M, int m, int n,  String target, String text, Boolean allowNsInText) {
        for (int row = 1; row <= m; row++) {
            for (int col = 1; col <= n; col++) {
                smithWatermanLastCellLogic(M, target,text, row, col, allowNsInText);
            }
        }
    }

    private static Alignment computeAlignment(
        String target, String text, boolean allowNsInText, String targetPamSuffix, String textPamSuffix) {
        int targetLen = target.length();
        int textLen = text.length();
        int[][] dp = sharedDpMatrix.get();
        // Auto-resize if sequences are unexpectedly large - should not happen
        if (dp.length < targetLen + 1 || dp[0].length < textLen + 1) {
            dp = new int[Math.max(dp.length, targetLen + 10)][Math.max(dp[0].length, textLen + 10)];
            sharedDpMatrix.set(dp);
        }
        
        // fill the first col
        for (int row = 1; row <= targetLen; row++) {
            dp[row][0] = row;
        }
        smithWatermanRowsFill(dp, targetLen, textLen, target, text, allowNsInText);

        // Find best alignment ending at target[targetLen] and text[textLen] (enitre target aligned and last text character aligned)
        int bestScore = dp[targetLen][textLen];
        if (bestScore > MAX_EDITS) {
            return null; // No valid alignment found
            // System.out.println(target + " " + text + " " + bestScore + "  " + EndPosition + " " + Arrays.toString(dp[targetLen]));
            // System.exit(1);
        }
        
        // TODO: The problem here is that in case of not allowPamEdits, I can filter the last end position and therefore miss this site
        // if (dp[targetLen][textLen - 1] < bestScore) {
        //     return null; // position textLen - 1 should have been aligned in other part
        // }

        // We use the naive post-processing alignment finder
        return naive_find_alignment(
            dp, allowNsInText, target, text, targetLen, textLen,
            new StringBuilder[] { new StringBuilder(), new StringBuilder() }, 0, 0,
            targetPamSuffix, textPamSuffix);
    }


    public static Alignment naive_find_alignment(
            int[][] M, Boolean allowNsInText, String target, String text, int target_i, int text_i,
            StringBuilder[] targetTextAlign, int mismatches, int bulges, String targetPamSuffix, String textPamSuffix) {
        int targetLen = target.length();

        if (bulges > MAX_BULGES) {
            return null;
        }

        if (bulges > 0 && mismatches > MAX_MISMATCHES_WITH_BULGES) {
            return null;
        }
        if (mismatches > MAX_MISMATCHES_WITHOUT_BULGES) {
            return null;
        }
        if (target_i == 0) {
            // Reverse the alignment strings
            StringBuilder finalTargetAlign = new StringBuilder(targetTextAlign[0]).reverse();
            StringBuilder finalTextAlign = new StringBuilder(targetTextAlign[1]).reverse();
            // Handle PAM suffix correctly
            if (!targetPamSuffix.isEmpty()) {
                finalTargetAlign.append(targetPamSuffix);
                finalTextAlign.append(textPamSuffix);
            }
            return new Alignment(mismatches, bulges, finalTargetAlign, finalTextAlign);
        }

        if (text_i <= 0) {
            return null;
        }
        int effective_max_edit = bulges == 0?  MAX_EDITS : EFFECTIVE_MAX_EDIT_WITH_BULGE;
        if (M[target_i][text_i] > (effective_max_edit - mismatches - bulges)) {
            return null;
        }
        
        // store the current length of the alignment strings to backtrack after recursive calls
        int len0 = targetTextAlign[0].length();
        int len1 = targetTextAlign[1].length();
        
        // mismatch
        targetTextAlign[0].append(target.charAt(target_i - 1));
        targetTextAlign[1].append(text.charAt(text_i - 1));
        // compute the mismatch score
        int mismatchSocre;
        if (allowNsInText) {
            mismatchSocre = (target.charAt(target_i - 1) != text.charAt(text_i - 1) &&
                    target.charAt(target_i - 1) != 'N' && text.charAt(text_i - 1) != 'N') ? 1 : 0;
        } else {
            mismatchSocre = (target.charAt(target_i - 1) != text.charAt(text_i - 1) &&
                    target.charAt(target_i - 1) != 'N') ? 1 : 0;
        }
    
        Alignment targetTextMissOrMatchAlignment = naive_find_alignment(M, allowNsInText, target, text, target_i - 1,
                text_i - 1, targetTextAlign, mismatches + mismatchSocre, bulges, targetPamSuffix, textPamSuffix);
        targetTextAlign[0].setLength(len0);
        targetTextAlign[1].setLength(len1);

        // target bulge
        targetTextAlign[0].append(target.charAt(target_i - 1));
        targetTextAlign[1].append('-');
        Alignment targetTextTargetBugleAlignment = naive_find_alignment(M, allowNsInText, target, text, target_i - 1,
                text_i, targetTextAlign, mismatches, bulges + 1, targetPamSuffix, textPamSuffix);
        targetTextAlign[0].setLength(len0);
        targetTextAlign[1].setLength(len1);

        // text bulge
        Alignment targetTextnTextBugleAlignment = null;
        if (target_i != targetLen || !targetPamSuffix.isEmpty()) {
            // we do to have text bulge in the begining/start of the alignment
            // note that it impossible to put text bulge when target_i == 0, unless we have a PAM suffix
            targetTextAlign[0].append('-');
            targetTextAlign[1].append(text.charAt(text_i - 1));
            targetTextnTextBugleAlignment = naive_find_alignment(M, allowNsInText, target, text, target_i,
                    text_i - 1, targetTextAlign, mismatches, bulges + 1, targetPamSuffix, textPamSuffix);
            targetTextAlign[0].setLength(len0);
            targetTextAlign[1].setLength(len1);
        }

        // choose the best alignment, prefering mismatches over bulges
        if (targetTextMissOrMatchAlignment != null || targetTextTargetBugleAlignment != null
                || targetTextnTextBugleAlignment != null) {
            int missOrMatchEdit = targetTextMissOrMatchAlignment != null
                    ? targetTextMissOrMatchAlignment.getBulges() + targetTextMissOrMatchAlignment.getMismatches()
                    : MAX_EDITS;
            int targetBugleEdit = targetTextTargetBugleAlignment != null
                    ? targetTextTargetBugleAlignment.getBulges() + targetTextTargetBugleAlignment.getMismatches()
                    : MAX_EDITS;
            int textBugleEdit = targetTextnTextBugleAlignment != null
                    ? targetTextnTextBugleAlignment.getBulges() + targetTextnTextBugleAlignment.getMismatches()
                    : MAX_EDITS;
            int minEdit = Math.min(missOrMatchEdit, Math.min(targetBugleEdit, textBugleEdit));
            if (targetTextMissOrMatchAlignment != null && missOrMatchEdit == minEdit) {
                return targetTextMissOrMatchAlignment;
            }
            if (targetTextTargetBugleAlignment != null && targetBugleEdit == minEdit) {
                return targetTextTargetBugleAlignment;
            }
            return targetTextnTextBugleAlignment;
        }
        return null;
    }

    private static void automataProcessAlignment(
            String target, Alignment targetTextAlignment, int targetTextAlignmentPos, int targetTextAlignmentEdit, String strand,
            List<Alignment> alignmentList, List<Integer> endPosList, List<Integer> editNumList, List<String> targetList) {
        // TODO: Maybe here instead on each one: rereverse the aligned target and text as the builded sequence is reversed
        targetTextAlignment.setAlignedTarget(targetTextAlignment.getAlignedTargetBuilder().toString());
        targetTextAlignment.setAlignedTargetBulider(null);
        targetTextAlignment.setAlignedText(targetTextAlignment.getAlignedTextBuilder().toString());
        targetTextAlignment.setAlignedTextBulider(null);
        // add to alignment list
        alignmentList.add(targetTextAlignment);
        endPosList.add(targetTextAlignmentPos);
        editNumList.add(targetTextAlignmentEdit);
        // add to target list
        targetList.add(target);
    }

    private static Alignment automataPostprocessing(
            int textEndPosition, String strand, String target, String pam, String text, Boolean allowNsInText, Boolean chooseBestInWindow,
            List<Alignment> alignmentList, List<Integer> endPosList, List<Integer> editNumList, List<String> targetList,
            StringBuilder[] targetTextEmptyAlign, Alignment targetTextAlignment, int targetTextAlignmentPos, int targetTextAlignmentEdit) {
        Alignment targetTextAlignmentTemp = find_alignment(
            pam, allowNsInText, ALLOW_PAM_EDITS,
            target,
            text.substring(Math.max(textEndPosition - target.length() - MAX_EDITS, 0), textEndPosition));

        if (targetTextAlignmentTemp != null) {
            if (targetTextAlignment == null) {
                return targetTextAlignmentTemp;
            }
            if (((textEndPosition - targetTextAlignmentPos) > SITE_WINDOW_SIZE) || !chooseBestInWindow) {
                automataProcessAlignment(
                        target, targetTextAlignment, targetTextAlignmentPos, targetTextAlignmentEdit, strand, alignmentList,
                        endPosList, editNumList, targetList);

                // return the next possible alignment
                return targetTextAlignmentTemp;
            }
            int targetTextAlignmentTempMisBulge = targetTextAlignmentTemp.getMismatches()
                    + targetTextAlignmentTemp.getBulges();
            int targetTextAlignmentMisBulge = targetTextAlignment.getMismatches() + targetTextAlignment.getBulges();

            if (targetTextAlignmentTempMisBulge < targetTextAlignmentMisBulge) {
                return targetTextAlignmentTemp;
            }
            if ((targetTextAlignmentTempMisBulge == targetTextAlignmentMisBulge) &&
                    (targetTextAlignmentTemp.getBulges() < targetTextAlignment.getBulges())) {
                return targetTextAlignmentTemp;
            }
        }
        return null;
    }

    public static class AutomataResultsProcessFileHandler implements Runnable {
        private FileWriter writer;
        private String text;
        private List<String> targets;
        private String pam;
        private String strand;
        private Boolean allowNsInText;
        private Boolean chooseBestInWindow;
        private List<Integer> endPositions;
        private List<Integer> targetIds;
        private String chr;

        public AutomataResultsProcessFileHandler(
                FileWriter writer, String text, List<String> targets, String pam,
                String strand, Boolean allowNsInText, Boolean chooseBestInWindow,
                PositionsRes endPositions, String chr) {
            this.writer = writer;
            this.text = text;
            this.targets = targets;
            this.pam = pam;
            this.strand = strand;
            this.allowNsInText = allowNsInText;
            this.chooseBestInWindow = chooseBestInWindow;
            this.endPositions = endPositions.positions;
            this.targetIds = endPositions.ids;
            this.chr = chr;
        }

        /**
         * Container class to hold alignment tracking state for each target
         */
        private static class AlignmentTracker {
            List<Alignment> alignments;
            List<Integer> positions;
            List<Integer> edits;
            
            AlignmentTracker(int targetCount) {
                alignments = new ArrayList<>(targetCount);
                positions = new ArrayList<>(targetCount);
                edits = new ArrayList<>(targetCount);
                
                for (int i = 0; i < targetCount; i++) {
                    alignments.add(null);
                    positions.add(0);
                    edits.add(0);
                }
            }
            
            Alignment getAlignment(int targetId) {
                return alignments.get(targetId);
            }
            
            Integer getPosition(int targetId) {
                return positions.get(targetId);
            }
            
            Integer getEdit(int targetId) {
                return edits.get(targetId);
            }
            
            void update(int targetId, Alignment alignment, int position, int edit) {
                alignments.set(targetId, alignment);
                positions.set(targetId, position);
                edits.set(targetId, edit);
            }
            
            void reset(int targetId) {
                update(targetId, null, 0, 0);
            }
        }
        
        /**
         * Container class to hold collected off-target results
         */
        private static class OffTargetResults {
            List<Integer> endPosList = new ArrayList<>();
            List<Integer> editNumList = new ArrayList<>();
            List<Alignment> alignmentList = new ArrayList<>();
            List<String> targetList = new ArrayList<>();
        }
        
        /**
         * Processes a single position and updates the alignment tracker
         */
        private void processPosition(
                int targetId, String target, int endPos, 
                AlignmentTracker tracker, OffTargetResults results) {
            
            Alignment currentAlignment = tracker.getAlignment(targetId);
            Integer currentPos = tracker.getPosition(targetId);
            Integer currentEdit = tracker.getEdit(targetId);
            
            StringBuilder[] emptyAlign = new StringBuilder[2];
            emptyAlign[0] = new StringBuilder();
            emptyAlign[1] = new StringBuilder();
            
            Alignment newAlignment = automataPostprocessing(
                    endPos, strand, target, pam, text, allowNsInText, chooseBestInWindow,
                    results.alignmentList, results.endPosList, results.editNumList, results.targetList,
                    emptyAlign, currentAlignment, currentPos, currentEdit);
            
            if (newAlignment != null) {
                int totalEdits = newAlignment.getBulges() + newAlignment.getMismatches();
                tracker.update(targetId, newAlignment, endPos, totalEdits);
            } else if (currentAlignment != null) {
                automataProcessAlignment(
                        target, currentAlignment, currentPos, currentEdit, strand,
                        results.alignmentList, results.endPosList, results.editNumList, results.targetList);
                tracker.reset(targetId);
            }
        }
        
        /**
         * Processes all end positions and collects alignments
         */
        private OffTargetResults processAllPositions() {
            OffTargetResults results = new OffTargetResults();
            AlignmentTracker tracker = new AlignmentTracker(targets.size());
            
            // Process each candidate position
            for (int i = 0; i < endPositions.size(); i++) {
                int targetId = targetIds.get(i);
                String target = targets.get(targetId);
                int endPos = endPositions.get(i);
                
                processPosition(targetId, target, endPos, tracker, results);
            }
            
            // Finalize remaining alignments
            finalizeRemainingAlignments(tracker, results);
            
            return results;
        }
        
        /**
         * Processes any remaining alignments after all positions have been checked
         */
        private void finalizeRemainingAlignments(AlignmentTracker tracker, OffTargetResults results) {
            for (int targetId = 0; targetId < targets.size(); targetId++) {
                Alignment alignment = tracker.getAlignment(targetId);
                if (alignment != null) {
                    automataProcessAlignment(
                            targets.get(targetId), alignment, 
                            tracker.getPosition(targetId), tracker.getEdit(targetId), 
                            strand, results.alignmentList, results.endPosList, 
                            results.editNumList, results.targetList);
                }
            }
        }
        
        /**
         * Converts results to OffTargetData structure
         */
        private OffTargetData createOffTargetData(OffTargetResults results) {
            int[] endPosArr = results.endPosList.stream().mapToInt(i -> i).toArray();
            int[] editNumArr = results.editNumList.stream().mapToInt(i -> i).toArray();
            Alignment[] alignmentArr = results.alignmentList.toArray(new Alignment[0]);
            
            return new OffTargetData(endPosArr.length, endPosArr, editNumArr, alignmentArr);
        }
        
        /**
         * Calculates the genomic position based on strand orientation
         */
        private int calculateGenomicPosition(int endPos, Alignment alignment) {
            if (strand.equals("+")) {
                return endPos;
            } else {
                return text.length() - endPos + alignment.getAlignedText().replace("-", "").length();
            }
        }
        
        /**
         * Builds a CSV row for a single off-target site
         */
        private void appendCsvRow(StringBuilder builder, OffTargetData data, int index, String target) {
            Alignment alignment = data.getAlignmentArr()[index];
            int endPos = data.getEndPosArr()[index];
            int genomicPos = calculateGenomicPosition(endPos, alignment);
            
            // Chromosome and strand
            builder.append(chr).append(",").append(strand).append(",");
            
            // Genomic position
            builder.append(genomicPos).append(",");
            
            // Target
            builder.append(target).append(",");
            
            // Site sequence with flanking region
            int seqStart = Math.max(endPos - target.length() - MAX_EDITS, 0);
            builder.append(text.substring(seqStart, endPos)).append(",");
            
            // Edit distance
            builder.append(data.getEditNumArr()[index]).append(",");
            
            // Aligned sequences
            builder.append(alignment.getAlignedTarget()).append(",");
            builder.append(alignment.getAlignedText()).append(",");
            
            // Mismatches and bulges
            builder.append(alignment.getMismatches()).append(",");
            builder.append(alignment.getBulges()).append("\n");
        }
        
        /**
         * Converts OffTargetData to CSV format
         */
        private String buildCsvOutput(OffTargetData data, List<String> targetList) {
            StringBuilder builder = new StringBuilder();
            int size = data.getSize();
            
            for (int i = 0; i < size; i++) {
                appendCsvRow(builder, data, i, targetList.get(i));
            }
            
            return builder.toString();
        }
        
        /**
         * Writes CSV output to file in thread-safe manner
         */
        private void writeOutput(String csvOutput) {
            synchronized (writer) {
                try {
                    writer.append(csvOutput);
                    writer.flush();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

        @Override
        public void run() {
            // Process all positions and collect alignments
            OffTargetResults results = processAllPositions();
            
            // Convert to structured data
            OffTargetData offTargetData = createOffTargetData(results);
            
            // Build CSV output
            String csvOutput = buildCsvOutput(offTargetData, results.targetList);
            
            // Write to file
            writeOutput(csvOutput);
        }
    }

    public static class PositionsRes {
        public final List<Integer> ids;
        public final List<Integer> positions;

        public PositionsRes(List<Integer> ids, List<Integer> positions) {
            this.ids = ids;
            this.positions = positions;
        }

        public int size() {
            return positions.size();
        }
    }

    /**
     * Processes the results of an automata search from a file and writes the
     * results to a FileWriter.
     * The processing is done in parallel using multiple threads.
     * @param endPositions       List of end positions to process.
     * @param text               The genome chromosome sequence (text).
     * @param chrName            The genome chromosome name.
     * @param target             The target sequence.
     * @param pam                The PAM sequence.
     * @param writer             The FileWriter to write the results to.
     * @param allowNsInText      Boolean flag to allow 'N' characters in the text.
     * @param chooseBestInWindow Boolean flag to choose the best result in a window.
     * @param strand             The strand information.
     * @throws IOException If an I/O error occurs.
     */
    public static void automataResultsProcessFile(
            PositionsRes endPositions, String text, String chrName, List<String> targets, String pam,
            FileWriter writer, Boolean allowNsInText, Boolean chooseBestInWindow,
            String strand) throws IOException {
        ExecutorService executorService = Executors.newFixedThreadPool(NUM_THREADS);
        List<Future<?>> futures = new ArrayList<>();
        // devide endPositions to NUM_THREADS parts
        int partSize = endPositions.size() / NUM_THREADS; // turncate the decimal part
        for (int i = 0; i < NUM_THREADS; i++) {
            int start = i * partSize;
            int end = (i == NUM_THREADS - 1) ? endPositions.size() : (i + 1) * partSize;
            PositionsRes subEndPositions = new PositionsRes(
                endPositions.ids.subList(start, end), endPositions.positions.subList(start, end));
            AutomataResultsProcessFileHandler handler = new AutomataResultsProcessFileHandler(
                writer, text, targets, pam, strand, allowNsInText,
                chooseBestInWindow, subEndPositions, chrName);
            Future<?> future = executorService.submit(handler);
            futures.add(future);
        }
        try {
            for (Future<?> future : futures) {
                future.get(); // wait for all tasks to complete
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        executorService.shutdown();

    }

    /**
     * Parses a text file to extract integer positions from lines
     *
     * @param filePath the path to the text file to be parsed
     * @return a list of lists (for each target) of integers representing the extracted positions
     */
    public static PositionsRes parseAutomataResultsFile(String filePath, int targetNum) {
        // Each cell in the positions list will contain the positions for a specific target
        // Positions are sorted for each target 
        List<Integer> positions = new ArrayList<>();
        List<Integer> ids = new ArrayList<>();
        int [] lastPositions = new int[targetNum];
        for (int i = 0; i < targetNum; i++) {
            lastPositions[i] = -1;
        }

        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] parts = line.trim().split(":");
                if (parts.length != 2) {
                     continue; // skip malformed lines
                }
                int targetID; int position;
                try {
                    position = Integer.parseInt(parts[0]);
                    targetID = Integer.parseInt(parts[1]);
                } catch (NumberFormatException e) {
                    // If the line does not start with a number, skip it
                    continue;
                }
                if (position == lastPositions[targetID]) {
                    continue;
                }
                lastPositions[targetID] = position;
                // Add the position to the corresponding target's list
                positions.add(position);
                ids.add(targetID);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        return new PositionsRes(ids, positions);
    }

    /**
     * Configuration class to hold parsed command-line arguments
     */
    private static class Config {
        String fastaFilePath;
        String targetsFilePath;
        String pam;
        String autoOutputFolder;
        boolean chooseBestInWindow;
        boolean allowNsInText;
        
        Config(String[] args) {
            this.fastaFilePath = args[0];
            this.targetsFilePath = args[1];
            this.pam = args[10];
            this.autoOutputFolder = args[12];
            this.chooseBestInWindow = args[8].equals("true");
            this.allowNsInText = false;
            
            // Set global configuration
            AutoOffTargetSearchAlign.setOutputPath(args[2]);
            AutoOffTargetSearchAlign.setMaxEdits(Integer.parseInt(args[3]));
            AutoOffTargetSearchAlign.setMaxMismatchesWithoutBulges(Integer.parseInt(args[4]));
            AutoOffTargetSearchAlign.setMaxMismatchesWithBulges(Integer.parseInt(args[5]));
            AutoOffTargetSearchAlign.setMaxBulges(Integer.parseInt(args[6]));
            AutoOffTargetSearchAlign.setNumThreads(Integer.parseInt(args[7]));
            AutoOffTargetSearchAlign.setSiteWindowSize(Integer.parseInt(args[9]));
            AutoOffTargetSearchAlign.setAllowPamEdits(args[11].equals("true"));
        }
    }
    
    /**
     * Prepares the targets file path. If the input is a sequence string (no dot),
     * creates a temporary file with the sequence.
     */
    private static String prepareTargetsFile(String targetsFilePath) throws IOException {
        if (targetsFilePath.indexOf('.') == -1) {
            FileWriter writerTemp = new FileWriter("target_temp.txt");
            writerTemp.append(targetsFilePath + "\n");
            writerTemp.close();
            return "target_temp.txt";
        }
        return targetsFilePath;
    }
    
    /**
     * Reads target sequences from a file
     */
    private static List<String> readTargets(String targetsFilePath) throws IOException {
        List<String> targets = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(targetsFilePath))) {
            String target;
            while ((target = reader.readLine()) != null) {
                targets.add(target);
            }
        }
        return targets;
    }
    
    /**
     * Prepares chromosome files for both forward and reverse complement strands
     */
    private static Map<String, File[]> prepareChromosomeFiles(String fastaFilePath) throws IOException {
        FastaReader.splitFastaToFiles(fastaFilePath);
        
        Map<String, File[]> strandToFiles = new HashMap<>();
        strandToFiles.put("+", FastaReader.getChrFiles(fastaFilePath, false));
        strandToFiles.put("-", FastaReader.getChrFiles(fastaFilePath, true));
        
        return strandToFiles;
    }
    
    /**
     * Constructs the automata output file path based on chromosome name and strand
     */
    private static String getAutomataOutputPath(File chrFile, String strand, String autoOutputFolder) {
        String name = chrFile.getName();
        int dotIndex = name.lastIndexOf('.');
        String baseName = (dotIndex == -1) ? name : name.substring(0, dotIndex);
        String extension = (dotIndex == -1) ? "" : name.substring(dotIndex);
        String suffix = strand.equals("-") ? "_rc" : "_fw";
        
        return java.nio.file.Paths.get(autoOutputFolder, baseName + suffix + extension).toString();
    }
    
    /**
     * Extracts chromosome name from file name
     */
    private static String extractChromosomeName(File file) {
        String chrName = file.getName();
        int dotIndex = chrName.lastIndexOf('.');
        return (dotIndex == -1) ? chrName : chrName.substring(0, dotIndex);
    }
    
    /**
     * Processes a single chromosome on a specific strand
     */
    private static void processChromosome(
            File chrFile, String strand, List<String> targets, Config config, FileWriter writer) 
            throws IOException {
        
        String text = Files.readString(Path.of(chrFile.getAbsolutePath()));
        String chrName = extractChromosomeName(chrFile);
        String autoFileOutput = getAutomataOutputPath(chrFile, strand, config.autoOutputFolder);
        
        if (new File(autoFileOutput).exists()) {
            PositionsRes endPositions = parseAutomataResultsFile(autoFileOutput, targets.size());
            System.out.println("Processing " + endPositions.size() + 
                             " end positions for " + chrName + " on strand " + strand);
            
            automataResultsProcessFile(
                endPositions, text, chrName, targets, config.pam, 
                writer, config.allowNsInText, config.chooseBestInWindow, strand);
        }
    }
    
    /**
     * Processes all chromosomes and strands
     */
    private static void processAllChromosomes(
            Map<String, File[]> strandToFiles, List<String> targets, Config config, FileWriter writer) 
            throws IOException {
        
        String[] strands = new String[]{"+", "-"};
        File[] chromosomeFiles = strandToFiles.get("+");
        
        for (int fileIndex = 0; fileIndex < chromosomeFiles.length; fileIndex++) {
            for (String strand : strands) {
                File chrFile = strandToFiles.get(strand)[fileIndex];
                processChromosome(chrFile, strand, targets, config, writer);
            }
        }
    }
    
    public static void main(String[] args) {
        Instant start = Instant.now();
        
        try {
            // Parse configuration
            Config config = new Config(args);
            
            // Prepare targets file
            String targetsFilePath = prepareTargetsFile(config.targetsFilePath);
            
            // Read targets
            List<String> targets = readTargets(targetsFilePath);
            
            // Prepare chromosome files
            Map<String, File[]> strandToFiles = prepareChromosomeFiles(config.fastaFilePath);
            
            // Process all chromosomes
            FileWriter writer = writeOutputFile();
            processAllChromosomes(strandToFiles, targets, config, writer);
            writer.close();
            
            // Cleanup
            new File("target_temp.txt").delete();
            System.gc();
            
        } catch (IOException e) {
            e.printStackTrace();
        }
        
        // Report execution time
        Instant end = Instant.now();
        Duration timeElapsed = Duration.between(start, end);
        System.out.println("total Execution time: " + timeElapsed.toMillis() + " milliseconds");
    }
}
