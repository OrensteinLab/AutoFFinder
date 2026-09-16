package PostAutoFFinder;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.util.AbstractList;
import java.util.Arrays;
import java.nio.file.Path;
import java.nio.file.Files;
import java.io.File;
import java.io.IOException;
import java.io.FileWriter;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.UncheckedIOException;

import java.time.Instant;
import java.time.Duration;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ExecutionException;


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
    // when true, skips both the 0-bulge mismatch-only scan and the 1-bulge fast scan,
    // always falling through to the DP + memoized traceback regardless of MAX_BULGES
    private static boolean DISABLE_FAST_PATH = false;
    private static final int TASKS_PER_THREAD = 4;
    private static final ThreadLocal<int[][]> sharedDpMatrix =
            ThreadLocal.withInitial(() -> new int[60][60]);
    private static final ThreadLocal<ReconstructionWorkspace> reconstructionWorkspace =
            ThreadLocal.withInitial(ReconstructionWorkspace::new);
    

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

    public static void setDisableFastPath(boolean value) {
        DISABLE_FAST_PATH = value;
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
        if (!DISABLE_FAST_PATH && MAX_BULGES <= 1) {
            SingleBulgeResult result = findSingleBulgeAlignment(
                    target, text, allowNsInText, targetPamSuffix, textPamSuffix);
            if (!result.ambiguous) {
                return result.alignment;
            }
        }

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

    private static SingleBulgeResult findSingleBulgeAlignment(
            String target, String text, boolean allowNsInText, String targetPamSuffix, String textPamSuffix) {
        int targetLen = target.length();
        int textLen = text.length();
        SingleBulgeResult result = new SingleBulgeResult();

        if (textLen >= targetLen) {
            int textStart = textLen - targetLen;
            int mismatches = countMismatches(target, text, textStart, 0, targetLen, allowNsInText);
            if (mismatches <= MAX_EDITS && mismatches <= MAX_MISMATCHES_WITHOUT_BULGES) {
                result.alignment = createAlignment(target, text.substring(textStart), mismatches, 0);
            }
        }

        if (MAX_BULGES == 0) {
            appendPamSuffix(result.alignment, targetPamSuffix, textPamSuffix);
            return result;
        }

        int minimumBulgedEdits = Math.min(
                findTargetBulgeMinimumEdits(target, text, allowNsInText),
                findTextBulgeMinimumEdits(target, text, allowNsInText, !targetPamSuffix.isEmpty()));
        if (result.alignment == null) {
            if (minimumBulgedEdits != Integer.MAX_VALUE) {
                result.ambiguous = true;
            }
            return result;
        }
        if (minimumBulgedEdits < editCount(result.alignment)) {
            result.ambiguous = true;
            return result;
        }

        appendPamSuffix(result.alignment, targetPamSuffix, textPamSuffix);
        return result;
    }

    private static void appendPamSuffix(
            Alignment alignment, String targetPamSuffix, String textPamSuffix) {
        if (alignment == null) {
            return;
        }
        alignment.getAlignedTargetBuilder().append(targetPamSuffix);
        alignment.getAlignedTextBuilder().append(textPamSuffix);
    }

    private static int editCount(Alignment alignment) {
        return alignment.getMismatches() + alignment.getBulges();
    }

    private static int findTargetBulgeMinimumEdits(String target, String text, boolean allowNsInText) {
        int targetLen = target.length();
        if (text.length() < targetLen - 1) {
            return Integer.MAX_VALUE;
        }
        int textStart = text.length() - targetLen + 1;
        int[] prefix = new int[targetLen + 1];
        int[] suffix = new int[targetLen + 1];
        for (int index = 0; index < targetLen - 1; index++) {
            prefix[index + 1] = prefix[index]
                    + mismatchScore(target.charAt(index), text.charAt(textStart + index), allowNsInText);
        }
        for (int index = targetLen - 1; index > 0; index--) {
            suffix[index] = suffix[index + 1]
                    + mismatchScore(target.charAt(index), text.charAt(textStart + index - 1), allowNsInText);
        }

        int minimumEdits = Integer.MAX_VALUE;
        for (int gap = 0; gap < targetLen; gap++) {
            int mismatches = prefix[gap] + suffix[gap + 1];
            if (validBulgeMismatchCount(mismatches)) {
                minimumEdits = Math.min(minimumEdits, mismatches + 1);
            }
        }
        return minimumEdits;
    }

    private static int findTextBulgeMinimumEdits(
            String target, String text, boolean allowNsInText, boolean allowRightmostGap) {
        int targetLen = target.length();
        if (text.length() < targetLen + 1) {
            return Integer.MAX_VALUE;
        }
        int textStart = text.length() - targetLen - 1;
        int[] prefix = new int[targetLen + 1];
        int[] suffix = new int[targetLen + 1];
        for (int index = 0; index < targetLen; index++) {
            prefix[index + 1] = prefix[index]
                    + mismatchScore(target.charAt(index), text.charAt(textStart + index), allowNsInText);
        }
        for (int index = targetLen - 1; index >= 0; index--) {
            suffix[index] = suffix[index + 1]
                    + mismatchScore(target.charAt(index), text.charAt(textStart + index + 1), allowNsInText);
        }

        int minimumEdits = Integer.MAX_VALUE;
        int lastGap = allowRightmostGap ? targetLen : targetLen - 1;
        for (int gap = 1; gap <= lastGap; gap++) {
            int mismatches = prefix[gap] + suffix[gap];
            if (validBulgeMismatchCount(mismatches)) {
                minimumEdits = Math.min(minimumEdits, mismatches + 1);
            }
        }
        return minimumEdits;
    }

    private static int countMismatches(
            String target, String text, int textStart, int targetStart, int length, boolean allowNsInText) {
        int mismatches = 0;
        for (int index = 0; index < length; index++) {
            mismatches += mismatchScore(
                    target.charAt(targetStart + index), text.charAt(textStart + index), allowNsInText);
        }
        return mismatches;
    }

    private static int mismatchScore(char target, char text, boolean allowNsInText) {
        return charactersMatch(target, text, allowNsInText) ? 0 : 1;
    }

    private static boolean validBulgeMismatchCount(int mismatches) {
        return mismatches + 1 <= EFFECTIVE_MAX_EDIT_WITH_BULGE
                && mismatches <= MAX_MISMATCHES_WITH_BULGES
                && mismatches <= MAX_MISMATCHES_WITHOUT_BULGES;
    }

    private static Alignment createAlignment(String target, String text, int mismatches, int bulges) {
        return new Alignment(mismatches, bulges, new StringBuilder(target), new StringBuilder(text));
    }

    private static final class SingleBulgeResult {
        private Alignment alignment;
        private boolean ambiguous;
    }


    public static Alignment naive_find_alignment(
            int[][] M, Boolean allowNsInText, String target, String text, int target_i, int text_i,
            StringBuilder[] targetTextAlign, int mismatches, int bulges, String targetPamSuffix, String textPamSuffix) {
        ReconstructionWorkspace workspace = reconstructionWorkspace.get();
        workspace.prepare(target.length() + 1, text.length() + 1,
            Math.max(MAX_MISMATCHES_WITHOUT_BULGES, MAX_MISMATCHES_WITH_BULGES) + 1, MAX_BULGES + 1);
        boolean hasPamSuffix = !targetPamSuffix.isEmpty();
        int result = findAlignmentMemoized(
            M, allowNsInText, hasPamSuffix, target, text,
            target_i, text_i, mismatches, bulges, workspace);
        if (result == -1) {
            return null;
        }

        StringBuilder finalTargetAlign = new StringBuilder(target.length() + MAX_BULGES);
        StringBuilder finalTextAlign = new StringBuilder(target.length() + MAX_BULGES);
        reconstructMemoizedAlignment(
            workspace, allowNsInText, target, text, target_i, text_i,
            mismatches, bulges, finalTargetAlign, finalTextAlign);
        finalTargetAlign.append(new StringBuilder(targetTextAlign[0]).reverse());
        finalTextAlign.append(new StringBuilder(targetTextAlign[1]).reverse());
        if (!targetPamSuffix.isEmpty()) {
            finalTargetAlign.append(targetPamSuffix);
            finalTextAlign.append(textPamSuffix);
        }
        return new Alignment(resultMismatches(result), resultBulges(result), finalTargetAlign, finalTextAlign);
    }

    private static int findAlignmentMemoized(
            int[][] M, boolean allowNsInText, boolean hasPamSuffix, String target, String text, int target_i, int text_i,
            int mismatches, int bulges, ReconstructionWorkspace workspace) {
        int targetLen = target.length();

        if (bulges > MAX_BULGES) {
            return -1;
        }

        if (bulges > 0 && mismatches > MAX_MISMATCHES_WITH_BULGES) {
            return -1;
        }
        if (mismatches > MAX_MISMATCHES_WITHOUT_BULGES) {
            return -1;
        }
        if (target_i == 0) {
            return packResult(mismatches, bulges);
        }

        if (text_i <= 0) {
            return -1;
        }
        int effective_max_edit = bulges == 0?  MAX_EDITS : EFFECTIVE_MAX_EDIT_WITH_BULGE;
        if (M[target_i][text_i] > (effective_max_edit - mismatches - bulges)) {
            return -1;
        }

        int stateKey = workspace.key(target_i, text_i, mismatches, bulges);
        if (workspace.isComputed(stateKey)) {
            return workspace.get(stateKey);
        }

        char targetChar = target.charAt(target_i - 1);
        char textChar = text.charAt(text_i - 1);
        int mismatchScore = charactersMatch(targetChar, textChar, allowNsInText) ? 0 : 1;

        int diagonal = findAlignmentMemoized(
            M, allowNsInText, hasPamSuffix, target, text, target_i - 1, text_i - 1,
                mismatches + mismatchScore, bulges, workspace);
        int targetBulge = findAlignmentMemoized(
            M, allowNsInText, hasPamSuffix, target, text, target_i - 1, text_i,
                mismatches, bulges + 1, workspace);
        int textBulge = -1;
        if (target_i != targetLen || hasPamSuffix) {
            textBulge = findAlignmentMemoized(
                M, allowNsInText, hasPamSuffix, target, text, target_i, text_i - 1,
                    mismatches, bulges + 1, workspace);
        }

        int best = diagonal;
        byte move = diagonal == -1 ? ReconstructionWorkspace.NO_MOVE : ReconstructionWorkspace.DIAGONAL;
        if (isBetterResult(targetBulge, best)) {
            best = targetBulge;
            move = ReconstructionWorkspace.TARGET_BULGE;
        }
        if (isBetterResult(textBulge, best)) {
            best = textBulge;
            move = ReconstructionWorkspace.TEXT_BULGE;
        }
        workspace.put(stateKey, best, move);
        return best;
    }

    private static boolean isBetterResult(int candidate, int current) {
        return candidate != -1 && (current == -1 || resultEdits(candidate) < resultEdits(current));
    }

    private static int packResult(int mismatches, int bulges) {
        return (mismatches << 16) | bulges;
    }

    private static int resultMismatches(int result) {
        return result >>> 16;
    }

    private static int resultBulges(int result) {
        return result & 0xffff;
    }

    private static int resultEdits(int result) {
        return resultMismatches(result) + resultBulges(result);
    }

    private static void reconstructMemoizedAlignment(
            ReconstructionWorkspace workspace, boolean allowNsInText,
            String target, String text, int target_i, int text_i, int mismatches, int bulges,
            StringBuilder alignedTarget, StringBuilder alignedText) {
        while (target_i > 0) {
            int stateKey = workspace.key(target_i, text_i, mismatches, bulges);
            byte move = workspace.getMove(stateKey);
            if (move == ReconstructionWorkspace.DIAGONAL) {
                char targetChar = target.charAt(target_i - 1);
                char textChar = text.charAt(text_i - 1);
                alignedTarget.append(targetChar);
                alignedText.append(textChar);
                mismatches += charactersMatch(targetChar, textChar, allowNsInText) ? 0 : 1;
                target_i--;
                text_i--;
            } else if (move == ReconstructionWorkspace.TARGET_BULGE) {
                alignedTarget.append(target.charAt(target_i - 1));
                alignedText.append('-');
                bulges++;
                target_i--;
            } else if (move == ReconstructionWorkspace.TEXT_BULGE) {
                alignedTarget.append('-');
                alignedText.append(text.charAt(text_i - 1));
                bulges++;
                text_i--;
            } else {
                throw new IllegalStateException("Missing traceback move for successful alignment");
            }
        }
        alignedTarget.reverse();
        alignedText.reverse();
    }

    private static final class ReconstructionWorkspace {
        private static final byte NO_MOVE = 0;
        private static final byte DIAGONAL = 1;
        private static final byte TARGET_BULGE = 2;
        private static final byte TEXT_BULGE = 3;

        private int[] values = new int[0];
        private byte[] moves = new byte[0];
        private int[] generations = new int[0];
        private int generation;
        private int textCapacity;
        private int mismatchCapacity;
        private int bulgeCapacity;

        private void prepare(int targetCapacity, int textCapacity, int mismatchCapacity, int bulgeCapacity) {
            int requiredSize = targetCapacity * textCapacity * mismatchCapacity * bulgeCapacity;
            if (values.length < requiredSize
                    || this.textCapacity != textCapacity
                    || this.mismatchCapacity != mismatchCapacity
                    || this.bulgeCapacity != bulgeCapacity) {
                values = new int[requiredSize];
                moves = new byte[requiredSize];
                generations = new int[requiredSize];
                generation = 0;
            }
            this.textCapacity = textCapacity;
            this.mismatchCapacity = mismatchCapacity;
            this.bulgeCapacity = bulgeCapacity;
            generation++;
        }

        private int key(int targetIndex, int textIndex, int mismatches, int bulges) {
            return (((targetIndex * textCapacity) + textIndex) * mismatchCapacity + mismatches)
                    * bulgeCapacity + bulges;
        }

        private boolean isComputed(int key) {
            return generations[key] == generation;
        }

        private int get(int key) {
            return values[key];
        }

        private byte getMove(int key) {
            return moves[key];
        }

        private void put(int key, int value, byte move) {
            values[key] = value;
            moves[key] = move;
            generations[key] = generation;
        }
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
            Alignment targetTextAlignment, int targetTextAlignmentPos, int targetTextAlignmentEdit) {
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
        private PositionsRes endPositions;
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
            this.endPositions = endPositions;
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
            
            Alignment newAlignment = automataPostprocessing(
                    endPos, strand, target, pam, text, allowNsInText, chooseBestInWindow,
                    results.alignmentList, results.endPosList, results.editNumList, results.targetList,
                    currentAlignment, currentPos, currentEdit);
            
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
                int targetId = endPositions.idAt(i);
                String target = targets.get(targetId);
                int endPos = endPositions.positionAt(i);
                
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
                    throw new UncheckedIOException("Failed to write post-processing output", e);
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
        private final int[] idValues;
        private final int[] positionValues;
        private final int offset;
        private final int size;

        public PositionsRes(int[] idValues, int[] positionValues) {
            this(idValues, positionValues, 0, positionValues.length);
        }

        public PositionsRes(List<Integer> ids, List<Integer> positions) {
            this(
                    ids.stream().mapToInt(Integer::intValue).toArray(),
                    positions.stream().mapToInt(Integer::intValue).toArray());
        }

        private PositionsRes(
                int[] idValues, int[] positionValues, int offset, int size) {
            if (idValues.length != positionValues.length) {
                throw new IllegalArgumentException("IDs and positions must have the same length");
            }
            this.idValues = idValues;
            this.positionValues = positionValues;
            this.offset = offset;
            this.size = size;
            this.ids = arrayView(idValues, offset, size);
            this.positions = arrayView(positionValues, offset, size);
        }

        private static List<Integer> arrayView(int[] values, int offset, int size) {
            return new AbstractList<>() {
                @Override
                public Integer get(int index) {
                    if (index < 0 || index >= size) {
                        throw new IndexOutOfBoundsException(index);
                    }
                    return values[offset + index];
                }

                @Override
                public Integer set(int index, Integer value) {
                    if (index < 0 || index >= size) {
                        throw new IndexOutOfBoundsException(index);
                    }
                    int valueIndex = offset + index;
                    int previous = values[valueIndex];
                    values[valueIndex] = value;
                    return previous;
                }

                @Override
                public int size() {
                    return size;
                }
            };
        }

        public int size() {
            return size;
        }

        public int idAt(int index) {
            return idValues[offset + index];
        }

        public int positionAt(int index) {
            return positionValues[offset + index];
        }

        private PositionsRes slice(int start, int end) {
            return new PositionsRes(
                    idValues, positionValues, offset + start, end - start);
        }
    }

    private static final class IntArrayBuilder {
        private int[] values = new int[1_024];
        private int size;

        private void add(int value) {
            if (size == values.length) {
                values = Arrays.copyOf(values, values.length * 2);
            }
            values[size++] = value;
        }

        private int[] toArray() {
            return Arrays.copyOf(values, size);
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
        if (endPositions.size() == 0) {
            return;
        }

        ExecutorService executorService = Executors.newFixedThreadPool(NUM_THREADS);
        List<Future<?>> futures = new ArrayList<>();
        int taskCount = Boolean.TRUE.equals(chooseBestInWindow)
                ? 1
                : (int) Math.min(
                        endPositions.size(), (long) NUM_THREADS * TASKS_PER_THREAD);
        int partSize = endPositions.size() / taskCount;
        for (int i = 0; i < taskCount; i++) {
            int start = i * partSize;
            int end = (i == taskCount - 1) ? endPositions.size() : (i + 1) * partSize;
            PositionsRes subEndPositions = endPositions.slice(start, end);
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
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new IOException("Interrupted while processing candidate alignments", e);
        } catch (ExecutionException e) {
            throw new IOException("Candidate alignment worker failed", e.getCause());
        } finally {
            executorService.shutdown();
        }

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
        IntArrayBuilder positions = new IntArrayBuilder();
        IntArrayBuilder ids = new IntArrayBuilder();
        int [] lastPositions = new int[targetNum];
        for (int i = 0; i < targetNum; i++) {
            lastPositions[i] = -1;
        }

        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = br.readLine()) != null) {
                int lineStart = 0;
                int lineEnd = line.length();
                while (lineStart < lineEnd && line.charAt(lineStart) <= ' ') {
                    lineStart++;
                }
                while (lineEnd > lineStart && line.charAt(lineEnd - 1) <= ' ') {
                    lineEnd--;
                }
                int separator = line.indexOf(':', lineStart);
                int extraSeparator = line.indexOf(':', separator + 1);
                if (separator <= lineStart || separator >= lineEnd - 1
                    || (extraSeparator != -1 && extraSeparator < lineEnd)) {
                    continue;
                }
                int position = parseNonNegativeInt(line, lineStart, separator);
                int targetID = parseNonNegativeInt(line, separator + 1, lineEnd);
                if (position < 0 || targetID < 0 || targetID >= targetNum) {
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
            throw new UncheckedIOException("Failed to read candidate file " + filePath, e);
        }

        return new PositionsRes(ids.toArray(), positions.toArray());
    }

    private static int parseNonNegativeInt(String value, int start, int end) {
        int result = 0;
        for (int index = start; index < end; index++) {
            char character = value.charAt(index);
            if (character < '0' || character > '9') {
                return -1;
            }
            int digit = character - '0';
            if (result > (Integer.MAX_VALUE - digit) / 10) {
                return -1;
            }
            result = result * 10 + digit;
        }
        return result;
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
            // optional trailing arg (defaults to false) to disable the single-bulge fast path
            AutoOffTargetSearchAlign.setDisableFastPath(args.length > 13 && args[13].equals("true"));
        }
    }
    
    /**
     * Prepares the targets file path. If the input is not an existing file,
     * treats it as a literal DNA sequence and creates a temporary guide file.
     */
    private static String prepareTargetsFile(String targetsFilePath) throws IOException {
        if (Files.isRegularFile(Path.of(targetsFilePath))) {
            return targetsFilePath;
        }
        if (targetsFilePath.matches("[ACGTNacgtn]+")) {
            FileWriter writerTemp = new FileWriter("target_temp.txt");
            writerTemp.append(targetsFilePath.toUpperCase()).append("\n");
            writerTemp.close();
            return "target_temp.txt";
        }
        throw new IOException(
                "Guide file does not exist and the argument is not a DNA sequence: "
                        + targetsFilePath);
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
    private static long processChromosome(
            File chrFile, String strand, List<String> targets, Config config,
            CandidateSource candidateSource, FileWriter writer)
            throws IOException {
        
        String text = Files.readString(Path.of(chrFile.getAbsolutePath()));
        String chrName = extractChromosomeName(chrFile);
        PositionsRes endPositions = candidateSource.load(chrFile, strand, targets.size());
        
        if (endPositions != null) {
            System.out.println("Processing " + endPositions.size() + 
                             " end positions for " + chrName + " on strand " + strand);
            
            long postprocessingStart = System.nanoTime();
            automataResultsProcessFile(
                endPositions, text, chrName, targets, config.pam, 
                writer, config.allowNsInText, config.chooseBestInWindow, strand);
            return System.nanoTime() - postprocessingStart;
        }
        return 0;
    }
    
    /**
     * Processes all chromosomes and strands
     */
    private static long processAllChromosomes(
            Map<String, File[]> strandToFiles, List<String> targets, Config config,
            CandidateSource candidateSource, FileWriter writer)
            throws IOException {
        
        String[] strands = new String[]{"+", "-"};
        File[] chromosomeFiles = strandToFiles.get("+");
        long postprocessingNanos = 0;
        
        for (int fileIndex = 0; fileIndex < chromosomeFiles.length; fileIndex++) {
            for (String strand : strands) {
                File chrFile = strandToFiles.get(strand)[fileIndex];
                postprocessingNanos += processChromosome(
                        chrFile, strand, targets, config, candidateSource, writer);
            }
        }
        return postprocessingNanos;
    }

    /**
     * Creates the candidate source used to obtain end positions for each chromosome/strand.
     * Defaults to reading previously generated automata text output files (the original
     * behavior). Set the "autoffinder.candidateSource" system property to "binary" to read
     * raw FPGA output files captured to disk, or "fpga" to run a real FPGA device live via a
     * native JNI library.
     */
    private static CandidateSource createCandidateSource(Config config, String targetsFilePath) {
        String mode = System.getProperty("autoffinder.candidateSource", "file");
        if (mode.equals("file")) {
            return new TextFileCandidateSource(config.autoOutputFolder);
        }
        if (mode.equals("fpga")) {
            String xclbinPath = System.getProperty("relev.xclbin");
            if (xclbinPath == null || xclbinPath.isEmpty()) {
                throw new IllegalArgumentException(
                        "FPGA mode requires -Drelev.xclbin=<path-to-xclbin>");
            }
            return new FpgaCandidateSource(
                    System.getProperty("relev.nativeLibrary"),
                    xclbinPath,
                    targetsFilePath,
                    MAX_EDITS);
        }
        if (mode.equals("binary")) {
            Path root = Path.of(config.autoOutputFolder);
            Path forwardDirectory = Path.of(System.getProperty(
                    "relev.binary.forwardDir",
                    root.resolve("hg38_only_chrs_split").toString()));
            Path reverseDirectory = Path.of(System.getProperty(
                    "relev.binary.reverseDir",
                    root.resolve("hg38_only_chrs_split_rc").toString()));
            int editDistance = Integer.parseInt(System.getProperty(
                    "relev.binary.editDistance", Integer.toString(MAX_EDITS)));
            return new RawBinaryCandidateSource(
                    forwardDirectory, reverseDirectory, editDistance);
        }
        throw new IllegalArgumentException(
                "Unknown candidate source '" + mode + "'; expected file, binary, or fpga");
    }

    /**
     * Prints raw buffer preparation/decoding timings when the candidate source reports them
     * (i.e. when reading raw FPGA output, either captured to disk or run live).
     */
    private static void printCandidateTimings(
            CandidateSource candidateSource, long postprocessingNanos) {
        long preparationNanos = candidateSource.preparationNanos();
        long decodingNanos = candidateSource.decodingNanos();
        if (preparationNanos == 0 && decodingNanos == 0) {
            return;
        }
        System.out.printf("Raw buffer preparation time: %.3f ms%n", preparationNanos / 1_000_000.0);
        System.out.printf("Java binary decoding time: %.3f ms%n", decodingNanos / 1_000_000.0);
        System.out.printf("Java candidate post-processing time: %.3f ms%n", postprocessingNanos / 1_000_000.0);
        System.out.printf(
                "Java stage 2 time (decode + post-processing): %.3f ms%n",
                (decodingNanos + postprocessingNanos) / 1_000_000.0);
    }
    
    public static void main(String[] args) {
        Instant start = Instant.now();

        if (args.length < 13 || args.length > 14) {
            System.err.println(
                    "Usage: AutoOffTargetSearchAlign <genome-fasta> <guide-file-or-sequence> "
                            + "<output-prefix> <maxE> <maxM> <maxMB> <maxB> <threads> "
                            + "<best-in-window> <best-window-size> <PAM> <allow-PAM-edits> "
                            + "<candidate-source-path> [disable-fast-path]");
            System.exit(2);
        }
        
        try {
            // Parse configuration
            Config config = new Config(args);
            
            // Prepare targets file
            String targetsFilePath = prepareTargetsFile(config.targetsFilePath);
            
            // Read targets
            List<String> targets = readTargets(targetsFilePath);
            
            // Prepare chromosome files
            Map<String, File[]> strandToFiles = prepareChromosomeFiles(config.fastaFilePath);
            
            // Prepare candidate source (defaults to reading automata text output files;
            // can be switched to raw FPGA binary files or a live FPGA device)
            CandidateSource candidateSource = createCandidateSource(config, targetsFilePath);

            // Process all chromosomes
            FileWriter writer = writeOutputFile();
            long postprocessingNanos = processAllChromosomes(
                    strandToFiles, targets, config, candidateSource, writer);
            writer.close();
            printCandidateTimings(candidateSource, postprocessingNanos);
            
            // Cleanup
            new File("target_temp.txt").delete();
            System.gc();
            
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }
        
        // Report execution time
        Instant end = Instant.now();
        Duration timeElapsed = Duration.between(start, end);
        System.out.println("total Execution time: " + timeElapsed.toMillis() + " milliseconds");
    }
}
