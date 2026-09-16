package PostAutoFFinder;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Stream;

public final class RawBinaryBatchBenchmark {
    private RawBinaryBatchBenchmark() {
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 2) {
            System.err.println("Usage: RawBinaryBatchBenchmark <raw-root> <edit-distance>");
            System.exit(2);
        }

        Path root = Path.of(args[0]);
        int editDistance = Integer.parseInt(args[1]);
        Path forwardDirectory = root.resolve("hg38_only_chrs_split");
        Path reverseDirectory = root.resolve("hg38_only_chrs_split_rc");
        RawBinaryCandidateSource source = new RawBinaryCandidateSource(
                forwardDirectory, reverseDirectory, editDistance);

        List<String> chromosomeNames = chromosomeNames(forwardDirectory, editDistance);
        if (chromosomeNames.isEmpty()) {
            throw new IllegalArgumentException(
                    "No edit-distance-" + editDistance + " captures found in " + forwardDirectory);
        }

        long totalCandidates = 0;
        int fileCount = 0;
        for (String chromosomeName : chromosomeNames) {
            for (String strand : new String[] {"+", "-"}) {
                long preparationBefore = source.preparationNanos();
                long decodingBefore = source.decodingNanos();
                AutoOffTargetSearchAlign.PositionsRes candidates =
                        source.load(new File(chromosomeName), strand, 128);
                long preparation = source.preparationNanos() - preparationBefore;
                long decoding = source.decodingNanos() - decodingBefore;

                totalCandidates += candidates.size();
                fileCount++;
                System.out.printf(
                        "%s strand %s: %,d candidates, preparation %.3f ms, decoding %.3f ms%n",
                        chromosomeName, strand, candidates.size(),
                        preparation / 1_000_000.0, decoding / 1_000_000.0);
            }
        }

        System.out.println();
        System.out.printf("Files decoded: %,d%n", fileCount);
        System.out.printf("Candidates decoded: %,d%n", totalCandidates);
        System.out.printf(
                "Total raw buffer preparation time: %.3f ms%n",
                source.preparationNanos() / 1_000_000.0);
        System.out.printf(
                "Total Java binary decoding time: %.3f ms%n",
                source.decodingNanos() / 1_000_000.0);
        System.out.printf(
                "Average preparation time per file: %.3f ms%n",
                source.preparationNanos() / 1_000_000.0 / fileCount);
        System.out.printf(
                "Average decoding time per file: %.3f ms%n",
                source.decodingNanos() / 1_000_000.0 / fileCount);
    }

    private static List<String> chromosomeNames(Path directory, int editDistance) throws Exception {
        String suffix = "_ed" + editDistance + ".out";
        try (Stream<Path> files = Files.list(directory)) {
            return files
                    .map(path -> path.getFileName().toString())
                    .filter(name -> name.endsWith(suffix))
                    .map(name -> name.substring(0, name.length() - suffix.length()))
                    .sorted(Comparator.comparingInt(RawBinaryBatchBenchmark::chromosomeOrder))
                    .toList();
        }
    }

    private static int chromosomeOrder(String chromosomeName) {
        String name = chromosomeName;
        if (name.startsWith("chr")) {
            name = name.substring(3);
        }
        int dotIndex = name.indexOf('.');
        if (dotIndex >= 0) {
            name = name.substring(0, dotIndex);
        }
        if (name.equals("X")) {
            return 23;
        }
        if (name.equals("Y")) {
            return 24;
        }
        if (name.equals("M")) {
            return 25;
        }
        try {
            return Integer.parseInt(name);
        } catch (NumberFormatException exception) {
            return Integer.MAX_VALUE;
        }
    }
}
