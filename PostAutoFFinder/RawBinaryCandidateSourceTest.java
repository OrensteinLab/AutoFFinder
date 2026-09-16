package PostAutoFFinder;

import java.io.File;
import java.nio.file.Path;

public final class RawBinaryCandidateSourceTest {
    private RawBinaryCandidateSourceTest() {
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 4) {
            System.err.println(
                    "Usage: RawBinaryCandidateSourceTest <raw-root> <chromosome-file-name> <+|-> <edit-distance>");
            System.exit(2);
        }
        if (!args[2].equals("+") && !args[2].equals("-")) {
            throw new IllegalArgumentException("Strand must be + or -");
        }

        Path root = Path.of(args[0]);
        RawBinaryCandidateSource source = new RawBinaryCandidateSource(
                root.resolve("hg38_only_chrs_split"),
                root.resolve("hg38_only_chrs_split_rc"),
                Integer.parseInt(args[3]));
        AutoOffTargetSearchAlign.PositionsRes candidates = source.load(
                new File(args[1]), args[2], 128);

        System.out.printf("Raw buffer preparation time: %.3f ms%n", source.preparationNanos() / 1_000_000.0);
        System.out.printf("Java binary decoding time: %.3f ms%n", source.decodingNanos() / 1_000_000.0);
        System.out.printf("Decoded candidates: %,d%n", candidates.size());
    }
}