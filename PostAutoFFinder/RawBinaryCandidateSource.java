package PostAutoFFinder;

import java.io.File;
import java.io.IOException;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;

public final class RawBinaryCandidateSource implements CandidateSource {
    private static final int MATCH_RECORD_BYTES = 8;
    // The FPGA image always operates on 128 lanes; captures for fewer guides are derived
    // offline by bit-masking the 128-lane records down to the first N target lanes, so any
    // target count from 1 up to the full lane width is a valid capture to decode.
    private static final int MAX_FPGA_TARGET_COUNT = 128;

    private final Path forwardDirectory;
    private final Path reverseDirectory;
    private final int editDistance;
    private long preparationNanos;
    private long decodingNanos;

    public RawBinaryCandidateSource(
            Path forwardDirectory, Path reverseDirectory, int editDistance) {
        if (editDistance < 0 || editDistance > 6) {
            throw new IllegalArgumentException("FPGA edit distance must be between 0 and 6");
        }
        this.forwardDirectory = forwardDirectory;
        this.reverseDirectory = reverseDirectory;
        this.editDistance = editDistance;
    }

    @Override
    public AutoOffTargetSearchAlign.PositionsRes load(
            File chromosomeFile, String strand, int targetCount) throws IOException {
        if (targetCount <= 0 || targetCount > MAX_FPGA_TARGET_COUNT) {
            throw new IllegalArgumentException(
                    "Raw FPGA captures support at most " + MAX_FPGA_TARGET_COUNT
                            + " targets, but received " + targetCount);
        }

        Path directory = strand.equals("-") ? reverseDirectory : forwardDirectory;
        Path capture = directory.resolve(
                chromosomeFile.getName() + "_ed" + editDistance + ".out");
        if (!capture.toFile().isFile()) {
            throw new IOException("Raw FPGA output does not exist: " + capture);
        }

        long preparationStart = System.nanoTime();
        MappedByteBuffer mapped;
        try (FileChannel channel = FileChannel.open(capture, StandardOpenOption.READ)) {
            long size = channel.size();
            if (size == 0 || size > Integer.MAX_VALUE || size % MATCH_RECORD_BYTES != 0) {
                throw new IOException("Invalid raw FPGA output size " + size + ": " + capture);
            }
            mapped = channel.map(FileChannel.MapMode.READ_ONLY, 0, size);
        }
        mapped.load();
        mapped.position(0);
        preparationNanos += System.nanoTime() - preparationStart;

        long decodingStart = System.nanoTime();
        AutoOffTargetSearchAlign.PositionsRes result =
                RelevBinaryDecoder.decode(mapped, targetCount);
        decodingNanos += System.nanoTime() - decodingStart;
        return result;
    }

    @Override
    public long preparationNanos() {
        return preparationNanos;
    }

    @Override
    public long decodingNanos() {
        return decodingNanos;
    }
}