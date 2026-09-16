package PostAutoFFinder;

import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;

public final class FpgaCandidateSource implements CandidateSource {
    private static final int IO_BURST_BYTES = 64;
    private static final int CONFIG_BYTES = 32 * 128;
    private static final int FPGA_TARGET_COUNT = 128;

    private final String xclbinPath;
    private final String patternsPath;
    private final int editDistance;

    public FpgaCandidateSource(
            String nativeLibraryPath, String xclbinPath, String patternsPath, int editDistance) {
        if (editDistance < 0 || editDistance > 6) {
            throw new IllegalArgumentException("FPGA edit distance must be between 0 and 6");
        }
        FpgaNative.loadLibrary(nativeLibraryPath);
        this.xclbinPath = xclbinPath;
        this.patternsPath = patternsPath;
        this.editDistance = editDistance;
    }

    @Override
    public AutoOffTargetSearchAlign.PositionsRes load(
            File chromosomeFile, String strand, int targetCount) throws IOException {
        if (targetCount != FPGA_TARGET_COUNT) {
            throw new IllegalArgumentException(
                    "This FPGA image requires exactly " + FPGA_TARGET_COUNT
                            + " targets, but received " + targetCount);
        }

        int outputCapacity = calculateOutputCapacity(Files.size(chromosomeFile.toPath()));
        ByteBuffer output = ByteBuffer.allocateDirect(outputCapacity).order(ByteOrder.LITTLE_ENDIAN);
        int bytesWritten = FpgaNative.run(
                xclbinPath,
                chromosomeFile.getAbsolutePath(),
                patternsPath,
                editDistance,
                output);
        if (bytesWritten <= 0 || bytesWritten > output.capacity()) {
            throw new IOException("Native FPGA runner returned invalid byte count " + bytesWritten);
        }

        output.limit(bytesWritten);
        return RelevBinaryDecoder.decode(output, targetCount);
    }

    private static int calculateOutputCapacity(long chromosomeBytes) throws IOException {
        long inputBytes = chromosomeBytes + CONFIG_BYTES;
        long capacity = ((inputBytes / IO_BURST_BYTES) + 1) * IO_BURST_BYTES;
        if (capacity > Integer.MAX_VALUE) {
            throw new IOException("Chromosome is too large for one FPGA buffer");
        }
        return (int) capacity;
    }
}