package PostAutoFFinder;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;

public final class RelevBinaryDecoderTest {
    private RelevBinaryDecoderTest() {
    }

    public static void main(String[] args) throws Exception {
        decodesInterleavedGroupsWithoutSorting();
        matchesTextParserOutput();
        rejectsDecreasingPositionsForOneTarget();
        rejectsMissingTerminator();
        System.out.println("ReLev binary decoder checks passed");
    }

    private static void decodesInterleavedGroupsWithoutSorting() {
        ByteBuffer records = buffer(7);
        putMatch(records, 100, 0b0000_0000_0000_0101, 0);
        putMatch(records, 300, 0b0000_0000_0000_0001, 0);
        putMatch(records, 150, 0b0000_0000_0000_0010, 1);
        putMatch(records, 200, 0b0000_0000_0000_0010, 1);
        putMatch(records, 300, 0b0000_0000_0000_0001, 0);
        putMatch(records, 500, 0b0000_0000_0000_0001, 0);
        putTerminator(records);

        AutoOffTargetSearchAlign.PositionsRes result =
                RelevBinaryDecoder.decode(flip(records), 32);

        assertCandidate(result, 0, 0, 100);
        assertCandidate(result, 1, 2, 100);
        assertCandidate(result, 2, 0, 300);
        assertCandidate(result, 3, 17, 150);
        assertCandidate(result, 4, 17, 200);
        assertCandidate(result, 5, 0, 500);
        if (result.size() != 6) {
            throw new AssertionError("Expected 6 decoded candidates, got " + result.size());
        }
    }

    private static void matchesTextParserOutput() throws Exception {
        ByteBuffer records = buffer(5);
        putMatch(records, 100, 0b0101, 0);
        putMatch(records, 300, 0b0001, 0);
        putMatch(records, 150, 0b0010, 1);
        putMatch(records, 200, 0b0010, 1);
        putTerminator(records);
        AutoOffTargetSearchAlign.PositionsRes binary =
                RelevBinaryDecoder.decode(flip(records), 32);

        String temporaryDirectory = System.getenv("TMPDIR");
        Path textFile = temporaryDirectory == null
            ? Files.createTempFile("relev-decoder", ".txt")
            : Files.createTempFile(Path.of(temporaryDirectory), "relev-decoder", ".txt");
        try {
            Files.writeString(textFile, String.join("\n",
                    "100:0", "100:2", "300:0", "150:17", "200:17") + "\n");
            AutoOffTargetSearchAlign.PositionsRes text =
                    AutoOffTargetSearchAlign.parseAutomataResultsFile(
                            textFile.toString(), 32);
            if (binary.size() != text.size()) {
                throw new AssertionError("Binary and text candidate counts differ");
            }
            for (int index = 0; index < binary.size(); index++) {
                assertCandidate(binary, index, text.idAt(index), text.positionAt(index));
            }
        } finally {
            Files.deleteIfExists(textFile);
        }
    }

    private static void rejectsDecreasingPositionsForOneTarget() {
        ByteBuffer records = buffer(3);
        putMatch(records, 300, 1, 0);
        putMatch(records, 200, 1, 0);
        putTerminator(records);
        expectFailure(IllegalStateException.class,
                () -> RelevBinaryDecoder.decode(flip(records), 16));
    }

    private static void rejectsMissingTerminator() {
        ByteBuffer records = buffer(1);
        putMatch(records, 100, 1, 0);
        expectFailure(IllegalArgumentException.class,
                () -> RelevBinaryDecoder.decode(flip(records), 16));
    }

    private static ByteBuffer buffer(int recordCount) {
        return ByteBuffer.allocate(recordCount * 8).order(ByteOrder.LITTLE_ENDIAN);
    }

    private static ByteBuffer flip(ByteBuffer buffer) {
        return buffer.flip().order(ByteOrder.LITTLE_ENDIAN);
    }

    private static void putMatch(ByteBuffer buffer, int position, int mask, int groupId) {
        buffer.putInt(position);
        buffer.putShort((short) mask);
        buffer.put((byte) groupId);
        buffer.put((byte) 1);
    }

    private static void putTerminator(ByteBuffer buffer) {
        buffer.putLong(0L);
    }

    private static void assertCandidate(
            AutoOffTargetSearchAlign.PositionsRes result,
            int index, int expectedId, int expectedPosition) {
        if (result.idAt(index) != expectedId || result.positionAt(index) != expectedPosition) {
            throw new AssertionError(String.format(
                    "Candidate %d: expected %d:%d, got %d:%d",
                    index, expectedPosition, expectedId,
                    result.positionAt(index), result.idAt(index)));
        }
    }

    private static void expectFailure(Class<? extends Throwable> type, Runnable action) {
        try {
            action.run();
        } catch (Throwable throwable) {
            if (type.isInstance(throwable)) {
                return;
            }
            throw new AssertionError("Expected " + type.getSimpleName(), throwable);
        }
        throw new AssertionError("Expected " + type.getSimpleName());
    }
}