package PostAutoFFinder;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Arrays;

public final class RelevBinaryDecoder {
    private static final int MATCH_RECORD_BYTES = 8;
    private static final int TARGETS_PER_GROUP = 16;

    private RelevBinaryDecoder() {
    }

    public static AutoOffTargetSearchAlign.PositionsRes decode(
            ByteBuffer source, int targetCount) {
        if (targetCount <= 0) {
            throw new IllegalArgumentException("targetCount must be positive");
        }

        ByteBuffer records = source.duplicate().order(ByteOrder.LITTLE_ENDIAN);
        IntArrayBuilder positions = new IntArrayBuilder();
        IntArrayBuilder targetIds = new IntArrayBuilder();
        int[] lastPosition = new int[targetCount];
        Arrays.fill(lastPosition, -1);

        boolean foundTerminator = false;
        while (records.remaining() >= MATCH_RECORD_BYTES) {
            int position = records.getInt();
            int matchMask = Short.toUnsignedInt(records.getShort());
            int groupId = Byte.toUnsignedInt(records.get());
            int valid = Byte.toUnsignedInt(records.get());

            if (valid == 0) {
                foundTerminator = true;
                break;
            }
            if (valid != 1 || position < 0 || matchMask == 0) {
                throw new IllegalArgumentException("Invalid FPGA match record");
            }

            int groupOffset = groupId * TARGETS_PER_GROUP;
            for (int bit = 0; bit < TARGETS_PER_GROUP; bit++) {
                if ((matchMask & (1 << bit)) == 0) {
                    continue;
                }

                int targetId = groupOffset + bit;
                if (targetId >= targetCount) {
                    throw new IllegalArgumentException(
                            "FPGA match references target " + targetId
                                    + " but only " + targetCount + " targets were provided");
                }
                if (position < lastPosition[targetId]) {
                    throw new IllegalStateException(
                            "FPGA positions decreased for target " + targetId + ": "
                                    + lastPosition[targetId] + " to " + position);
                }
                if (position == lastPosition[targetId]) {
                    continue;
                }

                positions.add(position);
                targetIds.add(targetId);
                lastPosition[targetId] = position;
            }
        }

        if (!foundTerminator) {
            throw new IllegalArgumentException("FPGA output does not contain a termination record");
        }
        return new AutoOffTargetSearchAlign.PositionsRes(
                targetIds.toArray(), positions.toArray());
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
}