package PostAutoFFinder;

import java.io.File;
import java.io.IOException;

public interface CandidateSource {
    AutoOffTargetSearchAlign.PositionsRes load(
            File chromosomeFile, String strand, int targetCount) throws IOException;

    default long preparationNanos() {
        return 0;
    }

    default long decodingNanos() {
        return 0;
    }
}