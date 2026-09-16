package PostAutoFFinder;

import java.io.File;
import java.nio.file.Path;

public final class TextFileCandidateSource implements CandidateSource {
    private final String outputFolder;

    public TextFileCandidateSource(String outputFolder) {
        this.outputFolder = outputFolder;
    }

    @Override
    public AutoOffTargetSearchAlign.PositionsRes load(
            File chromosomeFile, String strand, int targetCount) {
        Path outputPath = getOutputPath(chromosomeFile, strand);
        if (!outputPath.toFile().exists()) {
            return null;
        }
        return AutoOffTargetSearchAlign.parseAutomataResultsFile(
                outputPath.toString(), targetCount);
    }

    private Path getOutputPath(File chromosomeFile, String strand) {
        String name = chromosomeFile.getName();
        int dotIndex = name.lastIndexOf('.');
        String baseName = dotIndex == -1 ? name : name.substring(0, dotIndex);
        String extension = dotIndex == -1 ? "" : name.substring(dotIndex);
        String suffix = strand.equals("-") ? "_rc" : "_fw";
        return Path.of(outputFolder, baseName + suffix + extension);
    }
}