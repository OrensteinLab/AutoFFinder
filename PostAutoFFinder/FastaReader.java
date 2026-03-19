package PostAutoFFinder;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Comparator;
import  java.nio.file.Files;
import java.io.*;

/**
 * The FastaReader class provides methods to read a FASTA file, split it into individual sequence files,
 * and generate reverse complement sequences for each sequence.
 */
public class FastaReader {

    /**
     * Splits a FASTA file into multiple files, each containing a single chromosome sequence.
     * The sequences are saved in a target folder created based on the input file path.
     * Additionally, generates reverse complement files for each chromosome sequence file.
     *
     * @param fastaFilePath the path to the input FASTA file
     * @throws IOException if an I/O error occurs while reading the input file or writing the output files
     */
    public static void splitFastaToFiles(String fastaFilePath) throws IOException {
        
        String splitFolderPathStr = createTargetFolder(fastaFilePath, false);
        if (splitFolderPathStr == null) {
            return ;
        }

        // Open the file for reading
        BufferedReader reader = new BufferedReader(new FileReader(fastaFilePath));
        
        String line = reader.readLine();
        String header = null;
        // we will use this to append all the short sequences after the header into one string
        StringBuilder sequence = new StringBuilder();
        
        while (line != null) {
            if (line.startsWith(">")) {
                // If the line starts with ">", it is the header line
                if (header != null) {
                    // If we have already processed a sequence, save it to a file
                    saveSequenceToFile(header, sequence.toString(), splitFolderPathStr);
                    sequence = new StringBuilder();
                }
                header = line;
            } else {
                // Otherwise, it is part of the sequence
                sequence.append(line.toUpperCase());
            }
            line = reader.readLine();
        }

        // Save the last sequence to a file
        saveSequenceToFile(header, sequence.toString(), splitFolderPathStr);
        
        // Close the file
        reader.close();
        // Generate reverse complement of the sequences
        System.out.println("Saving reverse complement sequences.");
        splitFolderPathStr = createTargetFolder(fastaFilePath, true);
        File[] files = FastaReader.getChrFiles(fastaFilePath, false);
            // Process each file to generate reverse complement and save it
            for (File file : files) {
                FastaReader.generateReverseComplementFile(file);
                
            }
    }
    
    /**
     * Creates a target folder for the split files of the given FASTA file.
     * If the folder already exists, the method will not create a new one and will return null.
     * Otherwise, it will create the folder and return its path as a string.
     *
     * @param fastaFilePath the path to the FASTA file
     * @return the path to the created split folder as a string, or null if the folder already exists
     * @throws IOException if an I/O error occurs
     */
    private static String createTargetFolder(String fastaFilePath, Boolean reverse) throws IOException {
        // If there is already a split of the fasta file, then do not split
        Path filePath = Path.of(fastaFilePath);
        String fileName = filePath.getFileName().toString();
        String splitFolderPathStr = null;
        if (reverse) {
            splitFolderPathStr = filePath.toAbsolutePath().getParent().toString() +
            "/" + fileName.substring(0, fileName.lastIndexOf('.')) + "_split_rc/";
        }
        else{
            splitFolderPathStr = filePath.toAbsolutePath().getParent().toString() +
            "/" + fileName.substring(0, fileName.lastIndexOf('.')) + "_split/";
        }
        Path splitFolderPath = Path.of(splitFolderPathStr);

        if (Files.exists(splitFolderPath)) {
            System.out.println(fastaFilePath + " split is alredy exists, exits without splitting.");
            return null;
        }
        else {
            // Create the folder of the split
            new File(splitFolderPathStr).mkdirs();
        }

        return splitFolderPathStr;
    }

    /**
     * Saves a given sequence to a file.
     *
     * @param header the header of the sequence, from which the sequence ID is extracted
     * @param sequence the sequence to be saved
     * @param splitFolderPathStr the path to the folder where the file will be saved
     * @throws IOException if an I/O error occurs
     */
    private static void saveSequenceToFile(String header, String sequence, String splitFolderPathStr) throws IOException {
        // Extract the sequence ID from the header
        String sequenceID = header.substring(1).trim().replaceAll("\\s", "_");
        
        // Create a new file for the sequence
        String filename = splitFolderPathStr + sequenceID + ".txt";
        FileWriter writer = new FileWriter(filename);
        
        // Write the sequence to the file
        writer.write(sequence + "\n");
        
        // Close the file
        writer.close();
        
        System.out.println("Saved sequence " + sequenceID + " to file " + filename);
    }

    /**
     * Retrieves an array of File objects representing the chromosome files
     * located in a directory derived from the given FASTA file path.
     *
     * @param fastaFilePath the path to the FASTA file
     * @return an array of File objects representing the chromosome files
     * @throws IOException if an I/O error occurs
     */
    public static File[] getChrFiles(
        String fastaFilePath, Boolean reverse) throws IOException {
            Path filePath = Path.of(fastaFilePath);
            String fileName = filePath.getFileName().toString();
            String splitFolderPathStr = null;
            if (reverse) {
                splitFolderPathStr = filePath.toAbsolutePath().getParent().toString() +
                    "/" + fileName.substring(0, fileName.lastIndexOf('.')) + "_split_rc/";
                
            }
            else{
                splitFolderPathStr = filePath.toAbsolutePath().getParent().toString() +
                    "/" + fileName.substring(0, fileName.lastIndexOf('.')) + "_split/";
            }
            File dir = new File(splitFolderPathStr);
            File[] files = dir.listFiles();
            // Sort the files by name
            if (files != null) {
                Arrays.sort(files, Comparator.comparing(File::getName));
            }

            return files;
    }

    /**
     * Retrieves an array of File objects representing the chromosome files
     * located in a directory derived from the given FASTA file path.
     *
     * @param fastaFilePath the path to the FASTA file
     * @return an array of File objects representing the chromosome files
     * @throws IOException if an I/O error occurs
     */
    public static File[] getRcChrFiles(String fastaFilePath) throws IOException {
        Path filePath = Path.of(fastaFilePath);
        String fileName = filePath.getFileName().toString();
        String splitFolderPathStr = filePath.toAbsolutePath().getParent().toString() +
                "/" + fileName.substring(0, fileName.lastIndexOf('.')) + "_split/";
        File dir = new File(splitFolderPathStr);
        File[] files = dir.listFiles();

        return files;
    }

    /**
     * Generates the reverse complement of the sequence in the given file and saves it as another file.
     *
     * @param file The input file containing the sequence.
     * @throws IOException If an I/O error occurs.
     */
    public static void generateReverseComplementFile(File file) throws IOException {
        // Read the sequence from the file
        String sequence = Files.readString(file.toPath());

        // Generate the reverse complement of the sequence
        String reverseComplement = reverseComplement(sequence);

        // Define the output file path
        String outputFilePath = file.getParent() + "_rc/" + file.getName();

        // Write the reverse complement to the output file
        try (FileWriter writer = new FileWriter(outputFilePath)) {
            writer.write(reverseComplement);
        }
    }

    /**
     * Generates the reverse complement of a given DNA sequence.
     *
     * This method takes a DNA sequence as input and returns its reverse complement.
     * The reverse complement is formed by reversing the input sequence and then
     * replacing each nucleotide with its complement: 'A' with 'T', 'C' with 'G',
     * 'G' with 'C', and 'T' with 'A'. The method also handles 'N' and '-' characters
     * by keeping them unchanged.
     *
     * @param seq the DNA sequence to be reversed and complemented
     * @return the reverse complement of the input DNA sequence
     */
    public static String reverseComplement(String seq) {
        StringBuilder rcSeq = new StringBuilder();

        for (int i = seq.length() - 1; i >= 0; i--) {
            char c = seq.charAt(i);

            switch (c) {
                case 'A':
                    rcSeq.append('T');
                    break;
                case 'C':
                    rcSeq.append('G');
                    break;
                case 'G':
                    rcSeq.append('C');
                    break;
                case 'T':
                    rcSeq.append('A');
                    break;
                case 'N':
                    rcSeq.append('N');
                    break;
                case '-':
                    rcSeq.append('-');
                    break;
                default:
                    break;
            }
        }

        return rcSeq.toString();
    }

    /**
     * Generates the complementary DNA sequence for a given DNA sequence.
     * 
     * This method takes a DNA sequence as input and returns its complementary
     * sequence. The complementary bases are as follows:
     * - 'A' is complemented by 'T'
     * - 'C' is complemented by 'G'
     * - 'G' is complemented by 'C'
     * - 'T' is complemented by 'A'
     * - 'N' remains 'N'
     * - '-' remains '-'
     * 
     * @param seq the input DNA sequence
     * @return the complementary DNA sequence
     */
    public static String complement(String seq) {
        StringBuilder cSeq = new StringBuilder();

        for (int i = 0; i < seq.length(); i++) {
            char c = seq.charAt(i);

            switch (c) {
                case 'A':
                    cSeq.append('T');
                    break;
                case 'C':
                    cSeq.append('G');
                    break;
                case 'G':
                    cSeq.append('C');
                    break;
                case 'T':
                    cSeq.append('A');
                    break;
                case 'N':
                    cSeq.append('N');
                    break;
                case '-':
                    cSeq.append('-');
                    break;
                default:
                    break;
            }
        }

        return cSeq.toString();
    }

    public static void main(String[] args) throws IOException {
        splitFastaToFiles("genomes/hg38_only_chrs.fa");
    }
}