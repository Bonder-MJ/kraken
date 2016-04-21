/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package systemsgenetics.krakengenomesizenormalization;

import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

/**
 *
 * @author MarcJan
 */
public class CalculateNormalizationFactor {

    static final Pattern TAB_PATTERN = Pattern.compile("\\t");
    static final Pattern NCBI_PATTERN = Pattern.compile("\\t\\|\\t");
    static final Pattern PIPE_PATTERN = Pattern.compile("\\|");

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        CommandLineParser parser = new GnuParser();
        Options options = new Options();

        Option refferenceGenomeSizesIn = OptionBuilder.withArgName("path(s)").hasArg().withDescription("List of NCBI genome information or information directly calculated on fasta's (can be comma seperated). Alternatively the folder with Fastas.").withLongOpt("genomeInformation").create("g");
        Option secondRefferenceGenomeSizesIn = OptionBuilder.withArgName("path(s)").hasArg().withDescription("List a secondary set of genome size information (can be comma seperated), this is only used when entries are missing from the first file (cannot be comma seperated!).").withLongOpt("secondaryGenomeInformation").create("g2");
        Option ncbiTaxanomicTree = OptionBuilder.withArgName("path").hasArg().withDescription("NCBI taxonomy, as taken from nodes.dmp").withLongOpt("taxonomyInformation").create("t");
        Option ncbiGiToTaxId = OptionBuilder.withArgName("path").hasArg().withDescription("NCBI linking file from Gi to TaxId, as taken from gi_taxid_nucl.dmp").withLongOpt("linkingInformation").create("l");
        Option fileOut = OptionBuilder.withArgName("path").hasArg().withDescription("Location and name of the output file, which can be used for normalization.").withLongOpt("OutputFile").create("o");
        Option fasta = OptionBuilder.withArgName("boolean").withDescription("If set will run the program in fasta mode, which requires the -l option.").create("f");
        Option fasta2 = OptionBuilder.withArgName("boolean").withDescription("If set will run the program in fasta2 mode, which requires the -l option and secondary refference genome size in.").create("f2");
        options.addOption(fileOut).addOption(refferenceGenomeSizesIn).addOption(ncbiTaxanomicTree).addOption(secondRefferenceGenomeSizesIn).addOption(fasta).addOption(fasta2).addOption(ncbiGiToTaxId);

        String[] referenceFiles = null;
        String[] secondReferenceFiles = null;
        File ncbiTaxanomicTreeFile = null;
        File ncbiGiToTaxIdFile = null;
        File outputFile = null;
        boolean fastaMode = false;
        boolean fastaMode2 = false;
        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
            HelpFormatter formatter = new HelpFormatter();

            if (cmd.hasOption("OutputFile") || cmd.hasOption("o")) {
                // initialise the member variable
                outputFile = new File(cmd.getOptionValue("OutputFile"));
            } else {
                System.out.println("Missing necesarray information");
                formatter.printHelp("ant", options);
                System.exit(0);
            }
            if (cmd.hasOption("genomeInformation") || cmd.hasOption("g")) {
                // initialise the member variable
                referenceFiles = cmd.getOptionValue("genomeInformation").split(",");
            } else {
                System.out.println("Missing necesarray information");
                formatter.printHelp("ant", options);
                System.exit(0);
            }
            
            if (cmd.hasOption("secondaryGenomeInformation") || cmd.hasOption("g2")) {
                // initialise the member variable
                secondReferenceFiles = cmd.getOptionValue("secondaryGenomeInformation").split(",");
            }
            
            if (cmd.hasOption("taxonomyInformation") || cmd.hasOption("t")) {
                // initialise the member variable
                ncbiTaxanomicTreeFile = new File(cmd.getOptionValue("taxonomyInformation"));
            } else {
                System.out.println("Missing necesarray information");
                formatter.printHelp("ant", options);
                System.exit(0);
            }
            if (cmd.hasOption("linkingInformation") || cmd.hasOption("l")) {
                // initialise the member variable
                ncbiGiToTaxIdFile = new File(cmd.getOptionValue("linkingInformation"));
            }
            fastaMode = cmd.hasOption("f");
            fastaMode2 = cmd.hasOption("f2");
            if ((fastaMode ^ fastaMode2) && (ncbiGiToTaxId == null || referenceFiles.length > 1)) {
                if(fastaMode2 && secondReferenceFiles!=null){
                    System.out.println("Check input requirements, either missing file to link gi to taxId or trying to use multiple folders for fasta parsing.");
                    formatter.printHelp("ant", options);
                    System.exit(0);
                } else {
                    System.out.println("Check input requirements, either missing file to link gi to taxId or trying to use multiple folders for fasta parsing.");
                    formatter.printHelp("ant", options);
                    System.exit(0);
                }
            }

        } catch (org.apache.commons.cli.ParseException ex) {
            Logger.getLogger(CalculateNormalizationFactor.class.getName()).log(Level.SEVERE, null, ex);
        }

        TIntObjectHashMap<taxonInformation> taxInfo = readTaxaInformation(ncbiTaxanomicTreeFile);
        System.out.println("Loaded: " + taxInfo.size() + " nodes from the NCBI nodes.dmp");

        if (!fastaMode && !fastaMode2) {
            readGenomeSizeInformation(referenceFiles, taxInfo);
            if (secondReferenceFiles != null) {
                addSecondGenomeSizeInformation(secondReferenceFiles, taxInfo);
                System.out.println("Loaded: " + referenceFiles.length + secondReferenceFiles.length + " files with taxonomy sizes");
            } else {
                System.out.println("Loaded: " + referenceFiles.length + " files with taxonomy sizes");
            }
        } else if(fastaMode) {
            readFastaAndCalculteGenomeSizeInformation(referenceFiles, ncbiGiToTaxIdFile, taxInfo);
            if (secondReferenceFiles != null) {
                addSecondGenomeSizeInformation(secondReferenceFiles, taxInfo);
                System.out.println("Loaded sizes from fasta files, and NCBI information files.");
            } else {
                System.out.println("Loaded sizes from fasta files.");
            }
        } else {
            readGenomeSizeInformation(referenceFiles, taxInfo);
            readSecondFastaAndCalculteGenomeSizeInformation(secondReferenceFiles, ncbiGiToTaxIdFile, taxInfo);
            System.out.println("Loaded sizes from fasta files, and NCBI information files.");
        }

        returnGenomeSizeInformation(taxInfo, outputFile);

    }

    private static TIntObjectHashMap<taxonInformation> readTaxaInformation(File ncbiTaxanomicTreeFile) {
        TIntObjectHashMap<taxonInformation> taxaInformation = new TIntObjectHashMap();

        if (ncbiTaxanomicTreeFile != null && ncbiTaxanomicTreeFile.exists()) {
            try {
                BufferedReader r = new BufferedReader(new FileReader(ncbiTaxanomicTreeFile.getAbsoluteFile()));
                String row;
                while ((row = r.readLine()) != null) {
                    String[] parts = NCBI_PATTERN.split(row);
                    taxaInformation.put(Integer.parseInt(parts[0]), new taxonInformation(new String(parts[1]), new String(parts[0])));
                }
                r.close();
            } catch (Exception ex) {
                Logger.getLogger(CalculateNormalizationFactor.class.getName()).log(Level.SEVERE, null, ex);
            }
        } else {
            System.out.println("Not a valid ncbi tree information.");
            System.exit(0);
        }

        return taxaInformation;
    }

    private static void readGenomeSizeInformation(String[] referenceFiles, TIntObjectMap<taxonInformation> taxInfo) {
        for (String fileName : referenceFiles) {
            if (fileName != null && new File(fileName).exists()) {
                try {
                    BufferedReader r = new BufferedReader(new FileReader(fileName));
                    String row;

                    int taxaId = 0;
                    int genomeSize = 1;
                    int factor = 1;
                    boolean firstRow = true;
                    while ((row = r.readLine()) != null) {
                        String[] parts = TAB_PATTERN.split(row);
                        if (firstRow) {
                            if (parts.length == 2) {
                                if (parts[taxaId] != null || parts[taxaId].equals("")) {
                                    Integer taxId = Integer.parseInt(parts[taxaId]);
                                    double gs = Double.parseDouble(parts[genomeSize]) * factor;
                                    if (taxInfo.containsKey(taxId)) {
                                        while (taxId != null) {
                                            if (taxInfo.containsKey(taxId)) {
                                                taxId = taxInfo.get(taxId).insertAverageGenomeSizeRecursive(gs);
                                            } else {
                                                System.out.println("TaxId not in map: " + taxId);
                                            }
                                        }
                                    }
                                }
                            } else {
                                int i = 0;
                                boolean foundTxId = false;
                                boolean foundSize = false;
                                for (String s : parts) {
                                    if (s.equals("TaxID")) {
                                        taxaId = i;
                                        foundTxId = true;
                                    } else if (s.startsWith("Size ")) {
                                        foundSize = true;
                                        genomeSize = i;
                                        if (s.equals("Size (Mb)")) {
                                            factor = 1000000;
                                        } else if (s.equals("Size (Kb)")) {
                                            factor = 1000;
                                        } else {
                                            System.out.println("Unrecognized factor for genome size: " + s);
                                            System.exit(0);
                                        }
                                    }
                                    i++;
                                }
                                if (!foundSize || !foundTxId) {
                                    System.out.println("Error, no TaxID or Size  column identified.");
                                    System.exit(-1);
                                }
                            }
                            firstRow = false;
                        } else if (parts[taxaId] != null || parts[taxaId].equals("")) {
                            Integer taxId = Integer.parseInt(parts[taxaId]);
                            double gs = Double.parseDouble(parts[genomeSize]) * factor;
                            if (taxInfo.containsKey(taxId)) {
                                while (taxId != null) {
                                    if (taxInfo.containsKey(taxId)) {
                                        taxId = taxInfo.get(taxId).insertAverageGenomeSizeRecursive(gs);
                                    } else {
                                        System.out.println("TaxId not in map: " + taxId);
                                    }
                                }
                            }
                        }
                    }
                    r.close();
                } catch (Exception ex) {
                    Logger.getLogger(CalculateNormalizationFactor.class.getName()).log(Level.SEVERE, null, ex);
                }
            } else {
                System.out.println("Not a valid ncbi tree information.");
                System.exit(0);
            }
        }
    }

    private static void readFastaAndCalculteGenomeSizeInformation(String[] referenceFiles, File ncbiGiToTaxIdFile, TIntObjectHashMap<taxonInformation> taxInfo) {
        TObjectIntHashMap<String> tempMap = new TObjectIntHashMap<String>();
        File folder = new File(referenceFiles[0]);
        System.out.print("\tNumber of fastas to read: "+folder.listFiles().length);
        for (File fileInFolder : folder.listFiles()) {
            try {
                BufferedReader r = new BufferedReader(new FileReader(fileInFolder.getAbsolutePath()));
                String row;
                String fastaHeader = null;
                StringBuilder body = new StringBuilder();
                while ((row = r.readLine()) != null) {
                    if (row.startsWith(">")) {
                        if(body.length()!=0){
                            tempMap.put(fastaHeader, body.length());
                            body = new StringBuilder();
                        }
                        fastaHeader = PIPE_PATTERN.split(row)[1];
                    } else {
                        body.append(row);
                    }
                }
                //Process
                tempMap.put(fastaHeader, body.length());
                
                r.close();
            } catch (Exception ex) {
                Logger.getLogger(CalculateNormalizationFactor.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        System.out.print("\tDone. \n\t Matching GI numbers to TaxIds now.");
        HashSet<String> visited = new HashSet<>();
        boolean done = false;
        try {
            BufferedReader r = new BufferedReader(new FileReader(ncbiGiToTaxIdFile.getAbsolutePath()));
            String row;
            while ((row = r.readLine()) != null && !done) {
                String[] parts = TAB_PATTERN.split(row);
                if(!visited.contains(parts[0]) && tempMap.keySet().contains(parts[0])){
                    Integer id = Integer.parseInt(parts[1]);
                    if(taxInfo.containsKey(id)){
                        while(id!=null){
                            id = taxInfo.get(id).insertAverageGenomeSizeRecursive(tempMap.get(parts[0]));
                        }
                    }
                    visited.add(parts[0]);
                    if(visited.size() == tempMap.size()){
                        done = true;
                    }
                }
                
            }
            r.close();
        } catch (Exception ex) {
            Logger.getLogger(CalculateNormalizationFactor.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        System.out.println("\tDone.");

    }

    private static void addSecondGenomeSizeInformation(String[] referenceFiles, TIntObjectMap<taxonInformation> taxInfo) {
        for (String fileName : referenceFiles) {
            if (fileName != null && new File(fileName).exists()) {
                try {
                    BufferedReader r = new BufferedReader(new FileReader(fileName));
                    String row;

                    int taxaId = 0;
                    int genomeSize = 1;
                    int factor = 1;
                    boolean firstRow = true;
                    while ((row = r.readLine()) != null) {
                        String[] parts = TAB_PATTERN.split(row);
                        if (firstRow) {
                            if (parts.length == 2) {
                                if (parts[taxaId] != null || parts[taxaId].equals("")) {
                                    Integer taxId = Integer.parseInt(parts[taxaId]);
                                    double gs = Double.parseDouble(parts[genomeSize]) * factor;
                                    if (taxInfo.containsKey(taxId) && taxInfo.get(taxId).getAverageGenomeSize().isNaN()) {
                                        while (taxId != null) {
                                            if (taxInfo.containsKey(taxId)) {
                                                taxId = taxInfo.get(taxId).insertAverageGenomeSizeRecursive(gs);
                                            } else {
                                                System.out.println("TaxId not in map: " + taxId);
                                            }
                                        }
                                    }
                                }
                            } else {
                                int i = 0;
                                boolean foundTxId = false;
                                boolean foundSize = false;
                                for (String s : parts) {
                                    if (s.equals("TaxID")) {
                                        taxaId = i;
                                        foundTxId = true;
                                    } else if (s.startsWith("Size ")) {
                                        foundSize = true;
                                        genomeSize = i;
                                        if (s.equals("Size (Mb)")) {
                                            factor = 1000000;
                                        } else if (s.equals("Size (Kb)")) {
                                            factor = 1000;
                                        } else {
                                            System.out.println("Unrecognized factor for genome size: " + s);
                                            System.exit(0);
                                        }
                                    }
                                    i++;
                                }
                                if (!foundSize || !foundTxId) {
                                    System.out.println("Error, no TaxID or Size  column identified.");
                                    System.exit(-1);
                                }
                            }
                            firstRow = false;
                        } else if (parts[taxaId] != null || parts[taxaId].equals("")) {
                            Integer taxId = Integer.parseInt(parts[taxaId]);
                            double gs = Double.parseDouble(parts[genomeSize]) * factor;
                            if (taxInfo.containsKey(taxId)) {
                                while (taxId != null) {
                                    if (taxInfo.containsKey(taxId) && taxInfo.get(taxId).getAverageGenomeSize().isNaN()) {
                                        taxId = taxInfo.get(taxId).insertAverageGenomeSizeRecursive(gs);
                                    } else {
                                        System.out.println("TaxId not in map: " + taxId);
                                    }
                                }
                            }
                        }
                    }
                } catch (Exception ex) {
                    Logger.getLogger(CalculateNormalizationFactor.class.getName()).log(Level.SEVERE, null, ex);
                }
            } else {
                System.out.println("Not a valid ncbi tree information.");
                System.exit(0);
            }
        }

    }

    private static void readSecondFastaAndCalculteGenomeSizeInformation(String[] referenceFiles, File ncbiGiToTaxIdFile, TIntObjectHashMap<taxonInformation> taxInfo) {
        TObjectIntHashMap<String> tempMap = new TObjectIntHashMap<String>();
        File folder = new File(referenceFiles[0]);
        System.out.print("\tNumber of fastas to read: "+folder.listFiles().length);
        for (File fileInFolder : folder.listFiles()) {
            try {
                BufferedReader r = new BufferedReader(new FileReader(fileInFolder.getAbsolutePath()));
                String row;
                String fastaHeader = null;
                StringBuilder body = new StringBuilder();
                while ((row = r.readLine()) != null) {
                    if (row.startsWith(">")) {
                        if(body.length()!=0){
                            tempMap.put(fastaHeader, body.length());
                            body = new StringBuilder();
                        }
                        fastaHeader = PIPE_PATTERN.split(row)[1];
                    } else {
                        body.append(row);
                    }
                }
                //Process
                tempMap.put(fastaHeader, body.length());
                
                r.close();
            } catch (Exception ex) {
                Logger.getLogger(CalculateNormalizationFactor.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        System.out.print("\tDone. \n\t Matching GI numbers to TaxIds now.");
        HashSet<String> visited = new HashSet<>();
        boolean done = false;
        try {
            BufferedReader r = new BufferedReader(new FileReader(ncbiGiToTaxIdFile.getAbsolutePath()));
            String row;
            while ((row = r.readLine()) != null && !done) {
                String[] parts = TAB_PATTERN.split(row);
                if(!visited.contains(parts[0]) && tempMap.keySet().contains(parts[0])){
                    Integer id = Integer.parseInt(parts[1]);
                    if(taxInfo.containsKey(id) && taxInfo.get(id).getAverageGenomeSize().isNaN()){
                        while(id!=null){
                            id = taxInfo.get(id).insertAverageGenomeSizeRecursive(tempMap.get(parts[0]));
                        }
                    }
                    visited.add(parts[0]);
                    if(visited.size() == tempMap.size()){
                        done = true;
                    }
                }
                
            }
            r.close();
        } catch (Exception ex) {
            Logger.getLogger(CalculateNormalizationFactor.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        System.out.println("\tDone.");

    }
    
    private static void returnGenomeSizeInformation(TIntObjectMap<taxonInformation> taxInfo, File outputFile) {
        
        try {
            BufferedWriter outSnp = new BufferedWriter(new FileWriter(outputFile));
            for (int t : taxInfo.keySet().toArray()) {
                if (!taxInfo.get(t).getAverageGenomeSize().isNaN()) {
                    long gs = Math.round(taxInfo.get(t).getAverageGenomeSize());
                    outSnp.write(t + "\t" + gs);
                    outSnp.write('\n');
                }
            }
            outSnp.close();
        } catch (IOException ex) {
            Logger.getLogger(CalculateNormalizationFactor.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

}
