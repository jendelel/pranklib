package cz.siret.prank.lib.utils;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.io.PDBFileReader;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

public class BioUtils {

    public static String pdbToFasta(File pdbFile) throws IOException, StructureException {
        return pdbToFasta(pdbFile, null);
    }

    public static String chainToFasta(Chain chain, String header) {
        StringBuilder result = new StringBuilder();
        String seq = chain.getAtomSequence().trim();
        if (seq.length() == 0) return null;
        // Print the header like this: >4X09:A
        result.append(header).append(chain.getChainID()).append('\n');
        // Print the chain sequence and wrap lines at 80 characters
        for (int i = 0; i < seq.length(); i++) {
            if (i != 0 && i % 80 == 0) result.append('\n');
            result.append(seq.charAt(i));
        }
        return result.toString();
    }

    public static String pdbToFasta(File pdbFile, String chainId) throws IOException,
            StructureException {
        PDBFileReader pdbReader = new PDBFileReader();
        boolean isGzipped = pdbFile.getName().endsWith(".gz");
        try (InputStream inputStream = isGzipped ?
                new GZIPInputStream(new FileInputStream(pdbFile)) :
                new FileInputStream(pdbFile)) {
            Structure structure = pdbReader.getStructure(inputStream);
            StringBuilder output = new StringBuilder();
            String header = ">" + structure.getPDBHeader().getIdCode() + ":";
            if (chainId == null) {
                for (Chain chain : structure.getChains()) {
                    String chainFasta = chainToFasta(chain, header);
                    if (chainFasta != null) {
                        output.append(chainFasta);
                        output.append('\n');
                    }
                }
            } else {
                Chain chain = chainId == null ?
                        structure.getChains().get(0) :
                        structure.getChainByPDB(chainId);
                String chainFasta = chainToFasta(chain, header);
                if (chainFasta != null) {
                    output.append(chainFasta);
                    output.append('\n');
                }
            }
            return output.toString();
        }
    }

    public static void dirToFasta(File dir) throws IOException, StructureException {
        File[] listOfFiles = dir.listFiles();

        for (File f : listOfFiles) {
            if (f.isFile()) {
                if (f.getName().endsWith(".pdb")) {
                    String newFileName = f.getName().concat(".fasta");
                    String fastaContent = pdbToFasta(f);
                    File newFile = new File(f.getParent(), newFileName);
                    Utils.stringToFile(fastaContent, newFile);
                }
            }
        }
    }


}
