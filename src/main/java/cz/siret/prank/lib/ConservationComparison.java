package cz.siret.prank.lib;

import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.PDBFileReader;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.Serializable;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;

public class ConservationComparison implements Serializable {
    private String fileName;
    private String scoreOrigin;
    private double[] proteinScores;
    private double[] ligandScores;

    public ConservationComparison(String fileName, String scoreOrigin, double[] proteinScores,
                                  double[] ligandScores) {
        this.fileName = fileName;
        this.scoreOrigin = scoreOrigin;
        this.proteinScores = proteinScores;
        this.ligandScores = ligandScores;
    }

    public String getFileName() {
        return fileName;
    }

    public void setFileName(String fileName) {
        this.fileName = fileName;
    }

    public String getScoreOrigin() {
        return scoreOrigin;
    }

    public void setScoreOrigin(String scoreOrigin) {
        this.scoreOrigin = scoreOrigin;
    }

    public double[] getProteinScores() {
        return proteinScores;
    }

    public void setProteinScores(double[] proteinScores) {
        this.proteinScores = proteinScores;
    }

    public double[] getLigandScores() {
        return ligandScores;
    }

    public void setLigandScores(double[] ligandScores) {
        this.ligandScores = ligandScores;
    }

    public double KSTest() {
        return new KolmogorovSmirnovTest().kolmogorovSmirnovTest(proteinScores, ligandScores);
    }

    public double proteinAvg() {
        return Arrays.stream(proteinScores).average().getAsDouble();
    }

    public double ligandAvg() {
        return Arrays.stream(ligandScores).average().getAsDouble();
    }

    public double avgDifference() {
        return this.ligandAvg() - this.proteinAvg();
    }

    @Override
    public String toString() {
        return String.format("%s_%s:\t\t%e\t\t%f\t\t%f\t\t%f",
                getFileName(), getScoreOrigin(), KSTest(), proteinAvg(), ligandAvg(),
                avgDifference());
    }


    public static ConservationComparison fromFiles(File pdbFile, File conservationScore, String org,
                                                   File lingadIndices) throws IOException {
        try (InputStream scoreIn = new GZIPInputStream(new FileInputStream(conservationScore))) {
            String[] tokens = Utils.convertStreamToString(scoreIn).split(",");
            double[] proteinScores = Arrays.stream(tokens)
                    .mapToDouble((s) -> Double.parseDouble(s))
                    .map((i) -> i < 0 ? 0 : i).toArray();

            try (InputStream pdbIn = new FileInputStream(pdbFile)) {
                PDBFileReader reader = new PDBFileReader();
                Structure s = reader.getStructure(pdbIn);
                Chain ch = s.getChains().get(0);
                Sequence res = Sequence.fromChain(ch);

                try (InputStream indicesIn = new FileInputStream(lingadIndices)) {
                    String content = Utils.convertStreamToString(indicesIn);
                    tokens = content.split(",");
                    double[] ligandScores = Arrays.stream(tokens)
                            .mapToInt((str) -> Integer.parseInt(str))
                            .map((i)-> Arrays.binarySearch(res.getIndices(), i))
                            .mapToDouble((i) -> proteinScores[i]).toArray();

                    return new ConservationComparison(pdbFile.getName(), org,
                            proteinScores, ligandScores);
                }
            }
        }
    }

}
