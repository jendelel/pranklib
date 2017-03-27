package cz.siret.prank.lib;

import org.apache.commons.math3.exception.InsufficientDataException;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.PDBFileReader;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.Serializable;
import java.util.Arrays;
import java.util.stream.DoubleStream;

import cz.siret.prank.lib.utils.Utils;

public class ConservationComparison implements Serializable {
    private String fileName;
    private String scoreOrigin;
    private ConservationScore conservationScore;
    private ResidueNumber[] ligandIndices;
    private Structure protein;

    public ConservationComparison(String fileName,
                                  String scoreOrigin,
                                  ConservationScore conservationScore,
                                  ResidueNumber[] ligandIndices,
                                  Structure protein) {
        this.fileName = fileName;
        this.scoreOrigin = scoreOrigin;
        this.conservationScore = conservationScore;
        this.ligandIndices = ligandIndices;
        this.protein = protein;
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

    public ConservationScore getConservationScore() {
        return conservationScore;
    }

    public void setConservationScore(ConservationScore conservationScore) {
        this.conservationScore = conservationScore;
    }

    public double[] getLigandScores() {
        return Arrays.stream(ligandIndices)
                .mapToDouble((i) -> conservationScore.getScoreForResidue(i)).toArray();
    }

    public double[] getNonLigandScores() {
        return protein.getChains().stream()
                .flatMap(chain -> chain.getAtomGroups(GroupType.AMINOACID).stream()
                        .map(aa -> aa.getResidueNumber()))
                .filter(a -> Arrays.stream(ligandIndices).anyMatch(
                        b -> ResidueNumberWrapper.equalsPositional(a,b)))
                .mapToDouble(resNum -> conservationScore.getScoreForResidue(resNum))
                .toArray();
    }

    public double[] getProteinScores() {
        return protein.getChains().stream()
                .flatMap(chain -> chain.getAtomGroups(GroupType.AMINOACID).stream())
                .mapToDouble(aa -> conservationScore.getScoreForResidue(aa.getResidueNumber()))
                .toArray();
    }

    public ResidueNumber[] getLigandIndices() {
        return ligandIndices;
    }

    public void setLigandIndices(ResidueNumber[] ligandIndices) {
        this.ligandIndices = ligandIndices;
    }


    public double KSTest(boolean nonLigand) {
        try {
            if (nonLigand) {
                return new KolmogorovSmirnovTest().kolmogorovSmirnovTest(getNonLigandScores(),
                        getLigandScores());
            }
            return new KolmogorovSmirnovTest().kolmogorovSmirnovTest(getProteinScores(),
                    getLigandScores());
        } catch (InsufficientDataException e) {
            return -1;
        }
    }

    public double proteinAvg() {
        return Arrays.stream(getProteinScores()).average().getAsDouble();
    }

    public double ligandAvg() {
        return Arrays.stream(getLigandScores()).average().getAsDouble();
    }

    public double avgDifference() {
        return this.ligandAvg() - this.proteinAvg();
    }

    @Override
    public String toString() {
        return String.format("%s_%s:\t\t%e\t\t%e\t\t%f\t\t%f\t\t%f\t\t%f",
                getFileName(), getScoreOrigin(), KSTest(false), KSTest(true), proteinAvg(),
                DoubleStream.of(getNonLigandScores()).average().getAsDouble(),
                ligandAvg(), avgDifference());
    }

//    public static ConservationComparison merge(List<ConservationComparison> list, String origin) {
//        List<Double> proteinScores = new ArrayList<>();
//        List<Integer> ligandIndices = new ArrayList<>();
//        List<Pocket> pockets = new ArrayList<>();
//        List<Integer> truePockets = new ArrayList<>();
//        for (ConservationComparison cmp : list) {
//            // Shift the indices, so that they point to correct scores.
//            ligandIndices.addAll(IntStream.of(cmp.ligandIndices)
//                    .mapToObj(index -> Integer.valueOf(index + proteinScores.size()))
//                    .collect(Collectors.toList()));
//
//            // Merge score lists.
//            proteinScores.addAll(DoubleStream.of(cmp.getProteinScores()).mapToObj(Double::valueOf)
//                    .collect(Collectors.toList()));
//        }
//        return new ConservationComparison("merged", origin,
//                proteinScores.stream().mapToDouble(Double::doubleValue).toArray(),
//                ligandIndices.stream().mapToInt(Integer::intValue).toArray());
//    }


    public static ConservationComparison fromFiles(File pdbFile,
                                                   ConservationScore conservationScore,
                                                   String org,
                                                   File lingadIndices) throws
            IOException {
        try (InputStream pdbIn = new FileInputStream(pdbFile)) {
            PDBFileReader reader = new PDBFileReader();
            Structure s = reader.getStructure(pdbIn);

            try (InputStream indicesIn = new FileInputStream(lingadIndices)) {
                String content = Utils.convertStreamToString(indicesIn);
                String[] tokens = content.split("\n");
                ResidueNumber[] ligandIndices = Arrays.stream(tokens)
                        .map((str) -> ResidueNumber.fromString(str)).toArray(ResidueNumber[]::new);

                return new ConservationComparison(pdbFile.getName(), org,
                        conservationScore, ligandIndices, s);
            }
        }
    }

}
