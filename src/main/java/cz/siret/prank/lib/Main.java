package cz.siret.prank.lib;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Main {

    public static void main(String[] main) {
        analyzeConservationScores();
    }


    public static void analyzeConservationScores() {
        String[] files = new String[]{"c.097.001.004_1g8ma", "c.134.001.001_1q7ta",
                "d.014.001.005_1h72c", "d.016.001.003_2gf3a", "d.032.001.002_1lqpa",
                "d.058.004.003_1n5qb", "d.058.009.001_1ir2a", "d.068.002.002_1ejda",
                "d.080.001.003_1gcza", "d.096.001.004_2pesa", "d.117.001.001_1kzia",
                "d.128.001.002_1bg0a", "d.142.001.009_2r84a", "d.157.001.001_2gfka",
                "d.176.001.001_1soxb", "d.207.001.001_1o26a", "e.005.001.001_1dgfa",
                "e.007.001.001_2f3da", "e.008.001.004_1gx6a"};

        String[] orgins = new String[]{"90", "50", "swiss", "swiss90"};

        Path dir = Paths.get("e:/School/MFF/Projects/Prank2Web/Experiments" +
                "/analyze_binding-residues_chen11/");
        for (String origin : orgins) {
            List<ConservationComparison> comparisons = new ArrayList<>(orgins.length);
            for (String file : files) {
                File indices = dir.resolve(file + ".pdb_binding-residues.txt")
                        .toFile();
                File score = dir.resolve(String.format("%s.pdb.%s.hom.gz", file, origin)).toFile();
                File pdb = dir.resolve(file + ".pdb").toFile();
                try {
                    ConservationComparison comparison =
                            ConservationComparison.fromFiles(pdb, score, origin, indices);
                    System.out.println(comparison.toString());
                    comparisons.add(comparison);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            double[] proteinScores = comparisons.stream()
                    .map(ConservationComparison::getProteinScores)
                    .flatMapToDouble(Arrays::stream).toArray();
            double[] ligandScores = comparisons.stream()
                    .map(ConservationComparison::getLigandScores)
                    .flatMapToDouble(Arrays::stream).toArray();
            ConservationComparison all = new ConservationComparison("all", origin,
                    proteinScores, ligandScores);
            System.out.println(all.toString());
        }
    }
}
