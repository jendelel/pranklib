package cz.siret.prank.lib;

import com.google.gson.Gson;

import com.univocity.parsers.tsv.TsvParser;
import com.univocity.parsers.tsv.TsvParserSettings;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.Structure;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;

public class ConservationScore implements Serializable {
    private Map<Integer, Double> scores;
    private final transient Logger logger = LoggerFactory.getLogger(getClass());

    private ConservationScore(Map<Integer, Double> scores) {
        this.scores = scores;
    }

    private static class AA {
        public String letter;
        public double score;
        public int index;

        public AA(String letter, double score, int index) {
            this.letter = letter;
            this.score = score;
            this.index = index;
        }
    }

    public double getScoreForResidue(int residueNum) {
        Double res = scores.get(residueNum);
        if (res == null) {
            return 0;
        } else {
            return res.doubleValue();
        }
    }

    public Map<Integer, Double> getScoreMap() {
        return scores;
    }

    public int size() {
        return this.scores.size();
    }

    public void toJson(File scoreFile) throws FileNotFoundException {
        Utils.stringToFile((new Gson()).toJson(this), scoreFile);
    }

    public static ConservationScore fromJson(File scoreFile) throws IOException {
        Gson gson = new Gson();
        try (BufferedReader reader = new BufferedReader(new FileReader(scoreFile))) {
            return gson.fromJson(reader, ConservationScore.class);
        }
    }

    public enum ScoreFormat {
        ConCavityFormat,
        JSDFormat
    }

    private static List<AA> loadScoreFile(File scoreFile, ScoreFormat format) {
        TsvParserSettings settings = new TsvParserSettings();
        settings.setLineSeparatorDetectionEnabled(true);
        TsvParser parser = new TsvParser(settings);
        List<String[]> lines = parser.parseAll(scoreFile);
        List<AA> result = new ArrayList<>(lines.size());
        for (String[] line : lines) {
            int index = -1;
            double score = 0;
            String letter = "-";
            switch (format) {
                case ConCavityFormat:
                    index = Integer.parseInt(line[0]);
                    letter = line[1];
                    score = Double.parseDouble(line[2]);
                    break;
                case JSDFormat:
                    index = Integer.parseInt(line[0]);
                    score = Double.parseDouble(line[1]);
                    letter = line[2].substring(0, 1);
                    break;
            }
            if (letter != "-") {
                result.add(new AA(letter, score, index));
            }
        }
        return result;
    }

    public static ConservationScore fromFiles(Structure structure,
                                              Function<String, File> scoresFiles)
            throws FileNotFoundException {
        return fromFiles(structure, scoresFiles, ScoreFormat.JSDFormat);
    }

    /**
     * Parses conservation scores created from HSSP database and Jensen-Shannon divergence.
     *
     * @param structure  Protein BioJava structure
     * @param scoreFiles Map from chain ids to files
     * @param format     Score format (JSD or ConCavity), default: JSD
     * @return new instance of ConservationScore (map from residual numbers to conservation scores)
     */
    public static ConservationScore fromFiles(Structure structure,
                                              Function<String, File> scoreFiles,
                                              ScoreFormat format) throws FileNotFoundException {
        Map<Integer, Double> scores = new HashMap<>();
        for (Chain chain : structure.getChains()) {
            if (chain.getAtomGroups(GroupType.AMINOACID).size() <= 0) {
                continue;
            }
            String chainId = chain.getChainID();
            chainId = chainId.trim().isEmpty() ? "A" : chainId;
            List<AA> chainScores = null;
            File scoreFile = scoreFiles.apply(chainId);
            if (scoreFile.exists()) {
                chainScores = ConservationScore.loadScoreFile(scoreFile, format);
            }

            int indexCounter = 0;
            int AACounter = 0;
            for (Group aa : chain.getAtomGroups(GroupType.AMINOACID)) {
                int resNum = aa.getResidueNumber().getSeqNum();
                if (chainScores == null) { // Score file not found.
                    scores.put(resNum, 0.0);
                    continue;
                }

                while (AACounter < chainScores.size() &&
                        chainScores.get(AACounter).index < indexCounter) {
                    AACounter++;
                }
                if (AACounter >= chainScores.size()) {
                    scores.put(resNum, 0.0);
                    indexCounter++;
                    continue;
                }

                if (chainScores.get(AACounter).index != indexCounter) {
                    scores.put(resNum, 0.0);
                } else {
                    if (!aa.getChemComp().getOne_letter_code().toUpperCase().equals(
                            chainScores.get(AACounter).letter)) {
                        // "Letters do not match."
                        scores.put(resNum, 0.0);
                        continue;
                    }
                    Double score = chainScores.get(AACounter).score;
                    score = score < 0 ? 0.0 : score;
                    scores.put(resNum, score);
                }
                indexCounter++;
            }
        }
        return new ConservationScore(scores);
    }

}
