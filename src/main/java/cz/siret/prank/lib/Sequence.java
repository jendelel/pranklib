package cz.siret.prank.lib;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

public class Sequence implements Serializable {
    int[] indices;
    String[] seq;
    double[] scores;

    public String[] getSeq() {
        return seq;
    }

    public void setSeq(String[] seq) {
        this.seq = seq;
    }

    public double[] getScores() {
        return scores;
    }

    public void setScores(double[] scores) {
        this.scores = scores;
    }

    public int[] getIndices() {
        return indices;
    }

    public void setIndices(int[] indices) {
        this.indices = indices;
    }

    public static Sequence fromChain(Chain chain) {
        List<Integer> indices = new ArrayList<>();
        List<String> seq = new ArrayList<>();

        for (Group group : chain.getAtomGroups(GroupType.AMINOACID)) {
            String c = group.getChemComp().getOne_letter_code();
            if (!c.equals("?")) {
                seq.add(c);
                indices.add(group.getResidueNumber().getSeqNum());
            }
        }
        Sequence res = new Sequence();
        res.indices = indices.stream().mapToInt((s)->s).toArray();
        res.seq = seq.stream().toArray(String[]::new);
        return res;
    }
}
