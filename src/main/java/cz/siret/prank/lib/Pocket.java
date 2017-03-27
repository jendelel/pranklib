package cz.siret.prank.lib;

import org.biojava.nbio.structure.ResidueNumber;

import java.io.InputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;
import java.util.function.Function;

public class Pocket implements Serializable {
    private String name;
    private int rank;
    private float score;
    private int numOfConnollyPoints;
    private int numOfSurfaceAtoms;
    private float centerX;
    private float centerY;
    private float centerZ;
    private ResidueNumber[] residueIds;
    private Integer[] surfAtomIds;
    private ConservationScore conservationScores;
    private boolean truePocket;

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public int getRank() {
        return rank;
    }

    public void setRank(int rank) {
        this.rank = rank;
    }

    public float getScore() {
        return score;
    }

    public void setScore(float score) {
        this.score = score;
    }

    public int getNumOfConnollyPoints() {
        return numOfConnollyPoints;
    }

    public void setNumOfConnollyPoints(int numOfConnollyPoints) {
        this.numOfConnollyPoints = numOfConnollyPoints;
    }

    public int getNumOfSurfaceAtoms() {
        return numOfSurfaceAtoms;
    }

    public void setNumOfSurfaceAtoms(int numOfSurfaceAtoms) {
        this.numOfSurfaceAtoms = numOfSurfaceAtoms;
    }

    public float getCenterX() {
        return centerX;
    }

    public void setCenterX(float centerX) {
        this.centerX = centerX;
    }

    public float getCenterY() {
        return centerY;
    }

    public void setCenterY(float centerY) {
        this.centerY = centerY;
    }

    public float getCenterZ() {
        return centerZ;
    }

    public void setCenterZ(float centerZ) {
        this.centerZ = centerZ;
    }

    public ResidueNumber[] getResidueIds() {
        return residueIds;
    }

    public void setResidueIds(ResidueNumber[] residueIds) {
        this.residueIds = residueIds;
    }

    public Integer[] getSurfAtomIds() {
        return surfAtomIds;
    }

    public void setSurfAtomIds(Integer[] surfAtomIds) {
        this.surfAtomIds = surfAtomIds;
    }

    public ConservationScore getConservationScores() {
        return conservationScores;
    }
    public void setConservationScores(ConservationScore conservationScores) {
        this.conservationScores = conservationScores;
    }

    public double getConservationAvg() {
        return Arrays.stream(residueIds).mapToDouble(i-> conservationScores.getScoreForResidue(i))
                .average().getAsDouble();
    }

    private Function<Double, Double> sigmoid = x -> 1.0/(1.0+Math.exp(-x));
    public double getCombiningScore() {
        double a = -0.118042;
        double b = 0.667817;
        double c= -0.9389;
        double d = -0.0680766;
//        return getScore() * getConservationAvg();
        return sigmoid.apply((c * getConservationAvg() + d) * (a *getScore() + b));
    }

    public boolean isTruePocket() {
        return truePocket;
    }

    public void setTruePocket(boolean truePocket) {
        this.truePocket = truePocket;
    }

    public static List<Pocket> parseCSVPrediction(InputStream inputStream,
                                                  List<Integer> truePockets,
                                                  ConservationScore conservationScores) {
        // name,rank,score,connolly_points,surf_atoms,center_x,center_y,center_z,residue_ids,
        // surf_atom_ids
        Scanner scanner = new Scanner(inputStream);
        scanner.nextLine(); // Skip the header line
        List<Pocket> res = new ArrayList<>();
        while (scanner.hasNextLine()) {
            Pocket p = new Pocket();
            String[] tokens = scanner.nextLine().split(",");
            p.setName(tokens[0]);
            p.setRank(Integer.parseInt(tokens[1]));
            p.setScore(Float.parseFloat(tokens[2]));
            p.setNumOfConnollyPoints(Integer.parseInt(tokens[3]));
            p.setNumOfSurfaceAtoms(Integer.parseInt(tokens[4]));
            p.setCenterX(Float.parseFloat(tokens[5]));
            p.setCenterY(Float.parseFloat(tokens[6]));
            p.setCenterZ(Float.parseFloat(tokens[7]));
            p.setResidueIds(Arrays.stream(tokens[8].split(" "))
                    .map((s) -> ResidueNumber.fromString(s)).toArray(ResidueNumber[]::new));
            p.setSurfAtomIds(Arrays.stream(tokens[9].split(" "))
                    .map((s) -> Integer.parseInt(s)).toArray(Integer[]::new));
            if (conservationScores != null) {
                p.conservationScores = conservationScores;
            }
            res.add(p);
        }

        if (truePockets != null) {
            for (int index : truePockets) {
                if (index != -1) {
                    res.get(index - 1).setTruePocket(true);
                }
            }
        }

        return res;
    }
}
