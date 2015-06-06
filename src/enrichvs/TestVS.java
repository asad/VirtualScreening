/*
 * Copyright (c) 2015, Asad
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
package enrichvs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Locale;
import java.util.logging.Level;
import java.util.logging.Logger;
import utility.DynamicArrayOfBoolean;
import utility.DynamicArrayOfDouble;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class TestVS {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        File file = new File("data/dud_egfr.csv");
        Locale locale = new Locale("en", "UK");
        String pattern = "#.###";

        DecimalFormat decimalFormat = (DecimalFormat) NumberFormat.getNumberInstance(locale);
        decimalFormat.applyPattern(pattern);

        DynamicArrayOfDouble energy = new DynamicArrayOfDouble();
        DynamicArrayOfBoolean flags = new DynamicArrayOfBoolean();

        int index = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String header = br.readLine();
            String line = "";
            String cvsSplitBy = ",";

            while ((line = br.readLine()) != null) {
                // use comma as separator
                String[] country = line.split(cvsSplitBy);
                // System.out.println(Arrays.toString(country));

                energy.put(index, Double.parseDouble(country[1]));
                boolean parseBoolean = false;
                if (Integer.parseInt(country[2]) == 1) {
                    parseBoolean = true;
                }
                flags.put(index, parseBoolean);
                index++;
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(EnrichmentAssessment.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(EnrichmentAssessment.class.getName()).log(Level.SEVERE, null, ex);
        }

        double[] energies = energy.getArray();
        boolean[] boolArrays = flags.getArray();

        System.out.println("Number of data points: " + index);

        EnrichmentAssessment virtualScreening = new EnrichmentAssessment(energies, boolArrays);
        // Expected Bedroc 0.591155 (decreasing=FALSE), 0.3848914 (decreasing=TRUE)
        System.out.println("Virtual Screening BEDROC: " + decimalFormat.format(virtualScreening.bedroc(0.20d, false)));
        // Expected Enrichment_factor 3.108108 (decreasing=FALSE), 1.891892 (decreasing=TRUE)
        System.out.println("Virtual Screening EF: " + decimalFormat.format(virtualScreening.enrichment_factor(0.05d, false)));
        // Expected Enrichment_factor 1.021421 (decreasing=FALSE), 0.980375 (decreasing=TRUE)
        System.out.println("Virtual Screening RIE: " + decimalFormat.format(virtualScreening.rie(0.20d, false)));
        // Expected Enrichment_factor 0.07586939 (decreasing=FALSE), 0.05835606 (decreasing=TRUE)
        System.out.println("Virtual Screening AUAC: " + decimalFormat.format(virtualScreening.auac(0.05d, false)));
        // Expected Enrichment_factor 0.07668251 (decreasing=FALSE), 0.05862467 (decreasing=TRUE)
        System.out.println("Virtual Screening AUC: " + decimalFormat.format(virtualScreening.auc(0.05d, false)));

        System.out.println("Done");
    }

}
