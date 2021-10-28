/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2020 Erik Garrison
    Copyright © 2020      Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "Variant.h"
#include "split.h"
#include <string>
#include <sstream>
#include <iostream>
#include <getopt.h>
#include <iomanip>

using namespace std;
using namespace vcflib;


int main(int argc, char** argv) {

  if (argc == 2) {
    string h_flag = argv[1];
    if (h_flag == "-h" || h_flag == "--help") {
      cerr << "usage: " << argv[0] << " <[input file] >[output vcf]" << endl << endl
             << "Adds variant allele frequency to each record." << endl
             << "Uses:" << endl
             << "AO (alternate observation count), DP (read depth at the locus)" << endl;
        cerr << endl << "Type: statistics" << endl << endl;
        return 1;
    }
  }

    VariantCallFile variantFile;
    if (argc == 1) {
        variantFile.open(std::cin);
    } else {
        string filename = argv[argc-1];
        variantFile.open(filename);
        if (!variantFile.is_open()) {
            cerr << "could not open " << filename << endl;
            return 1;
        }
    }

    Variant var(variantFile);

    variantFile.removeInfoHeaderLine("VAF");
    variantFile.addHeaderLine("##INFO=<ID=VAF,Number=A,Type=Float,Description=\"The percentage of ALT reads\">");
    variantFile.addHeaderLine("##FORMAT=<ID=VAF,Number=A,Type=Float,Description=\"The percentage of ALT reads\">");
    // write the new header
    cout << variantFile.header << endl;

    // print the records, filtering is done via the setting of varA's output sample names
    while (variantFile.getNextVariant(var)) {
        int refobs = 0;
        vector<int> altobs(var.alt.size(), 0);
        for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin();
             s != var.samples.end(); ++s) {
            map<string, vector<string> >& sample = s->second;
            int x;
            if (sample.find("DP") != sample.end()) {
                convert(sample["DP"].front(), x);
                refobs += x;
            }
            if (sample.find("AO") != sample.end()) {
                vector<string>& aos = sample["AO"];
                for (int i = 0; i != var.alt.size(); ++i) {
                    convert(aos[i], x);
                    altobs[i] += x;
                }
            }

            for (int i = 0; i < var.alt.size(); i++) {
                if (refobs == 0) {
                    var.info["VAF"].push_back(convert(0));
                } else {
                    std::stringstream ss;
                    ss << std::fixed << std::setprecision(3) << (double)altobs[i]/(double)refobs;
                    var.info["VAF"].push_back(ss.str());
                    var.addFormatField("VAF");
                    sample["VAF"].push_back(ss.str());
                }
            }
            cout << var << endl;
            sample["VAF"].clear();
        }
    }

    return 0;

}
