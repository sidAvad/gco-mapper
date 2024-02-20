import stdpopsim

def main():
    
    species = stdpopsim.get_species("HomSap")
    contig = species.get_contig("chr1", genetic_map="HapMapII_GRCh38")
    contig.gene_conversion_fraction = 0.02 # From Browning, 2023
    contig.gene_conversion_length = 300 # From Browning, 2023
    engine = stdpopsim.get_engine("msprime") # default model is Hudson
    model = species.get_demographic_model("OutOfAfrica_4J17")
    samples={"CEU":10000,"YRI":10000}# Ouput sample size & population
    ts = engine.simulate(model, contig, samples)
    ts.dump("simulation-source-stdpopsim-20k.trees")
    
    return()

if __name__ == "__main__":
    main()
