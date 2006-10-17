package utils;

import java.util.HashMap;

public class MyVars {

	// Les genomes a telecharger
	final static public HashMap<String, String> organisms = buildOrganismsTable();
	final static public HashMap<String, String[]> badChromosomes = buildFiltersTable();
	
	final static public String downloadFolder = "/users/ldog/muffato/work/data/";
	final static public String genesDownloadFolder = downloadFolder + "genes/";
	final static public String fullGenesDownloadFolder = genesDownloadFolder + "full/";
	final static public String orthosDownloadFolder = downloadFolder + "orthos/";
	
	private static HashMap<String, String> buildOrganismsTable() {
		HashMap<String, String> table = new HashMap<String, String>();
		
		table.put("Human", "Homo Sapiens");
		table.put("Chimp", "Pan Troglodytes");
		table.put("Macaque", "Macaca Mulatta");
		table.put("Mouse", "Mus Musculus");
		table.put("Rat", "Rattus Norvegicus");
		table.put("Rabbit", "Oryctolagus Cuniculus");
		table.put("Armadillo", "Dasypus novemcinctus");
		table.put("Elephant", "Loxodonta africana");
		table.put("Hedgehog", "Echinops telfairi");
		table.put("Cow", "Bos Taurus");
		table.put("Dog", "Canis Familiaris");
		table.put("Opossum", "Monodelphis Domestica");
		table.put("Chicken", "Gallus gallus");
		table.put("Frog", "Xenopus Tropicalis");
		table.put("Fugu", "Fugu Rubripes");
		table.put("Tetraodon", "Tetraodon Nigroviridis");
		table.put("Stickleback", "Gasterosteus Aculeatus");
		table.put("Medaka", "Oryzias Latipes");
		table.put("Zebrafish", "Danio Rerio");
		table.put("C_elegans", "Caenorhabditis elegans");
		table.put("C_intestinalis", "Ciona Intestinalis");
		table.put("Fruitfly", "Drosophila Melanogaster");
		table.put("Mosquito", "Anopheles Gambiae");
		table.put("Yeast", "Saccharomyces Cerevisiae");
		
		return table;
	}
	
	private static HashMap<String, String[]> buildFiltersTable() {
		HashMap<String, String[]> table = new HashMap<String, String[]>();

		table.put("Human", new String[] {"c22_H2", "c5_H2", "c6_COX", "c6_QBL", "MT", "NT_.*"});
		table.put("Chimp", new String[] {"MT", "Un", ".*_random"});
		table.put("Macaque", new String[] {"1099.*", "MT"});
		table.put("Mouse", new String[] {"MT", "NT_.*"});
		table.put("Rat", new String[] {"MT"});
		table.put("Cow", new String[] {"ChrUn.*"});
		table.put("Dog", new String[] {"MT", "Un"});

		table.put("Rabbit", new String[] {"GeneScaffold_.*", "scaffold.*"});
		table.put("Armadillo", new String[] {"GeneScaffold_.*", "scaffold.*"});
		table.put("Elephant", new String[] {"GeneScaffold_.*", "scaffold.*"});
		table.put("Hedgehog", new String[] {"GeneScaffold_.*", "scaffold.*"});

		table.put("Opossum", new String[] {"MT", "Un"});
		table.put("Chicken", new String[] {".*_random", "E22C19W28", "E26C13", "E50C23", "MT", "Un"});
		table.put("Frog", new String[] {"scaffold_.*"});

		table.put("Tetraodon", new String[] {".*_random"});
		table.put("Stickleback", new String[] {"MT", "scaffold_.*"});
		table.put("Medaka", new String[] {"scaffold.*", "ultracontig.*"});
		table.put("Zebrafish", new String[] {"Zv6_NA.*", "Zv6_scaffold.*"});
		table.put("Fugu", new String[] {"scaffold_.*"});

		table.put("C_elegans", new String[] {"MtDNA"});
		table.put("C_intestinalis", new String[] {"scaffold_.*"});
		table.put("Fruitfly", new String[] {"U", ".h"});
		table.put("Mosquito", new String[] {"UNKN"});
		table.put("Yeast", new String[] {"Mito"});
		
		return table;
	}
}
