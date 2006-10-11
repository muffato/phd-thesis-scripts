package data;

import java.io.File;

import org.ensembl.driver.AdaptorException;

import data.genes.EnsemblGenome;

public class DataDownloader {

	final static private String[][] organisms = {
		{"Human", "Homo Sapiens"},
		{"Stickleback", "Gasterosteus Aculeatus"},
		{"Frog", "Xenopus Tropicalis"},
		{"Chimp", "Pan Troglodytes"},
		{"Macaque", "Macaca Mulatta"},
		{"Mouse", "Mus Musculus"},
		{"Rat", "Rattus Norvegicus"},
		{"Cow", "Bos Taurus"},
		{"Dog", "Canis Familiaris"},
		{"Opossum", "Monodelphis Domestica"},
		{"Chicken", "Gallus gallus"},
		{"Fugu", "Fugu Rubripes"},
		{"Tetraodon", "Tetraodon Nigroviridis"},
		{"Zebrafish", "Danio Rerio"},
		{"C_elegans", "Caenorhabditis elegans"},
		{"C_savignyi", "Ciona Savignyi"},
		{"C_intestinalis", "Ciona Intestinalis"},
		{"Fruitfly", "Drosophila Melanogaster"},
		{"Mosquito", "Anopheles Gambiae"},
		{"Yeast", "Saccharomyces Cerevisiae"},
	};
	final static private String downloadFolder = "/users/ldog/muffato/work/data/";
	
	public DataDownloader() {
		
		for (String[] names: organisms) {
			String simpleName = names[0];
			String fullName = names[1];
			try {
				if ((new File(downloadFolder + "genes/genes" + simpleName + ".tsv.gz")).exists())
					continue;
				(new EnsemblGenome(fullName)).printData(downloadFolder + "genes/genes" + simpleName + ".tsv.gz");
			} catch (AdaptorException e) {
			}
		}
	}
	
	public static void main(String[] args) {
		new DataDownloader();
	}
	
}
