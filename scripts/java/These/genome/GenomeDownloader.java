package genome;

import java.io.File;
import java.util.Hashtable;
import java.util.Map;

import org.ensembl.driver.AdaptorException;


public class GenomeDownloader {

	// Les genomes a telecharger
	final static public Hashtable<String, String> organisms = buildOrganismsTable();
	
	final static public String downloadFolder = "/users/ldog/muffato/work/data/";
	
	private static Hashtable<String, String> buildOrganismsTable() {
		Hashtable<String, String> table = new Hashtable<String, String>();
		table.put("Human", "Homo Sapiens");
		table.put("Chimp", "Pan Troglodytes");
		table.put("Macaque", "Macaca Mulatta");
		table.put("Mouse", "Mus Musculus");
		table.put("Rat", "Rattus Norvegicus");
		table.put("Cow", "Bos Taurus");
		table.put("Dog", "Canis Familiaris");
		table.put("Opossum", "Monodelphis Domestica");
		table.put("Chicken", "Gallus gallus");
		table.put("Frog", "Xenopus Tropicalis");
		table.put("Fugu", "Fugu Rubripes");
		table.put("Tetraodon", "Tetraodon Nigroviridis");
		table.put("Stickleback", "Gasterosteus Aculeatus");
		table.put("Zebrafish", "Danio Rerio");
		table.put("C_elegans", "Caenorhabditis elegans");
		table.put("C_intestinalis", "Ciona Intestinalis");
		table.put("Fruitfly", "Drosophila Melanogaster");
		table.put("Mosquito", "Anopheles Gambiae");
		table.put("Yeast", "Saccharomyces Cerevisiae");
		return table;
	}
	
	public static void main(String[] args) {
		for (Map.Entry<String, String> names: organisms.entrySet()) {
			String simpleName = names.getKey();
			String fullName = names.getValue();
			try {
				if ((new File(downloadFolder + "genes/genes" + simpleName + ".tsv.gz")).exists())
					continue;
				(new EnsemblGenome(fullName)).printData(downloadFolder + "genes/genes" + simpleName + ".tsv.gz");
			} catch (AdaptorException e) {
			}
		}
	}
	
}
