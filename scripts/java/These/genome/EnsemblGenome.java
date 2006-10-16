package genome;

import java.util.Iterator;
import java.util.List;

import org.ensembl.datamodel.Gene;
import org.ensembl.datamodel.Location;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.CoreDriver;
import org.ensembl.driver.CoreDriverFactory;
import org.ensembl.driver.GeneAdaptor;


/**
 * Un genome d'Ensembl. On fournit au constructeur le nom de l'espece et le genome est telecharge
 */
public class EnsemblGenome {

	// La liste qui contiendra le genome
	private List genome = null;
	private String name = "";
	
	// Le constructeur
	public EnsemblGenome(String species) throws AdaptorException {
		
		System.out.print("Telechargement du genome de '" + species + "' ");

		try {
			CoreDriver coreDriver = CoreDriverFactory.createCoreDriverUsingDatabasePrefix("ensembldb.ensembl.org", 3306, species.toLowerCase().replaceAll(" ", "_") + "_core", "anonymous", null);
			GeneAdaptor ad = coreDriver.getGeneAdaptor();
			System.out.print(".");
			//Location loc = null;
			int nb = 0;
			Iterator it = ad.fetchIterator(true);
			while (it.hasNext()) {
				nb++;
				//it.next();
				System.err.println(species + " - " + ((Gene) it.next()).getLocation().getSeqRegionName());
			}
			System.out.println(nb);
			/*try {
				loc = new Location("chromosome");
				genome = ad.fetch(loc);
				
			} catch (NullPointerException e1) {
				try {
					loc = new Location("group");
					genome = ad.fetch(loc);
				} catch (NullPointerException e2) {
					loc = new Location("scaffold");
					genome = ad.fetch(loc);
				}
			}
			System.out.print(".");
		} catch (ParseException e) {
			System.out.println(" Echec");
			throw new AdaptorException();*/
		} catch (NullPointerException e3) {
			System.out.println(" Echec");
			throw new AdaptorException();
		} catch (AdaptorException e) {
			System.out.println(" Echec");
			throw e;
		}
		
		name = species;
		//System.out.println(" OK (" + genome.size() + " genes)");
	}
	
	// Export du genome dans un fichier tabulaire avec compression GZIP au vol
	public void printData(String fileName) {
		System.out.print("Enregistrement du genome de '" + name + "' ... ");
		/*PrintStream output = null;
		try {
			output = new PrintStream(new GZIPOutputStream(new FileOutputStream(fileName)), true);
		} catch (IOException e) {
			System.out.println("Echec");
			return;
		}*/
		int nbGenes = 0;
		Gene g = null;
		Location loc = null;
		for (int j = 0; j < genome.size(); j++) {
		  g = (Gene) genome.get(j);
		  loc = g.getLocation();
		  if (loc.getSeqRegionName() == null)
			  continue;
		  if (loc.getSeqRegionName().startsWith("M"))
			  continue;
		  if (loc.getSeqRegionName().startsWith("U"))
			  continue;
		  if (loc.getSeqRegionName().startsWith("c"))
			  continue;
		  if (loc.getSeqRegionName().startsWith("E"))
			  continue;
		  if (loc.getSeqRegionName().endsWith("h"))
			  continue;
		  if (loc.getSeqRegionName().endsWith("_random"))
			  continue;
		  //output.println(loc.getSeqRegionName() + "\t" + loc.getStart() + "\t" + loc.getEnd() + "\t" + loc.getStrand() + "\t" + g.getAccessionID());
		  nbGenes ++;
		}
		//output.close();
		System.out.println("OK (" + nbGenes + " genes gardes)");
	}
}
