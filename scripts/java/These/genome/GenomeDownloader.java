package genome;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.zip.GZIPOutputStream;

import org.ensembl.datamodel.Gene;
import org.ensembl.datamodel.Location;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.CoreDriver;
import org.ensembl.driver.CoreDriverFactory;
import org.ensembl.driver.GeneAdaptor;

import utils.MyVars;


public class GenomeDownloader {
	
	private static void downloadGenome(String species) {
		
		System.out.print("Telechargement du genome de '" + species + "' ");

		// On ouvre les fichiers pour ecrire les genomes (en GZIP)
		PrintStream output = null;
		PrintStream outputFull = null;
		try {
			// Le fichier avec les genes des vrais chromosomes
			String fileName = MyVars.genesDownloadFolder + "genes" + species + ".tsv.gz"; 
			output = new PrintStream(new GZIPOutputStream(new FileOutputStream(fileName)), true);
			// Le fichier avec le genome complet
			fileName = MyVars.fullGenesDownloadFolder + "genes" + species + ".tsv.gz"; 
			outputFull = new PrintStream(new GZIPOutputStream(new FileOutputStream(fileName)), true);
		} catch (IOException e) {
			System.out.println("Echec");
			return;
		}

		// On telecharge
		try {
			// On se connecte a la vase
			String databaseName = MyVars.organisms.get(species).toLowerCase().replaceAll(" ", "_") + "_core";
			CoreDriver coreDriver = CoreDriverFactory.createCoreDriverUsingDatabasePrefix("ensembldb.ensembl.org", 3306, databaseName, "anonymous", null);
			GeneAdaptor ad = coreDriver.getGeneAdaptor();
			System.out.print(".");
			
			// On initialise
			String[] filtres = MyVars.badChromosomes.get(species);
			int nbGenes = 0;
			int nbGenesFiltres = 0;
			Iterator it = ad.fetchIterator(true);
			System.out.print(".");
			
			// C'est parti pour le show
			while (it.hasNext()) {
				
				// On recupere le gene
				Gene g = (Gene) it.next();
				Location loc = g.getLocation();
				String ligne = loc.getSeqRegionName() + "\t" + loc.getStart() + "\t" + loc.getEnd() + "\t" + loc.getStrand() + "\t" + g.getAccessionID();
				
				// On l'imprime dans la liste complete
				nbGenes++;
				outputFull.println(ligne);
				
				// Fait-il partie de la liste des chanceux
				boolean ok = true;
				for(String filtre: filtres) {
					if (loc.getSeqRegionName().matches(filtre)) {
						ok = false;
						break;
					}
				}
				
				// Oui, on l'affiche encore
				if (ok) {
					nbGenesFiltres++;
					output.println(ligne);
				}
			}
			
			// C'est fini
			System.out.println(" OK (" + nbGenesFiltres + "/" + nbGenes + " genes)");
			output.close();
			outputFull.close();
			
		} catch (AdaptorException e) {
			System.out.println(" Echec [" + e + "]");
		}
		
	}
	
	public static void main(String[] args) {
		for(String simpleName: MyVars.organisms.keySet())
			downloadGenome(simpleName);
	}
	
}
