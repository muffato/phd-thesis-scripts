//package org.ensembl;

import org.ensembl.registry.Registry;
import org.ensembl.compara.driver.ComparaDriver;
import org.ensembl.compara.datamodel.GenomeDB;
import org.ensembl.datamodel.FeaturePair;
import org.ensembl.datamodel.Gene;
import org.ensembl.datamodel.Location;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.CoreDriver;

import java.text.ParseException;
import java.util.List;

class Test {


	public static void main(String[] args) throws AdaptorException, ParseException {
		System.out.println("kikooooo la compagnie !");
		Registry registry = Registry.createDefaultRegistry();
		CoreDriver coreDriver = registry.getGroup("human").getCoreDriver();
		List genes = coreDriver.getGeneAdaptor().fetch(new Location("chromosome:::"));
		System.out.println(genes.size());
		/*for (int j = 0; j < genes.size(); j++) {
		  Gene g = (Gene) genes.get(j);
		  System.out.println(g.getAccessionID());
		}*/
		System.out.println("Vous me voyez ?");
		
		ComparaDriver compara = registry.getGroup("compara").getComparaDriver();

		List gdbs = compara.getGenomeDBAdaptor().fetch();
		for (int i = 0; i < gdbs.size(); i++) {
			GenomeDB gdb =(GenomeDB) gdbs.get(i);
			String targetSpecies = (String) gdb.getName();
			// if we were only interested in hits in chicken...
			// List hits = compara.getMemberAdaptor().fetch("Homo sapiens", new String[]{"ENSG00000118961"}, "Gallus gallus");
			List hits = compara.getMemberAdaptor().fetch("Homo sapiens", new String[]{"ENSG00000118961"}, targetSpecies);
			System.out.print(targetSpecies);
			for (int j = 0; j < hits.size(); j++) {
				FeaturePair p = (FeaturePair) hits.get(j);
				// hitDisplayName is the orthologue gene name; begins with
				// "ENS" if ensembl predicted gene otherwise name is assigned
				// by another group.
				System.out.print("\t" +p.getHitDisplayName());
			}
			System.out.println();
		}
	}
}
