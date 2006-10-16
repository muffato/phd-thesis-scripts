package orthos;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Map;
import java.util.zip.GZIPInputStream;

import org.ensembl.compara.driver.ComparaDriver;
import org.ensembl.driver.AdaptorException;
import org.ensembl.registry.Registry;

import utils.MyCombinator;
import utils.MyVars;


public class OrthosDownloader {

	private Hashtable<String, ArrayList<String>> dicGenomes = new Hashtable<String, ArrayList<String>>();
	private Hashtable<String, MyCombinator<String>> dicGroups = new Hashtable<String, MyCombinator<String>>();
	private Hashtable<String, Hashtable<String, ArrayList<String[]>>> dicOrthos = new Hashtable<String, Hashtable<String, ArrayList<String[]>>>();
	private Hashtable<String, ArrayList<String>> dicConstituants = new Hashtable<String, ArrayList<String>>();

	public static String
	    FETCH_HOMOLOGIES_BY_QUERY_SPECIES_STABLE_ID_AND_HIT_SPECIES =
	      "select "+
	      "  member1.stable_id, member2.stable_id " +
	      /*"  member1.member_id, "+
	      "  member1.chr_start, member1.chr_end, member1.chr_name, member1.stable_id, "+
	      "  member2.chr_start, member2.chr_end, member2.chr_name, member2.stable_id "+*/
	      "from "+
	      "  genome_db genome_db1, "+
	      "  member member1, "+
	      "  homology_member homology_member1, "+
	      "  homology_member homology_member2, "+
	      "  member member2, "+
	      "  genome_db genome_db2 "+
	      "where  "+
	      "  genome_db1.name = ? and "+ //1 - query species name
	      "  genome_db1.genome_db_id = member1.genome_db_id and "+
	      "  member1.member_id = homology_member1.member_id and "+
	      "  homology_member1.homology_id = homology_member2.homology_id and "+
	      "  homology_member2.member_id = member2.member_id and "+
	      "  member2.member_id != member1.member_id and "+
	      "  member2.genome_db_id = genome_db2.genome_db_id and "+
	      "  genome_db2.name = ? ";

	
	public OrthosDownloader() {
		
		// On charge les genomes
		for (Map.Entry<String, String> names: MyVars.organisms.entrySet()) {
			String simpleName = names.getKey();
			//loadSpecies(simpleName);
			dicConstituants.put(simpleName, new ArrayList<String>(Arrays.asList(new String[]{simpleName})));
		}
		
		downLoadOrthos("Human", "Chimp");
		
		// On telecharge tous les orthologues entre chaque paire d'espece
		for (Map.Entry<String, String> names1: MyVars.organisms.entrySet()) {
			String simpleName1 = names1.getKey();
			Hashtable<String, ArrayList<String[]>> hash = new Hashtable<String, ArrayList<String[]>>();
			for (Map.Entry<String, String> names2: MyVars.organisms.entrySet()) {
				String simpleName2 = names2.getKey();
				if (simpleName1 == simpleName2)
					continue;
				hash.put(simpleName2, downLoadOrthos(simpleName1, simpleName2));
			}
			dicOrthos.put(simpleName1, hash);
		}
		
		// On n'a plus besoin des genomes
		//dicGenomes.clear();
		
		buildGroup("tmp1", "Human", "Chimp");
		buildGroup("Primates", "tmp1", "Macaque");
		
		buildGroup("Rodents", "Mouse", "Rat");

		buildGroup("tmp2", "Primates", "Rodents");
		
		buildGroup("tmp3", "Dog", "Cow");
		
		buildGroup("Eutherians", "tmp2", "tmp3");

		buildGroup("Mammals", "Eutherians", "Opossum");
		
		buildGroup("Amniotes", "Mammals", "Chicken");
		
		buildGroup("Tetrapodes", "Amniotes", "Frog");
	}
	

	private ArrayList<String[]> downLoadOrthos(String esp1, String esp2) {


		try {
			System.out.print("Telechargement des orthologues entre " + esp1 + " et " + esp2 + " ... ");
			ComparaDriver compara = Registry.createDefaultRegistry().getGroup("compara").getComparaDriver();
			/*List gdbs = compara.getGenomeDBAdaptor().fetch();
			for(Object obj: gdbs)
				System.out.println(obj);*/
			String fullName1 = MyVars.organisms.get(esp1);
			String fullName2 = MyVars.organisms.get(esp2);
			//ArrayList<String> genome = dicGenomes.get(esp1);
			ArrayList<String[]> orthos = new ArrayList<String[]>();

			System.out.print("Ready ? ");
			
			
			
			
			
		    Connection conn = null;
		    StringBuffer sql;
		    ResultSet resultSet;
		    PreparedStatement statement = null;
		    String statementString;

		    statementString = FETCH_HOMOLOGIES_BY_QUERY_SPECIES_STABLE_ID_AND_HIT_SPECIES;

		    try {

		      sql = new StringBuffer();
		      String sqlString;

		      //get a partial prepared statement from the static string at the top.
		      
		      sql.append(statementString);

		      sql.append(" group by member2.member_id ");
		        
		      sqlString = sql.toString();
		      conn = compara.getConnection();
				System.out.print("getConnection ");
		      statement = conn.prepareStatement(sqlString);
				System.out.print("prepareStatement ");

		      statement.setObject(1, fullName1);
		      statement.setObject(2, fullName2);
				System.out.print("statement OK {"+statement.toString()+"} ");
		      
		      resultSet = statement.executeQuery();
				System.out.print("query ");
		      
		      while ( resultSet.next() ) {
		    	    String name   = resultSet.getString(1);
		    	    String hid    = resultSet.getString(2);
					orthos.add(new String[]{name, hid});
					
		      }
				System.out.print("END ");
		    	conn.close();


		    } catch ( SQLException exception ) {
		      throw new AdaptorException(exception.getMessage(), exception);
		    } finally {
		    }

			
			
			
			
			
			
			
			
			
			
			//List hits = compara.getMemberAdaptor().fetch(fullName1, genome.toArray(new String[nbGenes]), fullName2);
			/*List hits = compara.getMemberAdaptor().fetch(fullName1, new Location("chromosome:"), fullName2);
			System.out.print(hits.size());	
			for (int j = 0; j < hits.size(); j++) {
				String[] tab = new String[2];
				FeaturePair p = (FeaturePair) hits.get(j);
				tab[0] = p.getDisplayName();
				tab[1] = p.getHitDisplayName();
				orthos[j] = tab;
			}*/

			/*for(int i=0; i<nbGenes; i++) {
				
				//if (i % 100 == 0)
				System.out.println(i);
				
				String gene = genome.get(i);
				List hits = compara.getMemberAdaptor().fetch(fullName1, new String[]{gene}, fullName2);
				
				String[] tab = new String[hits.size()+1];
				tab[0] = gene;
				for (int j = 0; j < hits.size(); j++) {
					FeaturePair p = (FeaturePair) hits.get(j);
					tab[j] = p.getHitDisplayName();
				}
				orthos[i] = tab;
				
			}*/
			System.out.println(orthos.size() + " OK");
			return orthos;
		} catch (AdaptorException e) {
			System.out.println("Echec");
			return null;
		}
	}

	private void buildGroup(String node, String branch1, String branch2) {
		MyCombinator<String> tab = new MyCombinator<String>();
		
		MyCombinator<String> lst = dicGroups.get(branch1);
		if (lst == null) {
			for (String gene : dicGenomes.get(branch1)) {
				ArrayList<String> tmp = new ArrayList<String>(1);
				tmp.add(gene);
				tab.addElement(tmp);
			}
		} else {
			for (ArrayList<String> elt : lst.liste)
				tab.addElement(elt);
		}
		
		lst = dicGroups.get(branch2);
		if (lst == null) {
			for (String gene : dicGenomes.get(branch2)) {
				ArrayList<String> tmp = new ArrayList<String>(1);
				tmp.add(gene);
				tab.addElement(tmp);
			}
		} else {
			for (ArrayList<String> elt : lst.liste)
				tab.addElement(elt);
		}
			
		for(String esp1: dicConstituants.get(branch1)) {
			for(String esp2: dicConstituants.get(branch1)) {
				
				ArrayList<String[]> orthos12 = dicOrthos.get(esp1).get(esp2);
				ArrayList<String[]> orthos21 = dicOrthos.get(esp2).get(esp1);

				for(String[] link: orthos12)
					tab.addLink(link);
				for(String[] link: orthos21)
					tab.addLink(link);
				
			}
			
			dicGroups.put(node, tab);
		}
		
		ArrayList<String> cst = new ArrayList<String>();
		cst.addAll(dicConstituants.get(branch1));
		cst.addAll(dicConstituants.get(branch2));
		dicConstituants.put(node, cst);
	}

	private void loadSpecies(String name) {
		ArrayList<String> tab = new ArrayList<String>();
		String fileName = MyVars.downloadFolder + "genes/genes" + name + ".tsv.gz";
		System.out.print("Chargement de " + name + " ... ");
		try {
			BufferedReader input = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fileName))));
			String ligne = null;
			while ((ligne = input.readLine()) != null) {
				String[] t = ligne.split("\t");
				tab.add(t[t.length-1]);
			}
			input.close();
			dicGenomes.put(name, tab);
			
			System.out.println("OK");
		} catch (IOException e) {
			System.out.println("Echec");
		}
	}
	
	public static void main(String[] args) {
		new OrthosDownloader();
	}
	
}
