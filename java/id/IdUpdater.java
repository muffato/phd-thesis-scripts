package id;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.List;

import org.ensembl.datamodel.StableIDEvent;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.CoreDriver;
import org.ensembl.driver.CoreDriverFactory;
import org.ensembl.driver.StableIDEventAdaptor;

public class IdUpdater {

	CoreDriver coreDriver;
	
	public IdUpdater() {
		try {
			coreDriver = CoreDriverFactory.createCoreDriverUsingDatabasePrefix("ensembldb.ensembl.org", 3306, "canis_familiaris_core", "anonymous", null);
			//getNewID("GSTENG00022540001");
		} catch (AdaptorException e) {
			return;
		};
		BufferedReader input = new BufferedReader(new InputStreamReader(System.in));
		String ligne = null;
		try {
			while ((ligne = input.readLine()) != null) {
				String[] tab = ligne.split("[^A-Z0-9]");
				for(String id: tab) {
					getNewID(id);
				}
			}
		} catch (IOException e) {
			return;
		}
	}

	private void getNewID(String ancientID) throws AdaptorException {
		StableIDEventAdaptor adaptor = coreDriver.getStableIDEventAdaptor();

		// Find stableIDs in the current release that relate to the geneStableID
		List relatedIDs = adaptor.fetch(ancientID);
		System.out.print(ancientID + ": " + relatedIDs.size() + " ");
		for (Iterator iter = relatedIDs.iterator(); iter.hasNext();) {
			String relatedID = ((StableIDEvent) iter.next()).getStableID();
			System.out.print(relatedID + " ");
		}
		System.out.println();
	}

	public static void main(String[] args) {
		new IdUpdater();
	}

}
